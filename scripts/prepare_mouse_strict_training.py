#!/usr/bin/env python3
"""严格合并 NCBI GRCm39、GENCODE M39 Basic 和 IMGT 小鼠剪接标签。"""

from __future__ import annotations

import argparse
from pathlib import Path
import shutil
from types import SimpleNamespace

from openspliceai.create_data.gencode_splice import read_gencode_basic_splice_sites
from openspliceai.create_data.imgt_genome_match import (
    filter_imgt_records_by_grcm39,
    remove_imgt_sites_on_chromosomes,
)
from openspliceai.create_data.immune_splice import (
    find_missing_imgt_sites,
    imgt_records_to_data,
    merge_data,
    read_datafile_h5,
    read_imgt_splice_records,
    read_ncbi_immune_splice_sites,
    split_imgt_records_by_gene,
    write_audit_report,
    write_datafile_h5,
)
from openspliceai.create_data.immune_training import (
    combine_data,
    coverage_to_report,
    datafile_positive_sites,
    find_sites_missing_from_baseline,
    splice_sites_to_data,
    split_splice_sites_by_gene,
)

AUDIT_FILENAME = "ncbi_gencode_imgt_splice_audit.json"


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--ncbi-gff", type=Path, required=True)
    p.add_argument("--ncbi-fasta", type=Path, required=True)
    p.add_argument("--gencode-gtf", type=Path, required=True)
    p.add_argument("--imgt-flat", type=Path, required=True)
    p.add_argument("--base-data-dir", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--training-context-radius", type=int, default=1000)
    p.add_argument("--match-radius", type=int, default=20)
    p.add_argument("--imgt-match-radius", type=int, default=30)
    p.add_argument("--validation-ratio", type=float, default=0.10)
    p.add_argument("--include-nonfunctional", action="store_true")
    p.add_argument("--allow-nonunique-imgt", action="store_true")
    p.add_argument("--skip-dataset", action="store_true")
    return p.parse_args()


def exclude_sites(sites, chromosomes):
    blocked = set(str(x) for x in chromosomes)
    kept, removed = [], []
    for site in sites:
        (removed if site.seqid in blocked else kept).append(site)
    return kept, removed


def rename_imgt_report(report):
    report = dict(report)
    report["prior_reference_site_count"] = report.pop("ncbi_site_count")
    report["covered_by_prior_sources"] = report.pop("covered_by_ncbi")
    for item in report["sites"]:
        item["covered_by_prior_sources"] = item.pop("covered_by_ncbi")
        item["best_prior_context_identity"] = item.pop(
            "best_ncbi_context_identity"
        )
    return report


def main():
    a = parse_args()
    required = [
        "datafile_train.h5",
        "datafile_validation.h5",
        "datafile_test.h5",
    ]
    missing = [x for x in required if not (a.base_data_dir / x).exists()]
    if missing:
        raise FileNotFoundError(f"基线缺少：{', '.join(missing)}")
    if a.output_dir.resolve() == a.base_data_dir.resolve():
        raise ValueError("输出目录不能覆盖基线目录")
    a.output_dir.mkdir(parents=True, exist_ok=True)

    print("[1/8] 读取基线标签")
    base_train = read_datafile_h5(a.base_data_dir / "datafile_train.h5")
    base_val = read_datafile_h5(a.base_data_dir / "datafile_validation.h5")
    base_test = read_datafile_h5(a.base_data_dir / "datafile_test.h5")
    test_chromosomes = sorted(set(str(x) for x in base_test[1]))
    baseline_sites = datafile_positive_sites(
        combine_data(base_train, base_val),
        context_radius=max(30, a.match_radius),
    )

    print("[2/8] NCBI GRCm39 免疫位点")
    ncbi_sites = read_ncbi_immune_splice_sites(
        a.ncbi_gff,
        a.ncbi_fasta,
        context_radius=a.training_context_radius,
    )
    ncbi_missing_all, ncbi_cov = find_sites_missing_from_baseline(
        ncbi_sites,
        baseline_sites,
        context_radius=a.match_radius,
        identity_threshold=1.0,
    )
    ncbi_missing, ncbi_test_removed = exclude_sites(
        ncbi_missing_all,
        test_chromosomes,
    )

    print("[3/8] GENCODE M39 Basic level 1/2 规范位点")
    gencode_sites, gencode_parse = read_gencode_basic_splice_sites(
        a.gencode_gtf,
        a.ncbi_fasta,
        context_radius=a.training_context_radius,
        expected_release="M39",
        expected_assembly="GRCm39",
        accepted_levels=(1, 2),
    )
    gencode_missing_all, gencode_cov = find_sites_missing_from_baseline(
        gencode_sites,
        baseline_sites + ncbi_missing_all,
        context_radius=a.match_radius,
        identity_threshold=1.0,
    )
    gencode_missing, gencode_test_removed = exclude_sites(
        gencode_missing_all,
        test_chromosomes,
    )

    print("[4/8] IMGT functional 明确位点")
    imgt_records = read_imgt_splice_records(
        a.imgt_flat,
        functional_only=not a.include_nonfunctional,
        mouse_only=True,
    )
    imgt_mapped, imgt_map_report = filter_imgt_records_by_grcm39(
        imgt_records,
        a.ncbi_fasta,
        context_radius=a.imgt_match_radius,
        require_unique=not a.allow_nonunique_imgt,
    )
    imgt_mapped, imgt_test_removed = remove_imgt_sites_on_chromosomes(
        imgt_mapped,
        test_chromosomes,
    )
    imgt_missing, imgt_dup_report = find_missing_imgt_sites(
        imgt_mapped,
        ncbi_sites + gencode_sites,
        context_radius=a.match_radius,
        identity_threshold=1.0,
    )
    imgt_dup_report = rename_imgt_report(imgt_dup_report)

    print("[5/8] 稳定拆分 train/validation")
    seed = "openspliceai-strict-m39-v1"
    ncbi_train, ncbi_val = split_splice_sites_by_gene(
        ncbi_missing,
        validation_ratio=a.validation_ratio,
        seed=seed,
    )
    gencode_train, gencode_val = split_splice_sites_by_gene(
        gencode_missing,
        validation_ratio=a.validation_ratio,
        seed=seed,
    )
    imgt_train, imgt_val = split_imgt_records_by_gene(
        imgt_missing,
        validation_ratio=a.validation_ratio,
        seed=seed,
    )

    report = {
        "policy": {
            "priority": "NCBI GRCm39 -> GENCODE M39 Basic -> IMGT",
            "assembly": "GRCm39 only; no liftOver",
            "gencode": "Basic, level 1/2, canonical motifs only",
            "imgt": "functional; exact full-length unique GRCm39 context",
            "test_set": "original test unchanged; mapped test-chromosome supplements removed",
            "test_chromosomes": test_chromosomes,
        },
        "baseline_positive_sites": len(baseline_sites),
        "ncbi": {
            "sites": len(ncbi_sites),
            "already_covered": len(ncbi_sites) - len(ncbi_missing_all),
            "removed_test_chromosome": len(ncbi_test_removed),
            "supplemented": len(ncbi_missing),
            "details": coverage_to_report(ncbi_cov),
        },
        "gencode": {
            "parse": gencode_parse,
            "sites": len(gencode_sites),
            "already_covered": len(gencode_sites) - len(gencode_missing_all),
            "removed_test_chromosome": len(gencode_test_removed),
            "supplemented": len(gencode_missing),
            "details": coverage_to_report(gencode_cov),
        },
        "imgt": {
            "genome_match": imgt_map_report,
            "removed_test_chromosome": imgt_test_removed,
            "duplicate_filter": imgt_dup_report,
        },
    }
    write_audit_report(a.output_dir / AUDIT_FILENAME, report)

    print("[6/8] 写入合并 datafile")
    for split, base, ns, gs, ir in (
        ("train", base_train, ncbi_train, gencode_train, imgt_train),
        ("validation", base_val, ncbi_val, gencode_val, imgt_val),
    ):
        nd = splice_sites_to_data(ns, name_prefix="NCBI-IMMUNE")
        gd = splice_sites_to_data(gs, name_prefix="GENCODE-M39-BASIC")
        idata = imgt_records_to_data(ir)
        supplement = combine_data(nd, gd, idata)
        write_datafile_h5(a.output_dir / f"ncbi_supplement_{split}.h5", nd)
        write_datafile_h5(a.output_dir / f"gencode_supplement_{split}.h5", gd)
        write_datafile_h5(a.output_dir / f"imgt_supplement_{split}.h5", idata)
        write_datafile_h5(
            a.output_dir / f"datafile_{split}.h5",
            merge_data(base, supplement),
        )
    shutil.copy2(
        a.base_data_dir / "datafile_test.h5",
        a.output_dir / "datafile_test.h5",
    )

    print("[7/8] 生成 dataset")
    if not a.skip_dataset:
        from openspliceai.create_data import create_dataset

        create_dataset.create_dataset(
            SimpleNamespace(
                output_dir=str(a.output_dir),
                chr_split="train-test",
                biotype="protein-coding",
            )
        )

    print("[8/8] 完成")
    print(f"审计报告：{a.output_dir / AUDIT_FILENAME}")
    print(f"NCBI 补充：{len(ncbi_missing)}")
    print(f"GENCODE 补充：{len(gencode_missing)}")
    print(f"IMGT 补充：{imgt_dup_report['supplemented_from_imgt']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
