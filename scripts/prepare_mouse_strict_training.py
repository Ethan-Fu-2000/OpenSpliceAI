#!/usr/bin/env python3
"""严格合并 NCBI GRCm39 与 IMGT 小鼠剪接标签。"""

from __future__ import annotations

import argparse
from pathlib import Path
import shutil
from types import SimpleNamespace

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

AUDIT_FILENAME = "ncbi_imgt_strict_splice_audit.json"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "以 NCBI GRCm39 为主体，只补充严格匹配 GRCm39 的 IMGT 明确位点。"
        )
    )
    parser.add_argument("--ncbi-gff", type=Path, required=True)
    parser.add_argument("--ncbi-fasta", type=Path, required=True)
    parser.add_argument("--imgt-flat", type=Path, required=True)
    parser.add_argument("--base-data-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--training-context-radius", type=int, default=1000)
    parser.add_argument("--match-radius", type=int, default=20)
    parser.add_argument("--imgt-match-radius", type=int, default=30)
    parser.add_argument("--validation-ratio", type=float, default=0.10)
    parser.add_argument("--include-nonfunctional", action="store_true")
    parser.add_argument("--allow-nonunique-imgt", action="store_true")
    parser.add_argument("--skip-dataset", action="store_true")
    return parser.parse_args()


def exclude_sites(sites, chromosomes):
    blocked = {str(value) for value in chromosomes}
    kept, removed = [], []
    for site in sites:
        (removed if site.seqid in blocked else kept).append(site)
    return kept, removed


def main() -> int:
    args = parse_args()
    required = [
        "datafile_train.h5",
        "datafile_validation.h5",
        "datafile_test.h5",
    ]
    missing = [
        name for name in required if not (args.base_data_dir / name).exists()
    ]
    if missing:
        raise FileNotFoundError(f"基线缺少：{', '.join(missing)}")
    if args.output_dir.resolve() == args.base_data_dir.resolve():
        raise ValueError("输出目录不能覆盖基线目录")
    if args.training_context_radius < args.imgt_match_radius:
        raise ValueError("training-context-radius 不能小于 imgt-match-radius")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    print("[1/7] 读取 NCBI 基线标签")
    base_train = read_datafile_h5(args.base_data_dir / "datafile_train.h5")
    base_validation = read_datafile_h5(
        args.base_data_dir / "datafile_validation.h5"
    )
    base_test = read_datafile_h5(args.base_data_dir / "datafile_test.h5")
    test_chromosomes = sorted(set(base_test[1]))
    baseline_sites = datafile_positive_sites(
        combine_data(base_train, base_validation),
        context_radius=max(30, args.match_radius),
    )

    print("[2/7] 读取 NCBI GRCm39 的 IG/TR 明确剪接位点")
    ncbi_sites = read_ncbi_immune_splice_sites(
        args.ncbi_gff,
        args.ncbi_fasta,
        context_radius=args.training_context_radius,
    )
    ncbi_missing_all, ncbi_coverage = find_sites_missing_from_baseline(
        ncbi_sites,
        baseline_sites,
        context_radius=args.match_radius,
        identity_threshold=1.0,
    )
    ncbi_missing, ncbi_test_removed = exclude_sites(
        ncbi_missing_all,
        test_chromosomes,
    )

    print("[3/7] 读取 IMGT functional 小鼠明确位点")
    imgt_records = read_imgt_splice_records(
        args.imgt_flat,
        functional_only=not args.include_nonfunctional,
        mouse_only=True,
    )

    print("[4/7] 将 IMGT 位点零错配定位到 GRCm39，并排除测试染色体")
    imgt_mapped, imgt_map_report = filter_imgt_records_by_grcm39(
        imgt_records,
        args.ncbi_fasta,
        context_radius=args.imgt_match_radius,
        require_unique=not args.allow_nonunique_imgt,
    )
    imgt_mapped, imgt_test_removed = remove_imgt_sites_on_chromosomes(
        imgt_mapped,
        test_chromosomes,
    )
    imgt_missing, imgt_duplicate_report = find_missing_imgt_sites(
        imgt_mapped,
        ncbi_sites,
        context_radius=args.match_radius,
        identity_threshold=1.0,
    )

    print("[5/7] 按基因稳定拆分 train/validation")
    seed = "openspliceai-strict-imgt-v2"
    ncbi_train, ncbi_validation = split_splice_sites_by_gene(
        ncbi_missing,
        validation_ratio=args.validation_ratio,
        seed=seed,
    )
    imgt_train, imgt_validation = split_imgt_records_by_gene(
        imgt_missing,
        validation_ratio=args.validation_ratio,
        seed=seed,
    )

    report = {
        "policy": {
            "priority": "NCBI GRCm39 baseline -> missing NCBI IG/TR -> IMGT",
            "assembly": "GCF_000001635.27 / GRCm39 only; no liftOver",
            "gencode": "not used",
            "imgt": (
                "mouse functional records; site-centred sequence must map "
                "full-length to GRCm39 with NM=0; unique mapping required by default"
            ),
            "test_set": (
                "原始 NCBI test 完全不变；映射到测试染色体的补充位点全部排除"
            ),
            "test_chromosomes": test_chromosomes,
        },
        "baseline_positive_sites": len(baseline_sites),
        "ncbi": {
            "sites": len(ncbi_sites),
            "already_covered": len(ncbi_sites) - len(ncbi_missing_all),
            "removed_test_chromosome": len(ncbi_test_removed),
            "supplemented": len(ncbi_missing),
            "details": coverage_to_report(ncbi_coverage),
        },
        "imgt": {
            "genome_match": imgt_map_report,
            "removed_test_chromosome": imgt_test_removed,
            "duplicate_filter": imgt_duplicate_report,
        },
    }
    write_audit_report(args.output_dir / AUDIT_FILENAME, report)

    print("[6/7] 写入 NCBI+IMGT 合并数据")
    for split, base, ncbi_split, imgt_split in (
        ("train", base_train, ncbi_train, imgt_train),
        (
            "validation",
            base_validation,
            ncbi_validation,
            imgt_validation,
        ),
    ):
        ncbi_data = splice_sites_to_data(
            ncbi_split,
            name_prefix="NCBI-IMMUNE",
        )
        imgt_data = imgt_records_to_data(imgt_split)
        supplement = combine_data(ncbi_data, imgt_data)
        write_datafile_h5(
            args.output_dir / f"ncbi_supplement_{split}.h5",
            ncbi_data,
        )
        write_datafile_h5(
            args.output_dir / f"imgt_supplement_{split}.h5",
            imgt_data,
        )
        write_datafile_h5(
            args.output_dir / f"immune_supplement_{split}.h5",
            supplement,
        )
        write_datafile_h5(
            args.output_dir / f"datafile_{split}.h5",
            merge_data(base, supplement),
        )

    shutil.copy2(
        args.base_data_dir / "datafile_test.h5",
        args.output_dir / "datafile_test.h5",
    )

    print("[7/7] 生成 dataset 并完成")
    if not args.skip_dataset:
        from openspliceai.create_data import create_dataset

        create_dataset.create_dataset(
            SimpleNamespace(
                output_dir=str(args.output_dir),
                chr_split="train-test",
                biotype="protein-coding",
            )
        )

    print(f"审计报告：{args.output_dir / AUDIT_FILENAME}")
    print(f"NCBI IG/TR 补充：{len(ncbi_missing)}")
    print(
        "IMGT 补充："
        f"{imgt_duplicate_report['supplemented_from_imgt']}"
    )
    print("GENCODE：未使用")
    print("测试集：保持原始 NCBI 基线不变")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
