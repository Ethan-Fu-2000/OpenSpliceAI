#!/usr/bin/env python3
"""把 NCBI 和 IMGT 的小鼠 IG/TR 明确剪接位点补入 OpenSpliceAI 训练数据。

正确的数据优先级是：

1. 保留作者使用 NCBI 普通转录本生成的基线数据；
2. 从 NCBI GFF3 读取 IG/TR 免疫区段明确或可由相邻 exon 得到的位点；
3. 检查这些 NCBI 位点是否已经真正存在于基线 HDF5；
4. 基线中没有的 NCBI 免疫位点加入 train / validation；
5. IMGT 中明确标注、且 NCBI 没有覆盖的位点再补入；
6. test 保持作者的 NCBI 基线不变，避免测试集污染。
"""

from __future__ import annotations

import argparse
from pathlib import Path
import shutil
from types import SimpleNamespace

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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="NCBI 优先、IMGT 补缺，生成小鼠 IG/TR 增强训练集。"
    )
    parser.add_argument("--ncbi-gff", type=Path, required=True)
    parser.add_argument("--ncbi-fasta", type=Path, required=True)
    parser.add_argument(
        "--imgt-flat",
        type=Path,
        required=True,
        help="IMGT/LIGM-DB 的 imgt.dat、imgt.dat.gz 或 imgt.dat.Z",
    )
    parser.add_argument(
        "--base-data-dir",
        type=Path,
        required=True,
        help="作者原始 openspliceai create-data 生成的目录",
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--training-context-radius",
        type=int,
        default=5000,
        help=(
            "NCBI 免疫位点两侧提取的真实基因组上下文长度；默认每侧 5000 bp，"
            "可支持 flanking-size=10000"
        ),
    )
    parser.add_argument(
        "--match-radius",
        type=int,
        default=20,
        help="NCBI/IMGT/基线去重时比较位点两侧多少 bp，默认 20",
    )
    parser.add_argument("--identity-threshold", type=float, default=0.90)
    parser.add_argument("--validation-ratio", type=float, default=0.10)
    parser.add_argument(
        "--include-nonfunctional",
        action="store_true",
        help="同时加入 ORF/假基因；默认只使用 IMGT functional 记录",
    )
    parser.add_argument(
        "--skip-dataset",
        action="store_true",
        help="只生成合并后的 datafile_*.h5，不生成 dataset_*.h5",
    )
    return parser.parse_args()


def require_base_files(base_dir: Path) -> None:
    missing = [
        name
        for name in (
            "datafile_train.h5",
            "datafile_validation.h5",
            "datafile_test.h5",
        )
        if not (base_dir / name).exists()
    ]
    if missing:
        raise FileNotFoundError(
            f"基线目录 {base_dir} 缺少文件：{', '.join(missing)}"
        )


def main() -> int:
    args = parse_args()
    require_base_files(args.base_data_dir)
    if args.output_dir.resolve() == args.base_data_dir.resolve():
        raise ValueError(
            "--output-dir 必须与 --base-data-dir 不同，以保留原始 NCBI 基线"
        )
    if args.training_context_radius < args.match_radius:
        raise ValueError("--training-context-radius 不能小于 --match-radius")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    print("[1/7] 读取作者 NCBI 基线中的已有剪接正例……")
    base_train = read_datafile_h5(args.base_data_dir / "datafile_train.h5")
    base_validation = read_datafile_h5(
        args.base_data_dir / "datafile_validation.h5"
    )
    baseline_sites = datafile_positive_sites(
        combine_data(base_train, base_validation),
        context_radius=max(30, args.match_radius),
    )

    print("[2/7] 读取 NCBI GRCm39 的 IG/TR 免疫剪接位点……")
    ncbi_sites = read_ncbi_immune_splice_sites(
        args.ncbi_gff,
        args.ncbi_fasta,
        context_radius=args.training_context_radius,
    )

    print("[3/7] 检查哪些 NCBI 免疫位点尚未真正进入基线 HDF5……")
    ncbi_missing, ncbi_coverage = find_sites_missing_from_baseline(
        ncbi_sites,
        baseline_sites,
        context_radius=args.match_radius,
        identity_threshold=1.0,
    )
    ncbi_train, ncbi_validation = split_splice_sites_by_gene(
        ncbi_missing,
        validation_ratio=args.validation_ratio,
        seed="openspliceai-imgt-v1",
    )

    print("[4/7] 读取 IMGT 明确标注的 functional 小鼠剪接位点……")
    imgt_records = read_imgt_splice_records(
        args.imgt_flat,
        functional_only=not args.include_nonfunctional,
        mouse_only=True,
    )

    print("[5/7] 仅保留 NCBI 未覆盖的 IMGT 位点……")
    imgt_missing, imgt_report = find_missing_imgt_sites(
        imgt_records,
        ncbi_sites,
        context_radius=args.match_radius,
        identity_threshold=args.identity_threshold,
    )
    imgt_train, imgt_validation = split_imgt_records_by_gene(
        imgt_missing,
        validation_ratio=args.validation_ratio,
        seed="openspliceai-imgt-v1",
    )

    report = {
        "policy": {
            "priority": "NCBI baseline -> missing NCBI immune sites -> missing IMGT sites",
            "test_set": "原始 NCBI test 原样复制，未加入 NCBI/IMGT 免疫补充记录",
            "imgt_functional_only": not args.include_nonfunctional,
        },
        "parameters": {
            "training_context_radius": args.training_context_radius,
            "match_radius": args.match_radius,
            "imgt_identity_threshold": args.identity_threshold,
            "validation_ratio": args.validation_ratio,
        },
        "baseline_positive_site_count": len(baseline_sites),
        "ncbi_immune_site_count": len(ncbi_sites),
        "ncbi_already_in_baseline": len(ncbi_sites) - len(ncbi_missing),
        "supplemented_from_ncbi": len(ncbi_missing),
        "ncbi_supplement_train": len(ncbi_train),
        "ncbi_supplement_validation": len(ncbi_validation),
        "ncbi_sites": coverage_to_report(ncbi_coverage),
        "imgt": imgt_report,
        "imgt_supplement_train_records": len(imgt_train),
        "imgt_supplement_validation_records": len(imgt_validation),
    }
    write_audit_report(
        args.output_dir / "ncbi_imgt_splice_audit.json",
        report,
    )

    print("[6/7] 合并 NCBI 与 IMGT 补充记录到 train / validation……")
    for split, base, ncbi_split, imgt_split in (
        ("train", base_train, ncbi_train, imgt_train),
        (
            "validation",
            base_validation,
            ncbi_validation,
            imgt_validation,
        ),
    ):
        ncbi_data = splice_sites_to_data(ncbi_split)
        imgt_data = imgt_records_to_data(imgt_split)
        immune_supplement = combine_data(ncbi_data, imgt_data)

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
            immune_supplement,
        )
        write_datafile_h5(
            args.output_dir / f"datafile_{split}.h5",
            merge_data(base, immune_supplement),
        )

    shutil.copy2(
        args.base_data_dir / "datafile_test.h5",
        args.output_dir / "datafile_test.h5",
    )

    if not args.skip_dataset:
        print("[7/7] 重新生成 dataset_train/validation/test.h5……")
        from openspliceai.create_data import create_dataset

        dataset_args = SimpleNamespace(
            output_dir=str(args.output_dir),
            chr_split="train-test",
            biotype="protein-coding",
        )
        create_dataset.create_dataset(dataset_args)
    else:
        print("[7/7] 已按参数跳过 dataset_*.h5 生成。")

    print()
    print("完成。")
    print(f"审计报告：{args.output_dir / 'ncbi_imgt_splice_audit.json'}")
    print(f"NCBI 免疫位点补充数：{len(ncbi_missing)}")
    print(f"IMGT 位点补充数：{imgt_report['supplemented_from_imgt']}")
    print("测试集保持原始 NCBI 基线不变。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
