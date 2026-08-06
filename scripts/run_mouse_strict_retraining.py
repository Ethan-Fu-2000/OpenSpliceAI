#!/usr/bin/env python3
"""一键下载、准备并训练严格的小鼠 NCBI+IMGT 模型。"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import random
import shutil
import subprocess
import sys
from types import SimpleNamespace

import h5py

ROOT = Path(__file__).resolve().parents[1]
PROFILE = "ncbi_imgt_strict_v2"
AUDIT = "ncbi_imgt_strict_splice_audit.json"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="一键准备 NCBI GRCm39 + 严格 IMGT 补缺数据并重新训练。"
    )
    parser.add_argument(
        "--data-root",
        type=Path,
        default=Path("data/mouse_immune"),
    )
    parser.add_argument(
        "--work-root",
        type=Path,
        default=Path("work/mouse_ncbi_imgt_strict_retraining"),
    )
    parser.add_argument(
        "--model-dir",
        type=Path,
        default=Path("models/mouse_ncbi_imgt_strict"),
    )
    parser.add_argument("--release-tag", default="")
    parser.add_argument(
        "--flanking-size",
        type=int,
        choices=[80, 400, 2000, 10000],
        default=2000,
    )
    parser.add_argument("--epochs", type=int, default=20)
    parser.add_argument("--batch-size", type=int, default=2)
    parser.add_argument("--num-gpus", type=int, default=1)
    parser.add_argument("--random-seed", type=int, default=42)
    parser.add_argument("--project-name", default="mouse_ncbi_imgt_strict")
    parser.add_argument("--exp-num", default="1")
    parser.add_argument("--skip-download", action="store_true")
    parser.add_argument("--skip-train", action="store_true")
    parser.add_argument("--force-baseline", action="store_true")
    parser.add_argument("--force-augment", action="store_true")
    parser.add_argument(
        "--include-nonfunctional",
        action="store_true",
        help="默认关闭；开启后会包含 IMGT ORF/假基因记录",
    )
    parser.add_argument(
        "--allow-nonunique-imgt",
        action="store_true",
        help="默认关闭；不建议开启",
    )
    parser.add_argument("--early-stopping", action="store_true")
    parser.add_argument("--patience", type=int, default=3)
    parser.add_argument(
        "--scheduler",
        choices=["MultiStepLR", "CosineAnnealingWarmRestarts"],
        default="MultiStepLR",
    )
    parser.add_argument(
        "--loss",
        choices=["cross_entropy_loss", "focal_loss"],
        default="cross_entropy_loss",
    )
    return parser.parse_args()


def absolute(path: Path) -> Path:
    return path if path.is_absolute() else (ROOT / path).resolve()


def run(command) -> None:
    print("+", " ".join(str(value) for value in command))
    subprocess.run(
        [str(value) for value in command],
        cwd=ROOT,
        check=True,
    )


def manifest_ready(data_root: Path) -> bool:
    manifest_path = data_root / "manifest.json"
    if not manifest_path.exists():
        return False
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        policy = manifest["policy"]
        ncbi = manifest["ncbi"]
        imgt = manifest["imgt"]
        return (
            policy.get("profile") == PROFILE
            and policy.get("gencode_included") is False
            and ncbi["assembly_accession"] == "GCF_000001635.27"
            and ncbi["assembly_name"] == "GRCm39"
            and (data_root / ncbi["genome_fasta"]).exists()
            and (data_root / ncbi["annotation_gff3"]).exists()
            and (data_root / imgt["flat_file"]).exists()
        )
    except (KeyError, TypeError, ValueError, json.JSONDecodeError):
        return False


def ensure_data(data_root: Path, tag: str, skip_download: bool) -> None:
    if manifest_ready(data_root):
        print(f"使用现有 NCBI+IMGT 严格数据：{data_root}")
        return
    if skip_download:
        raise FileNotFoundError(
            "本地数据不是 ncbi_imgt_strict_v2，或仍是包含 GENCODE 的旧数据。"
        )
    required = ["gh", "sha256sum", "tar", "zstd"]
    missing = [name for name in required if shutil.which(name) is None]
    if missing:
        raise RuntimeError("缺少命令：" + ", ".join(missing))
    command = [
        "bash",
        ROOT / "scripts/download_mouse_immune_release.sh",
        data_root.parent,
    ]
    if tag:
        command.append(tag)
    run(command)
    if not manifest_ready(data_root):
        raise RuntimeError(
            "下载完成后仍未检测到纯 NCBI+IMGT 严格数据，请检查 Release 标签。"
        )


def resolve_paths(data_root: Path, work_root: Path):
    manifest = json.loads(
        (data_root / "manifest.json").read_text(encoding="utf-8")
    )

    def resolve(relative: str) -> Path:
        path = (data_root / relative).resolve()
        if not path.exists():
            raise FileNotFoundError(path)
        return path

    return {
        "ncbi_fasta": resolve(manifest["ncbi"]["genome_fasta"]),
        "ncbi_gff": resolve(manifest["ncbi"]["annotation_gff3"]),
        "imgt_flat": resolve(manifest["imgt"]["flat_file"]),
        "baseline": work_root / "ncbi_baseline",
        "strict": work_root / "ncbi_imgt_strict",
    }


def complete(directory: Path) -> bool:
    return all(
        (directory / name).exists()
        for name in (
            "datafile_train.h5",
            "datafile_validation.h5",
            "datafile_test.h5",
            "dataset_train.h5",
            "dataset_validation.h5",
            "dataset_test.h5",
        )
    )


def build_baseline(paths, args: argparse.Namespace) -> None:
    output = paths["baseline"]
    if args.force_baseline and output.exists():
        shutil.rmtree(output)
    if complete(output):
        print(f"NCBI 基线已存在：{output}")
        return

    from openspliceai.create_data import create_datafile, create_dataset

    random.seed(args.random_seed)
    namespace = SimpleNamespace(
        annotation_gff=str(paths["ncbi_gff"]),
        genome_fasta=str(paths["ncbi_fasta"]),
        output_dir=str(output),
        parse_type="all_isoforms",
        biotype="protein-coding",
        chr_split="train-test",
        split_method="random",
        val_split_ratio=0.1,
        split_ratio=0.8,
        canonical_only=True,
        flanking_size=args.flanking_size,
        verify_h5=False,
        remove_paralogs=False,
        min_identity=0.8,
        min_coverage=0.5,
        write_fasta=False,
    )
    create_datafile.create_datafile(namespace)
    create_dataset.create_dataset(namespace)


def build_strict(paths, args: argparse.Namespace) -> None:
    output = paths["strict"]
    if args.force_augment and output.exists():
        shutil.rmtree(output)
    if complete(output) and (output / AUDIT).exists():
        print(f"严格增强数据已存在：{output}")
        return

    command = [
        sys.executable,
        ROOT / "scripts/prepare_mouse_strict_training.py",
        "--ncbi-gff",
        paths["ncbi_gff"],
        "--ncbi-fasta",
        paths["ncbi_fasta"],
        "--imgt-flat",
        paths["imgt_flat"],
        "--base-data-dir",
        paths["baseline"],
        "--output-dir",
        output,
        "--training-context-radius",
        str(args.flanking_size // 2),
    ]
    if args.include_nonfunctional:
        command.append("--include-nonfunctional")
    if args.allow_nonunique_imgt:
        command.append("--allow-nonunique-imgt")
    run(command)


def summarize(directory: Path) -> None:
    result = {}
    for split in ("train", "validation", "test"):
        with h5py.File(directory / f"datafile_{split}.h5", "r") as handle:
            labels = [
                value.decode() if isinstance(value, bytes) else str(value)
                for value in handle["LABEL"][:]
            ]
            result[split] = {
                "records": len(handle["NAME"]),
                "acceptor": sum(value.count("1") for value in labels),
                "donor": sum(value.count("2") for value in labels),
            }
    audit = json.loads((directory / AUDIT).read_text(encoding="utf-8"))
    result["supplement"] = {
        "ncbi": audit["ncbi"]["supplemented"],
        "imgt": audit["imgt"]["duplicate_filter"]["supplemented_from_imgt"],
        "gencode": 0,
    }
    output = directory / "training_data_summary.json"
    output.write_text(
        json.dumps(result, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(result, indent=2, ensure_ascii=False))
    print(f"数据摘要：{output}")


def train(paths, args: argparse.Namespace, model_dir: Path) -> None:
    if args.skip_train:
        print("已按 --skip-train 跳过训练。")
        return
    from openspliceai.train import train as train_module

    train_module.train(
        SimpleNamespace(
            epochs=args.epochs,
            scheduler=args.scheduler,
            early_stopping=args.early_stopping,
            patience=args.patience,
            output_dir=str(model_dir),
            project_name=args.project_name,
            exp_num=str(args.exp_num),
            flanking_size=args.flanking_size,
            random_seed=args.random_seed,
            train_dataset=str(paths["strict"] / "dataset_train.h5"),
            test_dataset=str(paths["strict"] / "dataset_test.h5"),
            loss=args.loss,
            model="SpliceAI",
            batch_size=args.batch_size,
            num_gpus=args.num_gpus,
        )
    )


def main() -> int:
    args = parse_args()
    data_root = absolute(args.data_root)
    work_root = absolute(args.work_root)
    model_dir = absolute(args.model_dir)

    ensure_data(data_root, args.release_tag, args.skip_download)
    paths = resolve_paths(data_root, work_root)

    print("=== 1. 生成 NCBI GRCm39 基线 ===")
    build_baseline(paths, args)
    print("=== 2. 严格补充 NCBI IG/TR 与 IMGT 位点 ===")
    build_strict(paths, args)
    print("=== 3. 检查最终数据 ===")
    summarize(paths["strict"])
    print("=== 4. 训练 ===")
    train(paths, args, model_dir)
    print(f"训练数据：{paths['strict']}")
    print(f"审计报告：{paths['strict'] / AUDIT}")
    print(f"模型目录：{model_dir}")
    print("GENCODE：未使用")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
