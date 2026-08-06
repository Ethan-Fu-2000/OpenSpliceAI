#!/usr/bin/env python3
"""一键下载、准备并训练严格的小鼠 NCBI+GENCODE M39+IMGT 模型。"""

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
AUDIT = "ncbi_gencode_imgt_splice_audit.json"


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--data-root", type=Path, default=Path("data/mouse_immune"))
    p.add_argument(
        "--work-root",
        type=Path,
        default=Path("work/mouse_gencode_imgt_retraining"),
    )
    p.add_argument(
        "--model-dir",
        type=Path,
        default=Path("models/mouse_ncbi_gencode_imgt"),
    )
    p.add_argument("--release-tag", default="")
    p.add_argument("--flanking-size", type=int, choices=[80, 400, 2000, 10000], default=2000)
    p.add_argument("--epochs", type=int, default=20)
    p.add_argument("--batch-size", type=int, default=2)
    p.add_argument("--num-gpus", type=int, default=1)
    p.add_argument("--random-seed", type=int, default=42)
    p.add_argument("--project-name", default="mouse_ncbi_gencode_imgt")
    p.add_argument("--exp-num", default="1")
    p.add_argument("--skip-download", action="store_true")
    p.add_argument("--skip-train", action="store_true")
    p.add_argument("--force-baseline", action="store_true")
    p.add_argument("--force-augment", action="store_true")
    p.add_argument("--include-nonfunctional", action="store_true")
    p.add_argument("--allow-nonunique-imgt", action="store_true")
    p.add_argument("--early-stopping", action="store_true")
    p.add_argument("--patience", type=int, default=3)
    p.add_argument(
        "--scheduler",
        choices=["MultiStepLR", "CosineAnnealingWarmRestarts"],
        default="MultiStepLR",
    )
    p.add_argument(
        "--loss",
        choices=["cross_entropy_loss", "focal_loss"],
        default="cross_entropy_loss",
    )
    return p.parse_args()


def absolute(path: Path) -> Path:
    return path if path.is_absolute() else (ROOT / path).resolve()


def run(command):
    print("+", " ".join(str(x) for x in command))
    subprocess.run([str(x) for x in command], cwd=ROOT, check=True)


def manifest_ready(data_root: Path) -> bool:
    path = data_root / "manifest.json"
    if not path.exists():
        return False
    try:
        m = json.loads(path.read_text(encoding="utf-8"))
        g = m["gencode"]
        return (
            g["release"] == "M39"
            and g["assembly"] == "GRCm39"
            and (data_root / g["basic_annotation_gtf"]).exists()
        )
    except Exception:
        return False


def ensure_data(data_root: Path, tag: str, skip: bool):
    if manifest_ready(data_root):
        return
    if skip:
        raise FileNotFoundError("本地数据缺少 GENCODE M39/GRCm39")
    required = ["gh", "sha256sum", "tar", "zstd"]
    missing = [x for x in required if shutil.which(x) is None]
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
        raise RuntimeError("下载完成后仍未检测到 GENCODE M39 数据")


def resolve_paths(data_root: Path, work_root: Path):
    m = json.loads((data_root / "manifest.json").read_text(encoding="utf-8"))

    def r(value):
        path = (data_root / value).resolve()
        if not path.exists():
            raise FileNotFoundError(path)
        return path

    return {
        "ncbi_fasta": r(m["ncbi"]["genome_fasta"]),
        "ncbi_gff": r(m["ncbi"]["annotation_gff3"]),
        "gencode_gtf": r(m["gencode"]["basic_annotation_gtf"]),
        "imgt_flat": r(m["imgt"]["flat_file"]),
        "baseline": work_root / "ncbi_baseline",
        "strict": work_root / "ncbi_gencode_imgt",
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


def build_baseline(paths, args):
    out = paths["baseline"]
    if args.force_baseline and out.exists():
        shutil.rmtree(out)
    if complete(out):
        print(f"NCBI 基线已存在：{out}")
        return

    from openspliceai.create_data import create_datafile, create_dataset

    random.seed(args.random_seed)
    ns = SimpleNamespace(
        annotation_gff=str(paths["ncbi_gff"]),
        genome_fasta=str(paths["ncbi_fasta"]),
        output_dir=str(out),
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
    create_datafile.create_datafile(ns)
    create_dataset.create_dataset(ns)


def build_strict(paths, args):
    out = paths["strict"]
    if args.force_augment and out.exists():
        shutil.rmtree(out)
    if complete(out) and (out / AUDIT).exists():
        print(f"严格增强数据已存在：{out}")
        return

    command = [
        sys.executable,
        ROOT / "scripts/prepare_mouse_strict_training.py",
        "--ncbi-gff", paths["ncbi_gff"],
        "--ncbi-fasta", paths["ncbi_fasta"],
        "--gencode-gtf", paths["gencode_gtf"],
        "--imgt-flat", paths["imgt_flat"],
        "--base-data-dir", paths["baseline"],
        "--output-dir", out,
        "--training-context-radius", str(args.flanking_size // 2),
    ]
    if args.include_nonfunctional:
        command.append("--include-nonfunctional")
    if args.allow_nonunique_imgt:
        command.append("--allow-nonunique-imgt")
    run(command)


def summarize(directory: Path):
    result = {}
    for split in ("train", "validation", "test"):
        with h5py.File(directory / f"datafile_{split}.h5", "r") as h:
            labels = [
                x.decode() if isinstance(x, bytes) else str(x)
                for x in h["LABEL"][:]
            ]
            result[split] = {
                "records": len(h["NAME"]),
                "acceptor": sum(x.count("1") for x in labels),
                "donor": sum(x.count("2") for x in labels),
            }
    audit = json.loads((directory / AUDIT).read_text(encoding="utf-8"))
    result["supplement"] = {
        "ncbi": audit["ncbi"]["supplemented"],
        "gencode": audit["gencode"]["supplemented"],
        "imgt": audit["imgt"]["duplicate_filter"]["supplemented_from_imgt"],
    }
    output = directory / "training_data_summary.json"
    output.write_text(
        json.dumps(result, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(result, indent=2, ensure_ascii=False))


def train(paths, args, model_dir):
    if args.skip_train:
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


def main():
    args = parse_args()
    data_root = absolute(args.data_root)
    work_root = absolute(args.work_root)
    model_dir = absolute(args.model_dir)

    ensure_data(data_root, args.release_tag, args.skip_download)
    paths = resolve_paths(data_root, work_root)

    print("=== 1. NCBI GRCm39 基线 ===")
    build_baseline(paths, args)
    print("=== 2. NCBI + GENCODE M39 Basic + IMGT 严格补缺 ===")
    build_strict(paths, args)
    print("=== 3. 检查数据 ===")
    summarize(paths["strict"])
    print("=== 4. 训练 ===")
    train(paths, args, model_dir)
    print(f"数据：{paths['strict']}")
    print(f"审计：{paths['strict'] / AUDIT}")
    print(f"模型：{model_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
