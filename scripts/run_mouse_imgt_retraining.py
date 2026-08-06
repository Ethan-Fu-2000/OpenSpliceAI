#!/usr/bin/env python3
"""一键完成小鼠 NCBI + IMGT 数据准备和 OpenSpliceAI 重新训练。

默认流程：
1. 若本地没有数据，则从本仓库最新 mouse-immune-data-* Release 下载并解压；
2. 读取 manifest.json，自动定位 GRCm39 FASTA、GFF3 和 IMGT flat file；
3. 生成作者原始 NCBI protein-coding 基线；
4. 将基线缺失的 NCBI IG/TR 位点补入；
5. 再补入 NCBI 未覆盖的 IMGT functional 位点；
6. 检查最终 HDF5；
7. 从头训练模型。

可以用 ``--skip-train`` 只准备数据，也可以在确认审计报告后再次运行开始训练。
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import random
import shutil
import subprocess
import sys
from types import SimpleNamespace
from typing import Iterable, Sequence

import h5py


REPO_ROOT = Path(__file__).resolve().parents[1]


@dataclass
class PipelinePaths:
    data_root: Path
    work_root: Path
    baseline_dir: Path
    augmented_dir: Path
    model_dir: Path
    manifest: Path
    ncbi_fasta: Path
    ncbi_gff: Path
    imgt_flat: Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="一键准备 NCBI+IMGT 小鼠训练集并重新训练 OpenSpliceAI。"
    )
    parser.add_argument(
        "--data-root",
        type=Path,
        default=Path("data/mouse_immune"),
        help="Release 解压后的数据目录",
    )
    parser.add_argument(
        "--work-root",
        type=Path,
        default=Path("work/mouse_imgt_retraining"),
        help="基线和增强训练集输出目录",
    )
    parser.add_argument(
        "--model-dir",
        type=Path,
        default=Path("models/mouse_ncbi_imgt"),
    )
    parser.add_argument("--release-tag", default="")
    parser.add_argument("--project-name", default="mouse_ncbi_imgt")
    parser.add_argument("--exp-num", default="1")
    parser.add_argument(
        "--flanking-size",
        type=int,
        choices=[80, 400, 2000, 10000],
        default=2000,
    )
    parser.add_argument("--epochs", type=int, default=20)
    parser.add_argument(
        "--batch-size",
        type=int,
        default=2,
        help="8 GB 显存建议先用 2；显存充足可逐步调大",
    )
    parser.add_argument("--num-gpus", type=int, default=1)
    parser.add_argument("--random-seed", type=int, default=42)
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
    parser.add_argument("--early-stopping", action="store_true")
    parser.add_argument("--patience", type=int, default=3)
    parser.add_argument(
        "--skip-download",
        action="store_true",
        help="不自动从 Release 下载；本地必须已有 manifest.json",
    )
    parser.add_argument(
        "--skip-train",
        action="store_true",
        help="只准备并检查数据，不启动模型训练",
    )
    parser.add_argument(
        "--force-baseline",
        action="store_true",
        help="删除并重新生成作者 NCBI 基线",
    )
    parser.add_argument(
        "--force-augment",
        action="store_true",
        help="删除并重新生成 NCBI+IMGT 增强数据",
    )
    parser.add_argument(
        "--include-nonfunctional",
        action="store_true",
        help="同时加入 IMGT ORF/假基因；默认不要开启",
    )
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def absolute(path: Path) -> Path:
    return path if path.is_absolute() else (REPO_ROOT / path).resolve()


def run_command(
    command: Sequence[str],
    *,
    dry_run: bool = False,
    env: dict[str, str] | None = None,
) -> None:
    print("+", " ".join(str(item) for item in command))
    if dry_run:
        return
    subprocess.run(
        [str(item) for item in command],
        cwd=REPO_ROOT,
        env=env,
        check=True,
    )


def require_commands(commands: Iterable[str]) -> None:
    missing = [command for command in commands if shutil.which(command) is None]
    if missing:
        raise RuntimeError(
            "缺少系统命令："
            + ", ".join(missing)
            + "。请先按照中文文档安装。"
        )


def maybe_download(
    data_root: Path,
    release_tag: str,
    *,
    skip_download: bool,
    dry_run: bool,
) -> None:
    manifest = data_root / "manifest.json"
    if manifest.exists():
        print(f"已存在数据清单，跳过下载：{manifest}")
        return
    if skip_download:
        raise FileNotFoundError(
            f"--skip-download 已开启，但不存在 {manifest}"
        )
    if data_root.name != "mouse_immune":
        raise ValueError(
            "自动下载时 --data-root 的最后一级目录必须是 mouse_immune；"
            "例如 data/mouse_immune"
        )
    require_commands(["bash", "gh", "sha256sum", "tar", "zstd"])
    command = [
        "bash",
        str(REPO_ROOT / "scripts/download_mouse_immune_release.sh"),
        str(data_root.parent),
    ]
    if release_tag:
        command.append(release_tag)
    run_command(command, dry_run=dry_run)


def resolve_manifest_paths(
    data_root: Path,
    work_root: Path,
    model_dir: Path,
) -> PipelinePaths:
    manifest_path = data_root / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    ncbi = manifest["ncbi"]
    imgt = manifest["imgt"]

    def resolve(relative: str) -> Path:
        path = (data_root / relative).resolve()
        if not path.exists():
            raise FileNotFoundError(f"清单中的数据文件不存在：{path}")
        return path

    return PipelinePaths(
        data_root=data_root,
        work_root=work_root,
        baseline_dir=work_root / "ncbi_baseline",
        augmented_dir=work_root / "ncbi_imgt",
        model_dir=model_dir,
        manifest=manifest_path,
        ncbi_fasta=resolve(ncbi["genome_fasta"]),
        ncbi_gff=resolve(ncbi["annotation_gff3"]),
        imgt_flat=resolve(imgt["flat_file"]),
    )


def baseline_complete(directory: Path) -> bool:
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


def augmented_complete(directory: Path) -> bool:
    return baseline_complete(directory) and (
        directory / "ncbi_imgt_splice_audit.json"
    ).exists()


def build_baseline(
    paths: PipelinePaths,
    args: argparse.Namespace,
) -> None:
    if args.force_baseline and paths.baseline_dir.exists() and not args.dry_run:
        shutil.rmtree(paths.baseline_dir)
    if baseline_complete(paths.baseline_dir):
        print(f"作者 NCBI 基线已存在：{paths.baseline_dir}")
        return
    if args.dry_run:
        print(f"[dry-run] 将生成 NCBI 基线：{paths.baseline_dir}")
        return

    from openspliceai.create_data import create_datafile, create_dataset

    random.seed(args.random_seed)
    namespace = SimpleNamespace(
        annotation_gff=str(paths.ncbi_gff),
        genome_fasta=str(paths.ncbi_fasta),
        output_dir=str(paths.baseline_dir),
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


def build_augmented(
    paths: PipelinePaths,
    args: argparse.Namespace,
) -> None:
    if args.force_augment and paths.augmented_dir.exists() and not args.dry_run:
        shutil.rmtree(paths.augmented_dir)
    if augmented_complete(paths.augmented_dir):
        print(f"NCBI+IMGT 增强数据已存在：{paths.augmented_dir}")
        return

    command = [
        sys.executable,
        str(REPO_ROOT / "scripts/prepare_mouse_immune_training.py"),
        "--ncbi-gff",
        str(paths.ncbi_gff),
        "--ncbi-fasta",
        str(paths.ncbi_fasta),
        "--imgt-flat",
        str(paths.imgt_flat),
        "--base-data-dir",
        str(paths.baseline_dir),
        "--output-dir",
        str(paths.augmented_dir),
        "--training-context-radius",
        "5000",
        "--match-radius",
        "20",
        "--identity-threshold",
        "0.90",
        "--validation-ratio",
        "0.10",
    ]
    if args.include_nonfunctional:
        command.append("--include-nonfunctional")
    run_command(command, dry_run=args.dry_run)


def inspect_datafile(path: Path) -> dict[str, int]:
    with h5py.File(path, "r") as handle:
        names = handle["NAME"][:]
        labels = handle["LABEL"][:]
        donor = 0
        acceptor = 0
        for raw in labels:
            text = raw.decode("utf-8") if isinstance(raw, bytes) else str(raw)
            donor += text.count("2")
            acceptor += text.count("1")
        return {
            "records": len(names),
            "acceptor_labels": acceptor,
            "donor_labels": donor,
        }


def inspect_dataset(path: Path) -> dict[str, int]:
    with h5py.File(path, "r") as handle:
        x_keys = sorted(key for key in handle if key.startswith("X"))
        y_keys = sorted(key for key in handle if key.startswith("Y"))
        if len(x_keys) != len(y_keys) or not x_keys:
            raise RuntimeError(f"dataset 结构异常：{path}")
        windows = sum(int(handle[key].shape[0]) for key in x_keys)
        return {"shards": len(x_keys), "windows": windows}


def validate_outputs(paths: PipelinePaths) -> dict[str, object]:
    summary: dict[str, object] = {
        "datafiles": {},
        "datasets": {},
    }
    for split in ("train", "validation", "test"):
        summary["datafiles"][split] = inspect_datafile(
            paths.augmented_dir / f"datafile_{split}.h5"
        )
        summary["datasets"][split] = inspect_dataset(
            paths.augmented_dir / f"dataset_{split}.h5"
        )

    audit_path = paths.augmented_dir / "ncbi_imgt_splice_audit.json"
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    summary["audit"] = {
        "supplemented_from_ncbi": audit["supplemented_from_ncbi"],
        "supplemented_from_imgt": audit["imgt"]["supplemented_from_imgt"],
        "test_policy": audit["policy"]["test_set"],
    }
    summary_path = paths.augmented_dir / "training_data_summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    print(f"数据检查报告：{summary_path}")
    return summary


def train_model(paths: PipelinePaths, args: argparse.Namespace) -> None:
    if args.skip_train:
        print("已按 --skip-train 跳过模型训练。")
        return
    if args.dry_run:
        print(f"[dry-run] 将训练模型并输出到：{paths.model_dir}")
        return

    from openspliceai.train import train as train_module

    namespace = SimpleNamespace(
        epochs=args.epochs,
        scheduler=args.scheduler,
        early_stopping=args.early_stopping,
        patience=args.patience,
        output_dir=str(paths.model_dir),
        project_name=args.project_name,
        exp_num=str(args.exp_num),
        flanking_size=args.flanking_size,
        random_seed=args.random_seed,
        train_dataset=str(paths.augmented_dir / "dataset_train.h5"),
        test_dataset=str(paths.augmented_dir / "dataset_test.h5"),
        loss=args.loss,
        model="SpliceAI",
        batch_size=args.batch_size,
        num_gpus=args.num_gpus,
    )
    train_module.train(namespace)


def write_run_config(paths: PipelinePaths, args: argparse.Namespace) -> None:
    paths.work_root.mkdir(parents=True, exist_ok=True)
    payload = {
        "arguments": vars(args),
        "paths": {key: str(value) for key, value in asdict(paths).items()},
        "python": sys.version,
    }
    (paths.work_root / "run_config.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False, default=str) + "\n",
        encoding="utf-8",
    )


def main() -> int:
    args = parse_args()
    data_root = absolute(args.data_root)
    work_root = absolute(args.work_root)
    model_dir = absolute(args.model_dir)

    maybe_download(
        data_root,
        args.release_tag,
        skip_download=args.skip_download,
        dry_run=args.dry_run,
    )
    if args.dry_run and not (data_root / "manifest.json").exists():
        print("[dry-run] 数据尚未解压，后续路径将在实际运行时从 manifest 解析。")
        return 0

    paths = resolve_manifest_paths(data_root, work_root, model_dir)
    write_run_config(paths, args)

    print("\n=== 阶段 1：生成作者 NCBI 基线 ===")
    build_baseline(paths, args)

    print("\n=== 阶段 2：NCBI 优先，IMGT 补缺 ===")
    build_augmented(paths, args)

    if not args.dry_run:
        print("\n=== 阶段 3：检查最终训练数据 ===")
        validate_outputs(paths)

    print("\n=== 阶段 4：重新训练 ===")
    train_model(paths, args)

    print("\n流程完成。")
    print(f"增强训练数据：{paths.augmented_dir}")
    print(f"审计报告：{paths.augmented_dir / 'ncbi_imgt_splice_audit.json'}")
    print(f"模型目录：{paths.model_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
