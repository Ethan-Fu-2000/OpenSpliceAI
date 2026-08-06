#!/usr/bin/env python3
"""Add IMGT-explicit mouse IG/TR splice sites to OpenSpliceAI datafiles.

Run the normal NCBI-based ``openspliceai create-data`` first. This script then
checks which functional mouse IMGT splice annotations are already represented
by the NCBI-derived labels and appends only the missing examples to train and
validation. The NCBI test set is copied unchanged.
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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Audit NCBI mouse IG/TR splice labels, supplement missing explicit "
            "sites from IMGT, and rebuild OpenSpliceAI datasets."
        )
    )
    parser.add_argument("--ncbi-gff", type=Path, required=True)
    parser.add_argument("--ncbi-fasta", type=Path, required=True)
    parser.add_argument(
        "--imgt-flat",
        type=Path,
        required=True,
        help="IMGT/LIGM-DB imgt.dat, imgt.dat.gz, or imgt.dat.Z",
    )
    parser.add_argument(
        "--base-data-dir",
        type=Path,
        required=True,
        help="Directory produced by the normal openspliceai create-data command",
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--context-radius", type=int, default=30)
    parser.add_argument("--match-radius", type=int, default=20)
    parser.add_argument("--identity-threshold", type=float, default=0.90)
    parser.add_argument("--validation-ratio", type=float, default=0.10)
    parser.add_argument(
        "--include-nonfunctional",
        action="store_true",
        help="Also include ORF/pseudogene IMGT records; functional only by default",
    )
    parser.add_argument(
        "--skip-dataset",
        action="store_true",
        help="Write merged datafile_*.h5 only; do not rebuild dataset_*.h5",
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
            f"Missing base OpenSpliceAI files in {base_dir}: {', '.join(missing)}"
        )


def main() -> int:
    args = parse_args()
    require_base_files(args.base_data_dir)
    if args.output_dir.resolve() == args.base_data_dir.resolve():
        raise ValueError(
            "--output-dir must differ from --base-data-dir to preserve the NCBI baseline"
        )
    args.output_dir.mkdir(parents=True, exist_ok=True)

    print("[1/5] Auditing NCBI IG/TR splice labels...")
    ncbi_sites = read_ncbi_immune_splice_sites(
        args.ncbi_gff,
        args.ncbi_fasta,
        context_radius=args.context_radius,
    )

    print("[2/5] Reading explicit functional mouse IMGT splice features...")
    imgt_records = read_imgt_splice_records(
        args.imgt_flat,
        functional_only=not args.include_nonfunctional,
        mouse_only=True,
    )

    print("[3/5] Removing sites already covered by NCBI context...")
    missing_records, report = find_missing_imgt_sites(
        imgt_records,
        ncbi_sites,
        context_radius=args.match_radius,
        identity_threshold=args.identity_threshold,
    )
    train_records, validation_records = split_imgt_records_by_gene(
        missing_records,
        validation_ratio=args.validation_ratio,
    )
    report.update(
        {
            "supplement_train_records": len(train_records),
            "supplement_validation_records": len(validation_records),
            "test_policy": "NCBI test data copied unchanged; no IMGT examples added",
        }
    )
    write_audit_report(args.output_dir / "ncbi_imgt_splice_audit.json", report)

    print("[4/5] Merging IMGT supplement into train/validation datafiles...")
    for split, records in (
        ("train", train_records),
        ("validation", validation_records),
    ):
        base = read_datafile_h5(args.base_data_dir / f"datafile_{split}.h5")
        supplement = imgt_records_to_data(records)
        write_datafile_h5(
            args.output_dir / f"imgt_supplement_{split}.h5",
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

    if not args.skip_dataset:
        print("[5/5] Rebuilding dataset_train/validation/test.h5...")
        from openspliceai.create_data import create_dataset

        dataset_args = SimpleNamespace(
            output_dir=str(args.output_dir),
            chr_split="train-test",
            biotype="protein-coding",
        )
        create_dataset.create_dataset(dataset_args)
    else:
        print("[5/5] Dataset rebuilding skipped by request.")

    print(f"Audit report: {args.output_dir / 'ncbi_imgt_splice_audit.json'}")
    print(f"IMGT missing sites added: {report['supplemented_from_imgt']}")
    print("NCBI test set was not changed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
