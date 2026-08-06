#!/usr/bin/env python3
"""Download current NCBI GRCm39 annotation and IMGT splice annotations.

Large biological data files are intentionally not committed to Git. This
script downloads the current NCBI annotation package for the GRCm39 reference
assembly and the current public IMGT/GENE-DB + IMGT/LIGM-DB files, then records
release metadata and SHA-256 checksums in a manifest.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import time
import urllib.request
import zipfile

DEFAULT_ASSEMBLY = "GCF_000001635.27"
IMGT_BASE = "https://www.imgt.org/download"
IMGT_FILES = {
    "GENE-DB/RELEASE": "imgt/gene_db/RELEASE",
    "GENE-DB/README.txt": "imgt/gene_db/README.txt",
    "GENE-DB/IMGTGENEDB-GeneList": "imgt/gene_db/IMGTGENEDB-GeneList",
    "LIGM-DB/currentRelease": "imgt/ligm_db/currentRelease",
    "LIGM-DB/README": "imgt/ligm_db/README",
    "LIGM-DB/accessionNumber.lst": "imgt/ligm_db/accessionNumber.lst",
    "LIGM-DB/imgt.dat.Z": "imgt/ligm_db/imgt.dat.Z",
}


def sha256(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def download_atomic(
    url: str,
    destination: Path,
    force: bool = False,
    attempts: int = 5,
) -> None:
    if destination.exists() and not force:
        print(f"Using existing {destination}")
        return

    destination.parent.mkdir(parents=True, exist_ok=True)
    last_error: Exception | None = None

    for attempt in range(1, attempts + 1):
        request = urllib.request.Request(
            url,
            headers={"User-Agent": "OpenSpliceAI-mouse-immune-data/1.0"},
        )
        with tempfile.NamedTemporaryFile(
            dir=destination.parent,
            prefix=f".{destination.name}.",
            delete=False,
        ) as temporary:
            temporary_path = Path(temporary.name)
            try:
                print(f"Downloading {url} (attempt {attempt}/{attempts})")
                with urllib.request.urlopen(request, timeout=600) as response:
                    shutil.copyfileobj(response, temporary)
                os.replace(temporary_path, destination)
                return
            except Exception as error:
                last_error = error
                temporary_path.unlink(missing_ok=True)

        if attempt < attempts:
            delay = 10 * attempt
            print(f"Download failed: {last_error}. Retrying in {delay} seconds.")
            time.sleep(delay)

    raise RuntimeError(f"Failed to download {url} after {attempts} attempts") from last_error


def run_with_retries(
    command: list[str],
    *,
    attempts: int = 5,
    cleanup_path: Path | None = None,
) -> None:
    last_error: subprocess.CalledProcessError | None = None

    for attempt in range(1, attempts + 1):
        if cleanup_path is not None:
            cleanup_path.unlink(missing_ok=True)

        print(f"Running (attempt {attempt}/{attempts}): {' '.join(command)}")
        result = subprocess.run(command, check=False)
        if result.returncode == 0:
            return

        last_error = subprocess.CalledProcessError(result.returncode, command)
        if attempt < attempts:
            delay = 15 * attempt
            print(
                f"Command failed with exit code {result.returncode}. "
                f"Retrying in {delay} seconds."
            )
            time.sleep(delay)

    assert last_error is not None
    raise last_error


def run_ncbi_datasets(
    output_dir: Path,
    assembly: str,
    force: bool,
) -> tuple[Path, str]:
    datasets = shutil.which("datasets")
    if not datasets:
        raise RuntimeError(
            "NCBI Datasets CLI was not found. Install it with: "
            "mamba install -c conda-forge ncbi-datasets-cli"
        )

    ncbi_dir = output_dir / "ncbi"
    package_zip = ncbi_dir / f"{assembly}.dehydrated.zip"
    extracted_dir = ncbi_dir / assembly
    ncbi_dir.mkdir(parents=True, exist_ok=True)

    if force:
        package_zip.unlink(missing_ok=True)
        if extracted_dir.exists():
            shutil.rmtree(extracted_dir)

    if not extracted_dir.exists():
        command = [
            datasets,
            "download",
            "genome",
            "accession",
            assembly,
            "--include",
            "genome,gff3,gbff,seq-report",
            "--dehydrated",
            "--filename",
            str(package_zip),
            "--no-progressbar",
        ]
        run_with_retries(command, cleanup_path=package_zip)

        if not zipfile.is_zipfile(package_zip):
            raise RuntimeError(f"NCBI dehydrated package is not a valid ZIP: {package_zip}")

        extracted_dir.mkdir(parents=True)
        with zipfile.ZipFile(package_zip) as archive:
            archive.extractall(extracted_dir)

    rehydrate_command = [
        datasets,
        "rehydrate",
        "--directory",
        str(extracted_dir),
        "--max-workers",
        "4",
        "--no-progressbar",
    ]
    run_with_retries(rehydrate_command)

    version = subprocess.run(
        [datasets, "version"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    return extracted_dir, version


def locate_one(root: Path, name: str) -> Path:
    matches = list(root.rglob(name))
    if len(matches) != 1:
        raise RuntimeError(f"Expected one {name} below {root}, found {len(matches)}")
    return matches[0]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data/mouse_immune"),
    )
    parser.add_argument("--assembly", default=DEFAULT_ASSEMBLY)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    root = args.output_dir.resolve()
    root.mkdir(parents=True, exist_ok=True)

    extracted, datasets_version = run_ncbi_datasets(
        root,
        args.assembly,
        args.force,
    )
    ncbi_fasta = locate_one(extracted, "genomic.fna")
    ncbi_gff = locate_one(extracted, "genomic.gff")
    ncbi_gbff = locate_one(extracted, "genomic.gbff")
    data_report_matches = list(extracted.rglob("assembly_data_report.jsonl"))
    data_report = data_report_matches[0] if data_report_matches else None

    for remote, relative in IMGT_FILES.items():
        download_atomic(
            f"{IMGT_BASE}/{remote}",
            root / relative,
            force=args.force,
        )

    gene_release = (root / "imgt/gene_db/RELEASE").read_text(
        encoding="utf-8",
        errors="replace",
    ).strip()
    ligm_release = (root / "imgt/ligm_db/currentRelease").read_text(
        encoding="utf-8",
        errors="replace",
    ).strip()

    files = [ncbi_fasta, ncbi_gff, ncbi_gbff]
    if data_report:
        files.append(data_report)
    files.extend(root / relative for relative in IMGT_FILES.values())

    manifest = {
        "downloaded_at_utc": datetime.now(timezone.utc).isoformat(),
        "ncbi": {
            "assembly_accession": args.assembly,
            "datasets_cli_version": datasets_version,
            "download_mode": "dehydrated_then_rehydrated",
            "package_root": str(extracted.relative_to(root)),
            "genome_fasta": str(ncbi_fasta.relative_to(root)),
            "annotation_gff3": str(ncbi_gff.relative_to(root)),
            "annotation_gbff": str(ncbi_gbff.relative_to(root)),
            "assembly_data_report": (
                str(data_report.relative_to(root)) if data_report else None
            ),
        },
        "imgt": {
            "gene_db_release": gene_release,
            "ligm_db_release": ligm_release,
            "flat_file": "imgt/ligm_db/imgt.dat.Z",
        },
        "files": [
            {
                "path": str(path.relative_to(root)),
                "size_bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for path in files
        ],
    }
    manifest_path = root / "manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    print("\nDownloaded data paths:")
    print(f"  NCBI FASTA: {ncbi_fasta}")
    print(f"  NCBI GFF3:  {ncbi_gff}")
    print(f"  IMGT flat:  {root / 'imgt/ligm_db/imgt.dat.Z'}")
    print(f"  Manifest:   {manifest_path}")
    print(f"  IMGT/GENE-DB release: {gene_release}")
    print(f"  IMGT/LIGM-DB release: {ligm_release}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
