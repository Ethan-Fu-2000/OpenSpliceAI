#!/usr/bin/env python3
"""Download the latest public IMGT datasets needed for IG/TR splice training.

The downloader discovers release identifiers from IMGT at runtime, downloads
files atomically, and writes a SHA-256 manifest so a training run can always be
reproduced. Large IMGT/LIGM-DB flat files are optional.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import tempfile
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

BASE = "https://www.imgt.org/download"

CORE_FILES = {
    "GENE-DB/RELEASE": "gene_db/RELEASE",
    "GENE-DB/README.txt": "gene_db/README.txt",
    "GENE-DB/IMGTGENEDB-GeneList": "gene_db/IMGTGENEDB-GeneList",
    "GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-nt-WithoutGaps-F+ORF+inframeP":
        "gene_db/IMGTGENEDB-ReferenceSequences.fasta-nt-WithoutGaps-F+ORF+inframeP",
    "LIGM-DB/currentRelease": "ligm_db/currentRelease",
    "LIGM-DB/README": "ligm_db/README",
}

LIGM_FLAT_FILES = {
    "LIGM-DB/imgt.dat.Z": "ligm_db/imgt.dat.Z",
    "LIGM-DB/imgt.fasta.Z": "ligm_db/imgt.fasta.Z",
    "LIGM-DB/accessionNumber.lst": "ligm_db/accessionNumber.lst",
}


def sha256(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def download_atomic(url: str, destination: Path, timeout: int = 120) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "OpenSpliceAI-IMGT-downloader/1.0"},
    )
    with tempfile.NamedTemporaryFile(
        dir=destination.parent, prefix=f".{destination.name}.", delete=False
    ) as tmp:
        tmp_path = Path(tmp.name)
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                shutil.copyfileobj(response, tmp)
            os.replace(tmp_path, destination)
        except Exception:
            tmp_path.unlink(missing_ok=True)
            raise


def read_release(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="replace").strip()


def iter_files(include_ligm_flat: bool) -> Iterable[tuple[str, str]]:
    yield from CORE_FILES.items()
    if include_ligm_flat:
        yield from LIGM_FLAT_FILES.items()


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Download current public IMGT data with release metadata and checksums."
    )
    parser.add_argument(
        "--output-dir", type=Path, default=Path("data/imgt"),
        help="Destination directory (default: data/imgt)",
    )
    parser.add_argument(
        "--include-ligm-flat", action="store_true",
        help="Also download the large weekly IMGT/LIGM-DB flat files.",
    )
    parser.add_argument(
        "--force", action="store_true", help="Redownload files that already exist."
    )
    args = parser.parse_args()

    root = args.output_dir.resolve()
    root.mkdir(parents=True, exist_ok=True)
    records = []

    for remote, relative in iter_files(args.include_ligm_flat):
        destination = root / relative
        url = f"{BASE}/{remote}"
        if args.force or not destination.exists():
            print(f"Downloading {url}")
            download_atomic(url, destination)
        else:
            print(f"Using existing {destination}")
        records.append({
            "source_url": url,
            "relative_path": str(destination.relative_to(root)),
            "size_bytes": destination.stat().st_size,
            "sha256": sha256(destination),
        })

    gene_release = read_release(root / "gene_db/RELEASE")
    ligm_release = read_release(root / "ligm_db/currentRelease")
    manifest = {
        "downloaded_at_utc": datetime.now(timezone.utc).isoformat(),
        "source": "IMGT, the international ImMunoGeneTics information system",
        "gene_db_release": gene_release,
        "ligm_db_release": ligm_release,
        "include_ligm_flat": args.include_ligm_flat,
        "files": records,
    }
    manifest_path = root / "manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(f"Wrote {manifest_path}")
    print(f"IMGT/GENE-DB release: {gene_release}")
    print(f"IMGT/LIGM-DB release: {ligm_release}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
