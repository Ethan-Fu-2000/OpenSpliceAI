"""Utilities for supplementing OpenSpliceAI mouse training data with IMGT sites.

The normal OpenSpliceAI data builder infers donor/acceptor labels from adjacent
exons in NCBI GFF3 transcripts. That misses biologically explicit IG/TR sites
such as the donor at the end of a single-exon J gene and the acceptor at the
first constant-region exon. IMGT/LIGM-DB annotates those sites directly as
DONOR-SPLICE, ACCEPTOR-SPLICE, INT-DONOR-SPLICE and
INT-ACCEPTOR-SPLICE.

This module:
* audits splice sites recoverable from NCBI GFF3;
* parses explicit splice features from an IMGT flat file;
* compares label-centred sequence contexts so NCBI-covered sites are not added
  twice; and
* writes/merges OpenSpliceAI ``datafile_*.h5`` records.
"""

from __future__ import annotations

from collections import defaultdict
from contextlib import contextmanager
from dataclasses import dataclass, field
from difflib import SequenceMatcher
import gzip
import hashlib
import io
import json
from pathlib import Path
import re
import shutil
import subprocess
from typing import Dict, Iterable, Iterator, List, Mapping, MutableMapping, Optional, Sequence, TextIO, Tuple

import h5py
import numpy as np

DONOR_LABEL = 2
ACCEPTOR_LABEL = 1

IMGT_DONOR_KEYS = {"DONOR-SPLICE", "INT-DONOR-SPLICE"}
IMGT_ACCEPTOR_KEYS = {"ACCEPTOR-SPLICE", "INT-ACCEPTOR-SPLICE"}

DIRECT_DONOR_TYPES = {
    "donor_splice_site",
    "splice_donor_site",
    "five_prime_splice_site",
    "5_prime_splice_site",
}
DIRECT_ACCEPTOR_TYPES = {
    "acceptor_splice_site",
    "splice_acceptor_site",
    "three_prime_splice_site",
    "3_prime_splice_site",
}
TRANSCRIPT_LIKE_TYPES = {
    "mRNA",
    "transcript",
    "V_gene_segment",
    "D_gene_segment",
    "J_gene_segment",
    "C_gene_segment",
}
IMMUNE_BIOTYPES = {
    "V_region",
    "V_segment",
    "D_segment",
    "J_segment",
    "C_region",
    "segment",
}
IMMUNE_PREFIXES = ("IGH", "IGK", "IGL", "TRA", "TRB", "TRG", "TRD")
DATASET_NAMES = ("NAME", "CHROM", "STRAND", "TX_START", "TX_END", "SEQ", "LABEL")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def normalize_gene_name(name: str) -> str:
    """Normalize NCBI/IMGT gene or allele names for cross-source matching."""
    value = (name or "").strip().strip('"').upper()
    value = value.split("*")[0]
    value = value.replace("-", "").replace("_", "").replace(" ", "")
    return value


def is_immune_gene(name: str, biotype: str = "") -> bool:
    normalized = normalize_gene_name(name)
    return biotype in IMMUNE_BIOTYPES or normalized.startswith(IMMUNE_PREFIXES)


@dataclass(frozen=True)
class SpliceSite:
    source: str
    accession: str
    gene: str
    site_type: str
    label: int
    index: int
    oriented_sequence: str
    feature_key: str
    functionality: str = ""
    seqid: str = ""
    genomic_position: Optional[int] = None
    strand: str = "+"

    def context(self, radius: int = 30) -> str:
        start = self.index - radius
        end = self.index + radius + 1
        left = "N" * max(0, -start)
        right = "N" * max(0, end - len(self.oriented_sequence))
        body = self.oriented_sequence[max(0, start):min(len(self.oriented_sequence), end)]
        return left + body + right

    @property
    def motif(self) -> str:
        if self.site_type == "donor":
            return self.oriented_sequence[self.index + 1:self.index + 3]
        return self.oriented_sequence[max(0, self.index - 2):self.index]


@dataclass
class ImgtRecord:
    accession: str
    species: str
    sequence: str
    gene: str
    functionality: str
    strand: str
    sites: List[SpliceSite] = field(default_factory=list)


@dataclass(frozen=True)
class FastaEntry:
    length: int
    offset: int
    line_bases: int
    line_width: int


class IndexedFasta:
    """Minimal random-access FASTA reader using a standard ``.fai`` index."""

    def __init__(self, fasta_path: Path | str, build_if_missing: bool = True):
        self.path = Path(fasta_path)
        self.index_path = Path(str(self.path) + ".fai")
        if not self.index_path.exists():
            if not build_if_missing:
                raise FileNotFoundError(f"Missing FASTA index: {self.index_path}")
            self._build_index()
        self.entries = self._read_index()
        self._handle = self.path.open("rb")

    def close(self) -> None:
        self._handle.close()

    def __enter__(self) -> "IndexedFasta":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self.close()

    def _build_index(self) -> None:
        rows: List[Tuple[str, FastaEntry]] = []
        with self.path.open("rb") as handle:
            name: Optional[str] = None
            length = 0
            seq_offset = 0
            line_bases = 0
            line_width = 0
            while True:
                line_offset = handle.tell()
                line = handle.readline()
                if not line:
                    break
                if line.startswith(b">"):
                    if name is not None:
                        rows.append((name, FastaEntry(length, seq_offset, line_bases, line_width)))
                    name = line[1:].strip().split(None, 1)[0].decode("utf-8")
                    length = 0
                    seq_offset = handle.tell()
                    line_bases = 0
                    line_width = 0
                    continue
                stripped = line.rstrip(b"\r\n")
                if name is None or not stripped:
                    continue
                if line_bases == 0:
                    line_bases = len(stripped)
                    line_width = len(line)
                    seq_offset = line_offset
                length += len(stripped)
            if name is not None:
                rows.append((name, FastaEntry(length, seq_offset, line_bases, line_width)))
        with self.index_path.open("w", encoding="utf-8") as out:
            for name, entry in rows:
                out.write(
                    f"{name}\t{entry.length}\t{entry.offset}\t"
                    f"{entry.line_bases}\t{entry.line_width}\n"
                )

    def _read_index(self) -> Dict[str, FastaEntry]:
        entries: Dict[str, FastaEntry] = {}
        with self.index_path.open(encoding="utf-8") as handle:
            for line in handle:
                if not line.strip():
                    continue
                name, length, offset, line_bases, line_width = line.rstrip().split("\t")[:5]
                entries[name] = FastaEntry(
                    int(length), int(offset), int(line_bases), int(line_width)
                )
        return entries

    def length(self, seqid: str) -> int:
        return self.entries[seqid].length

    def fetch(self, seqid: str, start: int, end: int) -> str:
        """Fetch a 0-based, half-open interval."""
        entry = self.entries[seqid]
        start = max(0, start)
        end = min(entry.length, end)
        if start >= end:
            return ""
        chunks: List[bytes] = []
        position = start
        while position < end:
            row = position // entry.line_bases
            col = position % entry.line_bases
            take = min(end - position, entry.line_bases - col)
            byte_offset = entry.offset + row * entry.line_width + col
            self._handle.seek(byte_offset)
            chunks.append(self._handle.read(take))
            position += take
        return b"".join(chunks).decode("ascii").upper()


def parse_gff_attributes(raw: str) -> Dict[str, List[str]]:
    result: Dict[str, List[str]] = defaultdict(list)
    for item in raw.strip().split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
        elif " " in item:
            key, value = item.split(" ", 1)
        else:
            continue
        for part in value.strip().strip('"').split(","):
            result[key].append(part)
    return dict(result)


def _first(attrs: Mapping[str, Sequence[str]], *keys: str, default: str = "") -> str:
    for key in keys:
        values = attrs.get(key)
        if values:
            return str(values[0])
    return default


def _extract_context(
    fasta: IndexedFasta,
    seqid: str,
    genomic_position: int,
    strand: str,
    radius: int,
) -> str:
    start0 = genomic_position - 1 - radius
    end0 = genomic_position + radius
    left_pad = "N" * max(0, -start0)
    right_pad = "N" * max(0, end0 - fasta.length(seqid))
    sequence = (
        left_pad
        + fasta.fetch(seqid, max(0, start0), min(fasta.length(seqid), end0))
        + right_pad
    )
    if strand == "-":
        sequence = reverse_complement(sequence)
    expected = 2 * radius + 1
    if len(sequence) < expected:
        sequence += "N" * (expected - len(sequence))
    return sequence


def _site_from_context(
    *,
    source: str,
    accession: str,
    gene: str,
    site_type: str,
    context: str,
    feature_key: str,
    functionality: str = "",
    seqid: str = "",
    genomic_position: Optional[int] = None,
    strand: str = "+",
) -> SpliceSite:
    label = DONOR_LABEL if site_type == "donor" else ACCEPTOR_LABEL
    return SpliceSite(
        source=source,
        accession=accession,
        gene=gene,
        site_type=site_type,
        label=label,
        index=len(context) // 2,
        oriented_sequence=context,
        feature_key=feature_key,
        functionality=functionality,
        seqid=seqid,
        genomic_position=genomic_position,
        strand=strand,
    )


def read_ncbi_immune_splice_sites(
    gff_path: Path | str,
    genome_fasta: Path | str,
    context_radius: int = 30,
) -> List[SpliceSite]:
    """Extract explicit and adjacent-exon IG/TR splice labels from NCBI GFF3.

    A single-exon J segment produces no inferred donor, and the first exon of a
    C segment produces no inferred acceptor. Those are intentionally left
    absent so IMGT can supplement them.
    """
    parent_gene: Dict[str, str] = {}
    parent_strand: Dict[str, str] = {}
    parent_seqid: Dict[str, str] = {}
    transcript_exons: MutableMapping[str, List[Tuple[int, int]]] = defaultdict(list)
    direct_features: List[Tuple[str, str, str, int, int, str, str]] = []

    with Path(gff_path).open(encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            if not raw_line or raw_line.startswith("#"):
                continue
            fields = raw_line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            seqid, _, feature_type, start_s, end_s, _, strand, _, attr_s = fields
            try:
                start, end = int(start_s), int(end_s)
            except ValueError:
                continue
            attrs = parse_gff_attributes(attr_s)
            feature_id = _first(attrs, "ID")
            parents = attrs.get("Parent", [])

            if feature_type == "gene":
                gene = _first(attrs, "gene", "Name", "gene_name", default=feature_id)
                biotype = _first(attrs, "gene_biotype", "gene_type")
                if feature_id and is_immune_gene(gene, biotype):
                    parent_gene[feature_id] = gene
                    parent_strand[feature_id] = strand
                    parent_seqid[feature_id] = seqid
                continue

            if feature_type in TRANSCRIPT_LIKE_TYPES:
                for parent in parents:
                    if parent in parent_gene and feature_id:
                        parent_gene[feature_id] = parent_gene[parent]
                        parent_strand[feature_id] = strand
                        parent_seqid[feature_id] = seqid
                        break
                continue

            if feature_type == "exon":
                for parent in parents:
                    if parent in parent_gene:
                        transcript_exons[parent].append((start, end))
                continue

            if feature_type in DIRECT_DONOR_TYPES | DIRECT_ACCEPTOR_TYPES:
                gene = ""
                for parent in parents:
                    if parent in parent_gene:
                        gene = parent_gene[parent]
                        break
                if gene:
                    direct_features.append(
                        (feature_type, gene, seqid, start, end, strand, feature_id)
                    )

    sites: List[SpliceSite] = []
    with IndexedFasta(genome_fasta) as fasta:
        for parent, exons in transcript_exons.items():
            if len(exons) < 2:
                continue
            strand = parent_strand[parent]
            seqid = parent_seqid[parent]
            gene = parent_gene[parent]
            ordered = sorted(exons, key=lambda value: value[0], reverse=(strand == "-"))
            for left, right in zip(ordered, ordered[1:]):
                if strand == "+":
                    donor_pos = left[1]
                    acceptor_pos = right[0]
                else:
                    donor_pos = left[0]
                    acceptor_pos = right[1]
                for site_type, position in (
                    ("donor", donor_pos),
                    ("acceptor", acceptor_pos),
                ):
                    context = _extract_context(
                        fasta, seqid, position, strand, context_radius
                    )
                    sites.append(
                        _site_from_context(
                            source="NCBI-inferred",
                            accession=parent,
                            gene=gene,
                            site_type=site_type,
                            context=context,
                            feature_key="adjacent-exon",
                            seqid=seqid,
                            genomic_position=position,
                            strand=strand,
                        )
                    )

        for feature_type, gene, seqid, start, end, strand, feature_id in direct_features:
            site_type = "donor" if feature_type in DIRECT_DONOR_TYPES else "acceptor"
            if strand == "+":
                if site_type == "donor":
                    position = start if end - start + 1 >= 3 else start - 1
                else:
                    position = start + 3 if end - start + 1 >= 5 else end + 1
            else:
                if site_type == "donor":
                    position = end if end - start + 1 >= 3 else end + 1
                else:
                    position = end - 3 if end - start + 1 >= 5 else start - 1
            if position < 1 or seqid not in fasta.entries:
                continue
            context = _extract_context(fasta, seqid, position, strand, context_radius)
            sites.append(
                _site_from_context(
                    source="NCBI-explicit",
                    accession=feature_id,
                    gene=gene,
                    site_type=site_type,
                    context=context,
                    feature_key=feature_type,
                    seqid=seqid,
                    genomic_position=position,
                    strand=strand,
                )
            )

    return deduplicate_sites(sites)


@contextmanager
def open_imgt_text(path: Path | str) -> Iterator[TextIO]:
    path = Path(path)
    if path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
            yield handle
        return
    if path.suffix == ".Z":
        gzip_exe = shutil.which("gzip")
        if not gzip_exe:
            raise RuntimeError("Reading IMGT .Z files requires the 'gzip' executable")
        process = subprocess.Popen(
            [gzip_exe, "-dc", str(path)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        assert process.stdout is not None
        wrapper = io.TextIOWrapper(process.stdout, encoding="utf-8", errors="replace")
        try:
            yield wrapper
        finally:
            wrapper.close()
            stderr = (
                process.stderr.read().decode("utf-8", errors="replace")
                if process.stderr
                else ""
            )
            return_code = process.wait()
            if return_code != 0:
                raise RuntimeError(f"Failed to decompress {path}: {stderr.strip()}")
        return
    with path.open(encoding="utf-8", errors="replace") as handle:
        yield handle


def _iter_imgt_entries(handle: Iterable[str]) -> Iterator[List[str]]:
    current: List[str] = []
    for line in handle:
        current.append(line.rstrip("\n"))
        if line.startswith("//"):
            yield current
            current = []
    if current:
        yield current


def _parse_location(location: str) -> Optional[Tuple[int, int, bool]]:
    complement = "complement" in location.lower()
    numbers = [int(value) for value in re.findall(r"\d+", location)]
    if not numbers:
        return None
    return min(numbers), max(numbers), complement


def _parse_qualifier(text: str) -> Optional[Tuple[str, str]]:
    match = re.match(r"/([^=\s]+)(?:=(.*))?$", text.strip())
    if not match:
        return None
    key = match.group(1)
    value = (match.group(2) or "true").strip().strip('"')
    return key, value


def _parse_imgt_entry(lines: Sequence[str]) -> Optional[ImgtRecord]:
    accession = ""
    species = ""
    definition = ""
    sequence_parts: List[str] = []
    features: List[Dict[str, object]] = []
    current_feature: Optional[Dict[str, object]] = None
    in_sequence = False

    for line in lines:
        if line.startswith("AC"):
            accession = line[2:].strip().split(";")[0].strip()
        elif line.startswith("OS"):
            species = line[2:].strip().rstrip(".")
        elif line.startswith("DE"):
            definition += " " + line[2:].strip()
        elif line.startswith("SQ"):
            in_sequence = True
        elif line.startswith("//"):
            in_sequence = False
        elif in_sequence:
            sequence_parts.extend(re.findall(r"[A-Za-z]+", line))
        elif line.startswith("FT"):
            payload = line[2:]
            match = re.match(r"\s{1,5}([A-Z0-9][A-Z0-9'_-]*)\s+(.+)$", payload)
            if match and not match.group(2).lstrip().startswith("/"):
                current_feature = {
                    "key": match.group(1),
                    "location": match.group(2).strip(),
                    "qualifiers": {},
                }
                features.append(current_feature)
            else:
                qualifier_text = payload.strip()
                if current_feature is not None and qualifier_text.startswith("/"):
                    parsed = _parse_qualifier(qualifier_text)
                    if parsed:
                        key, value = parsed
                        qualifiers = current_feature["qualifiers"]
                        assert isinstance(qualifiers, dict)
                        qualifiers.setdefault(key, []).append(value)

    sequence = "".join(sequence_parts).upper()
    sequence = re.sub(r"[^ACGTN]", "N", sequence)
    if not accession or not sequence:
        return None

    all_qualifiers: Dict[str, List[str]] = defaultdict(list)
    for feature in features:
        qualifiers = feature["qualifiers"]
        assert isinstance(qualifiers, dict)
        for key, values in qualifiers.items():
            all_qualifiers[key].extend(values)

    gene = ""
    for key in ("allele", "gene", "gene_name"):
        if all_qualifiers.get(key):
            gene = all_qualifiers[key][0]
            break
    if not gene:
        match = re.search(
            r"\b((?:IG[HKL]|TR[ABGD])[A-Z0-9-]+(?:\*\d+)?)\b",
            definition.upper(),
        )
        if match:
            gene = match.group(1)

    functionality = ""
    for key in ("functionality", "functional"):
        if all_qualifiers.get(key):
            functionality = all_qualifiers[key][0]
            break
    if not functionality and re.search(r"\bFUNCTIONAL\b", definition.upper()):
        functionality = "functional"

    splice_features = [
        feature
        for feature in features
        if feature["key"] in IMGT_DONOR_KEYS | IMGT_ACCEPTOR_KEYS
    ]
    if not splice_features:
        return ImgtRecord(accession, species, sequence, gene, functionality, "+", [])

    complement_count = sum(
        1
        for feature in splice_features
        if "complement" in str(feature["location"]).lower()
    )
    strand = "-" if complement_count > len(splice_features) / 2 else "+"
    oriented_sequence = reverse_complement(sequence) if strand == "-" else sequence
    sites: List[SpliceSite] = []

    for feature in splice_features:
        key = str(feature["key"])
        parsed_location = _parse_location(str(feature["location"]))
        if parsed_location is None:
            continue
        start, end, _ = parsed_location
        if strand == "-":
            oriented_start0 = len(sequence) - end
            oriented_end0 = len(sequence) - start + 1
        else:
            oriented_start0 = start - 1
            oriented_end0 = end
        length = oriented_end0 - oriented_start0
        if key in IMGT_DONOR_KEYS:
            site_type = "donor"
            index = oriented_start0 if length >= 3 else oriented_start0 - 1
            label = DONOR_LABEL
        else:
            site_type = "acceptor"
            index = oriented_start0 + 3 if length >= 5 else oriented_end0
            label = ACCEPTOR_LABEL
        if not 0 <= index < len(oriented_sequence):
            continue

        feature_qualifiers = feature["qualifiers"]
        assert isinstance(feature_qualifiers, dict)
        feature_gene = ""
        for gene_key in ("allele", "gene", "gene_name"):
            values = feature_qualifiers.get(gene_key)
            if values:
                feature_gene = values[0]
                break
        if not feature_gene:
            candidates = []
            for other in features:
                if other is feature:
                    continue
                other_location = _parse_location(str(other["location"]))
                if other_location is None:
                    continue
                other_start, other_end, _ = other_location
                other_qualifiers = other["qualifiers"]
                assert isinstance(other_qualifiers, dict)
                other_gene = ""
                for gene_key in ("allele", "gene", "gene_name"):
                    values = other_qualifiers.get(gene_key)
                    if values:
                        other_gene = values[0]
                        break
                if not other_gene:
                    continue
                if site_type == "donor":
                    distance = max(0, start - other_end)
                    direction_penalty = 0 if other_end <= end else 1_000_000
                else:
                    distance = max(0, other_start - end)
                    direction_penalty = 0 if other_start >= start else 1_000_000
                candidates.append((direction_penalty + distance, other_gene))
            if candidates:
                feature_gene = min(candidates, key=lambda value: value[0])[1]
        feature_gene = feature_gene or gene

        sites.append(
            SpliceSite(
                source="IMGT",
                accession=accession,
                gene=feature_gene,
                site_type=site_type,
                label=label,
                index=index,
                oriented_sequence=oriented_sequence,
                feature_key=key,
                functionality=functionality,
                strand=strand,
            )
        )

    return ImgtRecord(
        accession=accession,
        species=species,
        sequence=oriented_sequence,
        gene=gene,
        functionality=functionality,
        strand=strand,
        sites=deduplicate_sites(sites),
    )


def _is_mouse_species(species: str) -> bool:
    value = species.lower()
    return "mus musculus" in value or "house mouse" in value


def _is_functional(functionality: str) -> bool:
    value = functionality.strip().lower()
    return value in {"f", "functional", "true"} or value.startswith("functional")


def read_imgt_splice_records(
    path: Path | str,
    functional_only: bool = True,
    mouse_only: bool = True,
) -> List[ImgtRecord]:
    records: List[ImgtRecord] = []
    with open_imgt_text(path) as handle:
        for lines in _iter_imgt_entries(handle):
            record = _parse_imgt_entry(lines)
            if record is None or not record.sites:
                continue
            if mouse_only and not _is_mouse_species(record.species):
                continue
            if functional_only and not _is_functional(record.functionality):
                continue
            records.append(record)
    return records


def context_identity(first: str, second: str) -> float:
    if len(first) != len(second):
        return SequenceMatcher(None, first, second, autojunk=False).ratio()
    comparable = [(a, b) for a, b in zip(first, second) if a != "N" and b != "N"]
    if not comparable:
        return 0.0
    return sum(a == b for a, b in comparable) / len(comparable)


def deduplicate_sites(
    sites: Iterable[SpliceSite],
    context_radius: int = 20,
) -> List[SpliceSite]:
    seen = set()
    result: List[SpliceSite] = []
    for site in sites:
        key = (
            normalize_gene_name(site.gene),
            site.site_type,
            site.context(context_radius),
        )
        if key in seen:
            continue
        seen.add(key)
        result.append(site)
    return result


def find_missing_imgt_sites(
    imgt_records: Sequence[ImgtRecord],
    ncbi_sites: Sequence[SpliceSite],
    context_radius: int = 20,
    identity_threshold: float = 0.90,
) -> Tuple[List[ImgtRecord], Dict[str, object]]:
    """Keep IMGT sites not represented by an NCBI label-centred context."""
    ncbi_by_key: MutableMapping[Tuple[str, str], List[SpliceSite]] = defaultdict(list)
    for site in ncbi_sites:
        ncbi_by_key[(normalize_gene_name(site.gene), site.site_type)].append(site)

    missing_records: List[ImgtRecord] = []
    covered_count = 0
    missing_count = 0
    details: List[Dict[str, object]] = []

    for record in imgt_records:
        missing_sites: List[SpliceSite] = []
        for site in record.sites:
            candidates = ncbi_by_key[
                (normalize_gene_name(site.gene), site.site_type)
            ]
            query_context = site.context(context_radius)
            best_identity = max(
                (
                    context_identity(
                        query_context,
                        candidate.context(context_radius),
                    )
                    for candidate in candidates
                ),
                default=0.0,
            )
            is_covered = best_identity >= identity_threshold
            if is_covered:
                covered_count += 1
            else:
                missing_count += 1
                missing_sites.append(site)
            details.append(
                {
                    "accession": site.accession,
                    "gene": site.gene,
                    "site_type": site.site_type,
                    "feature_key": site.feature_key,
                    "motif": site.motif,
                    "best_ncbi_context_identity": round(best_identity, 4),
                    "covered_by_ncbi": is_covered,
                }
            )
        if missing_sites:
            missing_records.append(
                ImgtRecord(
                    accession=record.accession,
                    species=record.species,
                    sequence=record.sequence,
                    gene=record.gene,
                    functionality=record.functionality,
                    strand=record.strand,
                    sites=missing_sites,
                )
            )

    report: Dict[str, object] = {
        "ncbi_site_count": len(ncbi_sites),
        "imgt_record_count": len(imgt_records),
        "imgt_site_count": covered_count + missing_count,
        "covered_by_ncbi": covered_count,
        "supplemented_from_imgt": missing_count,
        "identity_threshold": identity_threshold,
        "context_radius": context_radius,
        "sites": details,
    }
    return missing_records, report


def imgt_records_to_data(records: Sequence[ImgtRecord]) -> List[List[str]]:
    data: List[List[str]] = [[] for _ in DATASET_NAMES]
    for record in records:
        labels = [0] * len(record.sequence)
        genes = set()
        for site in record.sites:
            labels[site.index] = site.label
            genes.add(site.gene)
        gene_text = ",".join(sorted(genes)) or record.gene or "unknown"
        values = (
            f"IMGT:{record.accession}:{gene_text}",
            f"IMGT_{record.accession}",
            "+",
            "1",
            str(len(record.sequence)),
            record.sequence,
            "".join(str(value) for value in labels),
        )
        for index, value in enumerate(values):
            data[index].append(value)
    return data


def read_datafile_h5(path: Path | str) -> List[List[str]]:
    result: List[List[str]] = []
    with h5py.File(path, "r") as handle:
        for name in DATASET_NAMES:
            values = []
            for value in handle[name][:]:
                values.append(
                    value.decode("utf-8") if isinstance(value, bytes) else str(value)
                )
            result.append(values)
    return result


def write_datafile_h5(path: Path | str, data: Sequence[Sequence[str]]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    text_dtype = h5py.string_dtype(encoding="utf-8")
    with h5py.File(path, "w") as handle:
        for name, values in zip(DATASET_NAMES, data):
            handle.create_dataset(
                name,
                data=np.asarray(list(values), dtype=text_dtype),
                dtype=text_dtype,
            )


def merge_data(
    base: Sequence[Sequence[str]],
    supplement: Sequence[Sequence[str]],
) -> List[List[str]]:
    if len(base) != len(DATASET_NAMES) or len(supplement) != len(DATASET_NAMES):
        raise ValueError("Unexpected OpenSpliceAI datafile structure")
    return [list(left) + list(right) for left, right in zip(base, supplement)]


def split_imgt_records_by_gene(
    records: Sequence[ImgtRecord],
    validation_ratio: float = 0.1,
    seed: str = "openspliceai-imgt-v1",
) -> Tuple[List[ImgtRecord], List[ImgtRecord]]:
    """Gene-grouped deterministic split; IMGT never enters test by default."""
    if not 0 <= validation_ratio < 1:
        raise ValueError("validation_ratio must be in [0, 1)")
    train: List[ImgtRecord] = []
    validation: List[ImgtRecord] = []
    threshold = int(validation_ratio * 10_000)
    for record in records:
        genes = sorted({normalize_gene_name(site.gene) for site in record.sites})
        group = "|".join(genes) or record.accession
        digest = hashlib.sha256(f"{seed}:{group}".encode("utf-8")).hexdigest()
        bucket = int(digest[:8], 16) % 10_000
        (validation if bucket < threshold else train).append(record)
    return train, validation


def write_audit_report(path: Path | str, report: Mapping[str, object]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
