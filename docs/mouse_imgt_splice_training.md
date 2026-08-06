# Mouse OpenSpliceAI training with NCBI + IMGT splice annotations

## Goal

Keep the author's normal NCBI GRCm39 training data, then add **only functional
mouse IG/TR donor and acceptor sites that NCBI-derived labels do not already
cover**.

This matters because the normal data builder labels splice sites only between
adjacent exons of one transcript-like feature. A single-exon `J_gene_segment`
has no downstream exon, so its biologically real J donor is not labelled. The
first exon of a `C_gene_segment` has no upstream exon inside that feature, so
its first acceptor is also not labelled. IMGT/LIGM-DB explicitly annotates such
positions as:

- `DONOR-SPLICE`
- `ACCEPTOR-SPLICE`
- `INT-DONOR-SPLICE`
- `INT-ACCEPTOR-SPLICE`

OpenSpliceAI learns a per-base class (`0` non-splice, `1` acceptor, `2` donor).
It does not learn a separate categorical "IGHJ1 pairs with IGHM" label. The
correct way to expose the J-C biology to this model is therefore to include the
explicit J donor and C acceptor bases, with their sequence context.

## Data policy

1. Download the current NCBI package for GRCm39 (`GCF_000001635.27`). The
   package content and annotation metadata are current at download time.
2. Download current IMGT/GENE-DB metadata and the current weekly
   IMGT/LIGM-DB `imgt.dat.Z` flat file.
3. Recover NCBI immune splice labels from:
   - explicit splice-site features, when present;
   - adjacent exons under `mRNA`, `transcript`, `V_gene_segment`,
     `D_gene_segment`, `J_gene_segment`, and `C_gene_segment`.
4. Parse only explicit IMGT splice labels from functional mouse records.
5. Match NCBI and IMGT by normalized gene name, donor/acceptor class, and
   label-centred sequence context. Add an IMGT site only when no NCBI site
   reaches the configured context identity threshold.
6. Add supplemental records to train/validation using a deterministic
   gene-grouped split. Keep the original NCBI test set unchanged.
7. Write a JSON audit containing every IMGT site, its motif, its best NCBI
   context identity, and whether it was supplemented.

No J-C pair is invented merely because two genes are close on the chromosome.
Only a site explicitly annotated by NCBI or IMGT becomes a positive label.

## 1. Download current data

```bash
mamba install -c conda-forge ncbi-datasets-cli

python scripts/download_mouse_immune_training_data.py \
  --output-dir data/mouse_immune
```

The command writes `data/mouse_immune/manifest.json` with:

- NCBI assembly accession and Datasets CLI version;
- NCBI FASTA/GFF3/GBFF paths;
- IMGT/GENE-DB and IMGT/LIGM-DB release identifiers;
- file sizes and SHA-256 checksums.

The large downloaded files are not committed to Git.

## 2. Build the normal NCBI baseline

Use the FASTA and GFF3 paths printed by the download script:

```bash
openspliceai create-data \
  --annotation-gff <path-from-manifest> \
  --genome-fasta <path-from-manifest> \
  --output-dir work/mouse_ncbi \
  --parse-type all_isoforms \
  --biotype protein-coding \
  --chr-split train-test \
  --split-method random \
  --canonical-only
```

Use the exact paths printed by the downloader or read them from
`manifest.json`; NCBI package basenames can change.

## 3. Audit NCBI and supplement from IMGT

```bash
python scripts/prepare_mouse_immune_training.py \
  --ncbi-gff <path-from-manifest> \
  --ncbi-fasta <path-from-manifest> \
  --imgt-flat data/mouse_immune/imgt/ligm_db/imgt.dat.Z \
  --base-data-dir work/mouse_ncbi \
  --output-dir work/mouse_ncbi_imgt
```

Important outputs:

```text
work/mouse_ncbi_imgt/
├── ncbi_imgt_splice_audit.json
├── imgt_supplement_train.h5
├── imgt_supplement_validation.h5
├── datafile_train.h5
├── datafile_validation.h5
├── datafile_test.h5          # unchanged NCBI test set
├── dataset_train.h5
├── dataset_validation.h5
└── dataset_test.h5
```

## 4. Train

Use the augmented datasets exactly as normal OpenSpliceAI inputs:

```bash
openspliceai train \
  --project-name mouse_ncbi_imgt \
  --exp-num 1 \
  --output-dir models/mouse_ncbi_imgt \
  --train-dataset work/mouse_ncbi_imgt/dataset_train.h5 \
  --test-dataset work/mouse_ncbi_imgt/dataset_test.h5 \
  --flanking-size 2000 \
  --epochs 20
```

## Quality checks before a long training run

Inspect `ncbi_imgt_splice_audit.json` and verify:

- `supplemented_from_imgt` is greater than zero;
- functional IGHJ donor and IGHC first-acceptor records are present;
- canonical donor motifs are usually `GT` and acceptor motifs usually `AG`;
- noncanonical sites are retained only when IMGT explicitly annotates them;
- the test-set policy says that no IMGT examples were added.

For a first run, use `--skip-dataset` to generate and inspect only the merged
`datafile_*.h5` files before one-hot encoding.
