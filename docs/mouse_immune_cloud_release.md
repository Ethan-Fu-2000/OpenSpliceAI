# Cloud download and GitHub Release storage

This repository keeps source code in Git and stores the large mouse NCBI + IMGT
reference package as GitHub Release assets.

## What the workflow does

`.github/workflows/build-mouse-immune-data.yml` runs on GitHub-hosted Ubuntu and:

1. installs the official NCBI Datasets CLI;
2. downloads the current annotation package for GRCm39
   (`GCF_000001635.27`), including genome FASTA, GFF3, GBFF and sequence report;
3. downloads current IMGT/GENE-DB metadata and the current IMGT/LIGM-DB flat file;
4. writes `manifest.json` with release identifiers, paths, sizes and SHA-256 values;
5. removes the redundant original NCBI package ZIP after extraction;
6. creates a Zstandard-compressed archive;
7. splits the archive into assets smaller than 2 GiB;
8. publishes the parts and checksum files to a versioned GitHub Release.

Release tags use this form:

```text
mouse-immune-data-YYYYMMDD-GITHUB_RUN_ID
```

The workflow has `contents: write` only because creating a Release and tag
requires repository content write permission. It does not open a pull request or
modify the upstream OpenSpliceAI repository.

## Trigger behavior on this feature branch

The workflow currently exists only on `feat/imgt-splice-annotations`.
`workflow_dispatch` can be used from the GitHub Actions interface only after the
workflow also exists on the repository default branch. Until then, the workflow
runs automatically when its workflow file or related download helpers are pushed
to this feature branch.

## Download the newest release locally

Install GitHub CLI and Zstandard, authenticate `gh`, then run:

```bash
bash scripts/download_mouse_immune_release.sh data
```

The helper will:

- select the newest release whose tag begins with `mouse-immune-data-`;
- download every numbered archive part;
- verify `PARTS_SHA256SUMS`;
- reconstruct and verify the complete archive;
- extract it to `data/mouse_immune`.

To download a specific release:

```bash
bash scripts/download_mouse_immune_release.sh \
  data \
  mouse-immune-data-YYYYMMDD-RUN_ID
```

The repository is public, so Release assets created here are public as well.
The package retains source and release metadata for both NCBI and IMGT.
