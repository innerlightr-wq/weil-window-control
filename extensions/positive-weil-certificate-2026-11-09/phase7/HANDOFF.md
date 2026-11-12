# Handoff — the DOI step is CLOSED

**Published 2026-11-12 — DOI [10.5281/zenodo.24171237](https://doi.org/10.5281/zenodo.24171237)**
([record](https://zenodo.org/records/24171237)), version 2.0.0, released at tag
`two-window-note-v2.0.0`.

## What was done

1. The note was finalized at 2.0.0 and the research package integrated under
   `phase7/`; the release gate passes (see `phase7/logs/acceptance_integrated.log`).
2. The manuscript PDF was deposited. Its MD5 on the public record,
   `5c2981418c515504312f2845e40baf85`, was compared against the local file and
   matches, as does the size (520469 bytes). **The deposit holds the PDF only** —
   the LaTeX source, the reproducibility archive and `SHA256SUMS.txt` were not
   uploaded and are distributed from this repository at the tag.
3. A post-deposit commit added the DOI to the citation metadata and the
   repository READMEs, and nothing else. The branch and the annotated tag were
   pushed.

**The publication snapshot is commit `575b591516dd98520b5236f9ca96bbec6a60f8cf`**,
not the tagged commit. Anything that must reproduce the deposited bytes must use
`575b591`. `note/RELEASE_PROVENANCE.md` records every hash and both manifests.

## Known, deliberately unreconciled

* The Zenodo record's version string reads **`v1`** (its default for a first
  upload) while the manuscript's version is **2.0.0**. The published record was
  not edited. Both refer to the same deposit.
* **No concept (all-versions) DOI** is displayed on the record, so none is
  claimed anywhere in this repository. If Zenodo later exposes one, it may be
  added to `note/zenodo.json` and the two `CITATION.*` files — nowhere else.
* Out of scope and unchanged: the irrational-window reparameterization and the
  interval enclosure of the prolate rungs `delta_{4j}`. Every ratio against
  those rungs remains labelled numerical.

## Repository state to be aware of

The main checkout at `~/GitHub/weil-window-control` is on branch
`dh-note-ziad-citation-2026-11-10` and has **unrelated** uncommitted work that
this task did not touch: a modified `README.md` and seven untracked paths
(`extensions/newtonian-weil-response-2026-10-08/`,
`extensions/odd-sector-residual-margin-2026-11-02/`,
`extensions/spectral-isometry-transfer-2026-11-02/`,
`extensions/suzuki-function-field-window-2026-11-02/`, two DH note PDFs and two
`paper/` latexmk files). Leave them alone or commit them separately; they are
not part of this release.

The release work was done in a separate worktree at
`~/scratch/rh-connes-continuation-2026-11-12/wt` on branch
`two-window-note-2026-11-12`, branched from
`7515396c7f2aa17b53d565a4f11f8e937ee164dc`.
