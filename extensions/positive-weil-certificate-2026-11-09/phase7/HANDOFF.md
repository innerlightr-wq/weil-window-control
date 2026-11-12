# Handoff — the DOI and tag step, NOT done here

**State at handoff.** The note is finalized at version 2.0.0, the research
package is integrated under `phase7/`, the release gate passes (61 checks, 0
failed, last run from a clean extraction of the archive with
`HOME=/nonexistent`), and the work is committed to the **local** branch
`two-window-note-2026-11-12`.

**Nothing was pushed. No tag was created, locally or remotely. No GitHub release
was made. Nothing was uploaded to Zenodo. No DOI was invented or reserved.**

## What is waiting on you

1. Upload to Zenodo the five files listed in `note/ZENODO_METADATA.txt`, using
   the metadata in that file (`note/zenodo.json` is the machine-readable copy).
2. Send back the minted DOI (both the version DOI and the concept DOI, if Zenodo
   issues both).

## What happens after you supply the DOI — a separate task, not started

1. Add the DOI to the citation metadata, in these places and no others:
   - `note/CITATION.cff` — add a top-level `doi:` field (it is deliberately
     absent now);
   - `note/CITATION.bib` — add a `doi = {...}` field to
     `@techreport{DeJesus2026WeilTwoWindow,...}` and remove the comment that
     forbids it;
   - `note/zenodo.json` — record the minted identifier;
   - `note/README.md` — replace the "no DOI has been assigned" status block;
   - `extensions/positive-weil-certificate-2026-11-09/README.md` — same.

   The repository **root** `CITATION.cff` refers to the separate sensitivity
   paper and must stay as it is.

2. **Preserve the identity of the deposited files.** The DOI edits above touch
   documentation only. If any of them changes a file that was uploaded, re-upload
   that file or publish the change as a new Zenodo version — do not leave the
   deposit and the repository disagreeing. The simplest path is: upload first,
   then make the DOI edits as a follow-up commit, and do not rebuild the archive
   or the PDF in that commit. `SHA256SUMS.txt` records the hashes of the exact
   uploaded files; check them against the deposit before editing anything.

3. Push and tag. Proposed tag name, **not created**:

   ```
   two-window-note-v2.0.0
   ```

   matching the repository's existing per-project convention (`dh-note-v1`,
   `lambda-note-v1`, `siegel-note-v1`, `dirichlet-windows-v1`) while carrying
   the manuscript version. Push `two-window-note-2026-11-12` first, then tag.

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
