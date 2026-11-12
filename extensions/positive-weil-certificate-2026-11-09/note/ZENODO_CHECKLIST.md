# Zenodo upload checklist — COMPLETED

> **Published 2026-11-12 — [doi:10.5281/zenodo.24171237](https://doi.org/10.5281/zenodo.24171237).** Version 2.0.0,
> first archived release. The deposit contains the manuscript PDF only
> (`md5 5c2981418c515504312f2845e40baf85`, 520469 bytes), verified against the
> public record. The archive, LaTeX source and checksum manifest are
> distributed from this repository at tag `two-window-note-v2.0.0`.
>
> Nothing below needs doing again. It is kept as the record of what was checked
> before the upload, and as the procedure for any future version. For the
> publication snapshot and the hashes, see `RELEASE_PROVENANCE.md`.

**Before uploading**

1. `cd extensions/positive-weil-certificate-2026-11-09 && python phase7/acceptance.py`
   must report **61 checks, 0 failed** and exit `0`.
2. `cd note && ./make_archive.sh`, then verify `SHA256SUMS.txt` matches the files
   you are about to upload (`sha256sum -c SHA256SUMS.txt`).
3. Confirm no DOI field has been pre-filled anywhere:
   `grep -rn "10\.5281/zenodo" note/CITATION.* note/zenodo.json` should only show
   the two *referenced* earlier notes, never a DOI for this one.

**Upload** the five files listed under FILES TO UPLOAD in `ZENODO_METADATA.txt`,
with the metadata from that file.

**After uploading** — see `../phase7/HANDOFF.md`. Do not edit the citation
metadata before the DOI is in hand, and preserve the identity of the deposited
files when you do.

---


> # ⛔ PUBLICATION ON HOLD — DO NOT UPLOAD
>
> **Status: DRAFT — certification correction in progress.** This package must
> not be deposited to Zenodo, released on GitHub, or designated an archival
> release, in its current state.
>
> Reason: the compact time-domain evaluator omitted an outer quadrature
> truncation allowance, so the improved upper bound of versions 1.0.x and every
> quantity derived from it are **not currently certified** (`../RETRACTIONS.md`
> R14). The lower-bound certificate is a separate result and is unaffected.
>
> **Two gates, both required before any upload:**
>
> 1. the outer quadrature remainder is derived, implemented and verified, or the
>    affected claims are permanently withdrawn rather than pending;
> 2. **the author has explicitly authorized publication.** The manuscript stays
>    marked *"Draft — awaiting author review"* until then, even after the repair
>    passes its checks.
>
> Until both gates are met, treat every box below as blocked.

**Nothing here has been uploaded.** No Zenodo record exists for this note, no
DOI has been reserved, and no GitHub Release or automatic archival has been
triggered. This file is the manual procedure, to be run by a person who decides
to deposit.

## Before uploading

- [ ] Run the fast check (`README.md` → *Reproducing*). All seven verifiers must
      exit `0`. Do not deposit on a failing tree.
- [ ] Rebuild the PDF and confirm it is 18 pages with no unresolved references
      and no overfull boxes, and that both Zenodo DOIs in the bibliography are
      live links.
- [ ] Re-run `sha256sum -c SHA256SUMS.txt` from this directory. Every line must
      report `OK`.
- [ ] Rebuild the archive with `./make_archive.sh` and extract it into an empty
      directory; run `consistency_check.py` there to confirm the package is
      self-contained.
- [ ] Confirm the archive contains no third-party material. In particular it
      must not contain any copy of arXiv:2608.24827, whose reproducibility
      materials were never obtained.

## Metadata

- [ ] Upload type: **Publication** → **Technical note**.
- [ ] Title, creator, ORCID, description, keywords and license exactly as in
      `zenodo.json`. That file is the source of truth for the form; copy from it
      rather than retyping.
- [ ] License: **CC BY 4.0** for the deposited manuscript. Note in the
      description that the accompanying code is MIT.
- [ ] Version: `2.0.0`.
- [ ] Related identifiers:
      - `isSupplementedBy` → `https://github.com/innerlightr-wq/weil-window-control`
      - `references` → `arXiv:2608.24827`
      - `references` → `10.5281/zenodo.23250932` (the sensitivity note)
      - `references` → `10.5281/zenodo.23624161` (the real-zero parity note, v2)
- [ ] Add the repository link (branch `positive-weil-certificate-2026-11-09`) so
      a reader can reach the certificates and verifiers.

## Files to deposit

- [ ] `weil_two_window_positivity.pdf` — the manuscript.
- [ ] `weil_two_window_positivity.tex` — its source.
- [ ] `weil-two-window-certificate-2.0.0.tar.gz` — the reproducibility archive.
- [ ] `SHA256SUMS.txt` — so a reader can check what they downloaded.

## DOI handling — read before clicking Publish

- [ ] **Do not reuse any existing DOI for this note.** In particular
      `10.5281/zenodo.23195177` belongs to the repository's separate sensitivity
      paper (*A δ²–L³ Sensitivity Law for Finite-Window Weil Positivity…*). That
      is a different work. The repository's root `CITATION.cff` refers to it and
      must be left as the repository's preferred citation.
- [ ] A DOI is minted by Zenodo **at publication**. Until then the note carries
      no identifier, and neither the PDF nor `CITATION.cff` prints one. Do not
      insert a placeholder or a "DOI pending" string that reads like a link.

## After uploading

- [ ] Record the **version DOI** and the **concept DOI** separately; they are
      different identifiers and the concept DOI is the one that tracks future
      versions.
- [ ] Add the version DOI to `CITATION.cff` and `CITATION.bib` in this
      directory, and to the `\date{}` line of the `.tex`, then rebuild the PDF.
      Do not change the repository's root `CITATION.cff`.
- [ ] Commit those edits as a separate, clearly labelled commit so the
      pre-deposit state stays legible in the history.
- [ ] Re-run `sha256sum` and update `SHA256SUMS.txt`, since the PDF changed.

## What not to claim, in the record or anywhere else

- Not a new arithmetic estimate, and not a priority claim.
- Not a new positivity framework: the reduction strategy is arXiv:2608.24827v2's
  and parity decoupling is likewise not new here.
- Not a verified range of zeta zeros, and not a result at any support larger
  than `L = 4/5`.
- Not peer-reviewed.
- No suggestion that any outside researcher uses, has used, or has endorsed this
  work.
