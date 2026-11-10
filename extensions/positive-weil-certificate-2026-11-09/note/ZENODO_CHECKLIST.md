# ZENODO_CHECKLIST.md — manual upload checklist

**Nothing here has been uploaded.** No Zenodo record exists for this note, no
DOI has been reserved, and no GitHub Release or automatic archival has been
triggered. This file is the manual procedure, to be run by a person who decides
to deposit.

## Before uploading

- [ ] Run the fast check (`README.md` → *Reproducing*). All seven verifiers must
      exit `0`. Do not deposit on a failing tree.
- [ ] Rebuild the PDF and confirm it is 14 pages with no unresolved references
      and no overfull boxes.
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
- [ ] Version: `1.0.0`.
- [ ] Related identifiers:
      - `isSupplementedBy` → `https://github.com/innerlightr-wq/weil-window-control`
      - `references` → `arXiv:2608.24827`
- [ ] Add the repository link (branch `positive-weil-certificate-2026-11-09`) so
      a reader can reach the certificates and verifiers.

## Files to deposit

- [ ] `weil_L45_certified_positivity.pdf` — the manuscript.
- [ ] `weil_L45_certified_positivity.tex` — its source.
- [ ] `weil-L45-certificate-v1.0.0.tar.gz` — the reproducibility archive.
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
