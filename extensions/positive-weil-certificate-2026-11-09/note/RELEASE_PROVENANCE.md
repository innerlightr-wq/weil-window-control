# Release provenance — version 2.0.0

**Published 2026-11-12.** DOI [10.5281/zenodo.24171237](https://doi.org/10.5281/zenodo.24171237) ·
[Zenodo record](https://zenodo.org/records/24171237) · release tag
`two-window-note-v2.0.0`.

## The immutable publication snapshot

The deposited manuscript was built from, and is byte-identical to, the file
committed at

```
    pre-DOI publication commit   575b591516dd98520b5236f9ca96bbec6a60f8cf
      "Quote the upper bound through a safe endpoint, not the enclosure midpoint"
    parent                       de2b8c155f68b8b0d5e2912c3b338d2cba4ad604
    integrated candidate         05601e2d89378b77d125544f62a4ceb199713331
    base                         7515396c7f2aa17b53d565a4f11f8e937ee164dc
                                   (= origin/positive-weil-certificate-2026-11-09)
    branch                       two-window-note-2026-11-12
```

**Use `575b591` for anything that must reproduce the deposited bytes**, not the
tagged commit. The tag points at a later commit that adds post-deposit citation
metadata only (see below).

## Artifact hashes at `575b591`

| file | sha256 | notes |
|---|---|---|
| `note/weil_two_window_positivity.pdf` | `22a6d35d3216cf528e1c1a252b7052fc2f735212fecec8b7bc3ec942e22f8159` | **deposited**; `md5 5c2981418c515504312f2845e40baf85`, 520469 bytes |
| `note/weil_two_window_positivity.tex` | `ee17ddb162f5d903bef1844fe5803ad3674fd996e5fc32088633c6034706ef34` | not deposited |
| `note/weil-two-window-certificate-2.0.0.tar.gz` | `322dc2617e25f222e46bc1648d6c7f52788b197d049827172d1c6bf733b187d5` | not deposited; gitignored build product |
| `note/SHA256SUMS.txt` | `06aa8cbae2056c808083aa076ecc413575e7b289f825ac8143604550dc787717` | not deposited; describes the tree **at `575b591`** |

### What the Zenodo deposit actually contains

**The manuscript PDF only** — `weil_two_window_positivity.pdf`. Its MD5 on the
public record (`5c2981418c515504312f2845e40baf85`) was compared against the
local file and **matches**, as does the size. The LaTeX source, the
reproducibility archive and the checksum manifest were **not** uploaded; they
are distributed from this repository at tag `two-window-note-v2.0.0`.

The record's own version string reads **`v1`** (Zenodo's default label for a
first upload). The manuscript's version is **2.0.0**, and the repository
metadata states 2.0.0 throughout. The two refer to the same deposit; no attempt
is made here to reconcile them by editing the published record.

No separate concept (all-versions) DOI is displayed on the record, so none is
claimed anywhere in this repository.

## What the tagged commit adds

The tag `two-window-note-v2.0.0` points at a commit **after** the deposit, which
changes **citation metadata only**:

```
    README.md                     (repository root: a findable entry + DOI link)
    CITATION.cff                  (repository root: a related-work entry only;
                                   the primary-work identity is unchanged)
    note/CITATION.cff             (doi, url, identifiers)
    note/CITATION.bib             (doi, url)
    note/zenodo.json              (doi, publication_date)
    note/README.md                (status block -> published)
    README.md (extension)         (status -> published)
    note/ZENODO_CHECKLIST.md      (marked done)
    note/RELEASE_PROVENANCE.md    (this file)
    note/SHA256SUMS-post-doi.txt  (manifest of the post-DOI tree)
    phase7/HANDOFF.md             (DOI step closed)
```

No manuscript source, no certificate, no verifier, no log and no frozen datum is
touched. The scientific content is frozen at `575b591`.

## Two manifests, two scopes

* **`note/SHA256SUMS.txt`** — unchanged, and describes the tree at `575b591`,
  i.e. the publication snapshot. Running `sha256sum -c` on it from the **tagged**
  commit reports **seven** mismatches, and they are exactly the seven
  citation-metadata files changed by the metadata commit: the extension
  `README.md`, `note/CITATION.bib`, `note/CITATION.cff`, `note/README.md`,
  `note/ZENODO_CHECKLIST.md`, `note/zenodo.json` and `phase7/HANDOFF.md`. That
  is expected, not drift: no manuscript, certificate, verifier, script or frozen
  datum differs.

  One inaccuracy is inherited rather than fixed. `SHA256SUMS.txt`'s header says
  it covers "every research log", but its generator excluded every `*.log` path
  in order to skip the LaTeX logs, so the logs under `*/logs/` are in fact
  **absent** from it. The file is left byte-identical because it is the record
  of the publication snapshot; the corrected coverage, including the research
  logs, is in `SHA256SUMS-post-doi.txt`.
* **`note/SHA256SUMS-post-doi.txt`** — the same coverage recomputed for the
  post-DOI tree. Use this one to check out the tag; use the other one to check
  the deposit.

## Rebuilding the deposited archive

`note/make_archive.sh` is deterministic given

```sh
export SOURCE_DATE_EPOCH=1794528000 FORCE_SOURCE_DATE=1
```

and reproduces `322dc2617e…` **from `575b591`** — verified by exporting that
commit with `git archive` and rebuilding. It does **not** reproduce that hash
from the tagged commit, because the tagged commit's metadata files differ and
are included in the archive. That is the reason the pre-DOI commit is recorded
above.
