#!/usr/bin/env python3
"""PROVENANCE ONLY -- NOT RUNNABLE IN THE INTEGRATED TREE.

This generated notes/proposed_2026-11-12.patch against the baseline commit
7515396c7f2aa17b53d565a4f11f8e937ee164dc before the phase-7 work was
integrated.  The patch is now applied, so re-running this would fail on its
anchors.  It is kept, with its section sources in ../patch_sections/, so the
documentation edits of this release can be audited edit by edit.
The release gate is ../acceptance.py.
"""
raise SystemExit(__doc__)

# ---- original generator, retained verbatim below for audit ----
# #!/usr/bin/env python3
# """Regenerate notes/proposed_2026-11-12.patch from the clean worktree.
# 
# Applies every proposed edit, exports `git diff`, then restores the worktree.
# Idempotent and reproducible: run it instead of hand-editing.
# """
# import os, subprocess, sys
# 
# WT = os.path.expanduser("~/scratch/rh-connes-continuation-2026-11-12/wt/"
#                         "extensions/positive-weil-certificate-2026-11-09")
# OUT = os.path.expanduser("~/scratch/rh-connes-continuation-2026-11-12/"
#                          "notes/proposed_2026-11-12.patch")
# SEC = os.path.expanduser("~/scratch/rh-connes-continuation-2026-11-12/work/patch_sections")
# os.chdir(WT)
# subprocess.run(["git", "checkout", "--", "."], check=True)
# 
# def sub(path, pairs):
#     s = open(path).read()
#     for a, b in pairs:
#         if a not in s:
#             sys.exit("ANCHOR MISSING in %s: %r" % (path, a[:60]))
#         s = s.replace(a, b, 1)
#     open(path, "w").write(s)
# 
# def append(path, section):
#     open(path, "a").write(open(os.path.join(SEC, section)).read())
# 
# # 1. DERIVATION.md -- the divergent archimedean series
# sub("DERIVATION.md", [
#  ("    Arch_{jk}  = -gamma C(0) + sum_{m>=0}[ C(0)/(m+1) - 2 int_0^{2L} C e^{-(m+1/2)u} du ]",
#   "    Arch_{jk}  = -(gamma + log 4) C(0)\n"
#   "                 + 2 sum_{m>=0}[ C(0)/(2m+1) - int_0^{2L} C(u) e^{-(2m+1/2)u} du ]"),
#  ("Two facts matter for error accounting and are easy to lose:",
#   open(os.path.join(SEC, "derivation_note.md")).read()
#   + "Two facts matter for error accounting and are easy to lose:"),
# ])
# 
# # 2. note/CLAIM_INDEX.md -- the discard ratios are measurements (R14)
# sub("note/CLAIM_INDEX.md", [
#  ("* the discard ratios (see Section 5) — they are certified *enclosures*, but no\n"
#   "  decay law is inferred from them and they enter no bound.",
#   open(os.path.join(SEC, "claim_index_ratios.md")).read().rstrip("\n")),
# ])
# 
# # 3. note/README.md -- the check counts
# sub("note/README.md", [
#  ("$V phase5/test_prec.py     # 12 checks", "$V phase5/test_prec.py     # 17 checks"),
#  ("$V note/consistency_check.py # 47 checks", "$V note/consistency_check.py # 52 checks"),
#  ("| `phase5/test_prec.py` | 12 | 0 s |", "| `phase5/test_prec.py` | 17 | 0 s |"),
#  ("| `note/consistency_check.py` | 47 | 0 s |", "| `note/consistency_check.py` | 52 | 0 s |"),
#  ("for **166 checks** in total, all passing.", "for **176 checks** in total, all passing."),
# ])
# 
# # 4. new research sections
# append("REPORT.md", "report_sec12.md")
# append("REPORT.md", "report_sec13.md")
# append("REPORT.md", "report_sec14.md")
# append("REPORT.md", "report_sec15.md")
# append("RETRACTIONS.md", "retractions_r15.md")
# 
# # 5. README.md headline banner
# sub("README.md", [
#  ("> # **CERTIFIED: `1.41e-17 <= lambda_Q(4/5) <= 1.674e-17`**",
#   open(os.path.join(SEC, "readme_banner.md")).read().rstrip("\n")),
# ])
# 
# diff = subprocess.run(["git", "diff"], capture_output=True, text=True, check=True).stdout
# open(OUT, "w").write(diff)
# stat = subprocess.run(["git", "diff", "--stat"], capture_output=True, text=True).stdout
# subprocess.run(["git", "checkout", "--", "."], check=True)
# left = subprocess.run(["git", "status", "--porcelain"], capture_output=True, text=True).stdout
# print(stat)
# print("patch -> %s (%d lines)" % (OUT, diff.count("\n")))
# print("worktree clean after restore: %s" % ("yes" if not left.strip() else "NO: " + left))
