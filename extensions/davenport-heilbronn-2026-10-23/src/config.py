r"""Single recorded configuration for the Davenport-Heilbronn extension.

Every script imports its conventions from here; nothing numeric is hard-coded twice."""

# ------------------------------------------------------------------------- precision
DPS_BASE = 60          # the working precision ("at least 60 decimal digits")
DPS_HIGH = 100         # the rerun precision for every reported value
DPS_XHIGH = 160        # third precision, used for the frozen witness only

# ------------------------------------------------------------------------- the function
CONDUCTOR = 5
PARITY = 1             # chi(-1) = -1, gamma factor Gamma((s+1)/2)

# the selected off-line quartet O_1 (search seed; refined in Phase A)
SEED_BETA = '0.8085171824566374'
SEED_GAMMA = '85.6993484853776'

# ------------------------------------------------------------------------- windows
WINDOW_SCAN = ['1.0', '1.2', '1.4', '1.6', '1.8', '2.0', '2.2', '2.4', '2.6']
GATE_WINDOWS = ['1.0', '1.5', '2.0']

# ------------------------------------------------------------------------- bases
# Both sectors use bases that VANISH at +-L, so g is continuous on R and G(it) = O(t^-2);
# this is what makes the horizontal contour segments in the explicit formula negligible
# and the zero-side sum absolutely convergent.
SECTOR_EVEN = 'even_d'      # cos((2k+1) pi x / 2L)
SECTOR_ODD = 'odd'          # sin((k+1) pi x / L)

LOW_BLOCK = 8               # modes 0..LOW_BLOCK-1  (the low-frequency block)
RES_HALFWIDTH = 4           # resonant block: k0 +- RES_HALFWIDTH around w_k ~ gamma_1

# ------------------------------------------------------------------------- numerics
GRAM_COND_WARN = 1e12       # report, never silently regularise
ZERO_SCAN_TMAX = 200.0
ZERO_SCAN_STEPS_PER_GAP = 10

SEED = 20261023             # no randomness is used; recorded for completeness

# ------------------------------------------------- witness selection (src/freeze.py)
# A detected negative eigenvalue of size 1e-21 is a real negative eigenvalue of the
# computed matrix, but it is not something a ball-arithmetic certificate can be built
# around at an affordable precision, and it is not what "the detector sees the off-line
# zero" means in any useful sense.  The certifiable witness is therefore required to be
# macroscopic, and the two are frozen separately and reported separately.
MACRO_THRESHOLD = '1e-3'
CERT_MAX_N = 40            # basis-size cap for the Arb re-computation
