* the discard ratios (see Section 5) — **measurements, not enclosures.** They
  are quotients of `e_T = q - r_T`, and `q` comes from the time-domain
  evaluator whose outer quadrature remainder is unbounded (`RETRACTIONS.md`
  R14); an arithmetic-only radius on such a quantity is not an enclosure. The
  checker's `D1`-`D3` verify the *arithmetic of the division*, which is correct;
  that is not a certification of the inputs. No decay law is inferred from them
  and they enter no bound.
