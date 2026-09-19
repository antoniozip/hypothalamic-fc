# Superseded reports — historical record only

**Do not cite any number in this folder.** These documents were written against output from a
C GLMCC port that inverted the sign of every coupling. They are kept because they are the audit
trail for why the previous body of results was discarded, and because a reviewer asking "why did
your connectivity change?" needs to be able to follow it.

The current record is `METHODS_AND_RESULTS.md` at the repository root.

| file | what it was | why it is superseded |
|---|---|---|
| `glmcc_positive_control.md` | first positive control | fixture was 1-indexed; then rerun through the defective port |
| `glmcc_python_vs_c.md` | Python vs C comparison | correctly identified that the port was at fault, wrong about the mechanism |
| `glmcc_c_sign_bug.md` | bug hunt | superseded by the root cause; two earlier fixes in it were wrong guesses |
| `glmcc_regeneration.md` | first regeneration | run with the defective port |
| `glmcc_regeneration_postfix.md` | regeneration after the fix | folded into `METHODS_AND_RESULTS.md` |
| `downstream_postfix.md` | downstream after the fix | folded into `METHODS_AND_RESULTS.md` |
| `corrected_lme_results.md` | duty-cycle-matched LME | built on defective matrices |
| `glmcc_sign_asymmetry_diagnosis.md` | duty-cycle diagnosis | its headline contrast compared a Python-derived arm against a defective-C arm |
| `glmcc_dutycycle_experiment/` | duty-cycle experiment bundle | same |

## The one thing worth carrying forward

`glmcc_sign_asymmetry_diagnosis.md` section 3 established, from spike times alone and
independently of any GLMCC output, that the `lightON` epoch is a ~10% duty-cycle concatenation
of ~108 stimulus windows while `ongoing` is ~89% continuous, and that per-observed-time firing
rates are comparable between the two (lightON slightly higher in every animal). That measurement
does not depend on the estimator and still holds.

Its *inference* from that fact — the mechanism for the sign asymmetry — does not, because the
asymmetry it was explaining was partly an artifact of the two arms having been computed by
different implementations.
