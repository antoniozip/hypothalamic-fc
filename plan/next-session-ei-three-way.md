# Next session: three-way E / I / unclassified assignment

**Status:** planned, not started · **Branch:** `analysis/glmcc-corrected`
**Context:** `reports/ei_differential/` (all prior E/I work), §5.5–5.8 of `METHODS_AND_RESULTS.md`

## Why the previous attempt was set up wrong

Every classification run so far forced a binary split: each unit had to be E or I. Real practice
does not do that. A waveform classifier assigns the confident tails and **leaves an explicit
unclassified middle** — commonly 30–60% of units — because the features genuinely do not separate
cells near the boundary.

Forcing a binary split on a continuum guarantees a result near the cluster-size prior. That is
visible in the region analysis: k-means with k = 4–5, taking the single fastest cluster, returned
23–46% "narrow" in **every** region regardless of neurochemistry (DMH 31.1%, VM-thalamus 33.3%),
which is roughly 1/k by construction.

So the binary result does not establish that the confident tails are meaningless. It only
establishes that the middle is a continuum. **Those are different claims and the previous work
conflated them.** This is worth re-testing.

## What to do

Work from `data/processed/waveform_features_harmonised.csv` — the filter-compensated features
(see §5.7). Do **not** use the raw features or the shipped `put_inh` labels.

1. **Fit a 2-component Gaussian mixture** per region (or pooled with region as a covariate,
   whichever the sample sizes allow) on (FW, PP), rather than k-means. GMM gives a posterior
   probability per unit, which is what an unclassified band needs.
2. **Assign three ways** on the posterior: `I` if p(narrow) ≥ 0.8, `E` if p(narrow) ≤ 0.2,
   `unclassified` otherwise. Report the coverage — what fraction each class takes. A coverage of
   40–70% classified is normal and is not a failure.
3. **Sensitivity:** repeat at 0.9/0.1 and 0.7/0.3. The regional pattern should be stable across
   thresholds if it is real.

## Success criterion — fix this before looking at the answer

The dataset contains its own controls. Define the test in advance:

| region | n | prior | required |
|---|---|---|---|
| **VM-thalamus** | 237 | rodent relay nuclei, local interneurons near-absent | confident-I **< 10%** |
| **DMH** | 122 | substantially GABAergic | confident-I materially **> VM-thalamus** |
| Mammillary Complex | 58 | glutamatergic projection | low |
| Ventromedial HT | 28 | predominantly SF1⁺ glutamatergic | low |
| ZI | 10 | almost entirely GABAergic | too small — cannot be used |

**The method passes only if VM-thalamus comes out low AND separates from DMH.** If confident-I is
~30% in VM-thalamus again, or if the regional ordering is flat, the features carry no E/I
information and the question is closed for this dataset — that would be the second independent
falsification and should be treated as conclusive.

Test the VM-thalamus vs DMH difference formally (proportion test or a mixed model with animal as
a random effect), not by eye.

## If it passes

4. Re-run the differential analyses (`scripts/ei_differential.py`) on **confident units only**,
   dropping the unclassified middle. Firing rate by class is the first check: a confident-narrow
   group that still shows no rate separation is still suspect.
5. Re-run the edge-type composition (E→E / E→I / I→E / I→I) on the matched arms in
   `results/glmcc_rm/`, paired per animal — the pooled version is Simpson's paradox (§5.5).
6. Report coverage alongside every number, so a reader can see how much of the data the claim
   rests on.

## Carry these constraints forward

- **Use the harmonised features.** Six recordings effectively lack the nominal 600 Hz high-pass;
  raw features contain a confirmed block artifact (§5.7).
- **Harmonisation is imperfect.** Compensating a mean, already-filtered waveform cannot undo phase
  distortion or averaging blur, and no pre-sorting data survives to do better. Whatever emerges
  carries that caveat.
- **No firing-rate separation exists** in the current labelling (I/E = 1.01, p = 0.33). If the
  confident-narrow group also shows none, that is strong evidence against, independent of region.
- **ARH (16), Anterior HT (12) and ZI (10) are too small** to cluster. ZI is the loss that hurts —
  almost entirely GABAergic, it would have been the positive control to pair with VM-thalamus.
- **Region labels and stereotaxic coordinates are sound** and independent of all of this.

## Estimated effort

Small — the features are computed and cached, the region merge is done, and
`scripts/ei_classify_waveforms.py` already has the clustering scaffolding. Expect under an hour,
most of it in the write-up and the sensitivity sweep.
