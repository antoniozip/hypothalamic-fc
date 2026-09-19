# Differential analysis by putative excitatory / inhibitory class

**Date:** 2026-09-19 · **Branch:** `analysis/glmcc-corrected`
**Artifacts:** `data/processed/neurons_ei_strict.csv`, `results/ei/`,
`figures/main/figE1_ei_differential.png`

## Summary

Putative E/I class, from an independent waveform criterion, **explains nothing** in this
connectivity. Node metrics, coupling sign and edge-type composition are all null. Two
methodological problems surfaced on the way, and both matter more than the null itself.

## Classification

`{animal}_put_inh.txt` comes from `tim_stuff.m:714`:

```matlab
putative_inh = intersect(find(half_w < 0.0006), find(trough_to_peak < 0.00065));
```

Narrow-spiking: half-peak width < 0.6 ms **and** trough-to-peak < 0.65 ms. MATLAB `find()` is
1-indexed, matching `neuron_id` in `neurons.csv` and in the metrics files. The criterion is
derived from spike waveforms and is **independent of GLMCC**, so using it to interrogate the
connectivity is not circular. That independence is what makes this analysis worth running.

### Only 7 of 13 animals were ever classified

| classified | % inhibitory | | not classified | % inhibitory |
|---|---|---|---|---|
| day3 | 51.1% | | night1 | 3.2% |
| night6 | 48.1% | | night4 | 2.9% |
| day6 | 46.4% | | day1 | 1.6% |
| day5 | 43.7% | | night2 | 0.0% |
| day4 | 42.6% | | night3 | 0.0% |
| night7 | 41.8% | | night5 | 0.0% |
| day2 | 38.7% | | | |

The gap runs 3.2% → 38.7% with nothing between. That is not a biological range; it is the
difference between the classification having been run and not. Analysable set: **399 units in 7
animals, 178 inhibitory (44.6%) and 221 excitatory.**

> **`src/py/ei_classification.py` is wrong and should not be used.** It assigns `"E"` to every
> unit not listed as inhibitory, so the six unclassified animals come out ~100% excitatory. Those
> animals hold 283 units, of which 4 appear in a stray `put_inh` entry, so **279 units are
> silently mislabelled excitatory** — 41% of the 682-unit dataset. `data/processed/neurons_ei.csv`
> carries those labels. This analysis uses `neurons_ei_strict.csv`, which marks them
> `unclassified`.

Matrices are the rate- and geometry-matched arms in `results/glmcc_rm/`, since the raw contrast
is a detection-power artifact (`reports/rate_matched/`).

## 1. The expected confound is absent — and that is itself a warning

Fast-spiking interneurons normally fire several times faster than pyramidal cells, and GLMCC's
threshold falls with spike count, so an E/I difference in connectivity could be firing rate in
disguise. It is not, because there is no rate difference to exploit:

| | median spikes per unit |
|---|---|
| excitatory | 4,095 |
| inhibitory | 4,500 |

Ratio I/E = **1.01×**, paired Wilcoxon across 14 animal × condition p = **0.33**.

Good for the analysis — no rate adjustment is needed. **Bad for the classification.** A
narrow-spiking population that does not fire faster than the broad-spiking one is not behaving
like a fast-spiking interneuron population. Together with the implausibly high 44.6% inhibitory
fraction, this suggests the waveform criterion is not cleanly separating cell types here. Every
null below inherits that doubt: a classification that does not track cell type cannot be
expected to predict anything.

## 2. Coupling sign cannot test the classification

If GLMCC's signs were meaningful, putative-inhibitory units should emit negative couplings more
often than putative-excitatory ones. They do not — but the test is **uninformative rather than
negative**, because the sign distribution is degenerate:

- **94.9% of classified units emit only negative couplings.**
- mean fraction negative: E **0.988** (n = 131), I **0.984** (n = 112), lightON, MWU p = 0.93
- ongoing: E 0.965 (n = 139), I 0.971 (n = 111), p = 0.95

With ~98% of all couplings negative there is no variance for the class to explain. This is the
same near-total inhibition already recorded in `METHODS_AND_RESULTS.md` §4.1 (97.4% of edges
negative), now shown to be independent of presynaptic cell class.

## 3. No node-metric differences

Paired per animal × condition (n = 14), median over units within each class:

| metric | E | I | Wilcoxon p |
|---|---|---|---|
| node_strength | −1.1646 | −1.1853 | 0.90 |
| clustering_coefficient | 0.3250 | 0.0208 | 0.37 |
| local_efficiency | 0.2591 | 0.0277 | 0.95 |
| hub_score | 0.0000 | 0.0343 | 0.53 |

Nothing approaches significance before correction. Note the large apparent gaps in clustering
and local efficiency with p ≈ 0.4–1.0: the between-animal spread swamps the class difference.

## 4. The edge-type gradient is a pooling artifact

Pooled over the seven animals, the composition looks orderly and interpretable:

| | lightON | ongoing |
|---|---|---|
| E→E | 9.378% | 9.720% |
| E→I | 7.841% | 8.479% |
| I→E | 7.858% | 8.294% |
| I→I | 7.294% | 7.662% |

E→E densest, I→I sparsest, same ordering in both conditions. **It does not survive the paired
test.** Taking the animal as the unit of analysis (n = 14 animal × condition):

| contrast | median A | median B | difference | p raw | p Holm | direction |
|---|---|---|---|---|---|---|
| E→E vs I→I | 4.455% | 5.000% | **−0.545** | 0.43 | 1.00 | 6/14 |
| E→E vs E→I | 4.455% | 4.167% | +0.288 | 0.54 | 1.00 | 5/14 |
| E→E vs I→E | 4.455% | 4.391% | +0.064 | 0.50 | 1.00 | 6/14 |
| E→I vs I→E | 4.167% | 4.391% | −0.224 | 0.81 | 1.00 | 5/14 |

The pooled E→E vs I→I difference even **reverses sign** in the per-animal medians. Pooling
counts across animals of very different size (24 to 126 units, 6 to 1,243 edges) lets the
largest animals dictate the aggregate — Simpson's paradox. The paired analysis is the correct
one and it is null: the direction holds in only 5–6 of 14.

Panel D of the figure shows both curves on the same axes.

## Conclusion

Across all four analyses there is **no differential effect of putative E/I class** on this
connectivity. The two findings worth carrying forward are methodological:

1. **`ei_classification.py` mislabels 279 units** in six animals as excitatory when they are
   simply unclassified — 41% of the dataset. Anything downstream of `neurons_ei.csv` inherits it.
2. **The waveform classification does not behave like a cell-type classification** — 44.6%
   inhibitory and no firing-rate separation. Before reading the nulls as biology, the criterion
   should be re-examined on the waveforms themselves.

A null from a classification that may not classify is weak evidence about circuitry. It is
strong evidence that this E/I labelling, as it stands, is not usable.

## Limitations

- **7 of 13 animals**, and the six excluded include day1 (126 units, the largest).
- **The classification may not be valid** (§1). This is the dominant limitation; it caps what any
  of these results can mean.
- Analysis runs on the rate-matched arms, which are thinned and therefore lower-powered; night3
  and similar small animals contribute very few edges.
- Sign validation is limited by the degenerate sign distribution, not by sample size.
- No spatial control: E and I units are not matched on region, and region composition differs
  between classes. A region-stratified version was not run.
- Edge-type analysis treats all classified pairs as independent within an animal; the paired test
  over animals is the conservative reading and is the one reported.

## Next actions

1. **Re-derive the waveform classification** from the raw waveforms, checking half-width and
   trough-to-peak distributions for genuine bimodality and verifying the firing-rate separation
   that should accompany it. Until that is done these nulls cannot be attributed to biology.
2. **Fix or delete `src/py/ei_classification.py`** and regenerate `neurons_ei.csv`, or point
   consumers at `neurons_ei_strict.csv`.
3. **Re-run the six unclassified animals** through the MATLAB criterion if the waveforms are
   still available — that would roughly double the analysable set and add day1.
4. If the classification is rehabilitated, repeat this analysis region-stratified.
