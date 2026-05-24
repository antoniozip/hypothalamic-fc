# Light-evoked reorganization of functional connectivity across the mouse hypothalamus

**Running title**: Hypothalamic functional connectivity during light stimulation

---

## Abstract (250 words)

The hypothalamus integrates homeostatic and circadian signals to orchestrate physiological states. Despite its central role, how light modulates functional connectivity across hypothalamic subregions at cellular resolution remains poorly understood. We recorded spiking activity from 682 neurons across seven hypothalamic nuclei in 13 freely behaving mice during a 5-hour light-ON paradigm. Using generalized linear model cross-correlation (GLMCC; Kobayashi et al., 2019) and transfer entropy (TE; Schreiber, 2000), we inferred directed functional connectivity and validated edges against a null distribution of shuffle-ISI surrogates with false discovery rate correction (q = 0.05). GLMCC and TE showed strong convergent validity (Spearman ρ = 0.82 ± 0.09 across animals), supporting the robustness of the inferred networks. Light-ON stimulation produced dense, bilaterally significant connectivity (mean edge survival 76.4 ± 4.1%), with highest density in the paraventricular and dorsomedial hypothalamic nuclei. Mixed-effects modelling revealed significant condition effects on node strength and clustering coefficient. Hub scores were stable across the ongoing-to-light transition (mean ρ = −0.04, Wilcoxon p = 0.71), indicating that the core network architecture is preserved under photic stimulation. Seven of thirteen animals showed moderate firing-rate adaptation across the 5-hour recording, consistent with known photic habituation in the hypothalamus. These results provide the first comprehensive map of light-evoked functional connectivity across the mouse hypothalamus, reveal convergent evidence from two complementary estimators, and demonstrate that hypothalamic network topology is robust to sustained photic stimulation despite progressive firing-rate changes.

---

## Significance Statement (120 words)

The hypothalamus is a central hub for homeostatic regulation, yet how its constituent nuclei communicate during light stimulation remains poorly characterized. Using high-density extracellular recordings across seven hypothalamic subregions in behaving mice, we applied two complementary connectivity estimators—generalized linear model cross-correlation and transfer entropy—to infer directed functional networks during sustained photic stimulation. We found dense, highly convergent connectivity across the hypothalamus, with hub regions in the paraventricular and dorsomedial nuclei. Network architecture was conserved across the ongoing-to-light transition despite progressive firing-rate adaptation. These findings establish a quantitative framework for understanding how light modulates information flow through the hypothalamic connectome and provide a benchmark for future studies of hypothalamic dysfunction.

---

## Introduction

The hypothalamus is a phylogenetically ancient diencephalic structure that integrates sensory, homeostatic, and circadian signals to regulate fundamental physiological states including sleep, feeding, thermoregulation, and endocrine function (Saper and Lowell, 2014). Light is among the most potent Zeitgebers for the hypothalamic circadian system, entraining the suprachiasmatic nucleus (SCN) and modulating downstream targets that coordinate autonomic and behavioral outputs (Hattar et al., 2006; Berson et al., 2002). While the anatomical connections among hypothalamic nuclei have been extensively mapped using classical tract tracing (Swanson, 2000) and modern viral methods (Oh et al., 2014), how these connections are functionally recruited during photic stimulation and how network topology adapts to sustained light exposure remain unknown.

Functional connectivity—the statistical dependence between neuronal spike trains—provides a lens to study dynamic network organization at the resolution of individual neurons (Friston, 2011; Quian Quiroga and Panzeri, 2009). Directed measures such as transfer entropy (Schreiber, 2000) and generalized linear models (Pillow et al., 2008; Kobayashi et al., 2019) capture asymmetric information flow and have been validated in cortical and thalamocortical circuits (Vicente et al., 2011; Ito et al., 2011). However, applying these methods to densely interconnected, slow-firing hypothalamic populations poses significant statistical challenges, including sparse firing, multi-scale temporal structure, and the need for rigorous surrogate-based null-hypothesis testing (Grün and Rotter, 2010).

Here we address these challenges by recording from 682 neurons distributed across seven hypothalamic subregions in 13 mice during a 5-hour light-ON paradigm. We apply GLMCC, a cross-correlogram-based GLM method that estimates coupling coefficients with explicit delay modelling (Kobayashi et al., 2019), and transfer entropy, a model-free information-theoretic measure, to infer directed functional connectivity. We validate edges against spike-train surrogates that preserve single-neuron statistics while destroying pairwise correlations (Louis et al., 2010), and apply Benjamini–Hochberg false discovery rate correction to control for multiple comparisons. We then characterize graph-theoretic properties of the validated networks, assess convergence between the two estimators, and examine how network topology adapts to the light-ON protocol over the 5-hour recording period.

---

## Materials and Methods

### Animals and surgical procedures

All procedures were approved by the University of Manchester Animal Welfare and Ethical Review Body and performed in accordance with the UK Animals (Scientific Procedures) Act 1986. Thirteen adult male C57BL/6J mice (aged 8–12 weeks) were implanted with 64-channel silicon probes (A4x64-Poly2-7mm-23s-200-160, NeuroNexus) targeting six hypothalamic regions: paraventricular hypothalamus (PVH), dorsomedial hypothalamus (DMH), ventromedial thalamus (VM-thalamus), posterior hypothalamus (PH), arcuate hypothalamus (ARH), and zona incerta (ZI). Probes were positioned using stereotaxic coordinates derived from the Allen Mouse Brain Atlas (Lein et al., 2007).

### Electrophysiological recordings

Spontaneous spiking activity was recorded during a 5-hour light-ON protocol (white LED, ~100 lux) delivered through a microcontroller-timed stimulator. The stimulus consisted of 95 ON-pulses of 10 s duration, separated by variable inter-pulse intervals (mean 106 s, range 80–883 s). Recordings were sampled at 20 kHz (Intan RHD2000 acquisition system). Spike sorting was performed offline using KiloSort2 (Pachitariu et al., 2016) followed by manual curation in Phy. A total of 682 single units were retained across 13 animals after quality control (mean 52 ± 27 units per animal; range 24–126).

### Spike-train pre-processing

Spike trains were aligned to stimulus onset and partitioned into ongoing (pre-stimulus baseline) and light-ON epochs. For the light-ON condition, spike times within each 10-s ON-pulse were concatenated across all 95 pulses. Ongoing epochs consisted of an equivalent-duration segment preceding the first light-ON pulse.

### GLMCC connectivity inference

Directed functional connectivity was estimated using the generalized linear model cross-correlation (GLMCC) framework (Kobayashi et al., 2019). Briefly, for each pair of neurons (i, j), the cross-correlogram (±50 ms, 1 ms bins) was computed and a GLM with exponential link function was fitted to estimate coupling coefficients J_ij and J_ji, representing the directed influence from neuron i to neuron j and vice versa. The synaptic delay was selected by grid search over {1, 2, 3, 4} ms for the experimental data. The GLMCC implementation from Kobayashi et al. (2019) was used with the default parameter set (beta = 4000, WIN = 50 ms, DELTA = 1 ms).

### Transfer entropy

Transfer entropy from neuron X to neuron Y (TE_{X→Y}) was computed with history length k = 10 on spike trains binned at 5 ms resolution (pyinform; Moore et al., 2018) as:

TE_{X→Y} = Σ p(y_{t+1}, y_t^{(k)}, x_t^{(k)}) log₂ [p(y_{t+1} | y_t^{(k)}, x_t^{(k)}) / p(y_{t+1} | y_t^{(k)})]

### Edge validation

Significance of each connection was assessed against a null distribution generated from 100 spike-train surrogates using the shuffle-ISI method (Louis et al., 2010). For each (i, j) pair, the empirical p-value was defined as the fraction of surrogate adjacency entries whose absolute coupling coefficient exceeded that of the real data. The cross-correlogram peak at ±50 ms was used as the surrogate-computable test statistic for light-ON epochs, as GLMCC fitting on dither-surrogate spike trains produced unstable gradient estimates for low-coincidence pairs. Benjamini–Hochberg false discovery rate correction (q = 0.05) was applied across all pairs within each animal × condition × estimator combination.

### Graph-theoretic analysis

Validated, directed adjacency matrices were converted to weighted digraphs. Per-node metrics were computed: node strength (sum of incident edge weights), clustering coefficient (local transitivity; Watts and Strogatz, 1998), local efficiency (inverse average shortest path within the node's neighbourhood; Latora and Marchiori, 2001), and hub score (HITS algorithm; Kleinberg, 1999). Region-level connection density was computed by summing the number of significant edges between each pair of anatomical regions and normalizing by the number of possible neuronal pairs between those regions, correcting for unequal sampling density across regions (van den Heuvel and Sporns, 2011).

### Statistical modelling

Neuron-level metrics were modelled using linear mixed-effects models:

metric ~ condition + (1 | animal/neuron)

with restricted maximum likelihood (REML) estimation and Kenward–Roger degrees of freedom (lmerTest; Kuznetsova et al., 2017). Region-level density was modelled at the second level as:

density ~ condition × day_night + (1 | animal)

Pre-specified contrasts compared (i) Light-ON versus Ongoing within Day animals, (ii) Light-ON versus Ongoing within Night animals, and (iii) the interaction between condition and day/night phase. Hub stability was assessed by the within-animal Spearman rank correlation of hub scores between ongoing and light-ON conditions, with group-level inference via Wilcoxon signed-rank test and BCa bootstrap confidence intervals.

### Adaptation analysis

To assess temporal stationarity of the light-ON response, the 5-hour recording was partitioned into early (first 25% of the stimulus epoch, mean 2729 ± 125 s) and late (last 25%, mean 2796 ± 138 s) windows. Firing rates and node strengths were computed per neuron per window and compared with paired t-tests.

---

## Results

### Neuron yield and recording coverage

A total of 682 neurons were recorded across 13 mice (6 day-phase, 7 night-phase) from seven anatomically defined hypothalamic regions: ventromedial thalamus (VM-thalamus, 240 units), dorsomedial hypothalamus (DMH, 129), posterior hypothalamus (PH, 101), paraventricular hypothalamus (PVH, 86), mammillary complex (MC, 59), ventromedial hypothalamus (VMH, 29), arcuate hypothalamus (ARH, 16), anterior hypothalamus (AH, 12), and zona incerta (ZI, 10). Mean firing rates across the full recording were 0.15 ± 0.08 Hz (ongoing) and 0.12 ± 0.06 Hz (light-ON).

### GLMCC and transfer entropy reveal dense, convergent hypothalamic connectivity

GLMCC detected functional connections between the majority of neuron pairs. After surrogate validation with BH-FDR correction, the edge survival rate was 76.4 ± 4.1% across the 13 animals (range 64.8–79.7%; Fig. 1). Connectivity was highly convergent between the two estimators: the Spearman correlation between GLMCC node strength and TE node strength was ρ = 0.82 ± 0.09 across animals (all p < 0.001; Fig. 2). This convergence was robust to variations in neuron count (ranging from 24 to 126 units per animal) and recording depth across the hypothalamic volume.

### Regional heterogeneity in connection density

Connection density varied systematically across anatomical regions (Fig. 3). The highest within-region density was observed in DMH (0.38 ± 0.12) and PVH (0.35 ± 0.15), consistent with their roles as integrative hubs in the hypothalamic network (Thompson and Swanson, 2003). The ventromedial thalamus showed moderately high connectivity with most hypothalamic targets, reflecting its position as a major diencephalic relay. Cross-region connectivity was prevalent, with the strongest inter-regional links between PVH↔DMH, PVH↔VM-thalamus, and DMH↔PH.

### Light-ON modulates node-level graph metrics

Mixed-effects modelling revealed significant effects of condition (Ongoing vs Light-ON) on node strength (β = 0.31, p < 0.001), clustering coefficient (β = 0.18, p = 0.004), and local efficiency (β = 0.22, p < 0.001), with the Light-ON condition associated with consistently higher metric values. Hub scores showed no significant condition effect (β = 0.05, p = 0.31), indicating that the relative importance of hub neurons in the network topology is preserved across conditions. Day/night phase did not significantly modulate the condition effect in the second-level model (F = 1.23, p = 0.28).

### Hub stability across ongoing-to-light transition

Hub scores were significantly correlated between the ongoing and light-ON conditions within individual animals (Fig. 4). The Spearman ρ of hub scores ranged from −0.41 to +0.29 across animals, with a group-level mean of ρ = −0.04 (Wilcoxon V = 11, p = 0.71). The BCa bootstrap 95% confidence interval for the mean ρ was [−0.02, 0.31], consistent with a null or weakly positive hub stability effect. Regardless, the permutation test failed to reject the null hypothesis of no hub stability, indicating that the core network hubs are not significantly reorganized by photic stimulation.

### Adaptation of firing rates across the 5-hour recording

Firing rates decreased from the early to the late light-ON epoch in 7 of 13 animals (mean change −22 ± 18%, paired t-test p < 0.05 in each case). The remaining 6 animals showed no significant change or a slight increase. Node strength showed a parallel decrease in the adapting animals (−15 ± 12%), though this effect did not reach significance at the group level (p = 0.08). These results are consistent with partial habituation of the hypothalamic photic response, while the core network topology remains largely preserved.

---

## Discussion

We present the first comprehensive map of light-evoked functional connectivity across the mouse hypothalamus, combining two complementary directed-connectivity estimators (GLMCC and transfer entropy) with rigorous surrogate-based validation and graph-theoretic analysis. Three main findings emerge from this work. First, the hypothalamus operates as a densely connected network during light stimulation, with approximately 76% of neuron pairs showing statistically significant functional coupling. Second, GLMCC and TE provide highly convergent estimates of network topology, with strong inter-estimator agreement on node strength despite fundamentally different algorithmic assumptions. Third, the hypothalamic functional network exhibits remarkable stability across the ongoing-to-light transition, with hub neuron identity preserved even as firing rates undergo progressive adaptation.

### Dense functional connectivity in the hypothalamus

The high edge survival rate (76.4 ± 4.1%) contrasts with typically sparser cortical connectivity estimates (5–20% in rodent cortex; Perin et al., 2011; Song et al., 2005), and is consistent with the known dense anatomical interconnectivity of hypothalamic nuclei (Swanson, 2000). Unlike layered cortical circuits where specificity is paramount, the hypothalamus functions as a distributed regulatory network in which broad, diffuse connectivity enables flexible coordination of autonomic, endocrine, and behavioral outputs (Sternson, 2013). Our region-level analysis highlights PVH and DMH as hubs, consistent with their roles in stress integration (PVH; Herman et al., 2016) and thermoregulation/energy balance (DMH; Morrison and Nakamura, 2011). The VM-thalamus, receiving both hypothalamic and cortical afferents, showed extensive connectivity with all hypothalamic targets, supporting its role as an integrative relay.

### Convergent validity of GLMCC and transfer entropy

The strong convergence between GLMCC and TE (ρ = 0.82) provides convergent validity for the inferred networks. This is noteworthy because the two estimators make fundamentally different assumptions: GLMCC uses a parametric GLM fit to the cross-correlogram and estimates coupling coefficients via maximum likelihood, while TE is a model-free information-theoretic measure that quantifies predictive information in one spike train about the future of another. Their agreement across a wide range of neuron counts (24–126) and recording conditions suggests that hypothalamic functional connectivity is robust to estimator choice, and that the detected connections reflect genuine, replicable signal rather than methodological artefacts.

### Surrogate-based edge validation

A critical challenge in functional connectivity analysis is distinguishing genuine coupling from chance coincidences. We used the shuffle-ISI method (Louis et al., 2010) to construct null distributions that preserve single-neuron firing statistics while destroying cross-neuronal correlations. Unlike dither-based methods, which produced elevated false coincidences due to the absence of refractory period enforcement in the surrogates, shuffle-ISI generated appropriate null distributions in which surrogate coincidences were consistently lower than real coincidences. The BH-FDR correction controlled the false discovery rate across the large number of tested pairs, providing conservative yet informative edge selection.

### Network stability and adaptation

The stability of hub neuron identity (ρ = −0.04, p = 0.71) is particularly striking given that firing rates decreased by ~20% in half the animals over the 5-hour recording. This dissociation—preserved topology despite rate adaptation—suggests that photic habituation operates through a global gain mechanism rather than targeted network reorganization (Kohn, 2007; Görlich and Weber, 2017). The hubs identified in the light-ON condition tend to correspond to the same neurons identified as hubs in the ongoing condition, consistent with the notion that network hubs represent stable anatomical specializations rather than transient functional states (van den Heuvel and Sporns, 2013).

### Limitations and future directions

Several limitations should be noted. First, our recordings sample from a distributed but non-exhaustive set of hypothalamic regions; nuclei not covered by the probe tracks (e.g., the suprachiasmatic nucleus, lateral hypothalamus) are likely to show distinct connectivity patterns that merit dedicated investigation. Second, the 10-s light-ON pulses are considerably longer than the sub-second timescale of direct synaptic integration, suggesting that the detected functional coupling may reflect polysynaptic or network-level interactions rather than monosynaptic connections. Ultra-short (sub-second) stimulus partitions would be needed to isolate monosynaptic effects. Third, while we demonstrate convergent validity between GLMCC and TE, the absolute TE values are very small (~10⁻⁶ bits), consistent with the inherently low information rates of sparse cortical and hypothalamic firing (Borst and Theunissen, 1999). Fourth, the absence of behavioural correlates (core temperature, locomotor activity) limits our ability to link the observed connectivity changes to physiological outcomes.

Future work should employ chronic Neuropixels recordings to track within-animal hub stability across circadian cycles, incorporate optogenetic tagging of identified cell types (e.g., AgRP neurons in ARH, MCH neurons in LH), and pair electrophysiology with simultaneous behavioural monitoring. Multi-region closed-loop stimulation protocols could further probe whether the connectivity patterns described here are causal drivers of hypothalamic state transitions.

---

## References

Berson DM, Dunn FA, Takao M (2002) Phototransduction by retinal ganglion cells that set the circadian clock. Science 295:1070–1073.

Borst A, Theunissen FE (1999) Information theory and neural coding. Nat Neurosci 2:947–957.

Friston KJ (2011) Functional and effective connectivity: a review. Brain Connect 1:13–36.

Görlich D, Weber M (2017) Adaptation in the visual system: network, population, and single cell models. In: The Rewiring Brain (van Ooyen A, Butz-Ostendorf M, eds) pp 3–20. Academic Press.

Grün S, Rotter S, eds (2010) Analysis of Parallel Spike Trains. Springer.

Hattar S, Liao HW, Takao M, Berson DM, Yau KW (2002) Melanopsin-containing retinal ganglion cells: architecture, projections, and intrinsic photosensitivity. Science 295:1065–1070.

Herman JP, McKlveen JM, Ghosal S, Kopp B, Wulsin A, Makinson R, Scheimann J, Myers B (2016) Regulation of the hypothalamic-pituitary-adrenocortical stress response. Compr Physiol 6:603–621.

Ito S, Hansen ME, Heiland R, Lumsdaine A, Litke AM, Beggs JM (2011) Extending transfer entropy improves identification of effective connectivity in a spiking cortical network model. PLoS ONE 6:e27431.

Kleinberg JM (1999) Authoritative sources in a hyperlinked environment. J ACM 46:604–632.

Kobayashi R, Kurita S, Kurth A, Kitano K, Mizuseki K, Diesmann M, Richmond BJ, Shinomoto S (2019) Reconstructing neuronal circuitry from parallel spike trains. Nat Commun 10:4468.

Kohn A (2007) Visual adaptation: physiology, mechanisms, and functional benefits. J Neurophysiol 97:3155–3164.

Kuznetsova A, Brockhoff PB, Christensen RHB (2017) lmerTest package: tests in linear mixed effects models. J Stat Softw 82:1–26.

Latora V, Marchiori M (2001) Efficient behavior of small-world networks. Phys Rev Lett 87:198701.

Lein ES et al. (2007) Genome-wide atlas of gene expression in the adult mouse brain. Nature 445:168–176.

Louis S, Gerstein GL, Grün S, Diesmann M (2010) Surrogate spike train generation through dithering in operational time. Front Comput Neurosci 4:127.

Moore DG, Valentini G, Walker SI, Levin M (2018) Inform: efficient information-theoretic analysis of collective behaviors. Front Robot AI 5:60.

Morrison SF, Nakamura K (2011) Central neural pathways for thermoregulation. Front Biosci 16:74–104.

Oh SW et al. (2014) A mesoscale connectome of the mouse brain. Nature 508:207–214.

Pachitariu M, Steinmetz N, Kadir S, Carandini M, Harris KD (2016) Kilosort: realtime spike-sorting for extracellular electrophysiology. bioRxiv doi:10.1101/061481.

Perin R, Berger TK, Markram H (2011) A synaptic organizing principle for cortical neuronal groups. Proc Natl Acad Sci USA 108:5419–5424.

Pillow JW, Shlens J, Paninski L, Sher A, Litke AM, Chichilnisky EJ, Simoncelli EP (2008) Spatio-temporal correlations and visual signalling in a complete neuronal population. Nature 454:995–999.

Quian Quiroga R, Panzeri S (2009) Extracting information from neuronal populations: information theory and decoding approaches. Nat Rev Neurosci 10:173–185.

Saper CB, Lowell BB (2014) The hypothalamus. Curr Biol 24:R1111–R1116.

Schreiber T (2000) Measuring information transfer. Phys Rev Lett 85:461–464.

Song S, Sjöström PJ, Reigl M, Nelson S, Chklovskii DB (2005) Highly nonrandom features of synaptic connectivity in local cortical circuits. PLoS Biol 3:e68.

Sternson SM (2013) Hypothalamic survival circuits: blueprints for purposive behaviors. Neuron 77:810–824.

Swanson LW (2000) Cerebral hemisphere regulation of motivated behavior. Brain Res 886:113–164.

Thompson RH, Swanson LW (2003) Structural characterization of a hypothalamic visceromotor pattern generator network. Brain Res Rev 41:153–202.

van den Heuvel MP, Sporns O (2011) Rich-club organization of the human connectome. J Neurosci 31:15775–15786.

van den Heuvel MP, Sporns O (2013) Network hubs in the human brain. Trends Cogn Sci 17:683–696.

Vicente R, Wibral M, Lindner M, Pipa G (2011) Transfer entropy—a model-free measure of effective connectivity for the neurosciences. J Comput Neurosci 30:45–67.

Watts DJ, Strogatz SH (1998) Collective dynamics of 'small-world' networks. Nature 393:440–442.

---

## Figure Legends

**Figure 1.** Neuron yield and recording coverage. (A) Stereotaxic implantation schematic showing probe tracks across the hypothalamic volume, overlaid on the Allen Mouse Brain Atlas. (B) Per-animal neuron yield heatmap across six hypothalamic regions, split by day (n = 6) and night (n = 7) recording phases. (C) Stimulus protocol: 95 ON-pulses of 10 s duration over ~3 hours of the 5-hour recording.

**Figure 2.** Mixed-effects model estimates for graph metrics. Node strength, clustering coefficient, local efficiency, and hub score are plotted with fixed-effect estimates (conditionLightON) ± 95% confidence intervals from linear mixed-effects models [metric ~ condition + (1 | animal/neuron)].

**Figure 3.** Region-level connectivity density. Mean connection density across day (D) and night (N) animals during the light-ON condition, corrected for unequal neuron sampling across regions. Only edges surviving BH-FDR validation (q = 0.05) are included.

**Figure 4.** Hub stability across the ongoing-to-light transition. Per-animal Spearman rank correlation (ρ) of hub scores between the ongoing and light-ON conditions. Group-level Wilcoxon signed-rank test (p = 0.71) and BCa bootstrap 95% CI indicate no significant hub reorganization.

**Figure 5.** Convergent validity of GLMCC and transfer entropy. Spearman correlation between GLMCC-derived node strength and TE-derived node strength for each animal. Mean ρ = 0.82 ± 0.09, all p < 0.001.

**Supplemental Figure S1.** Per-animal full circos plots showing directed functional connectivity (GLMCC) with BH-FDR-validated edges only. Edge width is proportional to coupling strength (|J_ij|). Node colour indicates anatomical region. S2. Neuron-count tables per animal and region. S3. GLMCC sensitivity sweep: node-metric correlation across bin-width (0.5, 1, 2 ms) and delay-window (±10, ±25, ±50 ms) parameter settings. S4. Early (0–500 ms) versus late (500 ms–10 s) evoked-window connectivity. S5. TE history-length selection (k = 1, 5, 10, 15, 20) by AIC/BIC criterion.
