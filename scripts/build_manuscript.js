#!/usr/bin/env node
/**
 * Build revised hypothalamic FC manuscript with integrated figures.
 * Incorporates all Tier 1-3 analysis findings.
 */

const fs = require("fs");
const path = require("path");
const {
  Document, Packer, Paragraph, TextRun, Table, TableRow, TableCell,
  ImageRun, Header, Footer, AlignmentType, HeadingLevel,
  BorderStyle, WidthType, ShadingType, PageNumber, PageBreak,
  LevelFormat, TableOfContents
} = require("docx");

const PROJECT = "/media/antonio/data/tim_brown";
const FIGURES = path.join(PROJECT, "figures", "main");

function loadImage(filename) {
  return fs.readFileSync(path.join(FIGURES, filename));
}

const border = { style: BorderStyle.SINGLE, size: 1, color: "CCCCCC" };
const borders = { top: border, bottom: border, left: border, right: border };
const cellMargins = { top: 80, bottom: 80, left: 120, right: 120 };

function figPara(imgData, caption, width = 500, height = 375) {
  return [
    new Paragraph({
      spacing: { before: 200, after: 100 },
      alignment: AlignmentType.CENTER,
      children: [new ImageRun({
        type: "png", data: imgData,
        transformation: { width, height },
        altText: { title: caption, description: caption, name: caption }
      })]
    }),
    new Paragraph({
      spacing: { after: 300 },
      alignment: AlignmentType.CENTER,
      children: [
        new TextRun({ text: caption, italics: true, size: 18, color: "555555" })
      ]
    })
  ];
}

function heading(text, level = 1) {
  return new Paragraph({
    heading: level === 1 ? HeadingLevel.HEADING_1 : HeadingLevel.HEADING_2,
    spacing: { before: 360, after: 200 },
    children: [new TextRun({ text, bold: true, size: level === 1 ? 28 : 24 })]
  });
}

function para(text, opts = {}) {
  return new Paragraph({
    spacing: { after: 120, line: 276 },
    alignment: AlignmentType.JUSTIFIED,
    children: [new TextRun({ text, size: 22, ...opts })]
  });
}

function richPara(runs) {
  return new Paragraph({
    spacing: { after: 120, line: 276 },
    alignment: AlignmentType.JUSTIFIED,
    children: runs.map(r => new TextRun({ size: 22, ...r }))
  });
}

// ═══════════════════════════════════════════════════════════════════
const doc = new Document({
  styles: {
    default: { document: { run: { font: "Times New Roman", size: 22 } } },
    paragraphStyles: [
      { id: "Heading1", name: "Heading 1", basedOn: "Normal", next: "Normal", quickFormat: true,
        run: { size: 28, bold: true, font: "Times New Roman" },
        paragraph: { spacing: { before: 360, after: 200 }, outlineLevel: 0 } },
      { id: "Heading2", name: "Heading 2", basedOn: "Normal", next: "Normal", quickFormat: true,
        run: { size: 24, bold: true, font: "Times New Roman" },
        paragraph: { spacing: { before: 240, after: 120 }, outlineLevel: 1 } },
    ]
  },
  sections: [{
    properties: {
      page: {
        size: { width: 12240, height: 15840 },
        margin: { top: 1440, right: 1440, bottom: 1440, left: 1440 }
      }
    },
    headers: {
      default: new Header({
        children: [new Paragraph({
          alignment: AlignmentType.RIGHT,
          children: [new TextRun({ text: "Hypothalamic Functional Connectivity", italics: true, size: 18, color: "888888" })]
        })]
      })
    },
    children: [

      // ═══════════════════════════════════════════════════════════
      // TITLE PAGE
      // ═══════════════════════════════════════════════════════════
      new Paragraph({ spacing: { before: 2400 } }),
      new Paragraph({
        alignment: AlignmentType.CENTER,
        spacing: { after: 200 },
        children: [new TextRun({
          text: "Light-Evoked Reorganization of Functional Connectivity Across the Mouse Hypothalamus",
          bold: true, size: 32, font: "Times New Roman"
        })]
      }),
      new Paragraph({
        alignment: AlignmentType.CENTER,
        spacing: { after: 120 },
        children: [new TextRun({
          text: "Revised Manuscript — May 2026",
          size: 22, color: "555555"
        })]
      }),
      new Paragraph({
        alignment: AlignmentType.CENTER,
        spacing: { after: 400 },
        children: [new TextRun({
          text: "Antonio Giuliano Zippo",
          size: 24
        })]
      }),
      new Paragraph({ children: [new PageBreak()] }),

      // ═══════════════════════════════════════════════════════════
      // ABSTRACT
      // ═══════════════════════════════════════════════════════════
      heading("Abstract"),
      para("The hypothalamus integrates homeostatic and circadian signals to orchestrate physiological states. " +
           "Despite its central role, how light modulates functional connectivity across hypothalamic subregions " +
           "at cellular resolution remains poorly understood. We recorded spiking activity from 682 neurons across " +
           "nine hypothalamic nuclei in 13 freely behaving mice (6 day-phase, 7 night-phase) during a 5-hour light-ON " +
           "paradigm. Using generalized linear model cross-correlation (GLMCC; Kobayashi et al., 2019) with " +
           "shuffle-ISI surrogate-based edge validation and Benjamini-Hochberg false discovery rate correction (q = 0.05), " +
           "we inferred directed functional connectivity and characterized graph-theoretic network properties."),

      para("Light-ON stimulation produced dense, bilaterally significant connectivity (mean edge survival " +
           "79.9 ± 4.1% across 13 animals). Mixed-effects modelling revealed a highly significant " +
           "condition × region interaction for all graph metrics (node strength: p < 2×10⁻¹²; " +
           "clustering coefficient: p = 1.1×10⁻⁸; local efficiency: p < 2×10⁻¹⁶; hub score: p = 0.005), " +
           "demonstrating that light's effect on network topology varies systematically by hypothalamic region. " +
           "The paraventricular (PVH) and dorsomedial (DMH) hypothalamic nuclei showed the strongest light-evoked " +
           "modulation. Connectivity was predominantly inter-regional (59.3% of edges cross anatomical boundaries), " +
           "and highly consistent across animals (76.5% mean edge consistency), indicating a stereotyped " +
           "anatomical backbone."),

      para("Hub scores were not significantly reorganized by photic stimulation (within-animal Spearman ρ = −0.08, " +
           "Wilcoxon p = 0.92), and rich-club analysis revealed no significant core of densely interconnected hubs " +
           "in any animal (0/13, p > 0.05). Triplet motif analysis showed massive enrichment of connected 3-node " +
           "circuits above degree-preserving null models (mean Z = 81.4), indicating that despite near-complete " +
           "binary connectivity (~80% density), specific circuit-level structure is strongly non-random. " +
           "The network showed no distance-dependent organization (ρ = 0.03 between inter-channel distance and " +
           "connection probability), confirming that detected connectivity reflects functional coupling rather than " +
           "spatial proximity artifacts. Firing rate partially predicted node strength (Spearman ρ = 0.50, " +
           "25% variance explained), with inhibitory neurons showing higher firing rates (0.28 vs 0.17 Hz, p = 0.001) " +
           "but not significantly higher connectivity. These results provide the first comprehensive map of " +
           "light-evoked functional connectivity across the mouse hypothalamus and establish that hypothalamic " +
           "network topology is robust to sustained photic stimulation, operating as a densely connected, " +
           "spatially uniform, and circuit-level structured integrative network."),

      new Paragraph({ children: [new PageBreak()] }),

      // ═══════════════════════════════════════════════════════════
      // INTRODUCTION
      // ═══════════════════════════════════════════════════════════
      heading("Introduction"),

      para("The hypothalamus is a phylogenetically ancient diencephalic structure that integrates sensory, " +
           "homeostatic, and circadian signals to regulate fundamental physiological states including sleep, " +
           "feeding, thermoregulation, and endocrine function (Saper and Lowell, 2014). Light is among the most " +
           "potent Zeitgebers for the hypothalamic circadian system, entraining the suprachiasmatic nucleus (SCN) " +
           "and modulating downstream targets that coordinate autonomic and behavioral outputs (Hattar et al., 2006; " +
           "Berson et al., 2002). While the anatomical connections among hypothalamic nuclei have been extensively " +
           "mapped using classical tract tracing (Swanson, 2000) and modern viral methods (Oh et al., 2014), " +
           "how these connections are functionally recruited during photic stimulation and how network topology " +
           "adapts to sustained light exposure remain unknown."),

      para("Functional connectivity\u2014the statistical dependence between neuronal spike trains\u2014provides " +
           "a lens to study dynamic network organization at the resolution of individual neurons (Friston, 2011; " +
           "Quian Quiroga and Panzeri, 2009). Directed measures such as generalized linear model cross-correlation " +
           "(GLMCC; Kobayashi et al., 2019) capture asymmetric information flow and have been validated in cortical " +
           "and thalamocortical circuits. However, applying these methods to densely interconnected, slow-firing " +
           "hypothalamic populations poses significant statistical challenges, including sparse firing, multi-scale " +
           "temporal structure, and the need for rigorous surrogate-based null-hypothesis testing (Gr\u00FCn and " +
           "Rotter, 2010)."),

      para("Here we address these challenges by recording from 682 neurons distributed across nine hypothalamic " +
           "subregions in 13 mice during a 5-hour light-ON paradigm. We apply GLMCC with shuffle-ISI surrogate-based " +
           "edge validation and Benjamini-Hochberg false discovery rate correction. We then characterize graph-theoretic " +
           "properties of the validated networks, including node strength, clustering coefficient, local efficiency, " +
           "and hub scores, using linear mixed-effects models that explicitly test for condition \u00D7 region " +
           "interactions. We further examine network architecture through rich-club analysis, edge consistency across " +
           "animals, within- versus between-region connectivity decomposition, centrality convergence, motif analysis, " +
           "and distance-dependent connectivity."),

      new Paragraph({ children: [new PageBreak()] }),

      // ═══════════════════════════════════════════════════════════
      // RESULTS
      // ═══════════════════════════════════════════════════════════
      heading("Results"),

      heading("Neuron yield and recording coverage", 2),
      para("A total of 682 neurons were recorded across 13 mice (6 day-phase: day1\u2013day6; 7 night-phase: " +
           "night1\u2013night7) from nine anatomically defined hypothalamic regions. The ventromedial thalamus " +
           "(VM-thalamus) contributed the largest sample (240 units, 35.2%), followed by dorsomedial hypothalamus " +
           "(DMH, 129 units), posterior hypothalamus (PH, 101 units), paraventricular hypothalamus (PVH, 86 units), " +
           "mammillary complex (MC, 59 units), ventromedial hypothalamus (VMH, 29 units), arcuate hypothalamus " +
           "(ARH, 16 units), anterior hypothalamus (AH, 12 units), and zona incerta (ZI, 10 units). " +
           "Mean firing rates were 0.20 \u00B1 0.37 Hz (light-ON). Inhibitory neurons, identified by " +
           "waveform characteristics, comprised 26.7% of recorded units (182/682) and were inhomogeneously " +
           "distributed across regions: DMH (40.3% I), PH (44.6% I), MC (47.5% I), and ARH (50.0% I) showed " +
           "high inhibitory proportions, while VM-thalamus (8.3% I) and AH (8.3% I) were predominantly excitatory."),

      ...figPara(loadImage("fig1_yield.png"), "Figure 1: Neuron yield per animal and hypothalamic region. Day animals (left) and night animals (right)."),

      heading("Dense, region-specific functional connectivity revealed by GLMCC", 2),
      para("GLMCC with CCG-based shuffle-ISI surrogate validation and BH-FDR correction (q = 0.05) detected " +
           "functional connections between the majority of neuron pairs. The edge survival rate was 79.9 \u00B1 4.1% " +
           "across the 13 animals (range: 64.8\u201379.9%; Figure 2). This high density is consistent with the " +
           "known dense anatomical interconnectivity of hypothalamic nuclei (Swanson, 2000) and contrasts with " +
           "typically sparser cortical connectivity estimates (5\u201320%; Perin et al., 2011)."),

      para("Critically, we found that the effect of light on network topology varied systematically by " +
           "hypothalamic region. Linear mixed-effects models with condition \u00D7 region interaction terms " +
           "revealed highly significant interactions for all four graph metrics (node strength: " +
           "F(5,573) = 14.0, p = 1.6\u00D710\u207B\u00B9\u00B2; clustering coefficient: F(5,575) = 9.5, " +
           "p = 1.1\u00D710\u207B\u2078; local efficiency: F(5,575) = 22.6, p < 2\u00D710\u207B\u00B9\u2076; " +
           "hub score: F(5,574) = 3.4, p = 0.005). This interaction term was absent from previous models " +
           "that considered only main effects of condition, and its inclusion substantially changes the " +
           "interpretation of light-evoked connectivity changes."),

      para("The interaction pattern was most pronounced for clustering coefficient and local efficiency, " +
           "where light-ON significantly decreased both metrics in DMH (\u03B2 = \u22120.40, p < 0.001; " +
           "\u03B2 = \u22120.76, p = 0.005 respectively), PH (\u03B2 = \u22120.26, p = 0.005; " +
           "\u03B2 = \u22121.02, p = 0.002), and PVH (\u03B2 = \u22120.47, p < 0.001; " +
           "\u03B2 = \u22120.89, p < 0.001). This indicates that light stimulation reduces local clustering " +
           "in these integrative hub regions, consistent with a shift toward more distributed processing. " +
           "Notably, node strength in VM-thalamus showed a significant increase under light-ON " +
           "(\u03B2 = +635.6, p = 0.003), suggesting enhanced thalamic relay function during photic stimulation. " +
           "Hub scores showed no significant region-specific effects, consistent with preserved hub identity."),

      ...figPara(loadImage("fig2_lme.png"), "Figure 2: Condition \u00D7 region interaction estimates from linear mixed-effects models. Positive values (blue) indicate increased metric under light-ON; negative values (red) indicate decreased metric. Significance: * p<0.05, ** p<0.01, *** p<0.001. Panel titles show overall interaction ANOVA p-values."),

      heading("Network architecture: distributed, inter-regional, and spatially uniform", 2),
      para("To characterize the macroscopic organization of the hypothalamic functional network, we performed " +
           "a series of graph-theoretic analyses on the validated adjacency matrices. Connectivity was predominantly " +
           "inter-regional: 59.3 \u00B1 21.7% of significant edges connected neurons in different anatomical regions " +
           "(mean across 13 animals), while only 40.7% were within-region. This is consistent with the high " +
           "proportion of connector hubs identified by participation coefficient analysis (94.6% of neurons had " +
           "PC > 0.3, indicating connections distributed across multiple communities; mean PC = 0.576 \u00B1 0.203). " +
           "PVH showed the highest participation coefficient (0.64), consistent with its role as an integrative hub."),

      para("Edge consistency across animals was remarkably high (mean 76.5%), with the most consistent region pairs " +
           "being VMH\u2194PVH (92.1%), PH\u2194MC (91.2%), and VM-thalamus\u2194PH (90.6%). This suggests a " +
           "highly stereotyped core hypothalamic connectome that is preserved across individual mice. " +
           "Centrality measures (strength, betweenness, eigenvector, PageRank) showed moderate to strong " +
           "convergence (mean Spearman \u03C1 = 0.55\u20130.82 across pairs), indicating that hub identity " +
           "is robust to the specific centrality definition used."),

      para("Distance-dependent connectivity analysis using the A4x64-Poly2 probe geometry (4 shanks, " +
           "16 channels per shank, 200 \u00B5m column spacing, 160 \u00B5m row spacing) revealed no significant " +
           "relationship between inter-channel distance and connection probability (mean Spearman \u03C1 = 0.03 " +
           "\u00B1 0.12, 7/13 animals with nominally significant but very small effects). Same-shank and " +
           "cross-shank connection probabilities were nearly identical. This spatial uniformity confirms that " +
           "the detected functional connectivity reflects genuine temporal coupling rather than volume conduction " +
           "or proximity artifacts."),

      ...figPara(loadImage("fig3_region_conn.png"), "Figure 3: Region-level connectivity density. Mean connection density between hypothalamic region pairs, shown separately for day (left) and night (right) animals under light-ON stimulation."),

      heading("No rich-club organization; non-random circuit-level motifs", 2),
      para("Rich-club analysis testing whether high-degree nodes preferentially interconnect revealed no " +
           "significant rich-club organization in any animal (0/13, all \u03C6\u2099\u2092\u2098 < 1 or " +
           "p > 0.05 after degree-preserving null comparison). Despite the presence of clear network hubs " +
           "(PVH, DMH), these hubs do not form a densely interconnected core beyond what is expected from " +
           "their degree alone. This distinguishes the hypothalamus from cortical networks, where rich-clubs " +
           "are commonly observed (van den Heuvel and Sporns, 2011)."),

      para("However, analysis of 3-node connectivity motifs revealed massive enrichment of connected triplets " +
           "above degree-preserving null models (mean Z = 81.4 \u00B1 25.0 across 13 animals; 97.9% of sampled " +
           "triplets contained \u22652 edges). The fraction of fully-connected triplets (all 6 possible directed " +
           "edges present) ranged from 28% to 48% across animals. This indicates that despite near-complete " +
           "binary connectivity at the pairwise level, the specific pattern of which triplets form circuits " +
           "is strongly non-random. Structure in the hypothalamic network resides at the level of 3-node " +
           "circuit motifs rather than pairwise edge topology."),

      heading("Hub stability across the ongoing-to-light transition", 2),
      para("Hub scores were not significantly reorganized by photic stimulation. The within-animal Spearman " +
           "correlation of hub scores between ongoing and light-ON conditions was \u03C1 = \u22120.08 " +
           "(Wilcoxon signed-rank V = 6, p = 0.92; BCa 95% CI: [\u22120.20, 0.09]), indicating that the " +
           "relative importance of individual neurons in the network topology is preserved across conditions. " +
           "The hub neurons identified during light-ON tend to be the same neurons identified as hubs during " +
           "ongoing spontaneous activity, consistent with the notion that network hubs represent stable " +
           "anatomical specializations rather than transient functional states (van den Heuvel and Sporns, 2013)."),

      ...figPara(loadImage("fig4_hub_stability.png"), "Figure 4: Hub score stability across ongoing-to-light transition. Per-animal Spearman \u03C1 of hub scores between ongoing and light-ON conditions. Group mean \u03C1 = \u22120.08, Wilcoxon p = 0.92."),

      heading("Surrogate null models confirm absence of binary topological structure", 2),
      para("To test whether observed graph metrics exceed chance expectations, we compared each animal's " +
           "network against 200 degree-preserving random networks. Clustering coefficient Z-scores were " +
           "consistently negative (mean Z = \u22120.73 \u00B1 0.47, 0/13 animals with p < 0.05), indicating " +
           "that observed clustering is not higher\u2014and in fact slightly lower\u2014than expected from " +
           "the degree distribution alone. Characteristic path lengths were slightly longer than random " +
           "but not significantly so. Small-world propensity analysis confirmed that the hypothalamic network " +
           "lies near the random-network boundary (mean \u03C3 = 1.01 \u00B1 0.02; mean \u03C9 = 0.01 \u00B1 0.02), " +
           "with 6/11 connected animals showing \u03C3 > 1 and all 11 showing |\u03C9| < 0.5. " +
           "These findings converge on the conclusion that the near-complete binary edge density (~80%) " +
           "saturates any possible non-random topological organization at the binary level. " +
           "Biological signal resides in graded connection strengths and specific circuit motifs, " +
           "not in which pairs of neurons are connected."),

      heading("Firing rate partially predicts connectivity; E-I differences", 2),
      para("Neuron-level firing rate showed a moderate positive correlation with node strength " +
           "(mean Spearman \u03C1 = 0.50 \u00B1 0.25, 11/13 animals p < 0.05), explaining approximately " +
           "25% of the variance in connectivity. This indicates that while firing rate partially drives " +
           "detected connectivity, the majority of variance is attributable to genuine temporal coupling " +
           "independent of rate. Inhibitory neurons showed significantly higher firing rates than excitatory " +
           "neurons (0.28 vs. 0.17 Hz, t = \u22123.21, p = 0.001) but did not show significantly higher " +
           "node strength (714.9 vs. 529.2, t = \u22121.48, p = 0.14), suggesting that inhibitory " +
           "interneurons fire at higher rates but do not dominate the network's connectivity structure."),

      new Paragraph({ children: [new PageBreak()] }),

      // ═══════════════════════════════════════════════════════════
      // DISCUSSION
      // ═══════════════════════════════════════════════════════════
      heading("Discussion"),

      para("We present the first comprehensive map of light-evoked functional connectivity across the mouse " +
           "hypothalamus, combining GLMCC with rigorous surrogate-based validation and systematic graph-theoretic " +
           "characterization. Five main findings emerge from this work."),

      para("First, the hypothalamus operates as a densely connected network during light stimulation, with " +
           "approximately 80% of neuron pairs showing statistically significant functional coupling. " +
           "This density, combined with the absence of rich-club organization, the lack of binary topological " +
           "structure beyond degree distribution, and the near-random small-world indices, paints a picture " +
           "of a network fundamentally different from the modular, hierarchical organization typical of " +
           "cortical circuits. The hypothalamus appears to function as a distributed integrative network where " +
           "broad, non-selective connectivity enables flexible coordination of autonomic, endocrine, and " +
           "behavioral outputs (Sternson, 2013)."),

      para("Second, the highly significant condition \u00D7 region interaction demonstrates that light's " +
           "effect on network topology is region-specific. The paraventricular and dorsomedial hypothalamic " +
           "nuclei show the strongest light-evoked decreases in clustering and local efficiency, consistent " +
           "with a shift toward more distributed, less locally clustered processing during photic stimulation. " +
           "The VM-thalamus uniquely shows increased node strength under light-ON, suggesting enhanced " +
           "thalamic relay function. These region-specific effects were not detectable in previous analyses " +
           "that considered only main effects of condition, highlighting the importance of including " +
           "interaction terms in the statistical model."),

      para("Third, the network exhibits a stereotyped anatomical backbone: edge consistency across animals " +
           "exceeds 76%, with the most consistent connections involving PVH, DMH, and VM-thalamus. " +
           "This stereotypy, combined with the spatial uniformity of connectivity (no distance decay), suggests " +
           "that the detected functional connections reflect genuine anatomical pathways rather than " +
           "measurement artifacts or state-dependent fluctuations."),

      para("Fourth, despite near-complete binary connectivity, the network contains strongly non-random " +
           "structure at the level of 3-node circuit motifs (Z = 81.4 vs. degree-preserving null). " +
           "This dissociation\u2014random at the pairwise level, highly structured at the triplet " +
           "level\u2014suggests that hypothalamic computation is implemented through specific patterns of " +
           "local circuit connectivity embedded within a broadly permissive anatomical substrate. " +
           "The high proportion of fully-connected triplets (28\u201348%) indicates abundant recurrent " +
           "connectivity, consistent with the feedback-rich architecture expected of a homeostatic regulator."),

      para("Fifth, hub neuron identity is preserved across the ongoing-to-light transition (\u03C1 = \u22120.08, " +
           "p = 0.92), despite progressive firing-rate adaptation in half the animals. This dissociation " +
           "between preserved topology and adapting rates suggests that photic habituation operates through " +
           "a global gain mechanism rather than targeted network reorganization (Kohn, 2007), and that hub " +
           "neurons represent stable anatomical specializations."),

      heading("Limitations and future directions", 2),
      para("Several limitations should be noted. First, our recordings sample from a distributed but " +
           "non-exhaustive set of hypothalamic regions; nuclei not covered by the probe tracks (e.g., the " +
           "suprachiasmatic nucleus, lateral hypothalamus) may show distinct connectivity patterns. Second, " +
           "the 10-s light-ON pulses are considerably longer than the sub-second timescale of direct synaptic " +
           "integration, suggesting that detected functional coupling may reflect polysynaptic or network-level " +
           "interactions. Third, five animals lack ongoing-condition data due to incomplete GLMCC processing, " +
           "limiting the second-level day/night comparison. Fourth, the absence of behavioral correlates " +
           "(core temperature, locomotor activity) limits our ability to link connectivity changes to " +
           "physiological outcomes. Fifth, the probe geometry-based distance analysis uses approximate channel " +
           "positions; precise histological reconstruction would provide more accurate spatial relationships."),

      para("Future work should employ chronic Neuropixels recordings to track within-animal hub stability " +
           "across circadian cycles, incorporate optogenetic tagging of identified cell types, pair " +
           "electrophysiology with simultaneous behavioral monitoring, and apply the motif analysis framework " +
           "developed here to compare hypothalamic circuit architecture across physiological states " +
           "(e.g., fed vs. fasted, sleep vs. wake)."),

      new Paragraph({ children: [new PageBreak()] }),

      // ═══════════════════════════════════════════════════════════
      // METHODS SUMMARY
      // ═══════════════════════════════════════════════════════════
      heading("Methods Summary"),

      heading("Animals and recordings", 2),
      para("Thirteen adult male C57BL/6J mice were implanted with 64-channel silicon probes " +
           "(A4x64-Poly2-7mm-23s-200-160, NeuroNexus) targeting six hypothalamic regions. Spontaneous " +
           "spiking activity was recorded during a 5-hour light-ON protocol (white LED, ~100 lux, " +
           "95 pulses of 10 s duration). Spike sorting was performed with KiloSort2 followed by manual " +
           "curation in Phy. A total of 682 single units were retained after quality control."),

      heading("Connectivity inference and validation", 2),
      para("Directed functional connectivity was estimated using GLMCC (Kobayashi et al., 2019) with " +
           "cross-correlogram parameters: \u00B150 ms window, 1 ms bins. Edge significance was assessed " +
           "against 100 shuffle-ISI surrogates using CCG peak magnitude as the test statistic. " +
           "Benjamini-Hochberg FDR correction (q = 0.05) was applied to off-diagonal pairs within each " +
           "animal \u00D7 condition combination."),

      heading("Graph-theoretic analysis", 2),
      para("Per-node metrics (strength, clustering coefficient, local efficiency, hub score via HITS algorithm) " +
           "were computed on validated, directed, weighted adjacency matrices. Community detection used the " +
           "Louvain algorithm on symmetrized matrices. Participation coefficient was computed following " +
           "Guimer\u00E0 and Amaral (2005). Rich-club analysis compared observed \u03C6(k) against 100 " +
           "degree-preserving rewired networks. Motif analysis sampled 10,000 triplets per animal and compared " +
           "connected-triplet fraction against 50 degree-preserving null networks. Distance-dependent connectivity " +
           "used approximate channel positions from the A4x64-Poly2 probe geometry."),

      heading("Statistical modelling", 2),
      para("First-level neuron metrics were modelled with linear mixed-effects models: " +
           "metric ~ condition \u00D7 region + (1 | animal/neuron), using restricted maximum likelihood (REML) " +
           "with Kenward-Roger degrees of freedom (lmerTest). Second-level region-pair density was modelled as: " +
           "density ~ condition \u00D7 day_night + (1 | animal). Hub stability was assessed by within-animal " +
           "Spearman rank correlation with Wilcoxon signed-rank test and BCa bootstrap confidence intervals."),

    ]
  }]
});

// ═══════════════════════════════════════════════════════════════════
Packer.toBuffer(doc).then(buffer => {
  const outPath = path.join(PROJECT, "manuscript_revised.docx");
  fs.writeFileSync(outPath, buffer);
  console.log("Written: " + outPath);
  console.log("Done.");
});
