# Snakemake workflow for hypothalamic FC analysis.
#
# Keyed by wildcards {animal}, {condition}, {estimator}.
#
# Full DAG targets:
#   data/processed/spiketrains_{animal}_{condition}.txt
#   data/surrogates/{animal}_{condition}/{surr_id}/cell_*.txt
#   results/{estimator}/adj_{animal}_{condition}.csv
#   results/{estimator}/validated_adj_{animal}_{condition}.csv
#   results/{estimator}/metrics_{animal}_{condition}.csv
#   results/{estimator}/density_{animal}_{condition}.csv
#   results/stats_first_level_{estimator}.csv
#   results/stats_second_level_{estimator}.csv
#   figures/main/fig1_yield.png
#   figures/main/fig2_metrics.png
#   figures/main/fig3_circos.png
#   figures/main/fig4_hub_stability.png
#   figures/main/fig5_convergence.png
#
# Usage: snakemake -n        (dry run)
#        snakemake all -j4    (full run with 4 cores)

ANIMALS = [
    "171019", "171207", "171208", "171213", "180110", "180111",
    "180131", "180221", "180228", "180302", "180419", "180420", "180423"
]
CONDITIONS = ["ongoing", "lightON"]
ESTIMATORS = ["glmcc", "te"]
N_SURROGATES = 100


# ---------------------------------------------------------------------------
# Top-level target
# ---------------------------------------------------------------------------

rule all:
    input:
        expand("results/{estimator}/density_{animal}_{condition}.csv",
               estimator=ESTIMATORS, animal=ANIMALS, condition=CONDITIONS),
        expand("results/stats_first_level_{estimator}.csv", estimator=ESTIMATORS),
        expand("results/stats_second_level_{estimator}.csv", estimator=ESTIMATORS),
        "figures/main/fig1_yield.png",
        "figures/main/fig2_lme.png",
        "figures/main/fig3_region_conn.png",
        "figures/main/fig4_hub_stability.png",
        "figures/main/fig5_convergence.png",


# ---------------------------------------------------------------------------
# Phase 0 — Spike train extraction
# ---------------------------------------------------------------------------

rule extract_spike_trains:
    """MAT → processed spike train txt (placeholder — runs existing extraction)."""
    input:
        "data/raw/DATA{animal}_{condition}.mat",
    output:
        "data/processed/spiketrains_DATA{animal}_{condition}.txt",
    shell:
        "echo 'Extracting spike trains for {wildcards.animal}/{wildcards.condition}' && "
        "touch {output}"  # placeholder; actual extraction uses existing txt files


rule build_neuron_table:
    """Build canonical neuron table from cell files."""
    input:
        expand("DATA{animal}/cell_{{n}}.txt", animal=ANIMALS),  # placeholder
    output:
        "data/processed/neurons.csv",
    shell:
        "python src/py/build_neuron_table.py"


# ---------------------------------------------------------------------------
# Phase 1 — Surrogate generation
# ---------------------------------------------------------------------------

rule generate_surrogates:
    """Generate 100 shuffle-ISI surrogates per animal/condition."""
    input:
        "data/processed/spiketrains_DATA{animal}_{condition}.txt",
    output:
        expand(
            "data/surrogates/{animal}_{condition}/{sid}/cell_{i}.txt",
            sid=range(N_SURROGATES), i=range(64)  # max 64 channels
        ),
    shell:
        "python src/py/generate_surrogates.py "
        "--animal {wildcards.animal} --condition {wildcards.condition} "
        "--n {N_SURROGATES}"


# ---------------------------------------------------------------------------
# Phase 2 — Adjacency computation
# ---------------------------------------------------------------------------

rule run_glmcc_adjacency:
    """Run GLMCC to compute raw adjacency matrix.
    
    Set GLMCC_USE_DOCKER=1 to use containerized GLMCC.
    Otherwise uses local GLMCC_HOME installation.
    """
    input:
        "data/processed/spiketrains_DATA{animal}_{condition}.txt",
    output:
        "results/glmcc/adj_{animal}_{condition}.csv",
    shell:
        "python src/py/run_glmcc.py "
        "--animal {wildcards.animal} --condition {wildcards.condition}"


rule run_te_adjacency:
    """Run Transfer Entropy to compute raw adjacency matrix."""
    input:
        "data/processed/spiketrains_DATA{animal}_{condition}.txt",
    output:
        "results/te/te_{animal}_{condition}_k10.csv",
    shell:
        "python src/py/run_transfer_entropy.py "
        "--animal {wildcards.animal} --condition {wildcards.condition} --k 10"


# ---------------------------------------------------------------------------
# Phase 3 — Edge validation
# ---------------------------------------------------------------------------

rule validate_edges_ccg:
    """Validate edges using CCG peaks + shuffle-ISI surrogates (lightON)."""
    input:
        "results/{estimator}/adj_{animal}_{condition}.csv",
    output:
        "results/{estimator}/validated_adj_{animal}_{condition}.csv",
    shell:
        "python src/py/validate_edges.py "
        "--animal {wildcards.animal} --condition {wildcards.condition} "
        "--estimator {wildcards.estimator} --method ccg"


# ---------------------------------------------------------------------------
# Phase 4 — Graph metrics
# ---------------------------------------------------------------------------

rule compute_graph_metrics:
    input:
        "results/{estimator}/validated_adj_{animal}_{condition}.csv",
    output:
        "results/{estimator}/metrics_{animal}_{condition}.csv",
    shell:
        "Rscript src/r/compute_graph_metrics.R --animal {wildcards.animal} "
        "--condition {wildcards.condition} --estimator {wildcards.estimator}"


# ---------------------------------------------------------------------------
# Phase 5 — Region-pair density
# ---------------------------------------------------------------------------

rule aggregate_density:
    input:
        "results/{estimator}/validated_adj_{animal}_{condition}.csv",
    output:
        "results/{estimator}/density_{animal}_{condition}.csv",
    shell:
        "Rscript src/r/aggregate_density.R --animal {wildcards.animal} "
        "--condition {wildcards.condition} --estimator {wildcards.estimator}"


# ---------------------------------------------------------------------------
# Phase 6 — Statistics
# ---------------------------------------------------------------------------

rule first_level_stats:
    input:
        expand("results/{estimator}/metrics_{animal}_{condition}.csv",
               animal=ANIMALS, condition=CONDITIONS),
    output:
        "results/stats_first_level_{estimator}.csv",
    shell:
        "Rscript src/r/stats_first_level.R --estimator {wildcards.estimator}"


rule second_level_stats:
    input:
        expand("results/{estimator}/density_{animal}_{condition}.csv",
               animal=ANIMALS, condition=CONDITIONS),
    output:
        "results/stats_second_level_{estimator}.csv",
    shell:
        "Rscript src/r/stats_second_level.R --estimator {wildcards.estimator}"


rule hub_stability:
    input:
        expand("results/{estimator}/metrics_{animal}_{condition}.csv",
               estimator=["glmcc"], animal=ANIMALS, condition=CONDITIONS),
    output:
        "results/hub_stability.csv",
    shell:
        "Rscript src/r/hub_stability.R"


# ---------------------------------------------------------------------------
# Phase 7 — Figures
# ---------------------------------------------------------------------------

rule fig1_yield:
    """Figure 1: Paradigm, neuron yield heatmap, stimulus protocol."""
    input:
        "data/processed/neurons.csv",
    output:
        "figures/main/fig1_yield.png",
    shell:
        "Rscript src/r/make_figures.R --figure fig1"


rule fig4_hub_stability:
    """Figure 4: Hub stability across conditions."""
    input:
        "results/hub_stability.csv",
    output:
        "figures/main/fig4_hub_stability.png",
    shell:
        "Rscript src/r/hub_stability.R"


rule figures:
    """Figures 2, 3, 5: metrics, circos, convergence."""
    input:
        expand("results/stats_first_level_{estimator}.csv", estimator=ESTIMATORS),
        expand("results/stats_second_level_{estimator}.csv", estimator=ESTIMATORS),
        expand("results/{estimator}/density_{animal}_{condition}.csv",
               estimator=ESTIMATORS, animal=ANIMALS, condition=CONDITIONS),
    output:
        "figures/main/fig2_lme.png",
        "figures/main/fig3_region_conn.png",
        "figures/main/fig5_convergence.png",
    shell:
        "Rscript src/r/make_figures.R"
