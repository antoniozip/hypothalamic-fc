# Snakemake workflow for hypothalamic FC analysis.
#
# Keyed by wildcards {animal}, {condition}, {estimator}.
# Targets:
#   data/processed/spiketrains_{animal}_{condition}.txt
#   data/surrogates/{animal}_{condition}/{surr_id}/cell_*.txt
#   results/{estimator}/adj_{animal}_{condition}.csv
#   results/{estimator}/validated_adj_{animal}_{condition}.csv
#   results/{estimator}/metrics_{animal}_{condition}.csv
#   results/{estimator}/density_{animal}_{condition}.csv
#
# Usage: snakemake -n   (dry run)
#        snakemake all  (full run)

ANIMALS = [
    "171019", "171207", "171208", "171213", "180110", "180111",
    "180131", "180221", "180228", "180302", "180419", "180420", "180423"
]
CONDITIONS = ["ongoing", "lightON"]
ESTIMATORS = ["glmcc", "te"]
N_SURROGATES = 100

rule all:
    input:
        expand("results/{estimator}/density_{animal}_{condition}.csv",
               estimator=ESTIMATORS, animal=ANIMALS, condition=CONDITIONS),

rule validate_edges:
    input:
        real = "results/{estimator}/adj_{animal}_{condition}.csv",
        surrogates = expand(
            "results/{estimator}/surrogate_adj/{animal}_{condition}/{surr_id}.csv",
            surr_id=range(N_SURROGATES)
        ),
    output:
        "results/{estimator}/validated_adj_{animal}_{condition}.csv",
    shell:
        "python src/py/validate_edges.py --animal {wildcards.animal} "
        "--condition {wildcards.condition} --estimator {wildcards.estimator}"

rule compute_graph_metrics:
    input:
        "results/{estimator}/validated_adj_{animal}_{condition}.csv",
    output:
        "results/{estimator}/metrics_{animal}_{condition}.csv",
    shell:
        "Rscript src/r/compute_graph_metrics.R --animal {wildcards.animal} "
        "--condition {wildcards.condition} --estimator {wildcards.estimator}"

rule aggregate_density:
    input:
        "results/{estimator}/validated_adj_{animal}_{condition}.csv",
    output:
        "results/{estimator}/density_{animal}_{condition}.csv",
    shell:
        "Rscript src/r/aggregate_density.R --animal {wildcards.animal} "
        "--condition {wildcards.condition} --estimator {wildcards.estimator}"

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

rule figures:
    input:
        expand("results/stats_first_level_{estimator}.csv", estimator=ESTIMATORS),
        expand("results/stats_second_level_{estimator}.csv", estimator=ESTIMATORS),
        expand("results/{estimator}/density_{animal}_{condition}.csv",
               estimator=ESTIMATORS, animal=ANIMALS, condition=CONDITIONS),
    output:
        "figures/main/fig2_metrics.png",
        "figures/main/fig3_circos.png",
        "figures/main/fig5_convergence.png",
    shell:
        "Rscript src/r/make_figures.R"
