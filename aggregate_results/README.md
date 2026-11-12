# Aggregate results for the public-data HNSCC study

These tables accompany *Interpretable B-cell state modelling and SPINK5 survival associations in head and neck cancer*. They contain aggregate estimates, diagnostics, uncertainty and model-feature summaries for the analyses in the manuscript.

- kbridge/ contains the model and ablation metrics, repeat and permutation summaries, model coefficients, gene-to-program map, fold-level information gates, paired comparisons and pooled external assessment.
- tcga_hnsc/ and gse65858/ contain cohort model estimates, covariate missingness, proportional-hazards checks and numerical validation.
- incremental_discrimination/ contains cohort and repetition summaries, the aggregate resampling distribution, model coefficients and model-fit comparisons.
- clinical_and_cellular_sensitivities/ contains cohort-level sensitivity, composition, cell-compartment and statistical-reporting summaries.

Each table can be regenerated from the public inputs using the companion analysis software and its instructions. The source repositories and analysis lineage are listed in SOURCE_DATA.md. MANIFEST.tsv lists table-level provenance and SHA-256 hashes. EXCLUSIONS.tsv documents which members from the full set of review materials were omitted.

This release contains no records keyed by individual participant, specimen, donor or cell. In particular, the inputs, split membership, individual model predictions, PCA coordinates and detection values by patient are omitted. The original public inputs remain available at their source repositories. Their applicable terms and attribution are independent of these derived outputs.
