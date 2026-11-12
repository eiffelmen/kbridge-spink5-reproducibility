# GSE65858 independent SPINK5 survival replication plan

Written 2026-09-25 before loading expression values or fitting survival models.

- Single externally specified cohort GSE65858, selected for public HNSCC survival/HPV data; no result-based cohort search.
- Primary sources: GEO series GSE65858 and platform GPL10558; retain processed log2/RSN matrix.
- SPINK5 score: select probes with exact official SPINK5 gene-symbol annotation, exclude ambiguous multiple-gene symbol annotations, take each sample's median across all retained SPINK5 probes. No survival-based probe selection. Report mapping.
- Endpoint: overall survival, exact GEO time/event definitions retained. Restrict to positive nonmissing survival time and known 0/1 death event; no survival cutoff optimization.
- Primary model: continuous SPINK5 standardized by sample SD, adjusted for age (per 10 years), stage categorical, sex, and available HPV status. Prefer transcriptionally active HPV (DNA+/RNA+) if clearly annotated; otherwise official HPV definition. Missing covariates complete-case, tabulate excluded samples. No outcome imputation.
- Additional descriptive models: univariate SPINK5 and age/stage/sex adjustment, each reported independently with n, events, HR per SPINK5 SD, 95% CI, P, concordance, PH checks. No stepwise selection. These are secondary sensitivity checks, not alternative primary findings.
- Check rank-transform Schoenfeld tests, exact sample uniqueness, missingness, expression distribution. A non-replication is retained.
- Download cap 300 MB. Stop if survival or unambiguous probe annotation is unavailable.

## Metadata-only clarification before any expression–survival fit

The downloaded header identifies 253 Primary, 14 Secondary and 3 Relapse tumors. The primary analysis is restricted to `tumor_type: Primary`, matching TCGA primary-tumor sampling and avoiding ambiguous timing of recurrent/secondary sampling. A fully labelled all-tumor sensitivity model will be retained. UICC IVA/IVB/IVC map to IV, yielding stage I/II/III/IV categories with I reference. HPV-active = `hpv16_dna_rna: DNA+RNA+`; inactive = `DNA-` or `DNA+RNA-`; NA excluded from HPV-adjusted models (no inference from other HPV DNA categories). This clarification uses metadata only, before gene-expression outcome associations.
