# K-BRIDGE and SPINK5 analysis code

This release accompanies *Interpretable B-cell state modelling and SPINK5 survival associations in head and neck cancer*. It contains the analysis scripts, fixed model settings and instructions for rebuilding the analysis inputs. The companion aggregate results include the gene-to-program definitions and model and figure summaries. See the published version metadata for the creators and citation.

The source datasets are [GSE214265](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE214265), [GSE137275](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE137275), [GSE65858](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65858), [GSE103322](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103322) and [TCGA-HNSC on UCSC Xena](https://tcga.xenahubs.net). The fetch scripts download public data from their original providers. Access and reuse are governed by the respective providers and depositors. The input manifest records the original sources and available SHA-256 hashes. Commands write inputs and results locally. Do not add these generated records to the public repository.

## K-BRIDGE internal and independent validation

Use Python 3.9 and the pinned dependencies in `kbridge/requirements.txt`. Run from the software directory.

```sh
python -m pip install -r kbridge/requirements.txt
python fetch_public_inputs.py kbridge
python kbridge/reproduce_internal.py --permutations
python kbridge/run_validation.py --data kbridge/data --output kbridge/results
python kbridge/run_external.py --data kbridge/data --external kbridge/external --output kbridge/results
```

The commands regenerate out-of-fold estimates, paired bootstrap contrasts, repeated partitions, 100 permutation fits and the frozen-model evaluation in the two human tonsil specimens. The file `kbridge/kbridge_model.py` defines the gene programs, information gate, graph and decoder. `kbridge/README.md` records the exact public inputs and transformations.

## Survival and incremental discrimination

Use a separate environment for the dependencies in `public_hnsc/requirements.txt`. Run from the software directory.

```sh
python -m pip install -r public_hnsc/requirements.txt
python public_hnsc/download_inputs.py
python fetch_public_inputs.py external_hnsc
python public_hnsc/reanalyse.py
python external_hnsc/run_external_analysis.py
python external_hnsc/check_ph_sensitivity.py
python fetch_public_inputs.py public_hnsc_sensitivity
python helpers/rebuild_hnsc_inputs.py
python incremental_prognosis/run_incremental.py --check-inputs
python incremental_prognosis/run_incremental.py
python incremental_prognosis/validate_outputs.py
```

The input-preparation command recreates the inputs used by the incremental and sensitivity analyses. It validates the upstream outputs before writing the files. That runner checks the input files against the SHA-256 hashes in `incremental_prognosis/source_manifest.json`. Figure 3 is generated in `incremental_prognosis/results`.

## Additional survival models and cell-compartment summaries

Keep the outputs from the previous steps. Use a separate Python environment with the dependencies in `public_hnsc_sensitivity/requirements.txt`.

```sh
python -m pip install -r public_hnsc_sensitivity/requirements.txt
python public_hnsc_sensitivity/run_survival_sensitivity.py
python public_hnsc_sensitivity/report_statistical_details.py
python public_hnsc_sensitivity/run_cellular_source.py
python plot_cell_detection.py public_hnsc_sensitivity/results/puram_patient_detection.tsv figure4_cell_resolved_spink5.pdf
```

This regenerates the survival-sensitivity models, composition scores, test statistics and intervals, cell detection summaries and Figure 4. The MCP-counter marker definitions are fetched from the versioned upstream commit recorded in `fetch_public_inputs.py`. The original sample-level public datasets are available at their source accessions. The downloadable companion results contain aggregate analysis outputs and model-feature summaries.
