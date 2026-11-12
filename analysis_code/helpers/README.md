# Rebuild HNSCC analysis inputs

Run from the release root with the Python environment pinned in `public_hnsc/requirements.txt` (including pandas 2.3.3). The helper rebuilds incremental and sensitivity inputs from the two upstream analyses and from the public files fetched by `fetch_public_inputs.py`.

```bash
python public_hnsc/download_inputs.py
python public_hnsc/reanalyse.py
python fetch_public_inputs.py external_hnsc
python external_hnsc/run_external_analysis.py
python fetch_public_inputs.py public_hnsc_sensitivity
python helpers/rebuild_hnsc_inputs.py
python incremental_prognosis/run_incremental.py --check-inputs
python incremental_prognosis/run_incremental.py
python incremental_prognosis/validate_outputs.py
python public_hnsc_sensitivity/run_survival_sensitivity.py
python public_hnsc_sensitivity/report_statistical_details.py
python public_hnsc_sensitivity/run_cellular_source.py
```

If using the separate dependency environments documented by each module, keep the helper and upstream TCGA/GSE65858 scripts in the pinned `public_hnsc` environment. Follow `public_hnsc_sensitivity/requirements.txt` for the last three downstream scripts.

The helper performs no downloads. Its defaults are the locations populated by the public-input fetcher. Use `--xena-expression`, `--puram-matrix` and `--mcp-signatures` to use previously obtained copies. Use `--root` if running the helper outside the release tree. It checks the upstream SHA-256 values in the incremental manifest, verifies the TCGA design matrix before adding sample identifiers, checks sample and event counts for both models, and emits SHA-256 values for all 12 downstream inputs. Re-running against identical files leaves them in place. If files differ, the helper stops before writing any outputs.

The public-source files, generated model inputs and intermediate results are local analysis products. The public code release only contains scripts, dependency specifications and source instructions.
