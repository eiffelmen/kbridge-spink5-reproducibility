# K-BRIDGE

The input files and results support the reported within-study classification and external tonsil evaluation. The numerical environment used Python 3.9 with the versions in requirements.txt.

## Public input files

- data/GSE214265_expression.txt.gz is the expression file GSE214265_R1_expression.matrix.txt.gz.
- data/GSE214265_metadata.csv.gz is the metadata file GSE214265_R1_metadata.csv.gz.
- external/ contains the two original GEO endogenous-UMI matrices and matching phenotype workbooks for GSM4074307 and GSM4074308. Only the two human tonsil specimens enter this analysis. The workbook column PhenoOriginal records the original sorting phenotype.

The two GSE214265 files come from the same GEO series. In this analysis, duplicate expression rows are collapsed by uppercase gene symbol, the fixed 64-gene panel is selected and the values are transformed with log(1+x). The matched files contain 159 B-cell and 117 ASC transcriptomic labels.

The external analysis compares GC and Mem B-cell phenotypes with PB/PC phenotypes from the original FACS sorting. The external counts are library-size normalized and log-transformed; all 64 source genes are measured in both specimens. Frozen source models are assessed separately for each specimen. The sensitivity analysis standardizes genes using only unlabelled target expression.

- [GEO GSE214265](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE214265)
- [GEO GSE137275](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE137275)

See ../README.md for execution commands. The original public input checksums are listed in input_checksums.sha256.
