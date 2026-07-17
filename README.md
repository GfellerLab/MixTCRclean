# MixTCRclean
R package to predict inconsistencies between V/J annotations and CDR3 sequences in TCR repertoire data.


# INTRODUCTION

MixTCRclean is a robust quality-control tool that detects inconsistencies between the V/J annotations and the CDR3 sequences in TCR data. 
It is based on the expected number of conserved germline residues of V and J genes that was derived by statistical modelling of a large dataset of bulk TCR repertoire data.

The package also provides standardization and pre-cleaning steps.

## Installation

```r
# install.packages("remotes")
remotes::install_github("GfellerLab/MixTCRclean")
```

This installs `MixTCRviz` and `ggseqlogoMOD` automatically.


### Manual installation

If the command above fails, install the dependencies in order:

* Install the `ggseqlogoMOD` package from https://github.com/GfellerLab/ggseqlogo (this is needed even if you already have the standard `ggseqlogo` package).
```r
# install.packages("remotes")
remotes::install_github("GfellerLab/ggseqlogo")
```
* Download the `MixTCRviz` directory from https://github.com/GfellerLab/MixTCRviz, open RStudio setting its working directory as the MixTCRviz folder, then:

```r
devtools::build()
install.packages("../MixTCRviz_1.0.tar.gz", repos = NULL)
```

* Download the `MixTCRclean` directory from https://github.com/GfellerLab/MixTCRclean, set it as your working directory, then:

```r
devtools::build()
install.packages("../MixTCRclean_0.0.0.9000.tar.gz", repos = NULL)
```


### Testing the INSTALLATION:

In the MixTCRclean directory run:

`Rscript test_MixTCRclean.R`

Alternatively, you can run the code test_MixTCRclean.R in Rstudio (or any R interface)

The output in test/out should be the same as in test/out_compare

# RUNNING

MixTCRclean should be primarily run in R, by loading the MixTCRclean library and calling the MixTCRclean function (e.g., MixTCRclean(input="test/test.csv, output.path=YOUR_OUTPUT_PATH)).


## Mandatory parameter:

 - input: Input can one of the different values:
   1) A .csv or .txt or .tsv file with the input TCRs
   2) A data.frame with the input TCRs
   3) A list in the MixTCRclean format.

   If using a filename or a data.frame:
    * Columns should ideally consist of "TRAV","TRAJ","cdr3_TRA","TRBV","TRBJ","cdr3_TRB".
    * The "TRAV", "TRAJ", "TRBV", "TRBJ" entries should follow the IMGT
   nomenclature, with or without allele (see below for potential name correction).
    * The "cdr3_TRA" and "cdr3_TRB" columns should provide CDR3A/CDR3B sequences, following the standard definition (e.g., CAVNSDGQKLLF).
   Cases with non-amino acid characters, or length < 7 or > 22 will be not be considered (i.e., put to NA).
   * Other formats are supported (see below)


## Some optional parameters:

Below are some of the most comonly used parameters. Full documentation about other parameters is available in the R package.
 - check.cdr3.mode (either 0, 1, or 2(default): Defines the level of checking done at the beginning or end of the CDR3 sequences.
   0 means no check at all.
   1 means we check for the conservation of the number of amino acids given by the parameters start.lg and end.lg (default 1 and 2) with the V/J germline sequence.
   2 means that we will use the expected number of conserved residues derived by a statistical modelling of a large pool of TCRs to check the beginning and end of each CDR3 based on these values.
 - output.path: name of the output directory (if not already existing, it
   will be created). If existing the files with the same name will be overwritten.
   It can be left empty if the output is assigned to a variable (e.g., m <- MixTCRclean(input1="test/test.csv")).

- filename.output (default=NULL). Provide a name for the output file.
  If NULL, the filename will be the original name + "_processed".


# OUTPUT

If the output of MixTCRclean is assigned to a variable (e.g., m <- MixTCRclean(input1="test/test.csv")), MixTCRclean returns processed data.

If output.path is given, MixTCRclean creates also a directory (output.path). 

# Data format

By default, MixTCRclean uses column names c("TRAV","TRAJ","cdr3_TRA","TRBV","TRBJ","cdr3_TRB") to define a TCR,
"species" to indicate the species and "model" to define groups of TCRs (e.g., binding to the same epitope).
- For single-chain data, only one chain can be provided. In those cases, it is recommended to define the chain in chain="A" or "B".
- "species" can be skipped, in which case all TCRs are assumed to come from the species.default (default="HomoSapiens"). 
- Other column names are supported, including "Va", "V_alpha", "CDR3a", "CDR3A", "CDR3_alpha", etc. for data with both chain.
 Or "V", "v_gene","V-region","aaSeqCDR3","CDR3",etc for single chain data, see list in data_raw/TidyVJ/mapping_colnames.csv for a full description

Other supported formats treating each chain in a different row include:

 - VDJdb with the columns: c("V", "J", "CDR3")
 - 10X Genomics format with columns: c("v_gene", "j_gene", "cdr3")
 - Qiagen with the columns: c("V-region", "J-region", "CDR3 amino acid seq")
 - Adaptive Biotech with the columns: c("vGeneName", "jGeneName", "aminoAcid")
 - Adaptive Biotech v4 with the columns: c("v_resolved", "j_resolved", "amino_acid")
 - AIRR with the columns: c("v_call", "j_call", "junction_aa")
 - MiXCR with the columns: c("allVHitsWithScore", "allJHitsWithScore", "aaSeqCDR3")
  or c("allVGenes", "allJGenes", "aaSeqCDR3")

By default, both chains are treated independently.


# OTHER INFORMATION

* V/J genes are key to run the prediction in MixTCRclean and only V/J names compatible with the IMGT nomenclature can be considered. Even if correct.gene.names==1 allows to correct several wrong V/J names, we strongly encourage the users to use only V/J gene names compatible with IMGT.

* Some V/J genes in IMGT give rise to truncated V segments (e.g., TRAV8-5). All of them are pseudogenes. These are not supported in MixTCRclean and will be put to NA. Other pseudogenes / ORF are shown in grey in the plots.

* TRBV6-2 and TRBV6-3 have exactly the same nucleotide sequence, and therefore cannot be distinguished at the sequecing level. In MixTCRclean, these entries are mapped into a single 'TRBV6-2/6-3' gene.
