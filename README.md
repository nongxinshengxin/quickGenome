# quickGenome
## Overview
quickGenome is an R package for the convenient accounting of genomic information and extraction of gene sequences. quickGenome can be used to calculate genome size, N50, contig count, GC content and other information, to detect chromosomal telomere positions, to calculate the number of gene exons, to calculate gene length distribution, to extract gene, mRNA, CDS and protein sequences based on GFF3 files.
The implementation of `quickGenome` function relies on `Biostrings` and `GenomicFeatures`.

## Installation
Before installing, you will need to download the `quickGenome` dependency package `Biostrings` and `GenomicFeatures` by `BiocManager`.

```{r}
if (!require("BiocManager"))
  install.packages("BiocManager")
library(BiocManager)
if (!require("Biostrings"))
  BiocManager::install("Biostrings")
if (!require("GenomicFeatures"))
  BiocManager::install("GenomicFeatures")
```
Install `devtools`, which is used to install R packages from GitHub.
```{r}
if (!require("devtools"))
  install.packages("devtools")
```
Once you have completed the above steps, start the installation.
```{r}
devtools::install_github("nongxinshengxin/quickGenome")
```
## Function
`quickGenome` has ten **core functions**.
- **genome_basic_Info()**  Calculating basic genomic information, including size, GC content, N50 and so on.
- **find_telo()**  Identification of the location of telomeric repeat sequences at the ends of chromosomes.
- **binGC()**  The genomic sequence was split into multiple sliding windows of the same length and the GC content of each sliding window was calculated.
- **calculate_exonNum()** The number of gene exons was counted according to gff/gtf.
- **calculate_geneLength()**  The distribution of gene lengths was calculated from gff/gtf.
- **extract_gene()**  Extraction of gene sequences (including introns).
- **extract_mRNA()**  Extraction of mRNA sequences (exons only, not introns).
- **extract_CDS()**  Extraction of CDS sequences and protein sequences.
- **extract_upstream()**  Extraction of upstream sequences of genes.
- **extract_downstream()**  Extraction of downstream sequences of genes.

## Application


## Reference
Biostrings https://doi.org/doi:10.18129/B9.bioc.Biostrings

GenomicFeatures https://doi.org/doi:10.18129/B9.bioc.GenomicFeatures

## Documentation
The English documentation is available in - https://github.com/nongxinshengxin/quickGenome

The Chinese documentation is available in - 微信公众号农心生信工作室

## Citation
Please, when using `quickGenome`, cite us using the reference: https://github.com/nongxinshengxin/quickGenome

## Contact us
- Email: nongxinshengxin@163.com
- Wechat Official Account：
