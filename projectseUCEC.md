---
output: html_document
---
<!---
The following chunk of code, which should not be shown in the resulting document (echo=FALSE)
sets up global processing options, such as forcing 'knitr' to stop when an error
in the R code is encountered, caching of the results in the 'cache'
directory and asking 'knitr' to figure out automatically the dependencies among
code chunks to re-calculate cached results (autodep=TRUE).

Other options could be changing the name of the directory where figures end up
('figure' by default), etc. For a full account of 'knitr' options please consult
http://yihui.name/knitr/options

At the end of the chunk a 'cat()' call is made to dump a CSS file that gives
a better look-and-feel than the knitr default one. See the source css/ieo.css
and the resulting projectTemplate.html to understand where this is being dumpted.
--->

<style type="text/css">
/* Based on the bioconductor.css file from the BiocStyle package */

body, td {
   font-family: sans-serif;
   background-color: white;
   font-size: 13px;
}

body {
  max-width: 800px;
  margin: 0 auto;
  padding: 1em 1em 2em;
  line-height: 20px;
}

/* Table of contents style */

div#TOC li {
    list-style:none;
    background-image:none;
    background-repeat:none;
    background-position:0;
}

/* element spacing */

p, pre { 
  margin: 0em 0em 1em;
}

/* center images and tables */
img, table {
  margin: 0em auto 1em;
}

p {
  text-align: justify;
}

tt, code, pre {
   font-family: 'DejaVu Sans Mono', 'Droid Sans Mono', 'Lucida Console', Consolas, Monaco, monospace;
}

h2, h3, h4, h5, h6 { 
  font-family: Helvetica, Arial, sans-serif;
  margin: 1.2em 0em 0.6em 0em;
  font-weight: bold;
}

h1.title {
  font-size: 250%;
  font-weight: bold;
  color: #1a81c2;
  line-height: 1.4em;
  margin-top: 20px;
  border-bottom: 1px #1a81c2 solid;
}


h1 {
  font-size: 250%;
  font-weight: bold;
  line-height: 1.4em;
  margin-top: 20px;
  border-bottom: 1px #1a81c2 solid;
}

h2 {
  font-size: 160%;  
}

h1, h2, h3 {
  color: #1a81c2;
}


h3, h4, h5, h6 {
  font-size:115%;
} /* not expecting to dive deeper than four levels on a single page */

/* links are simply blue, hovering slightly less blue */
a { color: #1a81c2; text-decoration: none; }
a:active { outline: none; }
a:visited { color: #1a81c2; }
a:hover { color: #4c94c2; }

pre, img {
  max-width: 100%;
  display: block;
}

pre {
  border: 0px none;
  background-color: #F8F8F8;
  white-space: pre;
  overflow-x: auto;
}

pre code {
  border: 1px #aaa dashed;
  background-color: white;
  display: block;
  padding: 1em;  
  color: #111;
  overflow-x: inherit;
}

/* markdown v1 */
pre code[class] {
  background-color: inherit;
}

/* markdown v2 */
pre[class] code {
  background-color: inherit;
}

/* formatting of inline code */
code { 
  background-color: transparent;
  color: #87b13f;
  font-size: 92%;
}

/* formatting of tables */

table, td, th {
  border: none;
  padding: 0 0.5em;
}

/* alternating row colors */
tbody tr:nth-child(odd) td {
  background-color: #F8F8F8;
}

blockquote {
   color:#666666;
   margin:0;
   padding-left: 1em;
   border-left: 0.5em #EEE solid;
}

hr {
   height: 0px;
   border-bottom: none;
   border-top-width: thin;
   border-top-style: dotted;
   border-top-color: #999999;
}

span.header-section-number {
  padding-right: 1em;
}

span.toc-section-number::after {
    content: "  ";
    white-space: pre;
}

.caption {
  font-style: italic;
  text-align: center;
  margin-top: 2px;
  margin-bottom: 25px;
  font-size: 11pt;
}

@media print {
   * {
      background: transparent !important;
      color: black !important;
      filter:none !important;
      -ms-filter: none !important;
   }

   body {
      font-size:12pt;
      max-width:100%;
   }

   hr {
      visibility: hidden;
      page-break-before: always;
   }

   pre, blockquote {
      padding-right: 1em;
      page-break-inside: avoid;
   }

   tr, img {
      page-break-inside: avoid;
   }

   img {
      max-width: 100% !important;
   }

   @page :left {
      margin: 15mm 20mm 15mm 10mm;
   }

   @page :right {
      margin: 15mm 10mm 15mm 20mm;
   }

   p, h2, h3 {
      orphans: 3; widows: 3;
   }

   h2, h3 {
      page-break-after: avoid;
   }

   .caption {
     font-style: italic;
     text-align: center;
     margin-top: 2px;
     margin-bottom: 25px;
     font-size: 11pt;
   }

}
</style>

# Supplementary Information: 
# Analysis of a The Cancer Genome Atlas (TCGA) RNA-seq data set on Uterine Corpus Endometrial Carcinoma (UCEC)

### Joaquim Aguirre (joaquim.aguirre01@estudiant.upf.edu)
### Gerard Funosas (gfunosas64@gmail.com)
### Cristina Prat (cristina.prat.ferrer@gmail.com)

## Introduction

Endometrial cancer develops in the cells that form the inner lining of the uterus, or the endometrium, and is one of the most common cancers of the female reproductive system. In 2010, approximately 43,000 women in the United States were estimated to have been diagnosed and almost 8,000 to have died of endometrial cancer. This cancer occurs most commonly in women aged 60 years or older. About 69 percent of endometrial cancers are diagnosed at an early stage, and as a result about 83 percent of women will survive five years following the time of diagnosis.

TCGA researchers have: 
- Identified four subtypes of endometrial cancer: POLE ultramutated, Microsatellite instability hypermutated, Copy number low and Copy number high
- Uncovered shared genomic features between endometrial cancer and serous ovarian cancer, the Basal-like subtype of breast cancer as well as colorectal cancer
- Characterized the marked differences between the two types of endometrial tumors (endometrioid and serous), and found that some endometrioid tumors have developed a strikingly similar pattern to serous tumors, suggesting they may benefit from a common treatment 
    - The serous and some of the endometrioid tumors are characterized by frequent mutations in TP53, extensive copy number alterations and few DNA methylation changes
    - The rest of the endometrioid tumors are characterized by few copy number alterations, scarce mutations in TP53 and frequent mutations in PTEN and KRAS


This document should be processed from R and you need to install the packages [knitr](http://cran.r-project.org/web/packages/knitr/index.html) and [markdown](http://cran.r-project.org/web/packages/markdown/index.html). Once they are installed, you have to type the following instructions that generate a HTML document that you can open with a web browser:

```
# Load libraries to process the document
library(knitr)     ## required for "knitting" from Rmd to md
library(markdown)  ## required for processing from md to HTML
knit2html("projectseUCEC.Rmd", force_v1=TRUE)  ## process Rmd to HTML
browseURL("projectseUCEC.html") ## open the resulting HTML file from R
```


## Data import

The project starts importing the raw table of counts. It contains RNA-seq counts for 20115 genes and 589 samples.


```r
library(SummarizedExperiment)
```

```
Loading required package: GenomicRanges
```

```
Loading required package: BiocGenerics
```

```
Loading required package: parallel
```

```

Attaching package: 'BiocGenerics'
```

```
The following objects are masked from 'package:parallel':

    clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    clusterExport, clusterMap, parApply, parCapply, parLapply,
    parLapplyLB, parRapply, parSapply, parSapplyLB
```

```
The following objects are masked from 'package:stats':

    IQR, mad, xtabs
```

```
The following objects are masked from 'package:base':

    anyDuplicated, append, as.data.frame, cbind, colnames,
    do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    grepl, intersect, is.unsorted, lapply, lengths, Map, mapply,
    match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
    Position, rank, rbind, Reduce, rownames, sapply, setdiff,
    sort, table, tapply, union, unique, unsplit
```

```
Loading required package: S4Vectors
```

```
Loading required package: stats4
```

```

Attaching package: 'S4Vectors'
```

```
The following objects are masked from 'package:base':

    colMeans, colSums, expand.grid, rowMeans, rowSums
```

```
Loading required package: IRanges
```

```
Loading required package: GenomeInfoDb
```

```
Loading required package: Biobase
```

```
Welcome to Bioconductor

    Vignettes contain introductory material; view with
    'browseVignettes()'. To cite Bioconductor, see
    'citation("Biobase")', and for packages 'citation("pkgname")'.
```

```r
se <- readRDS(file.path("data", "seUCEC.rds"))
se
```

```
class: RangedSummarizedExperiment 
dim: 20115 589 
metadata(5): experimentData annotation cancerTypeCode
  cancerTypeDescription objectCreationDate
assays(1): counts
rownames(20115): 1 2 ... 102724473 103091865
rowData names(3): symbol txlen txgc
colnames(589): TCGA.2E.A9G8.01A.11R.A40A.07
  TCGA.4E.A92E.01A.11R.A37O.07 ... TCGA.FL.A1YV.11A.12R.A32Y.07
  TCGA.FL.A3WE.11A.11R.A22K.07
colData names(549): type bcr_patient_uuid ...
  lymph_nodes_aortic_pos_by_ihc lymph_nodes_aortic_pos_total
```

Using the following command it is possible to explore the phenotypic data information associated to the samples, which is associated in the R S4 object as metadata:


```r
dim(colData(se))
```

```
[1] 589 549
```

```r
colData(se)[1:5, 1:5]
```

```
DataFrame with 5 rows and 5 columns
                                 type                     bcr_patient_uuid
                             <factor>                             <factor>
TCGA.2E.A9G8.01A.11R.A40A.07    tumor 9583C10B-B21A-4863-98FA-61E735E64EA5
TCGA.4E.A92E.01A.11R.A37O.07    tumor                                   NA
TCGA.5B.A90C.01A.11R.A37O.07    tumor 16AC4341-CF8F-45E2-B90B-2D12D5F74A59
TCGA.5S.A9Q8.01A.11R.A40A.07    tumor 8751429B-4A11-451E-B978-DC9E9DB0EB36
TCGA.A5.A0G1.01A.11R.A118.07    tumor 53707bb3-426a-43cb-830f-3eeed930295f
                             bcr_patient_barcode form_completion_date
                                        <factor>             <factor>
TCGA.2E.A9G8.01A.11R.A40A.07        TCGA-2E-A9G8            2014-5-21
TCGA.4E.A92E.01A.11R.A37O.07                  NA                   NA
TCGA.5B.A90C.01A.11R.A37O.07        TCGA-5B-A90C            2014-4-17
TCGA.5S.A9Q8.01A.11R.A40A.07        TCGA-5S-A9Q8           2013-12-23
TCGA.A5.A0G1.01A.11R.A118.07        TCGA-A5-A0G1             2011-1-7
                             prospective_collection
                                           <factor>
TCGA.2E.A9G8.01A.11R.A40A.07                     NO
TCGA.4E.A92E.01A.11R.A37O.07                     NA
TCGA.5B.A90C.01A.11R.A37O.07                     NO
TCGA.5S.A9Q8.01A.11R.A40A.07                    YES
TCGA.A5.A0G1.01A.11R.A118.07                     NO
```

```r
mcols(colData(se), use.names=TRUE)
```

```
DataFrame with 549 rows and 2 columns
                                                         labelDescription
                                                              <character>
type                                           sample type (tumor/normal)
bcr_patient_uuid                                         bcr patient uuid
bcr_patient_barcode                                   bcr patient barcode
form_completion_date                                 form completion date
prospective_collection            tissue prospective collection indicator
...                                                                   ...
lymph_nodes_pelvic_pos_total                               total pelv lnp
lymph_nodes_aortic_examined_count                           total aor lnr
lymph_nodes_aortic_pos_by_he                          aln pos light micro
lymph_nodes_aortic_pos_by_ihc                                 aln pos ihc
lymph_nodes_aortic_pos_total                                total aor-lnp
                                        CDEID
                                  <character>
type                                       NA
bcr_patient_uuid                           NA
bcr_patient_barcode                   2673794
form_completion_date                       NA
prospective_collection                3088492
...                                       ...
lymph_nodes_pelvic_pos_total          3151828
lymph_nodes_aortic_examined_count     3104460
lymph_nodes_aortic_pos_by_he          3151832
lymph_nodes_aortic_pos_by_ihc         3151831
lymph_nodes_aortic_pos_total          3151827
```

These metadata provides with useful information about the clinical variables of the samples. The 'CDEID' column corresponds to the so-called `Common Data Element (CDE)` identifier, which can be used in https://cdebrowser.nci.nih.gov to search for further information about the associated clinical variable using the `Advanced search` form and the `Public ID` attribute search.

Associated to the row (feature) data, there are 455 sequences (1 circular) from hg38 genome:


```r
rowRanges(se)
```

```
GRanges object with 20115 ranges and 3 metadata columns:
            seqnames               ranges strand |      symbol     txlen
               <Rle>            <IRanges>  <Rle> | <character> <integer>
          1    chr19 [58345178, 58362751]      - |        A1BG      3322
          2    chr12 [ 9067664,  9116229]      - |         A2M      4844
          9     chr8 [18170477, 18223689]      + |        NAT1      2280
         10     chr8 [18391245, 18401218]      + |        NAT2      1322
         12    chr14 [94592058, 94624646]      + |    SERPINA3      3067
        ...      ...                  ...    ... .         ...       ...
  100996331    chr15 [20835372, 21877298]      - |       POTEB      1706
  101340251    chr17 [40027542, 40027645]      - |    SNORD124       104
  101340252     chr9 [33934296, 33934376]      - |   SNORD121B        81
  102724473     chrX [49303669, 49319844]      + |      GAGE10       538
  103091865    chr21 [39313935, 39314962]      + |   BRWD1-IT2      1028
                 txgc
            <numeric>
          1 0.5644190
          2 0.4882329
          9 0.3942982
         10 0.3895613
         12 0.5249429
        ...       ...
  100996331 0.4308324
  101340251 0.4903846
  101340252 0.4074074
  102724473 0.5055762
  103091865 0.5924125
  -------
  seqinfo: 455 sequences (1 circular) from hg38 genome
```

From the S4 object, it is possible to extract information about the gender of the patients who donated the samples. As the study is focused on endometrial cancer, all the samples are from female patients. There are also 33 'NA' samples which were considered to be discarded, but finally they have been mantained as they provide the project with some normal samples, which are not abundant in the dataset.


```r
table(se$gender)
```

```

FEMALE 
   556 
```


## Quality assessment and normalization

The fact that each RNA-seq sample may have been ultimately sequenced at slightly different depth and that there may be sample-specific biase implies that it may be needed to consider two normalization steps:

-  Within-sample: adjustments to compare across features in a sample. 
    - Scaling: using counts per million reads (CPM mapped to the genome) 
-  Between-sample: adjustments to compare a feature across samples.
    - Sample-specific normalization factors: using the TMM algorithm (package edgeR)
    - Quantile normalization: using the CQN algorithm (package cqn)

To perform quality assessment and normalization it is necessary to load the [edgeR](http://bioconductor.org/packages/edgeR) R/Bioconductor package and create a `DGEList' object. We adjust for sample- and gene-specific factors, to make gene expression values comparable across samples. 


```r
library(edgeR)
```

```
Loading required package: limma
```

```

Attaching package: 'limma'
```

```
The following object is masked from 'package:BiocGenerics':

    plotMA
```

```r
dge <- DGEList(counts=assays(se)$counts, genes=mcols(se))
```

```
Warning in as.data.frame(x, row.names = NULL, optional = optional, ...):
Arguments in '...' ignored
```

```r
names(dge)
```

```
[1] "counts"  "samples" "genes"  
```


### CPM scaling

Now that the 'DGEList' object has been created, it is possible to perform the scaling to CPM values. Therefore, $\log_2$ CPM values of expression are calculated and used as an additional assay element to ease their manipulation. $\log_2$ CPM units separate better high and low expression, than raw counts or non-logged CPM units.


```r
assays(se)$logCPM <- cpm(dge, log=TRUE, prior.count=0.5)
assays(se)$logCPM[1:5, 1:5]
```

```
   TCGA.2E.A9G8.01A.11R.A40A.07 TCGA.4E.A92E.01A.11R.A37O.07
1                     0.4787484                   0.08427463
2                     8.6045526                   9.69917999
9                    -5.6232476                  -5.62324755
10                   -5.6232476                  -5.62324755
12                    5.2785238                   5.65345688
   TCGA.5B.A90C.01A.11R.A37O.07 TCGA.5S.A9Q8.01A.11R.A40A.07
1                     -1.299515                     1.486274
2                      7.998470                     8.158385
9                     -5.623248                    -5.623248
10                    -5.623248                    -5.623248
12                     4.835107                     5.988568
   TCGA.A5.A0G1.01A.11R.A118.07
1                      2.005532
2                      8.518400
9                     -5.623248
10                    -5.623248
12                     5.593132
```


### Library sizes

It is primordial to examine the library sizes in terms of total number of sequence read counts per sample. Figure S1.1 below shows library sizes per sample in increasing order:

<!---
you can control the height and width in pixels of the figure with 'out.height' and 'out.width'
--->

<img src="figure/libsizes1-1.png" title="Figure S1.1: Library sizes in increasing order." alt="Figure S1.1: Library sizes in increasing order." width="600px" style="display: block; margin: auto;" /><p class="caption">Figure S1.1: Library sizes in increasing order.</p>

As there is a large amount of samples, it is difficult to distinguish between tumor and normal samples in this plot. In further steps, it may be necessary to work with a subset in order to obtain clearer results.
However, the figure S1.1 reveals substantial differences in sequencing depth between samples and it may be considered discarding those samples whose depth is substantially lower than the rest. 


```r
sampledepth <- round(dge$sample$lib.size / 1e6, digits=1)
sampledepth <- sort(sampledepth)
sampledepth[1:50]
```

```
 [1]  3.3  5.4  6.1  6.1  6.1  6.2  6.8  7.0  7.1  7.7  8.2  8.7  8.8  9.0
[15]  9.0  9.8 10.0 10.0 10.1 10.2 10.3 10.3 10.4 10.7 10.9 10.9 10.9 11.0
[29] 11.3 11.3 11.5 11.5 11.8 11.8 12.3 12.5 12.7 12.8 12.8 12.9 12.9 13.4
[43] 13.4 13.7 13.8 13.9 14.0 14.0 14.2 14.2
```

It has been considered to discard those samples corresponding to the 10% quartile of the sample depth distribution, as the quality of the sequentiation of these samples is poorer.


```r
mask <- round(dge$sample$lib.size / 1e6, digits=1) > quantile(sampledepth, .1)
dge <- dge[ , mask]
dim(dge)
```

```
[1] 20115   527
```

```r
se <- se[ , mask]
dim(se)
```

```
[1] 20115   527
```

```r
sampledepth <- round(dge$sample$lib.size / 1e6, digits=1)
sampledepth <- sort(sampledepth)
sampledepth[1:50]
```

```
 [1] 14.7 14.7 14.8 14.8 14.8 14.8 14.9 14.9 14.9 14.9 15.0 15.0 15.0 15.0
[15] 15.0 15.0 15.0 15.1 15.1 15.2 15.2 15.2 15.2 15.2 15.2 15.4 15.4 15.5
[29] 15.5 15.5 15.6 15.6 15.6 15.6 15.6 15.6 15.7 15.8 15.8 15.8 15.9 16.0
[43] 16.0 16.1 16.1 16.1 16.1 16.2 16.2 16.3
```

The filtered set has now 527 samples. Before, there was a range of sample depth from 3.3 to 60.1 millions of reads, and now the range starts at 14.7 million reads. 


### Subsetting

As it has been stated previously, a subsetting step is necessary in order to manage easily the set of samples and obtain clearer results. But it is important to work with a subset which is as much representative as the initial set of samples and that contains the samples with higher quality. Three types of subsetting have been considered:
- Random subsetting: As there is a small number of normal samples, this subsetting takes all normal samples, and  then the same number of tumor samples randomly
- Paired subsetting: It uses only normal and tumor samples which are paired together, so that they are from the same patient
- Independent subsetting: It uses only normal and tumor samples which are independent, so that they are not from the same patient

The three options have been considered and tested. The paired subsetting offers the advantage that as samples are paired, the posterior analysis of batch effect identification will be performed with a perfectly balanced set, which avoids confusions for not having samples of one of the variables. However, in this dataset there are only 36 paired samples, which is a very small subset of samples.

The other two options permit to work with a higher number of samples, but there are many problems related with having an unbalanced set.

Finally, it has been preferred to work with the paired subset of samples in order to perform a clearer analysis of possible batch effects. The code used to perform random or indepedent subsetting has been left for possible future use.

#### Random subsetting


```r
#mask <- se$type == "normal"
#tumor <- se[, !mask]
#subtumor <- tumor[, sample(1:ncol(tumor), length(mask[ which(mask==TRUE)]))]
#se <- cbind(subtumor, se[, mask])
#se
#dge <- DGEList(counts=assays(se)$counts, genes=mcols(se))
```

#### Independent subsetting


```r
# mask <- se$type == "normal"
# tumor <- se[, !mask]
# normal <- se[, mask]
# participant <- substr(colnames(se), 9, 12)
# paired <- names(which(table(participant) == 2))
# length(paired)

# mask <- substr(colnames(se), 9, 12) %in% paired
# paired_samples <- se[, mask]
# mask <- substr(colnames(paired_samples), 9, 12) %in% substr(colnames(normal), 9, 12)
# paired_samples <- paired_samples[, mask]
# tumor_paired <- paired_samples[ , paired_samples$type == "tumor"]
# normal_paired <- paired_samples[ , paired_samples$type =="normal"]
# mask <- colnames(se) %in% colnames(tumor_paired)
# independent <- se[ , !mask]

# se <- independent
# dge <- DGEList(counts=assays(independent)$counts, genes=mcols(independent))
```

#### Paired subsetting


```r
mask <- se$type == "normal"
tumor <- se[, !mask]
normal <- se[, mask]
participant <- substr(colnames(se), 9, 12)
paired <- names(which(table(participant) == 2))
mask <- substr(colnames(se), 9, 12) %in% paired
paired_samples <- se[, mask]
mask <- substr(colnames(paired_samples), 9, 12) %in% substr(colnames(normal), 9, 12)
paired_samples <- paired_samples[, mask]
dim(paired_samples)
```

```
[1] 20115    36
```

```r
colnames(paired_samples)
```

```
 [1] "TCGA.AJ.A3NC.01A.11R.A22K.07" "TCGA.AJ.A3NE.01A.11R.A22K.07"
 [3] "TCGA.AJ.A3NH.01A.11R.A22K.07" "TCGA.AX.A05Y.01A.11R.A00V.07"
 [5] "TCGA.AX.A0IZ.01A.11R.A118.07" "TCGA.AX.A0J0.01A.11R.A109.07"
 [7] "TCGA.AX.A1CF.01A.11R.A137.07" "TCGA.AX.A2H8.01A.11R.A17B.07"
 [9] "TCGA.AX.A2HC.01A.11R.A17B.07" "TCGA.AX.A2HD.01A.21R.A17B.07"
[11] "TCGA.BG.A2AD.01A.21R.A16F.07" "TCGA.BG.A3EW.01A.11R.A22K.07"
[13] "TCGA.BG.A3PP.01A.11R.A22K.07" "TCGA.BK.A0CB.01A.32R.A104.07"
[15] "TCGA.BK.A13C.01A.11R.A118.07" "TCGA.BK.A4ZD.01A.11R.A27V.07"
[17] "TCGA.DI.A2QY.01A.12R.A19W.07" "TCGA.E6.A1M0.01A.11R.A144.07"
[19] "TCGA.AJ.A3NC.11A.11R.A22K.07" "TCGA.AJ.A3NE.11A.11R.A22K.07"
[21] "TCGA.AJ.A3NH.11A.11R.A22K.07" "TCGA.AX.A05Y.11A.11R.A27V.07"
[23] "TCGA.AX.A0IZ.11A.11R.A27V.07" "TCGA.AX.A0J0.11A.11R.A27V.07"
[25] "TCGA.AX.A1CF.11A.11R.A137.07" "TCGA.AX.A2H8.11A.11R.A17B.07"
[27] "TCGA.AX.A2HC.11A.11R.A32Y.07" "TCGA.AX.A2HD.11A.11R.A32Y.07"
[29] "TCGA.BG.A2AD.11A.11R.A16F.07" "TCGA.BG.A3EW.11A.22R.A22K.07"
[31] "TCGA.BG.A3PP.11A.11R.A22K.07" "TCGA.BK.A0CB.11A.33R.A104.07"
[33] "TCGA.BK.A13C.11A.11R.A118.07" "TCGA.BK.A4ZD.11A.12R.A27V.07"
[35] "TCGA.DI.A2QY.11A.11R.A19W.07" "TCGA.E6.A1M0.11A.11R.A144.07"
```

```r
table(paired_samples$type)
```

```

normal  tumor 
    18     18 
```

```r
se <- paired_samples
dge <- DGEList(counts=assays(paired_samples)$counts, genes=mcols(paired_samples))
```

```
Warning in as.data.frame(x, row.names = NULL, optional = optional, ...):
Arguments in '...' ignored
```


The plot showing the filtered library sizes after the random subsetting is represented in figure S1.2:

<img src="figure/libsizes2-1.png" title="Figure S1.2: Filtered library sizes in increasing order." alt="Figure S1.2: Filtered library sizes in increasing order." width="600px" style="display: block; margin: auto;" /><p class="caption">Figure S1.2: Filtered library sizes in increasing order.</p>


### Distribution of expression levels among samples

The following plot (Figure S2) shows the distribution of expression values per sample in terms of logarithmic CPM units. Due to the large number of samples, tumor and normal samples are displayed separately.

<!---
the option echo=FALSE hides the R code. When plotting in general one does not
want to see the code. Options fig.height and fig.width control height and width
of the plot in inches while out.height and out.width do it in the final output
file; see http://yihui.name/knitr/options for full details.
--->


```
Loading required package: lattice
```

```
Loading required package: annotate
```

```
Loading required package: AnnotationDbi
```

```
Loading required package: XML
```

<img src="figure/distRawExp-1.png" title="Figure S2: Non-parametric density distribution of expression profiles per sample." alt="Figure S2: Non-parametric density distribution of expression profiles per sample." width="800px" style="display: block; margin: auto;" /><p class="caption">Figure S2: Non-parametric density distribution of expression profiles per sample.</p>

There are no substantial differences appreciated between the samples in the distribution of expression values.


### Distribution of expression levels among genes

The average expression per gene through all the samples is calculated. $\log_2$ CPM values per gene are examined. Figure S3 shows the distribution of those values across genes.

<img src="figure/exprdist-1.png" title="Figure S3: Distribution of average expression level per gene." alt="Figure S3: Distribution of average expression level per gene." width="400px" style="display: block; margin: auto;" /><p class="caption">Figure S3: Distribution of average expression level per gene.</p>


### Filtering of lowly-expressed genes

In the light of the plot in figure S3, it may be considers a cutoff of 1 log CPM unit (as represented by the red line) as minimum value of expression to select genes being expressed across samples. Using this cutoff we proceed to filter out lowly-expressed genes.


```r
mask <- avgexp > 1
se <- se[mask, ]
dim(se)
```

```
[1] 11571    36
```

```r
dge <- dge[mask, ]
dim(dge)
```

```
[1] 11571    36
```

### Normalization

The normalization factors are calculated on the filtered expression data set. They have been calculated using the Trimmed Mean of M-values (TMM) method which addresses the issue of the different RNA composition of the samples by estimating a scaling factor for each library.


```r
dge <- calcNormFactors(dge)
```

The raw log2 CPM units have been replaced in the corresponding assay element of the `SummarizedExperiment` object, by the normalized ones.


```r
assays(se)$logCPM <- cpm(dge, log=TRUE, prior.count=0.5)
```


### MA-plots

The MA-plots of the normalized expression profiles are performed. First, the plots (Figure S4) corresponding to tumor samples are built:

<!---
Here we make a MA-plot for each sample. The options 'fig.height' and 'fig.width'
control the relative image size in *inches*. The final image size results from
'height'x'dpi' and 'width'x'dpi', where 'dpi' is the image resolution in
"dots per inch" (by default dpi=72). To scale the image to a desired size use
'out.width' and 'out.height'. More information at http://yihui.name/knitr/options
--->

<img src="figure/maPlotsTumor-1.png" title="Figure S4: MA-plots of the tumor samples." alt="Figure S4: MA-plots of the tumor samples." style="display: block; margin: auto;" /><p class="caption">Figure S4: MA-plots of the tumor samples.</p>

In general, we do not observe samples with major expression-level dependent biases, although some of them show variations in low-expressed values.

Now, the plots (Figure S5) of the normal samples are performed:

<img src="figure/maPlotsNormal-1.png" title="Figure S5: MA-plots of the normal samples." alt="Figure S5: MA-plots of the normal samples." style="display: block; margin: auto;" /><p class="caption">Figure S5: MA-plots of the normal samples.</p>

In this case, most of the normals have correct MA plots, but for some of them we see slightly expression-level dependent biases. The most suspicious cases are TCGA-AJ-A3NH, TCGA-AX-A2HC, TCGA-BK-A13C and TCGA-DI-A2QY, showing sizable dependency between M and A values. 
We should consider discarding those samples from the dataset if they present further signs of problematic features.  


### Batch identification

This step will be focused on finding potential surrogate of batch effect indicators. Given that each sample names corresponds to a TCGA barcode (see https://wiki.nci.nih.gov/display/TCGA/TCGA+barcode), following the strategy described in http://bioinformatics.mdanderson.org/main/TCGABatchEffects:Overview we are going to derive different elements of the TCGA barcode and examine their distribution across samples.

From this information we can make the following observations:

  * Samples were collected across different tissue source sites (TSS):
  

```r
tss <- substr(colnames(se), 6, 7)
table(tss)
```

```
tss
AJ AX BG BK DI E6 
 6 14  6  6  2  2 
```

  * All tumor samples belong to the same type (Code: 01, primary solid Tumor), while all normal samples belong to the "solid tissue normal" type (Code: 11):


```r
samplevial <- substr(colnames(se), 14, 16)
table(samplevial)
```

```
samplevial
01A 11A 
 18  18 
```

  * The majority of the samples were sequenced using the same analyte:


```r
portionanalyte <- substr(colnames(se), 18, 20)
table(portionanalyte)
```

```
portionanalyte
11R 12R 21R 22R 32R 33R 
 29   2   2   1   1   1 
```
  
  * Samples were sequenced within different plates:
  

```r
plate <- substr(colnames(se), 22, 25)
table(plate)
```

```
plate
A00V A104 A109 A118 A137 A144 A16F A17B A19W A22K A27V A32Y 
   1    2    1    3    2    2    2    4    2   10    5    2 
```

  * All samples were sequenced at the same center 07:
  

```r
center <- substr(colnames(se), 27, 28)
table(center)
```

```
center
07 
36 
```

We are going to use the TSS as surrogate of batch effect indicator. Considering our outcome of interest as molecular changes between sample types, tumor vs. normal, we will examine now the cross-classification of this outcome with TSS.


```r
table(data.frame(TYPE=se$type, TSS=tss))
```

```
        TSS
TYPE     AJ AX BG BK DI E6
  normal  3  7  3  3  1  1
  tumor   3  7  3  3  1  1
```

The experimental design using paired samples has the advantage which is perfectly balanced, as individuals from every outcome occur the same number of times through every batch.

We examine now how samples group together by hierarchical clustering and multidimensional scaling, annotating the outcome of interest and the the surrogate of batch indicator. 

We calculate again log CPM values with a higher prior count to moderate extreme fold-changes produced by low counts. 


```r
logCPM <- cpm(dge, log=TRUE, prior.count=3)
```

Next, calculate the distance between every pair of samples using a non-parametric association measure such as Spearman correlation:


```r
d <- as.dist(1-cor(logCPM, method="spearman"))
```

Then, perform a hierarchical clustering of the samples:


```r
sampleClustering <- hclust(d)
```

Now define the batch indicator variable as factor(tss), and extract the corresponding dendrogram annotating batch and outcome on its leafs. The resulting dendrogram is shown in Figure S6.


```r
par(mfrow=c(1, 1))
batch <- as.integer(factor(tss))
sampleDendrogram <- as.dendrogram(sampleClustering, hang=0.1)
names(batch) <- colnames(se)
outcome <- paste(substr(colnames(se), 9, 12), as.character(se$type), sep="-")
names(outcome) <- colnames(se)
sampleDendrogram <- dendrapply(sampleDendrogram,
                               function(x, batch, labels) {
                                 if (is.leaf(x)) {
                                   attr(x, "nodePar") <- list(lab.col=as.vector(batch[attr(x, "label")]))
                                   attr(x, "label") <- as.vector(labels[attr(x, "label")])
                                 }
                                 x
                               }, batch, outcome)
plot(sampleDendrogram, main="Hierarchical clustering of samples")
legend("topright", paste("Batch", sort(unique(batch)), levels(factor(tss))), fill=sort(unique(batch)))
```

<img src="figure/sampleClustering-1.png" title="Figure S6: Hierarchical clustering of the samples." alt="Figure S6: Hierarchical clustering of the samples." style="display: block; margin: auto;" /><p class="caption">Figure S6: Hierarchical clustering of the samples.</p>

In Figure S6, it can be observed that samples cluster primarily by sample type, tumor or normal. In general, there is a heterogeneous distribution of the different batches among the plot, and therefore there is no batch effect observed.


```r
plotMDS(dge, labels=outcome, col=batch)
legend("bottomright", paste("Batch", sort(unique(batch)), levels(factor(tss))),
       fill=sort(unique(batch)))
```

<img src="figure/mds-1.png" title="Figure S7: Multidimensional scaling plot of the samples." alt="Figure S7: Multidimensional scaling plot of the samples." style="display: block; margin: auto;" /><p class="caption">Figure S7: Multidimensional scaling plot of the samples.</p>

In Figure S7 we show the corresponding multidimensional plot (MDS). Here it can be seen more clearly that the first source of variation separates tumor from normal samples.
It can be observed that one tumor sample, corresponding to individual `A2HC-tumor` is separated from the rest, just as it happens in the hierchical clustering. A closer examination of its corresponding MA-plot also reveals a slight dependence of expression changes on average expression, overall in its paired normal sample. This turns to be one of the problematic samples we found previously, so at that point that pair of samples should be discarded to avoid undesired variation. 
Another similar case is the sample A2QY-normal, very clustered away within the normal group in the MDS plot, and also stated before as a problematic sample in the MA plot. For this reason, the A2QY samples will be removed.


```r
mask <- -grep("TCGA.AX.A2HC", colnames(se))
se <- se[,mask]
dge <- dge[,mask]
mask <- -grep("TCGA.DI.A2QY", colnames(se))
se <- se[,mask]
dge <- dge[,mask]
```

With 4 paired samples removed, we need to recalculate the variables and make the hierarchical clustering shown in Figure S8:


```r
tss <- substr(colnames(se), 6, 7)
logCPM <- cpm(dge, log=TRUE, prior.count=3)
```


```r
par(mfrow=c(1, 1))
batch <- as.integer(factor(tss))
sampleDendrogram <- as.dendrogram(sampleClustering, hang=0.1)
names(batch) <- colnames(se)
outcome <- paste(substr(colnames(se), 9, 12), as.character(se$type), sep="-")
names(outcome) <- colnames(se)
sampleDendrogram <- dendrapply(sampleDendrogram,
                               function(x, batch, labels) {
                                 if (is.leaf(x)) {
                                   attr(x, "nodePar") <- list(lab.col=as.vector(batch[attr(x, "label")]))
                                   attr(x, "label") <- as.vector(labels[attr(x, "label")])
                                 }
                                 x
                               }, batch, outcome)
plot(sampleDendrogram, main="Hierarchical clustering of samples")
legend("topright", paste("Batch", sort(unique(batch)), levels(factor(tss))), fill=sort(unique(batch)))
```

<img src="figure/sampleClustering_2-1.png" title="Figure S8: Hierarchical clustering of the samples." alt="Figure S8: Hierarchical clustering of the samples." style="display: block; margin: auto;" /><p class="caption">Figure S8: Hierarchical clustering of the samples.</p>

And get the new MDS plot as it shown in Figure S9: 


```r
plotMDS(dge, labels=outcome, col=batch)
legend("bottomright", paste("Batch", sort(unique(batch)), levels(factor(tss))),
       fill=sort(unique(batch)))
```

<img src="figure/mds3-1.png" title="Figure S9: Multidimensional scaling plot of the samples." alt="Figure S9: Multidimensional scaling plot of the samples." style="display: block; margin: auto;" /><p class="caption">Figure S9: Multidimensional scaling plot of the samples.</p>

### Removing batch effect: ComBat

One of these techniques to remove batch effect is ComBat (http://www.ncbi.nlm.nih.gov/pubmed/16632515), which is an empirical Bayes method robust to outliers in small sample sizes.
The sva (http://www.bioconductor.org/packages/release/bioc/html/sva.html) package provides a function called ComBat() to use this technique as follows:


```r
library(sva)
```

```
Loading required package: mgcv
```

```
Loading required package: nlme
```

```

Attaching package: 'nlme'
```

```
The following object is masked from 'package:IRanges':

    collapse
```

```
This is mgcv 1.8-12. For overview type 'help("mgcv-package")'.
```

```
Loading required package: genefilter
```

```r
mod <- model.matrix(~ se$type, colData(se))
combatexp <- ComBat(logCPM, batch, mod)
```

```
Found 5 batches
Adjusting for 1 covariate(s) or covariate level(s)
Standardizing Data across genes
Fitting L/S model and finding priors
Finding parametric adjustments
Adjusting the Data
```

```r
d <- as.dist(1-cor(combatexp, method="spearman"))

sampleClustering <- hclust(d)
```

Let's verify the extent to which batch effect has been removed with the hierarchical clustering of the samples:


```r
par(mfrow=c(1, 1))
batch <- as.integer(factor(tss))
sampleDendrogram <- as.dendrogram(sampleClustering, hang=0.1)
names(batch) <- colnames(se)
outcome <- paste(substr(colnames(se), 9, 12), as.character(se$type), sep="-")
names(outcome) <- colnames(se)
sampleDendrogram <- dendrapply(sampleDendrogram,
                               function(x, batch, labels) {
                                 if (is.leaf(x)) {
                                   attr(x, "nodePar") <- list(lab.col=as.vector(batch[attr(x, "label")]))
                                   attr(x, "label") <- as.vector(labels[attr(x, "label")])
                                 }
                                 x
                               }, batch, outcome)
plot(sampleDendrogram, main="Hierarchical clustering of samples")
legend("topright", paste("Batch", sort(unique(batch)), levels(factor(tss))), fill=sort(unique(batch)))
```

<img src="figure/sampleClustering_3-1.png" title="Figure S10: Hierarchical clustering of the samples." alt="Figure S10: Hierarchical clustering of the samples." style="display: block; margin: auto;" /><p class="caption">Figure S10: Hierarchical clustering of the samples.</p>

The plot in Figure S10 shows a better stratification of the tumor and normal samples than in the last plot (Figure S6) without removing the batch effect. 


## Differential expression: Simple analysis

We perform a simple examination of expression changes and their associated p-values using the R/Bioconductor package [sva](http://bioconductor.org/packages/sva).

Surrogate variable analysis (SVA) is a technique that tries to capture sources of heterogenity in high-throughput profiling data, such as non-biological variability introduced by batch effects.
The output of SVA is an estimation of the number of so-called "surrogate variables" and their continuous values, which can be used later on to adjust for these unmeasured and unwanted effects.


```r
library(sva)
mod <- model.matrix(~ se$type, colData(se))
mod0 <- model.matrix(~ 1, colData(se))
pv <- f.pvalue(assays(se)$logCPM, mod, mod0)
sum(p.adjust(pv, method="fdr") < 0.01)
```

```
[1] 5093
```

There are 5093 genes changing significantly their expression at FDR < 1%. In Figure S11 below we show the distribution of the resulting p-values.

<img src="figure/pdist-1.png" title="Figure S11: Distribution of raw p-values for an F-test on every gene between tumor and normal samples." alt="Figure S11: Distribution of raw p-values for an F-test on every gene between tumor and normal samples." width="500px" style="display: block; margin: auto;" /><p class="caption">Figure S11: Distribution of raw p-values for an F-test on every gene between tumor and normal samples.</p>

Now, let's estimate surrogate variables using the `sva()` function.


```r
sv <- sva(assays(se)$logCPM, mod, mod0)
```

```
Number of significant surrogate variables is:  10 
Iteration (out of 5 ):1  2  3  4  5  
```

```r
sv$n
```

```
[1] 10
```

The SVA algorithm has found 10 surrogate variables. Let's use them to assess again the extent of differential expression this time adjusting for these surrogate variables.


```r
modsv <- cbind(mod, sv$sv)
mod0sv <- cbind(mod0, sv$sv)
pvsv <- f.pvalue(assays(se)$logCPM, modsv, mod0sv)
sum(p.adjust(pvsv, method="fdr") < 0.01)
```

```
[1] 6722
```

We have increased the number of changing genes to 6722.
Figure S12 shows the resulting distribution of p-values.

<img src="figure/psvdist-1.png" title="Figure S12: Distribution of raw p-values for an F-test on every gene between tumor and normal samples, adjusting for surrogate variables estimated with SVA." alt="Figure S12: Distribution of raw p-values for an F-test on every gene between tumor and normal samples, adjusting for surrogate variables estimated with SVA." width="500px" style="display: block; margin: auto;" /><p class="caption">Figure S12: Distribution of raw p-values for an F-test on every gene between tumor and normal samples, adjusting for surrogate variables estimated with SVA.</p>


## Differential expression: Extended analysis

In this section, a more extense analysis of differential expression is performed. Different types of linear regression models are built in order to assess differential expression.
The conceptual purpose of a linear regression model is to represent, as accurately as possible, something complex, the data denoted by `y`, which is n-dimensional, in terms of something much simpler, the `model`, which is p -dimensional.
Thus, if our model is successful, the structure in the data should be captured in those `p` dimensions, leaving just random variation in the residuals which lie in an `(n-p)-dimensional space`.

In the context of DE analysis, linear regression models can be written in matrix form, in `design matrices`. The design matrix contains as many rows as samples and as many columns as coefficients to be estimated. 

A typical workflow to calculate DE analysis with [limma](http://www.bioconductor.org/packages/release/bioc/html/limma.html)
consists of the following steps:

1. Build design matrix: `model.matrix()`

2. Calculate observational-level weights (not necessary for microarray data): `voom()`

3. Estimate correlation in repeated measurements of a blocking factor (if needed): `duplicateCorrelation()`

4. Fit linear model: `lmFit()`

5. Build contrast matrix (if needed): `makeContrasts()`

6. Fit contrasts (if needed): `contrasts.fit()`

7. Calculate moderated $t$-statistics: `eBayes()`

8. Output results: `topTable()`

Now, we are going to observe the result of different approaches to perform DE analysis.

### 1. Fit directly the linear model to every gene

We create the desing matrix using the function `model.matrix()` as follows:


```r
design <- model.matrix(~factor(type), data = colData(se))
head(design)
```

```
                             (Intercept) factor(type)tumor
TCGA.AJ.A3NC.01A.11R.A22K.07           1                 1
TCGA.AJ.A3NE.01A.11R.A22K.07           1                 1
TCGA.AJ.A3NH.01A.11R.A22K.07           1                 1
TCGA.AX.A05Y.01A.11R.A00V.07           1                 1
TCGA.AX.A0IZ.01A.11R.A118.07           1                 1
TCGA.AX.A0J0.01A.11R.A109.07           1                 1
```

We are going to work as if the logCPM values were microarray expression data and fit directly the linear model specified in the design matrix to every gene using `lmFit()`:


```r
fit <- lmFit(assays(se)$logCPM, design)
```

We calculate moderated $t$-statistics using `eBayes`:


```r
fit <- eBayes(fit)
```

A quick overview of the results can be obtained with `decideTests()`:


```r
# Matrix of boolean: DE or not (criteria: p-value < 0.05)
FDRcutoff <- 0.05
res <- decideTests(fit, p.value = FDRcutoff)
summary(res) 
```

```
   (Intercept) factor(type)tumor
-1          30              3249
0          260              5114
1        11281              3208
```

It can be observed that:
- 3249 genes are underexpressed in tumor samples
- 3208 genes are overexpressed in tumor samples
- 5114 genes are not differentially expressed

To obtain a full table of all results we should use the function `topTable()` as follows:


```r
tt <- topTable(fit, coef = 2, n = Inf) # coef = 2 is for the second column. The colname can also be used. 
                                       # n = Inf is for obtaining all DE genes. By default we would only see 10. 
FDRcutoff <- 0.05
table(tt$adj.P.Val < FDRcutoff) # 6407 genes are differentially expressed
```

```

FALSE  TRUE 
 5114  6457 
```

```r
# B is a logODDS. Probability of being DE / probability of not being DE
# LogFC is the coefficient
# rownames of the table correspond to feature identifiers, Entrez Gene IDs in this case. We can add further gene metadata with the following command:
genesmd <- data.frame(chr = as.character(seqnames(rowRanges(se))), symbol = rowData(se)[,1], stringsAsFactors = FALSE)
fit$genes <- genesmd
tt <- topTable(fit, coef = 2, n = Inf)
head(tt)
```

```
            chr  symbol      logFC  AveExpr         t      P.Value
100423021 chr11 MIR4298  -7.770322 2.915817 -19.83072 8.213668e-20
136647     chr7  MPLKIP -10.924134 3.807787 -19.40805 1.576050e-19
125        chr4   ADH1B -10.196670 2.839204 -16.34769 2.601714e-17
10844     chr10 TUBGCP2   2.407317 3.986429  16.32964 2.687202e-17
81626      chr1 SHCBP1L   2.287176 7.656542  15.95795 5.262783e-17
144348    chr12  ZNF664   5.587508 3.240038  15.82544 6.707689e-17
             adj.P.Val        B
100423021 9.118237e-16 34.98370
136647    9.118237e-16 34.36003
125       7.773403e-14 29.42196
10844     7.773403e-14 29.39044
81626     1.073157e-13 28.73444
144348    1.073157e-13 28.49739
```

```r
# Fetch gene identifiers of DE genes with FDR < 5%
DEgenes <- rownames(tt)[tt$adj.P.Val < FDRcutoff]
DEgenes_symbol <- tt$symbol[tt$adj.P.Val < FDRcutoff]
length(DEgenes)
```

```
[1] 6457
```

```r
# Examine the chormosome distribution of genes called DE at 5% of FDR
sort(table(tt$chr[tt$adj.P.Val < FDRcutoff]), decreasing = TRUE)
```

```

                chr1                chr19                chr11 
                 669                  462                  447 
                chr2                chr17                 chr3 
                 439                  368                  364 
               chr12                 chr7                 chr6 
                 339                  314                  285 
                chr5                chr16                 chr9 
                 279                  263                  261 
                chrX                chr14                chr10 
                 243                  238                  230 
                chr4                 chr8                chr20 
                 214                  206                  204 
               chr15                chr22                chr13 
                 198                  139                  117 
               chr21                chr18                 chrY 
                  84                   82                    9 
chr17_GL000258v2_alt chr19_KI270921v1_alt  chr1_KI270766v1_alt 
                   1                    1                    1 
```

Using this simple model, the results obtained are **6457 DE genes**. As it can be seen in the chromosome distribution, the chromosome with more DE genes is chr1.


Two useful diagnostic plots (Figure S13) for DE analysis are the distributions of p-values and moderated t-statistics:

<img src="figure/fit1_diagnostic-1.png" title="Figure S13: Distribution of p-values and moderated t-statistics for the directly fitted model" alt="Figure S13: Distribution of p-values and moderated t-statistics for the directly fitted model" height="300px" style="display: block; margin: auto;" /><p class="caption">Figure S13: Distribution of p-values and moderated t-statistics for the directly fitted model</p>

### 2. Fit adjusting for the mean-variance relationship

RNA-seq counts may vary considerably from sample to sample for the same gene and different samples may be sequenced to different depths. This may lead to identical CPM values for very different count sizes and motivates the need to model the mean-variance trend of log-CPM values at the individual observation level.

What we will do is to calculate weights that estimate the mean-variance relationship at individual observation-by-gene level as Figure S14 shows:


```r
v <- voom(dge, design, plot=TRUE)
```

<img src="figure/fit2_voom-1.png" title="Figure S14: Mean-variance relationship at individual observation-by-gene level" alt="Figure S14: Mean-variance relationship at individual observation-by-gene level" height="500px" style="display: block; margin: auto;" /><p class="caption">Figure S14: Mean-variance relationship at individual observation-by-gene level</p>

Now, we will fit again the linear model this time using the voom weights:


```r
fit2 <- lmFit(v, design)
```

Calculate again moderated t-statistics:


```r
fit2 <- eBayes(fit2)
```

Examine the extent of differential expression at 5% FDR:


```r
res2 <- decideTests(fit2, p.value = FDRcutoff)
summary(res2)
```

```
   (Intercept) factor(type)tumor
-1          21              3289
0          294              5130
1        11256              3152
```

And then, we will fetch the table of results and the number of DE genes:


```r
fit2$genes <- genesmd
tt2 <- topTable(fit2, coef = 2, n = Inf)
table(tt$adj.P.Val < FDRcutoff)
```

```

FALSE  TRUE 
 5114  6457 
```

```r
DEgenes2 <- rownames(tt2)[tt2$adj.P.Val < FDRcutoff]
DEgenes_symbol2 <- tt2$symbol[tt2$adj.P.Val < FDRcutoff]
length(DEgenes2)
```

```
[1] 6441
```

We will examine again the chromosome distribution of genes called DE at 5% FDR:


```r
sort(table(tt2$chr[tt2$adj.P.Val < FDRcutoff]), decreasing = TRUE)
```

```

                chr1                chr19                chr11 
                 670                  461                  446 
                chr2                 chr3                chr17 
                 441                  369                  359 
               chr12                 chr7                 chr6 
                 344                  311                  283 
                chr5                 chr9                chr16 
                 277                  260                  259 
                chrX                chr14                chr10 
                 239                  234                  227 
                chr4                 chr8                chr20 
                 212                  207                  203 
               chr15                chr22                chr13 
                 201                  141                  117 
               chr21                chr18                 chrY 
                  86                   83                    8 
chr17_GL000258v2_alt chr19_KI270921v1_alt  chr1_KI270766v1_alt 
                   1                    1                    1 
```

Using this approach, the results obtained are **6441 DE genes**. As it can be seen in the chromosome distribution, the chromosome with more DE genes is chr1.


The diagnostic plots (Figure S15) for limma DE analysis with voom weights are the following:

<img src="figure/fit2_diagnostic-1.png" title="Figure S15: Distribution of p-values and moderated t-statistics for the model with voom weights" alt="Figure S15: Distribution of p-values and moderated t-statistics for the model with voom weights" height="300px" style="display: block; margin: auto;" /><p class="caption">Figure S15: Distribution of p-values and moderated t-statistics for the model with voom weights</p>


### 3. Adjust for known covariates

We will try to identify **potential sources of unwanted variation** in the phenotypic data, in order to adjust one of them to the model. Possible options could be:
- Histologic diagnosis
- Tumor invasion percent
- Age at diagnosis
- Race
- Tissue Source Site

However, the majority of them are not possible because the tumor samples contain the corresponding phenotypic data but the normal samples not. The only possible source of variation which is in both tumor and normal samples is the **Tissue Source Site**, which is obtained from the patient barcode.


```r
table(substr(se$bcr_patient_barcode, 6, 7))
```

```

AJ AX BG BK E6 
 6 12  6  6  2 
```

We have the following possible sites inside TSS of our samples:
- AJ: International Genomics Consortium 
- AX: Gynecologic Oncology Group 
- BG: University of Pittsburgh 
- BK: Christiana Healthcare 
- E6: Roswell Park

So, we adjust the model for the `tissue source site` factor. We need to start over building a new design matrix that includes that factor.


```r
design <- model.matrix(~type + substr(se$bcr_patient_barcode, 6, 7), data = colData(se))
```

Fit again the linear models for each gene with the updated design matrix:


```r
fit3 <- lmFit(v, design)
fit3 <- eBayes(fit3)
```

Examine the extent of differential expression at 5% FDR:


```r
res3 <- decideTests(fit3, p.value = FDRcutoff)
summary(res3)
```

```
   (Intercept) typetumor substr(se$bcr_patient_barcode, 6, 7)AX
-1           9      3307                                      0
0          339      5080                                  11571
1        11223      3184                                      0
   substr(se$bcr_patient_barcode, 6, 7)BG
-1                                      0
0                                   11571
1                                       0
   substr(se$bcr_patient_barcode, 6, 7)BK
-1                                      0
0                                   11571
1                                       0
   substr(se$bcr_patient_barcode, 6, 7)E6
-1                                      7
0                                   11563
1                                       1
```

Add gene metadata and fetch table of results and the number of DE genes


```r
fit3$genes <- genesmd
tt3 <- topTable(fit3, coef = 2, n = Inf)
DEgenes3 <- rownames(tt3)[tt3$adj.P.Val < FDRcutoff]
DEgenes_symbol3 <- tt3$symbol[tt3$adj.P.Val < FDRcutoff]
length(DEgenes3)
```

```
[1] 6491
```

Examine the chromosome distribution of genes called DE at 5% FDR:


```r
sort(table(tt3$chr[tt3$adj.P.Val < FDRcutoff]), decreasing = TRUE)
```

```

                chr1                chr19                chr11 
                 673                  462                  449 
                chr2                 chr3                chr17 
                 438                  369                  367 
               chr12                 chr7                 chr6 
                 345                  311                  285 
                chr5                 chr9                chr16 
                 280                  264                  259 
                chrX                chr14                chr10 
                 245                  235                  227 
                chr4                 chr8                chr15 
                 216                  214                  206 
               chr20                chr22                chr13 
                 204                  144                  116 
               chr21                chr18                 chrY 
                  86                   85                    8 
chr17_GL000258v2_alt chr19_KI270921v1_alt  chr1_KI270766v1_alt 
                   1                    1                    1 
```

Using this approach, the results obtained are **6491 DE genes**. As it can be seen in the chromosome distribution, the chromosome with more DE genes is chr1.


### 4. Adjust for unknown covariates

The previous model is only adjusted for the known covariate. It is also possible to adjust it for unknown covariates using **surrogate variable analysis (SVA)**. First, estimate surrogate variables (SVs) from the log-CPM values calculated by voom:


```r
mod0 <- model.matrix(~substr(se$bcr_patient_barcode, 6, 7), colData(se))
sv <- sva(v$E, mod = design, mod0 = mod0)
```

```
Number of significant surrogate variables is:  9 
Iteration (out of 5 ):1  2  3  4  5  
```

```r
sv$n # number of surrogate variables
```

```
[1] 9
```

There are **9 SVs**. Second, we add these SVs to the design matrix, which has the following form:


```r
design <- cbind(design, sv$sv)
colnames(design) <- c(colnames(design)[1:6], paste0("SV", 1:sv$n))
head(design)
```

```
                             (Intercept) typetumor
TCGA.AJ.A3NC.01A.11R.A22K.07           1         1
TCGA.AJ.A3NE.01A.11R.A22K.07           1         1
TCGA.AJ.A3NH.01A.11R.A22K.07           1         1
TCGA.AX.A05Y.01A.11R.A00V.07           1         1
TCGA.AX.A0IZ.01A.11R.A118.07           1         1
TCGA.AX.A0J0.01A.11R.A109.07           1         1
                             substr(se$bcr_patient_barcode, 6, 7)AX
TCGA.AJ.A3NC.01A.11R.A22K.07                                      0
TCGA.AJ.A3NE.01A.11R.A22K.07                                      0
TCGA.AJ.A3NH.01A.11R.A22K.07                                      0
TCGA.AX.A05Y.01A.11R.A00V.07                                      1
TCGA.AX.A0IZ.01A.11R.A118.07                                      1
TCGA.AX.A0J0.01A.11R.A109.07                                      1
                             substr(se$bcr_patient_barcode, 6, 7)BG
TCGA.AJ.A3NC.01A.11R.A22K.07                                      0
TCGA.AJ.A3NE.01A.11R.A22K.07                                      0
TCGA.AJ.A3NH.01A.11R.A22K.07                                      0
TCGA.AX.A05Y.01A.11R.A00V.07                                      0
TCGA.AX.A0IZ.01A.11R.A118.07                                      0
TCGA.AX.A0J0.01A.11R.A109.07                                      0
                             substr(se$bcr_patient_barcode, 6, 7)BK
TCGA.AJ.A3NC.01A.11R.A22K.07                                      0
TCGA.AJ.A3NE.01A.11R.A22K.07                                      0
TCGA.AJ.A3NH.01A.11R.A22K.07                                      0
TCGA.AX.A05Y.01A.11R.A00V.07                                      0
TCGA.AX.A0IZ.01A.11R.A118.07                                      0
TCGA.AX.A0J0.01A.11R.A109.07                                      0
                             substr(se$bcr_patient_barcode, 6, 7)E6
TCGA.AJ.A3NC.01A.11R.A22K.07                                      0
TCGA.AJ.A3NE.01A.11R.A22K.07                                      0
TCGA.AJ.A3NH.01A.11R.A22K.07                                      0
TCGA.AX.A05Y.01A.11R.A00V.07                                      0
TCGA.AX.A0IZ.01A.11R.A118.07                                      0
TCGA.AX.A0J0.01A.11R.A109.07                                      0
                                      SV1         SV2         SV3
TCGA.AJ.A3NC.01A.11R.A22K.07  0.102406869 -0.12182878  0.27746076
TCGA.AJ.A3NE.01A.11R.A22K.07  0.065448684 -0.28952300  0.15155429
TCGA.AJ.A3NH.01A.11R.A22K.07 -0.003567403  0.05556203  0.07984721
TCGA.AX.A05Y.01A.11R.A00V.07  0.223533483  0.12012572 -0.38638515
TCGA.AX.A0IZ.01A.11R.A118.07  0.156016392 -0.27855529 -0.01448107
TCGA.AX.A0J0.01A.11R.A109.07 -0.203823846  0.08789727 -0.08187704
                                     SV4         SV5         SV6
TCGA.AJ.A3NC.01A.11R.A22K.07 -0.17994665  0.08821434 -0.18239974
TCGA.AJ.A3NE.01A.11R.A22K.07 -0.03438035 -0.30359959 -0.13958625
TCGA.AJ.A3NH.01A.11R.A22K.07 -0.05125137  0.08454969  0.36926080
TCGA.AX.A05Y.01A.11R.A00V.07  0.19848949  0.08507325 -0.26517806
TCGA.AX.A0IZ.01A.11R.A118.07 -0.23306995 -0.34590747 -0.11788951
TCGA.AX.A0J0.01A.11R.A109.07 -0.25701563  0.20416750  0.07020241
                                      SV7         SV8         SV9
TCGA.AJ.A3NC.01A.11R.A22K.07  0.117946610 -0.03003028  0.03190465
TCGA.AJ.A3NE.01A.11R.A22K.07  0.285120986  0.09296381 -0.05701625
TCGA.AJ.A3NH.01A.11R.A22K.07  0.003551297 -0.45397797  0.11067228
TCGA.AX.A05Y.01A.11R.A00V.07 -0.389586836 -0.24884878 -0.32586872
TCGA.AX.A0IZ.01A.11R.A118.07  0.084247726  0.13700105  0.01980260
TCGA.AX.A0J0.01A.11R.A109.07 -0.134303513  0.27996279 -0.13762299
```

Third, we fit again the linear models for each gene with the updated design matrix, and calculate the moderated t-statistics:


```r
fit4 <- lmFit(v, design)
fit4 <- eBayes(fit4)
```

We examine the extent of differential expression at 5% FDR:


```r
res4 <- decideTests(fit4, p.value = FDRcutoff)
```

Finally, the metadata is added, the table of results is computed and the number of DE genes is calculated:


```r
fit4$genes <- genesmd
tt4 <- topTable(fit4, coef = 2, n = Inf)
DEgenes4 <- rownames(tt4)[tt4$adj.P.Val < FDRcutoff]
DEgenes_symbol4 <- tt4$symbol[tt4$adj.P.Val < FDRcutoff]
length(DEgenes4)
```

```
[1] 7845
```

The chromosome distribution of genes called DE at 5% FDR is calculated:


```r
sort(table(tt4$chr[tt4$adj.P.Val < FDRcutoff]), decreasing = TRUE)
```

```

                chr1                chr19                chr11 
                 813                  562                  542 
                chr2                chr17                 chr3 
                 533                  445                  444 
               chr12                 chr7                 chr6 
                 405                  371                  341 
                chr5                 chr9                chr16 
                 339                  325                  323 
                chrX                chr14                chr10 
                 298                  297                  277 
                chr4                 chr8                chr15 
                 270                  258                  236 
               chr20                chr22                chr13 
                 233                  176                  142 
               chr21                chr18                 chrY 
                 102                   98                   11 
chr11_JH159137v1_alt chr17_GL000258v2_alt chr19_KI270921v1_alt 
                   1                    1                    1 
 chr1_KI270766v1_alt 
                   1 
```

Using this approach, there are **7845 DE genes**. As it can be seen in the chromosome distribution, the chromosome with more DE genes is chr1.


Now, it is possible to examine the diagnostic plots (Figure S16) for the DE analysis using this approach:

<img src="figure/fit4_diagnostic-1.png" title="Figure S16: Distribution of p-values and moderated t-statistics for the model with known covariates and SVA" alt="Figure S16: Distribution of p-values and moderated t-statistics for the model with known covariates and SVA" height="300px" style="display: block; margin: auto;" /><p class="caption">Figure S16: Distribution of p-values and moderated t-statistics for the model with known covariates and SVA</p>


### Summary of the obtained results and choice of the model

The ideal would be to do a performance comparison of the different models, but it is not possible as it would be needed a dataset of documented differentially expressed genes in endometrial cancer in order to calculate the accuracy of each model. As there is no such benchmark, an alternative is to look at the statistical power in terms of DE genes. 

The following table summarizes the results of DE genes for the different approaches tried during the analysis. 

Approach used | DE genes | Chromosome with more DE
--------------|----------|------------------------
1. Fit directly the linear model  | 6457 | chr1
2. Adjust for mean-variance relationship | 6441 | chr1
3. Adjust for known covariates | 6491 | chr1
4. Adjust for known covariates + SVA | 7845 | chr1

The number of DE genes importantly increases in the fourth model with respect to the other models. If the outcome of interest is not confounded with other sources of variation, the more variables that we put in the linear model, the fewer degrees of freedom it has, and therefore the less statistical power that it has. So, the model that is going to be used for the following analysis is the fourth model, **adjusting for known covariates plus the SVA**.

### Volcano plot

A useful diagnostic plot for DE analysis is the so-called volcano plot (Figure S17). As it was said, we will use the model 4 corresponding to the adjust of known covariates plus the SVA, but we will compare it to the volcano plot of the model 2. The remarked genes in volcano plots change. For model 4, they are MIR4298, TSN, RNF122, CALHM2, TSG101, CACNA1S and SCARNA5.

<img src="figure/volcanoPlot-1.png" title="Figure S17: Volcano plots for models 2 and 4" alt="Figure S17: Volcano plots for models 2 and 4" height="500px" style="display: block; margin: auto;" /><p class="caption">Figure S17: Volcano plots for models 2 and 4</p>

### MA-plot

Another useful diagnostic plot (Figure S18) is the MA-plot:

<img src="figure/maPlot-1.png" title="Figure S18: MA-plot for model 4" alt="Figure S18: MA-plot for model 4" height="500px" style="display: block; margin: auto;" /><p class="caption">Figure S18: MA-plot for model 4</p>

### Factorial designs

Factorial designs are experimental designs where the outcome of interest consists of two or more factors.
Here we consider for illustrative purposes the analysis of the interaction between `type` and `tissue source site` factors.
There are several ways to approach this kind of analysis. One that facilitates the interpretation is building a single factor out of the combination of the factors under study and set the intercept term to zero.


```r
newfac <- factor(paste(se$type, substr(se$bcr_patient_barcode, 6, 7), sep = "."))
head(newfac)
```

```
[1] tumor.AJ tumor.AJ tumor.AJ tumor.AX tumor.AX tumor.AX
10 Levels: normal.AJ normal.AX normal.BG normal.BK normal.E6 ... tumor.E6
```

Fit a linear regression model using this new factor variable without intercept term.


```r
design <- model.matrix(~0 + newfac, colData(se))
head(design, n = 3)
```

```
                             newfacnormal.AJ newfacnormal.AX
TCGA.AJ.A3NC.01A.11R.A22K.07               0               0
TCGA.AJ.A3NE.01A.11R.A22K.07               0               0
TCGA.AJ.A3NH.01A.11R.A22K.07               0               0
                             newfacnormal.BG newfacnormal.BK
TCGA.AJ.A3NC.01A.11R.A22K.07               0               0
TCGA.AJ.A3NE.01A.11R.A22K.07               0               0
TCGA.AJ.A3NH.01A.11R.A22K.07               0               0
                             newfacnormal.E6 newfactumor.AJ newfactumor.AX
TCGA.AJ.A3NC.01A.11R.A22K.07               0              1              0
TCGA.AJ.A3NE.01A.11R.A22K.07               0              1              0
TCGA.AJ.A3NH.01A.11R.A22K.07               0              1              0
                             newfactumor.BG newfactumor.BK newfactumor.E6
TCGA.AJ.A3NC.01A.11R.A22K.07              0              0              0
TCGA.AJ.A3NE.01A.11R.A22K.07              0              0              0
TCGA.AJ.A3NH.01A.11R.A22K.07              0              0              0
```

```r
fit6 <- lmFit(v, design)
```

Build the contrast matrix that specifies the contrasts of interest between a set of coefficients estimated from a linear regression model. For instance, let's say we want to search for DE genes between BG and BK tissue source sites separately in tumor and normal samples:


```r
#table(substr(se$bcr_patient_barcode, 6, 7))
cont.matrix <- makeContrasts(tssnormal = newfacnormal.BG - newfacnormal.BK, tsstumor = newfactumor.BG - newfactumor.BK, levels = design)
```

Estimate coefficients for a given set of contrasts.


```r
fit6 <- contrasts.fit(fit6, cont.matrix)
```

Calculate moderated t -statistics.


```r
fit6 <- eBayes(fit6)
head(fit6$t)
```

```
    Contrasts
      tssnormal   tsstumor
  1  -0.2655023 -0.2593460
  2  -1.9442266 -0.5884159
  12 -0.5951888  0.8157689
  14  1.6885430  0.7683610
  16 -0.4830050  2.3325289
  18  0.7823703  0.5678817
```

Fetch table of results for each coefficient.

```r
tttssnormal <- topTable(fit6, coef = "tssnormal", n = Inf)
tttsstumor <- topTable(fit6, coef = "tsstumor", n = Inf)
```

Now, it is possible to explore graphically the overlap of DE genes between contrasts of interest (Figure S19):

<img src="figure/vennDiagram-1.png" title="Figure S19: VennDiagram: Overlap of DE genes between contrasts of interest" alt="Figure S19: VennDiagram: Overlap of DE genes between contrasts of interest" height="300px" style="display: block; margin: auto;" /><p class="caption">Figure S19: VennDiagram: Overlap of DE genes between contrasts of interest</p>

As it can be seen, there is not a clear interaction between TSS and type factors. Probably, there is not much unwanted variability in the TSS factor. It could be interest for future analyses to use `hystologic_analysis`, as there are several studies reporting differentially expressed genes in function of the type of histologic diagnosis.

## Functional Enrichment: The Gene Ontology analysis

A popular kind of functional enrichment analysis is what people call a _Gene Ontology (GO) analysis_ which often corresponds to applying the one-tailed Fisher's exact test to every gene set in the [GO database](http://www.geneontology.org). Functional enrichment analyses constitute a straightforward way to approach the question of what pathways may be differentially expressed (DE) in our data.

The [GO database project](http://www.geneontology.org) provides a controlled vocabulary to describe gene and gene product attributes in any organism. It consists of so-called _GO terms_, which are pairs of term identifier (GO ID) and description.

There are several R packages at CRAN/Bioconductor that facilitate performing a functional enrichment analysis on the entire collection of GO gene sets. We are going to illustrate this analysis with the Bioconductor package [GOstats](http://www.bioconductor.org/packages/release/bioc/html/GOstats.html).

Doing this analysis with [GOstats](http://www.bioconductor.org/packages/release/bioc/html/GOstats.html) consists of the following three steps:

1. Build a parameter object with information specifying the gene universe, the set of DE genes, the annotation package to use, etc. as follows:


```r
library(org.Hs.eg.db)
```

```

```

```r
library(GOstats)
```

```
Loading required package: Category
```

```
Loading required package: Matrix
```

```

Attaching package: 'Matrix'
```

```
The following object is masked from 'package:S4Vectors':

    expand
```

```
Loading required package: graph
```

```

Attaching package: 'graph'
```

```
The following object is masked from 'package:XML':

    addNode
```

```

```

```

Attaching package: 'GOstats'
```

```
The following object is masked from 'package:AnnotationDbi':

    makeGOGraph
```

```r
geneUniverse <- rownames(se)
params <- new("GOHyperGParams", geneIds=DEgenes4, universeGeneIds=geneUniverse,
                annotation="org.Hs.eg.db", ontology="BP",
                pvalueCutoff=0.05, testDirection="over")
```

These type of techniques are limited by the amount of DE genes we have in our data. The total size of the genes involved in the calculation (the gene universe) is critical to the significance level of the results.

2. Run the functional enrichment analysis.
  A problem in a GO analysis is that the hierarchy of GO terms and their overlap render highly dependent tests. If parent and a child GO term contain the same genes and both are significant, the child node is more relevant because is more specific.
  Compute the significance of a GO term conditional on the significance of its children ([Alexa et al.](http://bioinformatics.oxfordjournals.org/content/22/13/1600.abstract, 2006)).
  
  Proceed bottom-up removing all genes in a significant GO term from its parents. This is regarded as a __conditional test__.
  In the [GOstats](http://www.bioconductor.org/packages/release/bioc/html/GOstats.html) package, the conditional test can be used by simply setting the argument `conditional=TRUE` in the parameter object or modifying it as follows:
  In the case of a GO analysis is important to always perform a conditional test that takes into account the hierarchical structure of GO terms.

```r
conditional(params) <- TRUE

hgOverCond <- hyperGTest(params)
hgOverCond
```

```
Gene to GO BP Conditional test for over-representation 
11783 GO BP ids tested (131 have p < 0.05)
Selected gene set size: 6090 
    Gene universe size: 8939 
    Annotation package: org.Hs.eg 
```

3. Store and visualize the results.


```r
htmlReport(hgOverCond, file="goconditionaltests.html")
```

We can find out what methods are available to explore the results in detail. The most directly useful one is the `summary()` method that returns a `data.frame` object with results.


```
      GOBPID       Pvalue OddsRatio   ExpCount Count Size
1 GO:0019054 0.0008703101  5.168565  24.526233    33   36
2 GO:0010565 0.0027682020  2.007778  68.128426    81  100
3 GO:0051047 0.0033070683  1.605012 128.708767   146  189
4 GO:0044003 0.0041594172  2.819444  33.382929    42   49
5 GO:0033014 0.0066931710  8.442688  12.944401    18   19
6 GO:0010762 0.0067841237       Inf   8.856695    13   13
                                                       Term
1                       modulation by virus of host process
2           regulation of cellular ketone metabolic process
3                          positive regulation of secretion
4 modification by symbiont of host morphology or physiology
5                         tetrapyrrole biosynthetic process
6                        regulation of fibroblast migration
```

Storing the results in a `data.frame` object enables an automatic processing and filtering of the results.


```r
goresults <- summary(hgOverCond)
goresults[1:10,]
```

```
       GOBPID       Pvalue OddsRatio   ExpCount Count Size
1  GO:0019054 0.0008703101  5.168565  24.526233    33   36
2  GO:0010565 0.0027682020  2.007778  68.128426    81  100
3  GO:0051047 0.0033070683  1.605012 128.708767   146  189
4  GO:0044003 0.0041594172  2.819444  33.382929    42   49
5  GO:0033014 0.0066931710  8.442688  12.944401    18   19
6  GO:0010762 0.0067841237       Inf   8.856695    13   13
7  GO:0048168 0.0072030358  5.161009  16.350822    22   24
8  GO:0044774 0.0093654411  2.178838  42.239624    51   62
9  GO:1903543 0.0099641119       Inf   8.175411    12   12
10 GO:0071156 0.0112560898  2.024630  47.008614    56   69
                                                        Term
1                        modulation by virus of host process
2            regulation of cellular ketone metabolic process
3                           positive regulation of secretion
4  modification by symbiont of host morphology or physiology
5                          tetrapyrrole biosynthetic process
6                         regulation of fibroblast migration
7                 regulation of neuronal synaptic plasticity
8                           mitotic DNA integrity checkpoint
9                  positive regulation of exosomal secretion
10                           regulation of cell cycle arrest
```

GO terms involving a few genes (e.g., < 5) in their size ($m$) and in their enrichment by DE genes are likely to be less reliable than those that involve many genes.

The `OddsRatio` (OR) is a measure of association between an exposure and an outcome. The OR represents the odds that an outcome will occur given a particular exposure, compared to the odds of the outcome occurring in the absence of that exposure.

In order to try to spot the more reliable GO terms we can filter the previous results by a minimum value on the `Count` and `Size` columns and order them by the `OddsRatio` column:


```r
goresults <- goresults[goresults$Size >= 5 & goresults$Count >= 5, ]
goresults <- goresults[order(goresults$OddsRatio, decreasing=TRUE), ]
head(goresults)
```

```
       GOBPID      Pvalue OddsRatio ExpCount Count Size
6  GO:0010762 0.006784124       Inf 8.856695    13   13
9  GO:1903543 0.009964112       Inf 8.175411    12   12
15 GO:0007213 0.014633919       Inf 7.494127    11   11
16 GO:0032735 0.014633919       Inf 7.494127    11   11
17 GO:0046636 0.014633919       Inf 7.494127    11   11
18 GO:0071380 0.014633919       Inf 7.494127    11   11
                                                         Term
6                          regulation of fibroblast migration
9                   positive regulation of exosomal secretion
15 G-protein coupled acetylcholine receptor signaling pathway
16           positive regulation of interleukin-12 production
17        negative regulation of alpha-beta T cell activation
18              cellular response to prostaglandin E stimulus
```

We can extract the genes that _enrich_ each GO term and paste it to the result as follows:


```r
geneIDs <- geneIdsByCategory(hgOverCond)[goresults$GOBPID]
geneSYMs <- sapply(geneIDs, function(id) select(org.Hs.eg.db, columns="SYMBOL", key=id, keytype="ENTREZID")$SYMBOL)
```

```
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
'select()' returned 1:1 mapping between keys and columns
```

```r
geneSYMs <- sapply(geneSYMs, paste, collapse=", ")
goresults <- cbind(goresults, Genes=geneSYMs)
rownames(goresults) <- 1:nrow(goresults)
```

We can generate an HTML page from a `data.frame` object using the `xtable` package:


```r
library(xtable)
xtab <- xtable(goresults, align="l|c|r|r|r|r|r|p{3cm}|p{3cm}|")
print(xtab, file="goresults.html", type="html")
```


## Conclusions

The main source of variation of the dataset seems to be driven by the tumor and normal condition of the samples, as seen in the figure 1. The hierarchical clustering of the samples showed a better stratification when adjusting for the tss factor, without this batch effect. 

The SVA analysis helps identifying more differentially expressed genes, since adjusting for other sources of variation reduces the degrees of variation of the studied variable. That is translated to more statistical power in the model and hence to more DEgenes found. 

The functional enrichment analysis with GO terms identified nine out of the ten most differentially expressed groups that perfectly fit in the pattern of a tumorous development. 

One of those nine cancer related GO terms, defined as negative regulation of smooth muscle contraction, can be specifically related to the endometrial cancer. 

At least two found differentially expressed genes (PIK3C2A, CTNNB1) are known to be involved specifically to endometrial cancer. 


## Session information


```r
sessionInfo()
```

```
R version 3.3.0 (2016-05-03)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 14.04.4 LTS

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=es_ES.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=es_ES.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=es_ES.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=es_ES.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats4    parallel  stats     graphics  grDevices utils     datasets 
[8] methods   base     

other attached packages:
 [1] xtable_1.8-2               GO.db_3.3.0               
 [3] GOstats_2.38.0             graph_1.50.0              
 [5] Category_2.38.0            Matrix_1.2-6              
 [7] org.Hs.eg.db_3.3.0         sva_3.20.0                
 [9] genefilter_1.54.2          mgcv_1.8-12               
[11] nlme_3.1-128               geneplotter_1.50.0        
[13] annotate_1.50.0            XML_3.98-1.4              
[15] AnnotationDbi_1.34.3       lattice_0.20-33           
[17] edgeR_3.14.0               limma_3.28.5              
[19] SummarizedExperiment_1.2.2 Biobase_2.32.0            
[21] GenomicRanges_1.24.0       GenomeInfoDb_1.8.1        
[23] IRanges_2.6.0              S4Vectors_0.10.1          
[25] BiocGenerics_0.18.0        markdown_0.7.7            
[27] knitr_1.13                

loaded via a namespace (and not attached):
 [1] formatR_1.4            RColorBrewer_1.1-2     XVector_0.12.0        
 [4] tools_3.3.0            zlibbioc_1.18.0        digest_0.6.9          
 [7] RSQLite_1.0.0          evaluate_0.9           DBI_0.4-1             
[10] stringr_1.0.0          grid_3.3.0             GSEABase_1.34.0       
[13] RBGL_1.48.1            survival_2.39-4        magrittr_1.5          
[16] codetools_0.2-14       splines_3.3.0          AnnotationForge_1.14.2
[19] KernSmooth_2.23-15     stringi_1.1.1         
```
