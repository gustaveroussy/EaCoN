# **EaCoN**

_**Ea**sy **Co**py **N**umber !_

---

## **DESCRIPTION**

EaCoN aims to be a user-friendly solution to perform relative or absolute copy-number analysis for multiple sources of data.
It consists in a series of R packages that perform such type of analysis, from raw CEL files of Affymetrix microarrays (GenomeWide snp6, OncoScan, CytoScan 750K, CytoScan HD) or from aligned reads as BAMs for WES (whole exome sequencing).

---

## **FEATURES**

- Full R (no OS-, language-, nor shell-specific dependencies)
- Platform-independent (ie, the three most common ones : Linux, Windows and OSX)
- Support of multiple sources of data :
  - Affymetrix (homo sapiens)
    - Genomide SNP6
    - CytoScan 750k
    - CytoScan HD
    - OncoScan / OncoScan_CNV
  - WES (any species)
    - Any design for which a capture BED and known genome are available
- Full workflow analysis, from the raw data to an annotated HTML report :
  - L2R and BAF normalization from raw data
    - Using **Affymetrix Power Tools** for Affymetrix arrays and **rawcopy**
    - Using internal methods for WES data
  - L2R and BAF bivariate segmentation (using **ASCAT**)
  - L2R profile centralization
  - L2R profile CNA calling
  - Generation of results as tables and plots
  - Generation of a portable and interactive HTML report
  - Total (TCN) and allele-specific copy number (ASCN), ploidy and cellularity estimations (using **ASCAT**)
- Required annotations for Affymetrix data provided (via additional packages)
- Pre-computed normalization (GC%, wave) tracks provided (either)
- Few functions to call
- Few parameters
- Recommended guidelines provided for supported microarray designs, as well as for WES

---

## **REQUIREMENTS**

While the current EaCoN package is the core of the process, multiple others are needed, depending on the type of the source (WES, Affymetrix microarray design) and corresponding annotations (genome build, Affymetrix annotation databases). Others are required for the normalization, especially pre-computed GC% tracks.

- R > 3.0
- Public R dependencies :
  - [@ CRAN](https://cran.r-project.org/) :
    - [bedr]() _to handle capture bed_
    - [data.table]() _for tab-files / tables / blocks processing_
    - [devtools]() _to install ASCAT_
    - [dplyr]() _for tables / blocks processing_
    - [doParallel]() _for multithreading_
    - [DT]() _to generate nice tables in the HTML report_
    - [foreach]() _for loops & multithreading_
    - [iotools]() _for fast TSV writings_
    - [matrixStats]() _for column-wise statistics computing_
    - [mclust]() _for BAF regions clustering_
    - [rmarkdown]() _to generate the HTML report_
  - [@ BioConductor](https://www.bioconductor.org/) :
    - [affxparser]() _to read CEL files_
    - [BSgenome]() _to handle the reference genome and its sequences_
    - [BSgenome.Hsapiens.UCSC.hg19]() _as default genome_
    - [changepoint]() _to perform small events rescue segmentation using PELT_
    - [copynumber]() _to perform PCF segmentation_
    - [GenomeInfoDb]() _to use genome annotations_
    - [GenomicRanges]() _to handle genomic intervals_
    - [limma]() _to perform lowess normalization_
    - [rhdf5]() _to read HDF5 file format_
    - [Rsamtools]() _to compute counts and variants from BAM files_
  - [@ GitHub](https://github.com)
    - [ASCAT](https://github.com/Crick-CancerGenomics/ascat/ASCAT) _to perform the L2R & BAF bivariate segmentation, and absolute and allele-specific copy number estimation_
- Elements of the EaCoN workflow :
  - [@ GitHub](https://github.com/gustaveroussy) :
    - [apt.cytoscan.2.4.0](https://github.com/gustaveroussy/apt.cytoscan.2.4.0) _to perform the normalization of CytoScan (750k, HD) arrays using Affymetrix Power Tools_
    - [apt.oncoscan.2.4.0](https://github.com/gustaveroussy/apt.oncoscan.2.4.0) _to perform the normalization of OncoScan / OncoScan\_CNV arrays using Affymetrix Power Tools_
    - [apt.snp6.1.20.0](https://github.com/gustaveroussy/apt.snp6.1.20.0) _to perform the normalization of GenomeWide SNP6 arrays using Affymetrix Power Tools_
    - [chromosomes](https://github.com/gustaveroussy/chromosomes) _to provide the structure of chromsomes for homo sapiens (hg17 to hg38), mus musculus (mm7 to mm10) and rattus norvegicus (rn5 to rn6)_
  - [@ Google Drive](http://bit.ly/EaCoNpackages) : _Containing numerous data, these packages could not fit on GitHub_
    - Affymetrix design nnotations :
      - [CytoScanHD.Array.na33.r4](http://bit.ly/CSHDna33) _Affymetrix design annotations for the CytoScan HD array (build na33.r4 / hg19)_
      - [CytoScanHD.Array.na36.r1](http://bit.ly/CSHDna36) _Affymetrix design annotations for the CytoScan HD array (build na36.r1 / hg38)_
      - [CytoScan750K.Array.na33.r4](http://bit.ly/CS750na33) _Affymetrix design annotations for the CytoScan 750K array (build na33.r4 / hg19)_
      - [CytoScan750K.Array.na36.r1](http://bit.ly/CS750na36) _Affymetrix design annotations for the CytoScan 750K array (build na36.r1 / hg38)_
      - [OncoScan.na33.r4](http://bit.ly/OSna33r4) _Affymetrix design annotations for the OncoScan array (build na33.r4 / hg19)_
      - [OncoScan.na36.r1](http://bit.ly/OSna36r1) _Affymetrix design annotations for the OncoScan array (build na36.r1 / hg38)_
      - [OncoScanCNV.na33.r2](http://bit.ly/OSCNVna33) _Affymetrix design annotations for the OncoScan\_CNV array (build na33.r2 / hg19)_
      - [OncoScanCNV.na36.r1](http://bit.ly/OSCNVna36) _Affymetrix design annotations for the OncoScan\_CNV array (build na36.r1 / hg38)_
      - [GenomeWideSNP.6.na35.r1](http://bit.ly/SNP6na35) _Affymetrix design annotations for the GenomeWide\_SNP.6 array (build na35.r1 / hg19)_
    - Re-normalization :
      - GC% & wave-effect normalization :
        - [affy.CN.norm.data](http://bit.ly/AffyCNnorm) _Pre-computed tracks for all supported Affmetrix designs (both hg19 and hg38)_
        - [rcnorm](http://bit.ly/chromosomespackage) _Code and tracks to re-normalize BAF for the CytoScan family of designs and SNP6, using *rawcopy*_
- Raw data :
  - For Affymetrix microarrays : the CEL files, fresh out of the Affymetrix Scanner
  - For WES data :
    - The aligned reads of the test and reference samples, in two different BAM files
    - The capture BED file (the file containing the genomic regions targeted by the exome capture kit used in the wet experiment)

---

## **USAGE**

The full workflow is decomposed in a few different functions, which roughly correspond to these steps :  

> normalization -> segmentation |-> reporting |-> copy-number estimation  

EaCoN allows different ways to perform the full workflow : considering the analysis of a single sample, you can either perform each step independently and write, then load the intermediate results, or you can _**pipe**_ all steps in a single line of code. You can also perform the step-by-step approach on multiple samples in a row, even possibly at the same time using multithreading, using a batch mode.

### **Step by step mode**

- First, under R, load EaCoN and choose a directory in which results will be written, by exemple : **/home/me/my_project/EaCoN_results**

```R
require(EaCoN)
setwd("/home/me/my_project/EaCoN_results")
```

#### **Raw data processing**

##### **Affymetrix OncoScan / OncoScan_CNV**

- Let's say we have a pair of OncoScan_CNV CEL files to analyse in a **/home/me/my_project/CEL/** directory (Affymetrix OncoScan experiments have 2 arrays for a single experiment, thus a pair) :
```R
OS.Process(ATChannelCel = "/home/me/my_project/CEL/SAMPLE1_OncoScan_CNV_A.CEL", GCChannelCel = "/home/me/my_project/CEL/SAMPLE1_OncoScan_CNV_C.CEL", samplename = "SAMPLE1_OS")
```

- This will perform the normalization step, create a **/home/me/my_project/EaCoN_results/SAMPLE1_OS/** subdirectory and write 5 files in it :
  - _**SAMPLE1_OS_OncoScan_CNV_hg19_processed.RDS**_ : contains the normalized data which will be provided to the next step
  - _**SAMPLE1_OS_OncoScan_CNV_hg19_rawplot.png**_ : shows a graphical representation of the normalized L2R and BAF data
  - _**SAMPLE1_OS_2.4.0_na33.r2.paircheck.txt**_ : gives some statistics to evaluate the probability that the two "A" "C" CEL files effectively belong to the same individual (generated by APT)
  - _**SAMPLE1_OS_2.4.0_na33.r2.qc.txt**_ : some quality metrics of the arrays and profiles (generated by APT)
  - _**SAMPLE1_OS_2.4.0_na33.r2.log**_ : a simple log file (generated by APT)

##### **Affymetrix CytoScan 750k / CytoScan HD**

- This is identical to OncoScan, except that we have a single CEL file, and we just have to change the function name :
```R
CS.Process(CEL = "/home/me/my_project/CEL/SAMPLE2_CytoScanHD.CEL", samplename = "SAMPLE2_CSHD")
```
- The same output files will be generated (except for the "paircheck" file, obviously)

##### **Affymetrix GenomeWide SNP6**

- Here again, this is identical to OncoScan, except that we have a single CEL file, and we just have to change the function name :
```R
SNP6.Process(CEL = "/home/me/my_project/CEL/SAMPLE3_GenomeWide_snp.6.CEL", samplename = "SAMPLE3_SNP6")
```
- Again, the same output files will be generated (except for the "paircheck" file, obviously)

##### **WES data**

- This time it is quite different : the processing will be performed in three steps :
  - **First**, we will use the capture BED (A text file containing the positions of the captured regions, usualy provided by the capture kit manufacturer), choose a genome version corresponding to our aligned BAM files, and choose a window size for the future binning of the data. Thses will be used to generate what we call a "BINpack", a set of pre-computed tracks containing the bins position and corresponding GC% values. Several tracks will be computed corresponding to different levels of elongation of the bin positions. In the example below, we used the BED corresponding to Agilent SureSelect v5 capture kit, a bin size of 50 nt, and chose the human hg19 genome build.

```R
BINpack.maker(bed.file = "/home/me/my_project/WES/SureSelect_v5.bed", bin.size = 50, genome.pkg =
"BSgenome.Hsapiens.UCSC.hg19")
```

  - This will generate a "BINpack" (with a ".rda" extension) that will be used in the next normalization steps : **/home/me/my_project/EaCoN_results/SureSelect_v5_merged_sorted_hg19_b50.GC.rda**
  - **PLEASE NOTE THAT THIS STEP IS SAMPLE-INDEPENDENT, THUS NEEDS TO BE PERFORMED AGAIN ONLY IF YOU CHANGE EITHER THE CAPTURE BED, THE BIN SIZE OR THE GENOME BUILD.** Thus, the generated BINpack can be used for any other sample in the same conditions.

  - **Second**, the WES data will be binned using the generated BINpack. We need three files as input : 
    - The aligned reads for the test sample (usualy in cancer, the patient tumor), in BAM format
    - The aligned reads for the reference sample (patient normal), in BAM format too
    - The BINpack itself

```R
WES.Bin(testBAM = "/home/me/my_project/WES/SAMPLE4_WES_hg19_Tumor.BAM", refBAM = "/home/me/my_project/WES/SAMPLE4_WES_hg19_Tumor.BAM", BINpack = "/home/me/my_project/EaCoN_results/SureSelect_v5_merged_sorted_hg19_b50.GC.rda", samplename = "SAMPLE4_WES")
```

  - This will generate a **/home/me/my_project/EaCoN_results/SAMPLE4_WES/** subdirectory which contains :
    - _**SAMPLE4_WES_hg19_b50_binned.RDS**_ : contains the binned data which will be provided to the next step
    - _**SAMPLE4_WES_hg19_b50_coverage.png**_ : shows a graphical representation of the proportion of the capture bed regions covered at different coverage levels
    - _**SAMPLE4_WES_hg19_b50_coverage.txt**_ : contins the numerical values corresponding to the coverage plot

  - **Third**, now that the data have been binned, the normalization step can be performed :

```R
WES.Normalize.ff(BIN.RDS.file = "/home/me/my_project/EaCoN_results/SAMPLE4_WES/SAMPLE4_WES_hg19_b50_binned.RDS", BINpack = "/home/me/my_project/EaCoN_results/SureSelect_v5_merged_sorted_hg19_b50.GC.rda")
```

  - Analogously to what was described for Affymetrix microarrays, a **/home/me/my_project/EaCoN_results/SAMPLE4_WES/ASCAT/L2R/** subdirectory will be created, containing the same type of files described earlierly.


#### **L2R & BAF Segmentation**

- Now that we described the normalization process specific to each type of source data, we can segment it. The good news is that it's the very same step for each source type, one just have to pass the **\_processed.RDS** normalized file. Here is an example with the one obtained for OncoScan data :

```R
Segment.ASCAT.ff(RDS.file = "/home/me/my_project/EaCoN_results/SAMPLE1/SAMPLE1_OncoScan_CNV_hg19_processed.RDS")
```

- This will perform the segmentation, centralization and calling steps, create a **/home/me/my_project/EaCoN_results/SAMPLE1/ASCAT/L2R/** subdirectory and write multiple files in it :
  - _**SAMPLE1.EaCoN.ASPCF.RDS**_ : contains the segmented data which will be provided to the optional next step(s)
  - _**SAMPLE1.ASPCF.png**_ : shows a graphical representation of the segmented, centered and called L2R and BAF data
  - _**SAMPLE1.Rorschach.png**_ : shows a graphical representation of BAF vs L2R of probes, by chromosome
  - _**SAMPLE1.Cut.cbs**_ : contains the L2R segmentation results in the standard CBS format. "Cut" means that L2R value for the segments called as normal were set to a value of 0.0
  - _**SAMPLE1.NoCut.cbs**_ : Same as above, but the segments considered as normal have their L2R value intact
  - _**SAMPLE1.SegmentedBAF.txt**_ : contains the BAF segmentation results

#### **Copy-number estimation**

- Then an estimation of the total and allele-specific copy-number profiles, as well as global ploidy and sample cellularity can be estimated :

```R
ASCN.ASCAT.ff(RDS.file = "/home/me/my_project/EaCoN_results/SAMPLE1/ASCAT/L2R/SAMPLE1.EaCoN.ASPCF.RDS")
```

- This will perform these estimations for a range of values (default os 0.35 to 0.95, with a step of 0.05) of the "gamma" parameters (see more details in the **ASCAT** R package help pages), create a **/home/me/my_project/EaCoN_results/SAMPLE1/ASCAT/ASCN/** subdirectory, in which other subdirectories will be created, one for each gamma value **/home/me/my_project/EaCoN_results/SAMPLE1/ASCAT/ASCN/gamma_0.xx/**. In each of those will be written :
  - _**SAMPLE1_ascat.ASCN.RDS**_ : contains the TCN and ASCN segmentations, and cellularity and ploidy data
  - _**SAMPLE1_ASCATprofile.png**_ : shows a graphical representation of the final segmented TCN and ASCN profiles
  - _**SAMPLE1_rawprofile.png**_ : shows a graphical representation of the raw (uncorrected) segmented TCN and ASCN profiles
  - _**SAMPLE1_Rorschach.clown.png**_ : shows a graphical representation of BAF vs L2R of probes, by chromosome, with a coloration corresponding to CN levels.
  - _**SAMPLE1_TCNvsL2R.png**_ : shows a graphical representation of the comparison of the TCN and L2R values of each segment, clustered by TCN level. This is usefull to identify some mistakes in the TCN modelization.
  - _**SAMPLE1_gamma0.xx.cn**_ : contains the TCN and ASCN segmentation results in a non-standard format derived the CBS format
  - _**SAMPLE1_gamma0.xx_model.txt**_ : Contains the ploidy, cellularity and model statistics

#### **HTML reporting**

- Endly, an annotated HTML report can be rendered with :

```R
Annotate.ff(RDS.file = "/home/me/my_project/EaCoN_results/SAMPLE1/ASCAT/L2R/SAMPLE1.EaCoN.ASPCF.RDS", author.name = "Me!")
```

- This will write the HTML report in **/home/me/my_project/EaCoN_results/SAMPLE1/ASCAT/L2R/**.

### **Batch mode (with multithreadng)**

All the steps described above in single sample mode can be run in batch mode, that is for multiple samples, even possibly with multithreading to process multiple at the same time. It simply consists into using different function with the same name but an added ".Batch" suffix. Those are just wrappers to the single-sample version of these.

#### **Raw data processing**

##### **Affymetrix OncoScan / OncoScan_CNV**

The **OS.Process.Batch** function replaces the *ATChannelCel*, *GCChannelCel* and *samplename* parameters by the *pairs.file* parameters, which consists in a tab-separated file with made of three columns with a header, and one sample per line :
- **ATChannelCel** : the path to the "A" OncoScan CEL file
- **GCChannelCel** : the path to the "C" OncoScan CEL file
- **SampleName** : the sample name to use

By default, the function will run all samples one by one, but multithreading can be set using the _nthread_ parameter with a value greater than 1. Beware not setting a value higher than the current number of available threads on your machine ! Please also remember that each new thread will use its own amount of RAM...

Here is a synthetic example with 4 samples :
- The *pairs.file* (stored as **/home/me/my_project/CEL/OS_pairs.txt**) :
>**ATChannelCel** **GCChannelCEL**  **SampleName**
>/home/me/my_project/CEL/SAMPLE1_OncoScan_CNV_A.CEL /home/me/my_project/CEL/SAMPLE1_OncoScan_CNV_A.CEL  SAMPLE1_OS
>/home/me/my_project/CEL/SAMPLE5_OncoScan_CNV_A.CEL /home/me/my_project/CEL/SAMPLE5_OncoScan_CNV_A.CEL  SAMPLE5_OS
>/home/me/my_project/CEL/SAMPLE6_OncoScan_CNV_A.CEL /home/me/my_project/CEL/SAMPLE6_OncoScan_CNV_A.CEL  SAMPLE6_OS
>/home/me/my_project/CEL/SAMPLE7_OncoScan_CNV_A.CEL /home/me/my_project/CEL/SAMPLE7_OncoScan_CNV_A.CEL  SAMPLE7_OS

- The command line (using 2 threads)
```R
OS.Process.Batch(pairs.file = "/home/me/my_project/CEL/OS_pairs.txt", nthread = 2)
```

##### **Affymetrix CytoScan 750k / CytoScan HD**

Same principle, but this time we have one column less and header changes a bit :
- **CEL** : The path to the CEL file
- **SampleName** : the sample name to use

Here is a synthetic example with 4 samples :
- The *CEL.list.file* (stored as **/home/me/my_project/CEL/CSHD_list.txt**) :
>**CEL**  **SampleName**
>/home/me/my_project/CEL/SAMPLE8_CytoScanHD_A.CEL SAMPLE8_CSHD
>/home/me/my_project/CEL/SAMPLE9_CytoScanHD_A.CEL SAMPLE9_CSHD
>/home/me/my_project/CEL/SAMPLE10_CytoScanHD_A.CEL SAMPLE10_CSHD
>/home/me/my_project/CEL/SAMPLE11_CytoScanHD_A.CEL SAMPLE11_CSHD

- The command line (using 2 threads)
```R
CS.Process.Batch(pairs.file = "/home/me/my_project/CEL/CSHD_list.txt", nthread = 2)
```

##### **Affymetrix GenomeWide SNP6**

Identical to CytoScan 750k / HD, but the function is named **SNP6.Process.Batch**.

##### **WES data**

Still the same principle with an external list file, with column names :
- **testBAM** : the path to the test BAM file
- **refBAM** : the path to the reference BAM file
- **SampleName** : the sample name to use

Here is a synthetic example with 4 samples :
- The *BAM.list.file* (stored as **/home/me/my_project/WES/BAM_list.txt**) :
>**testBAM**  **refBAM**  **SampleName**
>/home/me/my_project/WES/SAMPLE4_WES_hg19_Tumor.BAM /home/me/my_project/WES/SAMPLE4_WES_hg19_Tumor.BAM SAMPLE4_WES
>/home/me/my_project/WES/SAMPLE12_WES_hg19_Tumor.BAM /home/me/my_project/WES/SAMPLE12_WES_hg19_Tumor.BAM SAMPLE12_WES
>/home/me/my_project/WES/SAMPLE13_WES_hg19_Tumor.BAM /home/me/my_project/WES/SAMPLE13_WES_hg19_Tumor.BAM SAMPLE13_WES
>/home/me/my_project/WES/SAMPLE14_WES_hg19_Tumor.BAM /home/me/my_project/WES/SAMPLE14_WES_hg19_Tumor.BAM SAMPLE14_WES

- The binning, then normalization command lines (using 2 threads)
```R
WES.Bin.Batch(BAM.list.file = "/home/me/my_project/WES/BAM_list.txt", BINpack = "/home/me/my_project/EaCoN_results/SureSelect_v5_merged_sorted_hg19_b50.GC.rda", nthread = 2)
WES.Normalize.ff.Batch(BINpack = "/home/me/my_project/EaCoN_results/SureSelect_v5_merged_sorted_hg19_b50.GC.rda", nthread = 2)
```

You can notice that here we did not specify any RDS or list file to **WES.Normalize.ff.Batch**. This is because this fonction needs as its first argument _BIN.RDS.files_, a **list** of _"\_binned.RDS"_ files (generated at the former command line), and by default it will recursively search downwards the current working directory for any of these RDS files. You can of course design your own list of RDS files to process, if you know a bit of R.

#### **L2R & BAF Segmentation**

As for the **WES.Normalize.ff.Batch** function, the **Segment.ASCAT.ff.Batch** function needs as its first argument _RDS.files_, a list of _"\_processed.RDS"_ files (generated at the raw data processing step). Likewise, it will by default recursively search downwards for any compatible RDS file.

Here is a synthetic example that will segment our CytoScan HD samples (as defined by the _pattern_ below) :

```R
Segment.ASCAT.ff.Batch(RDS.files = list.files(path = getwd(), pattern = "_CSHD.*_processed.RDS$", full.names = TRUE, recursive = TRUE), smooth.k = 5, SER.pen = 20, nrf = 1.0, nthread = 2)
```

#### **Copy-number estimation**

Still the same, with the **ASCN.ASCAT.ff.Batch** :

```R
ASCN.ASCAT.ff.Batch(RDS.files = list.files(path = getwd(), pattern = "_CSHD.*_EaCoN.ASPCF.RDS$", full.names = TRUE, recursive = TRUE), nthread = 2)
```

#### **HTML reporting**

And here again with the **Annotate.ff.Batch** :

```R
Annotate.ff.Batch(RDS.files = list.files(path = getwd(), pattern = "_CSHD.*_EaCoN.ASPCF.RDS$", full.names = TRUE, recursive = TRUE), author.name = "Me!")
```

### **Piped**

EaCoN has been implemented in a way that one can also choose to launch the full workflow in a single command line for a single sample, using pipes from the [magrittr](https://cran.r-project.org/web/packages/magrittr/vignettes/magrittr.html) package. However, this is not recommended as default use : even though EaCoN is provided with recommandations that should fit most case, the user may have to deal with particular profiles that would require parameter tweaking, which is not possible in piped mode...
Here is an example :

```R
setwd("/home/me/my_project/EaCoN_results")
require(EaCoN)
require(magrittr)
OS.Process(ATChannelCel = "/home/me/my_project/CEL/SAMPLE1_OncoScan_CNV_A.CEL", GCChannelCel = "/home/me/my_project/CEL/SAMPLE1_OncoScan_CNV_C.CEL", samplename = "SAMPLE1_OS", return.data = TRUE) %>% Segment.ASCAT(out.dir = "/home/me/my_project/EaCoN_results/SAMPLE1_OS", return.data = TRUE) %T>% Annotate(out.dir = "/home/me/my_project/EaCoN_results/SAMPLE1_OS/ASCAT/L2R") %>% ASCN.ASCAT(out.dir = "/home/me/my_project/EaCoN_results/SAMPLE1_OS")
```

### **Conclusion on usage**

- So, as you may have already noticed, most functions exist in EaCoN in three flavors :
  - The **"directly"** named : that accept a R _data_ object as first parameter, designed to be used in pipe mode.
  - The **".ff"** named : that require a file as first parameter, designed to be used in step-by-step mode.
  - The **".ff.Batch"** named : that require a list as first parameter, designed to be used in step-by-step mode, but running multiple samples for each step, possibly with multithreading.

---

## **GUIDELINES**

- For each step, default values for each data source already correspond to recommendations. However, for the common **segmentation** step (Segment.ASCAT), adaptation to the data source is recommended, by changing few parameters :

SOURCE | SER.pen | smooth.k | nrf | BAF.filter
--- | --- | --- | --- | ---
OncoScan | `40` *(default)* | `NULL` *(default)* | `0.5` *(default)* | `0.9`
CytoScan HD | `20` | `5` | `1.0` | `0.75` *(default)*
SNP6 | `60` | `5` | `0.25` | `0.75` *(default)*
WES | `NULL` | `5` | `0.5` *(default)* | `0.75` *(default)*

- All the functions depicted above have other non-described parameters. As the above recommandations should do the trick in most case, they certainly won't fit all. To adjust parameters more finely, I suggest to refer to the R help pages for corresponding functions.


---

## **AUTHORS & CONTACT**
 - Bastien Job (bastien.job@inserm.fr, bastien.job@gustaveroussy.fr) : developer and maintainer
 - Thibault Dayris (thibault.dayris@gustaveroussy.fr) : tester, user

