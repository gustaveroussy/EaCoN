# **EaCoN**

_**Ea**sy **Co**py **N**umber !_

---

## **DESCRIPTION**

EaCoN aims to be an all-packed in,  user-friendly solution to perform relative or absolute copy-number analysis for multiple sources of data, with three different segmenters available (and corresponding three copy-number modelization methods).
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
  - L2R and BAF normalization from raw data using
    - **[Affymetrix Power Tools](https://www.thermofisher.com/fr/fr/home/life-science/microarray-analysis/microarray-analysis-partners-programs/affymetrix-developers-network/affymetrix-power-tools.html)** (binaries, embedded) for Affymetrix arrays and the **[rawcopy](http://rawcopy.org/)** BAF normalization methods for CytoScan and SNP6 arrays
    - internal methods for WES data
  - L2R and BAF bivariate segmentation using **ASCAT** _(recommended)_, **FACETS** or **SEQUENZA**
  - L2R profile centralization
  - L2R profile CNA calling
  - Generation of results as tables and plots
  - Generation of a portable and interactive HTML report (for hs, mm and rn species at the moment)
  - Total (TCN) and allele-specific copy number (ASCN), ploidy and cellularity estimations using **ASCAT** _(recommended)_, **FACETS** or **SEQUENZA** for both WES and arrays (except for SNP6 arrays for the two latter)
- Required annotations for Affymetrix data provided (via additional packages)
- Pre-computed normalization (GC%, wave) tracks provided for arrays (either)
- Few functions to call
- Few parameters
- Recommended guidelines provided for supported microarray designs, as well as for WES

---

## **NOTES**

* Support for R v4.0.x has been added thanks to suggestions by [ShenWei-wei](https://github.com/ShenWei-wei), but has not been extensively tested yet (few successful tests using CytoScanHD CEL files only).
* Support for R v3.x should remain, now.
* Links for dependency packages (for Affymetrix microarrays) have been moved to [Zenodo](https://zenodo.org).

## **QUICK NEWS**

### **2021-06-18 : v0.3.6-1 _(SweetSummerSweat)_ is out !**
* CORR : Modified the correction to support R4 that er... broke compatibility with R3. Now the code evaluates the R major version (version$major) and applies ShenWei-wei's correction if >=4.

### **2021-05-23 : v0.3.6 _(Barolo)_ is out !**

* CORR : Added correction suggested by [ShenWei-wei](https://github.com/ShenWei-wei) to allow support for R v4.0.x
* MOD : Updated README (fixed some links)
* CORR : Updated outdated links to dependences hosted on nextcloud.gustaveroussy.fr

### **2020-08-17 : v0.3.5 _(CloudyMonday)_ is out !**

* CORR : Segment.*() : Added a patch to handle the NA behavior in copynumber::winsorize (error raised by new handling of NA values in runmed()). The patch consists in applying winsorization on non-NA values only (whereas all values were transmitted in earlier versions).
* CORR : WES.Bin() : Better handling of a possible desynch in chr names (when a canonical chr had no remaining values, its level was kept. This raised a rare error).
* MOD : Many funcs : Fixed calls to the "%do%" and "%dopar" operators without loading it.

### **2018-12-10 : v0.3.4-1 _(PostRoscovite)_ is out !**

* CORR : SNP6.Process(), CSHD.Process() : Edited code to handle changes in the rcnorm package, to discard the "chromosomes" package dependency.
* MOD : Removed some dependencies (already called by other dependencies, like 'copynumber' from 'sequenza') to make installation easier and more convenient.
* MOD : Edited the README.md (rewrote the INSTALL section)
* MOD : WES.Bin() : Added support for BAI files that have the exact same rootname as BAM files (instead of rootname.bam.bai only)
* MOD : Added "call. = FALSE" top all stop() calls
* ADD : WES.Bin() : Added raw read depth plots (to control putative TEST / REF inversion, or sex mismatch)

### **2018-10-30 : v0.3.4 _(Papy60)_ is out !**

- Now FACETS can be used on OncoScan, OncoScan_CNV, CytoScan750k and CytoScanHD arrays (still not on SNP6, though).
- Some bugs corrected, too.
- Sequenza CN output fixed.
- Changed the formula for the "most width" ploidy version, so that a value of 0 can't be returned.
- Paving the way to handle non-canonical genomes in the report (not active yet, though).
- Small code changes to adapt to "non-completely covered" genomes (ie, few chromosomes for toy datasets, by example).
- Harmonized the penalty parameter for all segmenters (now just "penalty" rather than "ASCAT.pen", "SEQUENZA.pen" or "FACETS.pen")
- Officially dropping support for sequenza with SNP6 arrays as it leads to a huge RAM consumption by copynumber::aspcf, due to the few probes covering both L2R and BAF.
- Included chromosomes objects in package for hs, mm and rn, to avoid dependency to another non-public sourceable package.

### **2018-10-02 : v0.3.3-1 _(LittleWomanNoCry)_ is out !**

- Multiple little bugs correction, tweaks.

### **2018-09-12 : v0.3.3 _(Trinity)_ is out !**

- Now EaCoN supports the [**SEQUENZA**](https://www.ncbi.nlm.nih.gov/pubmed/25319062) segmenter and copy number estimator _(except for Affymatrix SNP6)_ ! This makes 3 different segmenters available for quick. Sequenza uses the same bivariate segmentation algorithm as ASCAT, _(PCF : piecewise constant curve fitting)_, but in a different implementation (sequenza relies on the _copynumber_ package, ASCAT has its own built-in implementation), so expect similar but not identical results !
- A few bugs solved.

### **2018-08-08 : v0.3.2 _(PapeMamiePichine)_ is out !**

- Now EaCoN supports the [**FACETS**](https://www.ncbi.nlm.nih.gov/pubmed/27270079) segmenter and copy number estimator ! This segmenter extends the famous **CBS** _(Circular Binary Segmentation)_ algorithm by making it compatible with bi-variate segmentation (thus, using both the L2R and BAF signals). However, please note that **FACETS is only available for WES data**.
- A few bugs solved.

---

## **INSTALLATION**

### **CORE**
-  Please first install the **_devtools_** package that will allow installing packages from _github_ :

  ``` r
  install.packages('devtools')
  ```

- Then install **_ASCAT_** and **_FACETS_** from github :

  ``` r
  devtools::install_github("Crick-CancerGenomics/ascat/ASCAT")
  devtools::install_github("mskcc/facets")
  ```

- Then install required **_BioConductor_** packages :

  ``` r
  ## try using http:// if https:// URLs are not supported
  if(!installed.packages('BiocManager')) install.packages('BiocManager')
  BiocManager::install(c("affxparser", "Biostrings", "aroma.light", "BSgenome", "copynumber", "GenomicRanges", "limma", "rhdf5", "sequenza"))
  ```

- Then install **_EaCoN_** from github !

  ``` r
  ## Install the most recent STABLE version (@master)
  devtools::install_github("gustaveroussy/EaCoN")
  ```

### **MICROARRAY-SPECIFIC**

While the current EaCoN package is the core of the process and will straightly work for WES data, multiple other packages are needed to properly handle Affymetrix microarray : APT (affymetrix power tools), designs and corresponding annotations (genome build, Affymetrix annotation databases) ; others are required for the (re)normalization, especially pre-computed GC% or Wavetracks.

#### **ALL AFFYMETRIX MICROARRAYS**

- The **_affy.CN.norm_** package provides pre-computed GC% and wave-effect (re)normalization datasets for all compatible Affymetrix designs, for both NA33/NA35 (hg19) and NA36 (hg38) human genome builds. Install from remote URL :

  ``` r
  install.packages("https://zenodo.org/record/5494853/files/affy.CN.norm.data_0.1.2.tar.gz", repos = NULL, type = "source")
  ```

#### **ONCOSCAN FAMILY (OncoScan / OncoScan_CNV)**

- First, install embedded APT tool from github :

  ``` r
  devtools::install_github("gustaveroussy/apt.oncoscan.2.4.0")
  ```

- Then install annotations from remote URL :
  - For the **NA33 (hg19)** build :
    - For the **OncoScan** design :

      ``` r
      install.packages("https://zenodo.org/record/5494853/files/OncoScan.na33.r4_0.1.0.tar.gz", repos = NULL, type = "source")
      ```

    - For the **OncoScan_CNV** design :

      ``` r
      install.packages("https://zenodo.org/record/5494853/files/OncoScanCNV.na33.r2_0.1.0.tar.gz", repos = NULL, type = "source")
      ```

  - For the **NA36 (hg38)** build :
    - For the **OncoScan** design :

      ``` r
      install.packages("https://zenodo.org/record/5494853/files/OncoScan.na36.r1_0.1.0.tar.gz", repos = NULL, type = "source")
      ```

    - For the **OncoScan_CNV** design :

      ``` r
      install.packages("https://zenodo.org/record/5494853/files/OncoScanCNV.na36.r1_0.1.0.tar.gz", repos = NULL, type = "source")
      ```

#### **CYTOSCAN FAMILY (CytoScan 750k / CytoScan HD)**

- First, install the embedded APT tool from github :

  ``` r
  devtools::install_github("gustaveroussy/apt.cytoscan.2.4.0")
  ```

- Then install annotations from remote URL :
  - For the **NA33 (hg19)** build :
    - For the **CytoScan 750k** design :

      ``` r
      install.packages("https://zenodo.org/record/5494853/files/CytoScan750K.Array.na33.r4_0.1.0.tar.gz", repos = NULL, type = "source")
      ```
    - For the **CytoScan HD** design :

      ``` r
      install.packages("https://zenodo.org/record/5494853/files/CytoScanHD.Array.na33.r4_0.1.0.tar.gz", repos = NULL, type = "source")
      ```

  - For the **NA36 (hg38)** build :
    - For the **CytoScan 750k** design :

      ``` r
      install.packages(https://zenodo.org/record/5494853/files/CytoScan750K.Array.na36.r1_0.1.0.tar.gz", repos = NULL, type = "source")
      ```
    - For the **CytoScan HD** design :

      ``` r
      install.packages("https://zenodo.org/record/5494853/files/CytoScanHD.Array.na36.r1_0.1.0.tar.gz", repos = NULL, type = "source")
      ```

- Lastly, install the **_rcnorm_** package to perform BAF normalization for the CytoScan family of arrays :

  ``` r
  install.packages("https://zenodo.org/record/5494853/files/rcnorm_0.1.5.tar.gz", repos = NULL, type = "source")
  ```

#### **GENOMEWIDE SNP6**

- First, install the embedded APT tool from github :

  ``` r
  devtools::install_github("gustaveroussy/apt.snp6.1.20.0")
  ```

- Then install annotations from remote URL (There is no other build available than **NA35 (hg19)**) :

  ``` r
  install.packages("https://zenodo.org/record/5494853/files/GenomeWideSNP.6.na35.r1_0.1.0.tar.gz", repos = NULL, type = "source")
  ```

- Lastly, install the **_rcnorm_** package to perform BAF normalization for SNP6 arrays **(if not already installed at the CytoScan step!)** :

  ``` r
  install.packages("https://zenodo.org/record/5494853/files/rcnorm_0.1.5.tar.gz", repos = NULL, type = "source")
  ```

### **GENOMES**

- EaCoN requires a genome as available thanks to the [BSgenome](https://bioconductor.org/packages/release/bioc/html/BSgenome.html) package, available at BioConductor. To check which genomes are publicly availabe at BioConductor :

  ``` r
  if(!'BiocManager' %in% installed.packages()) install.packages('BiocManager')
  BiocManager::install('BSgenome')
  BSgenome::available.genomes()
  ```

- To check genome(s) installed in your R library :

  ``` r
  BSgenome::installed.genomes()
  ```

- For Affymetrix microarrays, you need to install these human genomes depending on which annotation package you want to use :

  ``` r
  if(!'BiocManager' %in% installed.packages()) install.packages('BiocManager')

  ## To support NA33 / NA35 annotations (hg19)
  BiocManager::install('BSgenome.Hsapiens.UCSC.hg19')

  ## To support NA36 annotations (hg38)
  BiocManager::install('BSgenome.Hsapiens.UCSC.hg38')
  ```

- For TCGA WES data, you will need the **hs37d5** genome _(a variation of the hg19 build used in the 1000 Genomes project)_

  ``` r
  if(!'BiocManager' %in% installed.packages()) install.packages('BiocManager')
  BiocManager::install("BSgenome.Hsapiens.1000genomes.hs37d5")
  ```

- If your favorite genome is not available, **it is possible to build your own** !
  - Simply download locally any public genome package (by example [BSgenome.Hsapiens.UCSC.hg19](https://bioconductor.org/packages/BSgenome.Hsapiens.UCSC.hg19/))
  - Uncompress it locally (MS Windows users can use [7-zip](https://www.7-zip.org/))
  - Rename the decompressed folder to your genome name
  - Edit all possible files (**_DESCRIPTION_**, **_NAMESPACE_**, **_R/zzz.R_**, **_man/package.Rd_**) to insert all required information relative to your own genome
  - Replace the **_inst/extdata/single_sequenes.2bit_** file by your genome sequence in the [2bit](http://genome.ucsc.edu/FAQ/FAQformat.html#format7) format (which can be converted to from a regular fasta file thanks to converters developed by the [UCSC](https://genome.ucsc.edu/goldenpath/help/twoBit.html))
  - Re-compress the modified and renamed directory to a .tar.gz (Windows users : 7-zip can do that, too)
  - Install !

---

## **INPUT**

- Raw data :
  - For Affymetrix microarrays : the **CEL** file(s), fresh out of the Affymetrix Scanner
  - For WES data :
    - The aligned reads of the test and reference samples, in two different **BAM** files
    - The **capture BED** file (the file containing the genomic regions targeted by the exome capture kit used in the wet experiment)

---

## **USAGE**

The full workflow is decomposed into a few different functions, which roughly correspond to these steps :  

```
normalization -> segmentation +-> reporting
                              |
                              +-> copy-number estimation
```

EaCoN allows different ways of running the full workflow : considering the analysis of a single sample, you can either run each step independently and write, then load the intermediate results, or you can _**pipe**_ all steps in a single line of code. You can also run the step-by-step approach on multiple samples in a row, even possibly at the same time using multithreading, using a batch mode.

### **Step by step mode**

First, under R, load EaCoN and choose a directory for writing results, for example : **/home/project/EaCoN_results**

  ``` r
  require(EaCoN)
  setwd("/home/project/EaCoN_results")
  ```

#### **Raw data processing**

##### **Affymetrix OncoScan / OncoScan_CNV**

- Let's say we have a pair of OncoScan_CNV CEL files to analyse in a **/home/project/CEL/** directory (Affymetrix OncoScan experiments have 2 arrays for a single experiment, thus a pair) :

  ``` r
  OS.Process(ATChannelCel = "/home/project/CEL/S1_OncoScan_CNV_A.CEL", GCChannelCel = "/home/project/CEL/S1_OncoScan_CNV_C.CEL", samplename = "S1_OS")
  ```

- This will perform the normalization step, create a **/home/project/EaCoN_results/S1_OS/** subdirectory and write 5 files in it :
  - _**S1_OS_OncoScan_CNV_hg19_processed.RDS**_ : contains the normalized data to be provided to the next step
  - _**S1_OS_OncoScan_CNV_hg19_rawplot.png**_ : shows a graphical representation of the normalized L2R and BAF data
  - _**S1_OS_2.4.0_na33.r2.paircheck.txt**_ : gives some statistics to evaluate the probability that the two "A" "C" CEL files effectively belong to the same individual (generated by APT)
  - _**S1_OS_2.4.0_na33.r2.qc.txt**_ : some quality metrics of the arrays and profiles (generated by APT)
  - _**S1_OS_2.4.0_na33.r2.log**_ : a simple log file (generated by APT)

##### **Affymetrix CytoScan 750k / CytoScan HD**

- This is identical to OncoScan, except that we have a single CEL file, and we just have to change the function name :

  ``` r
  CS.Process(CEL = "/home/project/CEL/S2_CytoScanHD.CEL", samplename = "S2_CSHD")
  ```
- The same output files will be generated (except for the "paircheck" file, obviously)

##### **Affymetrix GenomeWide SNP6**

- Here again, this is identical to OncoScan, except that we have a single CEL file, and we just have to change the function name :

  ```r
  SNP6.Process(CEL = "/home/project/CEL/S3_GenomeWide_snp.6.CEL", samplename = "S3_SNP6")
  ```
- Again, the same output files will be generated (except for the "paircheck" file, obviously)

##### **WES data**

- This time it is quite different : the processing will be performed in three steps :
  - **First**, we will use the capture BED (A text file containing the positions of the captured regions, usualy provided by the capture kit manufacturer), choose a genome version corresponding to our aligned BAM files, and choose a window size for the future binning of the data. These will be used to generate what we call a "BINpack", a set of pre-computed tracks containing the bin positions and corresponding GC% values. Several tracks will be computed corresponding to different levels of elongation of the bin positions. In the example below, we used the BED corresponding to Agilent SureSelect v5 capture kit, a bin size of 50 nt, and chose the human hg19 genome build.

    ``` r
    BINpack.Maker(bed.file = "/home/project/WES/SureSelect_v5.bed", bin.size = 50, genome.pkg = "BSgenome.Hsapiens.UCSC.hg19")
    ```

    - This will generate a "BINpack" (with a ".rda" extension) to be used in the next normalization steps : **/home/project/EaCoN_results/SureSelect_v5_merged_sorted_hg19_b50.GC.rda**

    - **PLEASE NOTE THAT THIS STEP IS SAMPLE-INDEPENDENT, THUS NEEDS TO BE PERFORMED AGAIN ONLY IF YOU CHANGE EITHER THE CAPTURE BED, THE BIN SIZE OR THE GENOME BUILD.** Thus, the generated BINpack can be used for any other sample in the same conditions.

  - **Second**, the WES data will be binned using the generated BINpack. We need three files as input :
    - The aligned reads for the test sample (usualy in cancer, the patient tumor), in BAM format
    - The aligned reads for the reference sample (patient normal), in BAM format too
    - The BINpack itself

      ``` r
      WES.Bin(testBAM = "/home/project/WES/S4_WES_hg19_Tumor.BAM", refBAM = "/home/project/WES/S4_WES_hg19_Tumor.BAM", BINpack = "/home/project/EaCoN_results/SureSelect_v5_merged_sorted_hg19_b50.GC.rda", samplename = "S4_WES")
      ```

    - This will generate a **/home/project/EaCoN_results/S4_WES/** subdirectory which contains :
      - _**S4_WES_hg19_b50_binned.RDS**_ : contains the binned data to be provided to the next step
      - _**S4_WES_hg19_b50_coverage.png**_ : shows a graphical representation of the proportion of the capture bed regions covered at different coverage levels
      - _**S4_WES_hg19_b50_coverage.txt**_ : contains the numerical values corresponding to the coverage plot

  - **Third**, now that the data have been binned, the normalization step can be performed :

    ``` r
    WES.Normalize.ff(BIN.RDS.file = "/home/project/EaCoN_results/S4_WES/S4_WES_hg19_b50_binned.RDS", BINpack = "/home/project/EaCoN_results/SureSelect_v5_merged_sorted_hg19_b50.GC.rda")
    ```

    - Analogously to what was described for Affymetrix microarrays, a **/home/project/EaCoN_results/S4_WES/ASCAT/L2R/** subdirectory will be created, containing the same type of files described earlier.


#### **L2R & BAF Segmentation**

- Now that we described the normalization process specific to each type of source data, we can segment it. The good news is that it's the very same step for each source type, one just have to pass the **\_processed.RDS** normalized file. Here is an example with the one obtained for OncoScan data, using the **ASCAT** segmenter :

  ``` r
  Segment.ff(RDS.file = "/home/me/my_project/EaCoN_results/SAMPLE1/S1_OncoScan_CNV_hg19_processed.RDS", segmenter = "ASCAT")
  ```

  - This will perform the segmentation, centralization and calling steps, create a **/home/project/EaCoN_results/S1/ASCAT/L2R/** subdirectory and write multiple files in it :
    - _**S1.SEG.ASCAT.RDS**_ : contains the segmented data to be provided to the optional next step(s)
    - _**S1.SEG.ASCAT.png**_ : shows a graphical representation of the segmented, centered and called L2R and BAF data
    - _**S1.Rorschach.png**_ : shows a graphical representation of BAF vs L2R of probes, by chromosome
    - _**S1.Cut.cbs**_ : contains the L2R segmentation results in the standard CBS format. "Cut" means that L2R value for the segments called as normal were set to a value of 0.0
    - _**S1.NoCut.cbs**_ : Same as above, but the segments considered as normal have their L2R value intact
    - _**S1.SegmentedBAF.txt**_ : contains the BAF segmentation results

- To perform the same using the **FACETS** segmenter, just change the value of the _segmenter_ parameter !

- I suppose you guessed how to do the same with **SEQUENZA**, right ? ;)

#### **Copy-number estimation**

- Then an estimation of the total and allele-specific copy-number profiles, as well as global ploidy and sample cellularity can be estimated. Here is an example using the **ASCAT** ASCN estimation, from a RDS generated by the Segment.ASCAT() function :

  ``` r
  ASCN.ff(RDS.file = "/home/me/my_project/EaCoN_results/SAMPLE1/ASCAT/L2R/SAMPLE1.ASCAT.RDS")
  ```

  - This will perform these estimations for a range of values (default os 0.35 to 0.95, with a step of 0.05) of the "gamma" parameters (see more details in the **ASCAT** R package help pages), create a **/home/project/EaCoN_results/S1/ASCAT/ASCN/** subdirectory, in which other subdirectories will be created, one for each gamma value **/home/project/EaCoN_results/S1/ASCAT/ASCN/gamma_0.xx/**. In each of those will be written :
    - _**S1_ascat.ASCN.RDS**_ : contains the TCN and ASCN segmentations, and cellularity and ploidy data
    - _**S1_ASCATprofile.png**_ : shows a graphical representation of the final segmented TCN and ASCN profiles
    - _**S1_rawprofile.png**_ : shows a graphical representation of the raw (uncorrected) segmented TCN and ASCN profiles
    - _**S1_Rorschach.clown.png**_ : shows a graphical representation of BAF vs L2R of probes, by chromosome, with a coloration corresponding to CN levels.
    - _**S1_TCNvsL2R.png**_ : shows a graphical representation of the comparison of the TCN and L2R values of each segment, clustered by TCN level. This is usefull to identify mistakes in the TCN modelization.
    - _**S1_gamma0.xx.cn**_ : contains the TCN and ASCN segmentation results in a non-standard format derived the CBS format
    - _**S1_gamma0.xx_model.txt**_ : Contains the ploidy, cellularity and model statistics

- To perform the same using the **FACETS** or **SEQUENZA** estimator, just use a RDS generated with Segment.FACETS() or Segment.SEQUENZA(), respectively (or their ".ff" equivalent).

#### **HTML reporting**

- Finally, an annotated HTML report can be rendered with :

  ``` r
  Annotate.ff(RDS.file = "/home/project/EaCoN_results/S1/ASCAT/L2R/S1.EaCoN.ASPCF.RDS", author.name = "Me!")
  ```

  - This will write the HTML report in **/home/project/EaCoN_results/S1/ASCAT/L2R/**.

### **Batch mode (with multithreadng)**

All the steps described above in single sample mode can be run in batch mode, that is for multiple samples, possibly combined with multithreading to process multiple samples in parallel. It simply consists into using different functions with the same name but an added ".Batch" suffix. Those are just wrappers to the single-sample version of the functions.

#### **Raw data processing**

##### **Affymetrix OncoScan / OncoScan_CNV**

The **OS.Process.Batch** function replaces the *ATChannelCel*, *GCChannelCel* and *samplename* parameters by the *pairs.file* parameters, which consists in a tab-separated file with made of three columns with a header, and one sample per line :
- **ATChannelCel** : the path to the "A" OncoScan CEL file
- **GCChannelCel** : the path to the "C" OncoScan CEL file
- **SampleName** : the sample name to use

By default, the function will run all samples one by one, but multithreading can be set using the _nthread_ parameter with a value greater than 1. Beware not setting a value higher than the current number of available threads on your machine ! Please also remember that each new thread will use its own amount of RAM...

Here is a synthetic example with 4 samples :
- The *pairs.file* (stored as **/home/project/CEL/OS_pairs.txt**) :

ATChannelCel | GCChannelCEL | SampleName
--- | --- | ---
/home/project/CEL/S1_OncoScan_CNV_A.CEL | /home/project/CEL/S1_OncoScan_CNV_C.CEL |  S1_OS
/home/project/CEL/S5_OncoScan_CNV_A.CEL | /home/project/CEL/S5_OncoScan_CNV_C.CEL |  S5_OS
/home/project/CEL/S6_OncoScan_CNV_A.CEL | /home/project/CEL/S6_OncoScan_CNV_C.CEL |  S6_OS
/home/project/CEL/S7_OncoScan_CNV_A.CEL | /home/project/CEL/S7_OncoScan_CNV_C.CEL |  S7_OS

- The command line (using 2 threads)

  ``` r
  OS.Process.Batch(pairs.file = "/home/project/CEL/OS_pairs.txt", nthread = 2)
  ```

##### **Affymetrix CytoScan 750k / CytoScan HD**

Same principle, but this time we have one column less and header changes a bit :
- **CEL** : The path to the CEL file
- **SampleName** : the sample name to use

Here is a synthetic example with 4 samples :
- The *CEL.list.file* (stored as **/home/project/CEL/CSHD_list.txt**) :

CEL | SampleName
--- | ---
/home/project/CEL/S8_CytoScanHD.CEL | S8_CSHD
/home/project/CEL/S9_CytoScanHD.CEL | S9_CSHD
/home/project/CEL/S10_CytoScanHD.CEL | S10_CSHD
/home/project/CEL/S11_CytoScanHD.CEL | S11_CSHD

- The command line (using 2 threads)

  ``` r
  CS.Process.Batch(pairs.file = "/home/project/CEL/CSHD_list.txt", nthread = 2)
  ```

##### **Affymetrix GenomeWide SNP6**

Identical to CytoScan 750k / HD, but the function is named **SNP6.Process.Batch**.

##### **WES data**

Still the same principle with an external list file, with column names :
- **testBAM** : the path to the test BAM file
- **refBAM** : the path to the reference BAM file
- **SampleName** : the sample name to use

Here is a synthetic example with 4 samples :
- The *BAM.list.file* (stored as **/home/project/WES/BAM_list.txt**) :

testBAM | refBAM | SampleName
--- | --- | ---
/home/project/WES/S4_WES_hg19_Tumor.BAM | /home/project/WES/S4_WES_hg19_Normal.BAM | S4_WES
/home/project/WES/S12_WES_hg19_Tumor.BAM | /home/project/WES/S12_WES_hg19_Normal.BAM | S12_WES
/home/project/WES/S13_WES_hg19_Tumor.BAM | /home/project/WES/S13_WES_hg19_Normal.BAM | S13_WES
/home/project/WES/S14_WES_hg19_Tumor.BAM | /home/project/WES/S14_WES_hg19_Normal.BAM | S14_WES

- Binning, then normalization command lines (using 2 threads)

  ``` r
  WES.Bin.Batch(BAM.list.file = "/home/project/WES/BAM_list.txt", BINpack = "/home/project/EaCoN_results/SureSelect_v5_merged_sorted_hg19_b50.GC.rda", nthread = 2)
  WES.Normalize.ff.Batch(BINpack = "/home/project/EaCoN_results/SureSelect_v5_merged_sorted_hg19_b50.GC.rda", nthread = 2)
  ```

Note that here we did not specify any RDS or list file to **WES.Normalize.ff.Batch**. This is because this fonction needs as its first argument _BIN.RDS.files_, a **list** of _"\_binned.RDS"_ files (generated at the former command line), and by default it will recursively search downwards the current working directory for any of these RDS files. You can of course design your own list of RDS files to process, if you know a bit of R.

#### **L2R & BAF Segmentation**

As for the **WES.Normalize.ff.Batch** function, the **Segment.ff.Batch** function needs as its first argument _RDS.files_, a list of _"\_processed.RDS"_ files (generated at the raw data processing step). Likewise, it will by default recursively search downwards for any compatible RDS file.

Here is a synthetic example that will segment our CytoScan HD samples (as defined by the _pattern_ below) using ASCAT :

  ```r
  Segment.ff.Batch(RDS.files = list.files(path = getwd(), pattern = ".*_processed.RDS$", full.names = TRUE, recursive = TRUE), segmenter = "ASCAT", smooth.k = 5, SER.pen = 20, nrf = 1.0, nthread = 2)
  ```

- To perform the same using the **FACETS** segmenter, just change the value of the _segmenter_ parameter.

- I suppose you guessed how to do the same with **SEQUENZA**, right ? ;)

#### **Copy-number estimation**

Still the same, with the **ASCN.ff.Batch** :

  ``` r
  ASCN.ff.Batch(RDS.files = list.files(path = getwd(), pattern = "SEG\\..*\\.RDS$", full.names = TRUE, recursive = TRUE), nthread = 2)
  ```

#### **HTML reporting**

And here again with the **Annotate.ff.Batch** :

  ``` r
  Annotate.ff.Batch(RDS.files = list.files(path = getwd(), pattern = "SEG\\..*\\.RDS$", full.names = TRUE, recursive = TRUE), author.name = "Me!")
  ```

### **Piped**

EaCoN has been implemented in such a way that one can also opt to launch the full workflow in a single command line for a single sample, using pipes from the [magrittr](https://cran.r-project.org/web/packages/magrittr/vignettes/magrittr.html) package. However, this is not recommended as a default use : even though EaCoN is provided with recommendations that should fit most cases, users may have to deal with particular profiles requiring parameter tweaking, which is not possible in piped mode...
Here is an example using ASCAT :

  ```r
  samplename <- "SAMPLE1_OS"
  workdir <- "/home/me/my_project/EaCoN_results"
  setwd(workdir)
  require(EaCoN)
  require(magrittr)

  OS.Process(ATChannelCel = "/home/me/my_project/CEL/SAMPLE1_OncoScan_CNV_A.CEL", GCChannelCel = "/home/me/my_project/CEL/SAMPLE1_OncoScan_CNV_C.CEL", samplename = samplename, return.data = TRUE) %>% Segment(out.dir = paste0(workdir, "/", samplename), segmenter = "ASCAT", return.data = TRUE) %T>% Annotate(out.dir = paste0(workdir, "/", samplename, "/ASCAT/L2R")) %>% ASCN.ASCAT(out.dir = paste0(workdir, "/", samplename))
  ```

### **Conclusion on usage**

- So, as you may have already noticed, most functions exist in EaCoN in three flavors :
  - The **"directly"** named : that accept a R _data_ object as first parameter, designed for usein pipe mode.
  - The **".ff"** named ("from file") : that require a file as first parameter, designed for use in step-by-step mode.
  - The **".ff.Batch"** named ("from file, in batch") : that require a list as first parameter, designed for use in step-by-step mode, but running multiple samples for each step, possibly with multithreading.

---

## **GUIDELINES**

### **Segmentation**

- For each step, default values for each data source already correspond to recommendations. However, for the common **segmentation** step, adaptation to the data source is recommended, by changing a few parameters :

SOURCE | SER.pen | smooth.k | nrf | BAF.filter
--- | --- | --- | --- | ---
OncoScan | `40` *(default)* | `NULL` *(default)* | `0.5` *(default)* | `0.9`
CytoScan HD | `20` | `5` | `1.0` | `0.75` *(default)*
SNP6 | `60` | `5` | `0.25` | `0.75` *(default)*
WES | `2` to `10` | `5` | `0.5` *(default)* to `1` | `0.75` *(default)*

- The FACETS segmenter cannot currently be used on SNP6 data (due to missing normalized A and B signals).

- The SEQUENZA segmenter SHOULD NOT be used with SNP6 microarrays (it theoretically can, but requires huge amounts of RAM, ie more than 32 GB). This may halt / swap your computer !

- For WES data, any sorted BAM should work, but we recommend using BAMs for which duplicates were marked/removed (samtools markdup, Picard MarkDuplicates, etc...), for higher quality results.

### NOTES

- All the functions depicted above have other parameters not described here. As the above recommendations should do the trick in most cases, they certainly won't fit all. To adjust parameters more finely, please refer to the R help pages of corresponding functions.

---

## **AUTHORS & CONTACT**

 - Bastien Job (bastien.job@inserm.fr, bastien.job@gustaveroussy.fr) : developer and maintainer
 - Thibault Dayris (thibault.dayris@gustaveroussy.fr) : tester, user
