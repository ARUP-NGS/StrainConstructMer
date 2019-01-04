---

<img src="docs/ARUP_horz_2c_pos.jpg" width="50%">

---

# ARUP StrainTypeMer
__ARUP StrainTypeMer__ is a rapid program created to compare genomes of two or more samples. The primary use case is for 
epidemiological analysis. The program analyzes NGS data created by Whole Genome Fragmentation protocol. This version 
of __StrainTypeMer__ works as two Ion Torrent Plugins. __ARUP StrainTypeMer__ is an NGS analysis designed to replace 
Pulsed Field Gel Electrophoresis (PFGE).  While most epidemiological NGS methods require reference alignment
__ARUP StrainTypeMer__ is designed to run reference free, thus creating a universal analysis method to compared strains
unhinder by the amount of identity a sample has with the reference sequence. Because __ARUP StrainTypeMer__ is
reference-free it creates a universal strain comparison method that will work across many species and genera of 
organisms.  

# StrainConstructMer
__StrainConstructMer__ is the first plugin in the StrainTypeMer analysis. This plugin processes a sample and 
transform the data so that comparisons can be performed with the __StrainCompareMer__ plugin.

### Features
* Calculates coverage
* Determines the genome size of the organism
    * Determines if genome size is expected value
* Identifies the bacteria based in [NCBI 16S rRNA RefSeq references](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA33175)
* Determines MLST type using [PUBMLST](https://pubmlst.org/])
* Identify antibiotic resistant genes present in [NCBI AMR references](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA313047)
* Performs Quality checks to determine is sample can be compared in __StrainCompareMer__

# StrainCompareMer
__StrainCompareMerMer__ is the second plugin in the StrainTypeMer analysis. This plugin processes requires a CSV input 
file that indicates the strains to be processed.  THe strains must exist in the database and the plugin must be able to 
locate the files associated with the strain.  The database location set on the global configuration page must be 
identical to the location set on __StrainConstructMer__ (default location `/results/plugins/scratch/`). 

### Features
* Includes multiple comparisons
    * Full Genome comparison
    * Full Genome comparison with select references
    * Core genome comparison for organism in table below
    * Non-core comparison for organism in table below
    * Interpretation based on cutoffs:
        * Indistinguishable \>99.9% Kmer Identity 
        * Closely Related 98.7-99.9% Kmer Identity 
        * Possibly Related 95.0-98.7% Kmer Identity 
        * Unrelated \<95.0% Kmer Identity 

___

## Installing Plugin

1. Click releases on the github page.
2. Download the `Source Code (zip)`
3. Install the zip file through the Torrent Server Plugin interface

___
## Usage

__StrainConstructMer__ will process all the sample on a project.  No input or configuration is required.


#### Summary report

The summary report appears after the plugin has finished processing samples.  The report contains a table with each row
cooresponding to a sample on the run.  


| Column position | column Name | column information |
|-----------------|-------------|--------------------|
| 1               | Status      | <ul><li>Icon indicates if the sample __passed QC__ and is suitable for comparison</li><li>Clicking the icon brings up data files including antibiotic genes detected</li></ul> |
| 2               | Barcode     | Barcode ID of sample  |
| 3               | Sample      | The sample name of the sample|
| 4               | Ver         | The sample version.  The sample version is incremented automatically by __StrainConstructMer__. If the `bamfile_path` and `read count` are unique and `Sample_ID` and `Sample` are already in the database the `version` is incremented by one. Otherwise previous data is overwritten.|
| 5               | Coverage    | The coverage of the sample.  25X is the minimum coverage allowed. |
| 6               | Genome Size | The estimated genome size based on the number of distinct kmers observed. The genome size must be greater that 1,000,000bp.|
| 7               | With Expected Range | Is the estimated genome size is within 10% of the minimum and maximum genome size for the top classifier hit.|
| 8               | Q20 Bases           | The percentage of Q20 base observed in the reads|
| 9               | Top Classifier Hit  | <ul><li>The top classifier hit</li><li>Clicking link brings up page showing all classifier hits</li></ul>|
| 10              | Sample Saved        | Status of the database backup. Green check mark if sample successfully added to the database.|
| 11              | Resistant Genes     | Links to table displaying the antibiotic genes found in the sample |
| 12              | Kmer Histogram      | The histogram of kmer frequency and kmer count |


#### Screen shot of Report

![screenshot](docs/summary_screenshot.png)

___
#### 16S Classifier table

Clicking the top hit in the classifier table opens a new tab showing the complete list of classifier hits from NCBI's 
reference set

| Column position | column Name | column information |
|-----------------|-------------|--------------------|
| 1               | Kmer identity                   | The % kmer identity the sample shares with reference |
| 2               | Species                         | Species of reference  |
| 3               | Accession                       | Genbank accession of reference |
| 4               | Predicted Size for Strain       | Predicted Genomes size of the sample |
| 5               | Median Genome Size Species      | The Median genome size of the reference species based on completed NCBI genomes |
| 6               | Min Genome Size Species         | The Min genome size of the reference species based on completed NCBI genomes |
| 7               | Max Genome Size Species         | The Max genome size of the reference species based on completed NCBI genomes |
| 8               | Median Genome Size Genus        | The Median genome size of the reference genus based on completed NCBI genomes |
| 9               | Min Genome Size Genus           | The Min genome size of the reference genus based on completed NCBI genomes |
| 10              | Max Genome Size Genus           | The Max genome size of the reference genus based on completed NCBI genomes |


#### Screen shot of Classifier Hits table

![screenshot](docs/classifier_results.png)

___

#### Antibiotic Resistance Genes

Clicking the drug icon for a sample opens a new tab showing the antibiotic resistance genes found in the strain based on
NCBI's AMR reference set [NCBI AMR references](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA313047)

#### Screen shot of Antibiotic Resistance Genes

![screenshot](docs/antibiotic_gene_results.png)

| Column position | column Name | column information |
|-----------------|-------------|--------------------|
| 1               | Allele                       | The allele name of the AMR reference |
| 2               | Genbank Accession            | Genbank accession of reference  |
| 3               | Gene                         | The gene name of the AMR reference |
| 4               | Description                  | The description of the AMR reference |
| 5               | Percent Kmer Identity        | Kmer identity the sample shares with the AMR reference |
| 6               | Coverage Change from Genome  | The coverage change observed for the gene. This is an indication of the number of copies of the gene |
| 7               | Originating Species          | The species of origin for the AMR reference |



## Limitations and notes

#### Sample database and backups
`StrainCompareMer` relies on samples to be run through `StrainConstructMer`. `StrainConstructMer` therefore writes information
to an SQLite database which by default is placed into `/results/plugins/scratch/`.  The database only holds sparse information
about the sample. If the plugin results for a sample are deleted then `StrainComparMer` may fail. This can be easily be rectify
by rerunning `StrainConstructMer` on the project that contains the sample of interest. Alternatively, you can specify a
backup location. The backup location will hold the required data need to perform comparision and will not be affected if
a plugin is run multiple times or deleted.  Backup data will need to be deleted manually using the commandline interface.

Both the backup directory and SQLite database can be configured using the global configuration page on the plugin settings page.

By default the backup directory is not set. If you are interested in setting up a backup directory and do not need to have
multiple instruments accessing the same data `/results/plugins/scratch/` is a good option.

If you desire to have multiple instruments accessing to the same `StrainConstructMer` results, this can be achieved by mounting a
network drive and specifying the database and backup location point to this network drive.
___

#### Comparing Strains
* We have found that at least 25X coverage is need to preform a high quality comparison. Above this coverage it is
easier to remove errors from the data set.

* Determining if estimated genome size is accurate is based on the classifier results. The top hit from the classifier
  is compared to a list of known genome sized for the species. The Genome sizes are based on completed NCBI genomes.
 
* Access to the raw data files can be achieved by clicking the status icon

* QC parameters are hard coded
    * \>25X coverage
    * Genome size \> 1,000,000bp
    * Genome size within 10% of min and max for species. If the species is not in NCBI genomes size of genus will be used.
    if the genus is not present the Genome size is ignored as a QC parameter. 
    
---

_When using this resource in publications, please cite the following:
Keith Simmon PhD, ARUP Laboratories at the University of Utah, ARUP StrainTypeMer software_

&copy; 2018, ARUP Laboratories



___
