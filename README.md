
<img src="docs/ARUP_horz_2c_pos.jpg" width="50%">

---

# ARUP StrainTypeMer
__ARUP StrainTypeMer__ is a program that rapidly compares the nucleotide content of one or more samples/strains. 
Comparisons are made using 31bp overlapping kmers. The primary use case is epidemiological analysis. This tool can 
replace Pulsed Field Gel Electrophoresis (PFGE). Results are presented and interpreted in a manner similar to PFGE. 
Samples/strains are grouped as _Indistinguishable_, or _Closely_, _Possibly_, or _Unrelated_.


The program analyzes NGS data created by a Whole Genome Fragmentation Protocol. This version of __ARUP StrainTypeMer__ 
works as two Ion Torrent Plugins: __ARUP StrainConstructMer__ and __ARUP StrainCompareMer__. While most epidemiological 
NGS methods require reference alignment, __ARUP StrainTypeMer__ is designed to run reference free. This creates a 
universal analysis method to compare strains unhindered by the amount of identity a sample has with the reference 
sequence or the availability of a suitable reference. Because __ARUP StrainTypeMer__ is reference-free, it creates a 
universal strain comparison method that will work across many species and genera. Analysis can be completed in under 
30 minutes for 5-20 samples.

# ARUP StrainConstructMer
__ARUP StrainConstructMer__ is the first plugin in the StrainTypeMer analysis. This plugin processes a sample and 
transforms the data so the __ARUP StrainCompareMer__ plugin can perform comparisons.

## Features
* Calculates coverage
* Determines the genome size of the organism
    * o	Determines if genome size is within expected size for organism
* Identifies the bacteria based in [NCBI 16S rRNA RefSeq references](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA33175)
* Determines MLST type using [PUBMLST](https://pubmlst.org/])
* Identifies antibiotic resistant genes present in the sample using [NCBI AMR references](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA313047)
* Performs quality checks to determine if samples can be compared in __ARUP StrainCompareMer__

# Guide
__ARUP StrainConstructMer__ will process all the samples on a project. No input or configuration is required. Once 
samples are complete, a summary report will be displayed. Samples processed by the plugin are placed into an embedded 
database that ARUP StrainCompareMer can access for comparisons. The database location is set on the global 
configuration page. It must be identical to the location set in ARUP StrainCompareMer  
(default location `/results/plugins/scratch/`)
[see below for notes regarding the StrainConstructMer database.](#sample-database-and-backups)

### Summary Report
The summary report appears after the plugin has finished processing samples.  The report contains a table with each row
corresponding to a sample on the run.  

| Column Position | Column Name         | Column Information |
|-----------------|---------------------|--------------------|
| 1               | Status              | <ul><li>Icon indicates if the sample __passed QC__ and is suitable for comparison</li><li>Clicking the icon brings up data files including antibiotic genes detected</li></ul> |
| 2               | Barcode             | Barcode ID of sample  |
| 3               | Sample              | The sample name of the sample|
| 4               | Ver                 | The sample version.  The sample version is assigned automatically by __StrainConstructMer__. If the `bamfile_path` and `read count` are unique and `Sample_ID` and `Sample` are already in the database the `version` is incremented by one. Otherwise previous data is overwritten.|
| 5               | Coverage            | The coverage of the genome. Calculated by dividing the distinct kmers by the sum of their frequency. 25X is the minimum coverage needed to remove errors for the dataset. If the coverage is less the 25X it fails QC the sample will a red X in the status column. |
| 6               | Genome Size         | The estimated genome size is based on the number of distinct kmers observed.|
| 7               | With Expected Range | Is the estimated genome size is within 10% of the minimum and maximum genome size for the top classifier hit.|
| 8               | Q20 Bases           | The percentage of Q20 base observed in the reads|
| 9               | Top Classifier Hit  | <ul><li>The top classifier hit</li><li>Clicking link brings up page showing all classifier hits</li></ul>|
| 10              | Sample Saved        | Status of the database backup. A green check mark appears if sample successfully added to the database.|
| 11              | Resistant Genes     | Links to table displaying the antibiotic genes found in the sample |
| 12              | Kmer Histogram      | The histogram of kmer frequency and kmer count |


#### Example of Report output

<kbd>
  <img src="docs/summary_screenshot.png" border="1">
</kbd>

---
### 16S rRNA Gene Classifier Table

Clicking the value in the top hit column for a sample in the summary table opens a new tab showing the complete list of classifier hits
from NCBI's reference set for the selected sample.

| Column position | column Name | column information |
|-----------------|-------------|--------------------|
| 1               | Kmer identity                   | The % kmer identity the sample shares with reference |
| 2               | Species                         | Species of reference  |
| 3               | Accession                       | GenBank accession of reference |
| 4               | Predicted Size for Strain       | Predicted Genomes size of the sample |
| 5               | Median Genome Size Species      | The Median genome size of the reference species based on completed NCBI genomes |
| 6               | Min Genome Size Species         | The Min genome size of the reference species based on completed NCBI genomes |
| 7               | Max Genome Size Species         | The Max genome size of the reference species based on completed NCBI genomes |
| 8               | Median Genome Size Genus        | The Median genome size of the reference genus based on completed NCBI genomes |
| 9               | Min Genome Size Genus           | The Min genome size of the reference genus based on completed NCBI genomes |
| 10              | Max Genome Size Genus           | The Max genome size of the reference genus based on completed NCBI genomes |


### Example of Classifier Hits Table

<kbd>
  <img src="docs/classifier_results.png" border="1">
</kbd>

---
### Antibiotic Resistance Genes

Clicking the drug icon for a sample opens a new tab showing the antibiotic resistance genes found in the strain based on
NCBI's AMR reference set [NCBI AMR references](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA313047)


| Column position | column Name | column information |
|-----------------|-------------|--------------------|
| 1               | Allele                       | The allele name of the AMR reference |
| 2               | GenBank Accession            | GenBank accession of reference  |
| 3               | Gene                         | The gene name of the AMR reference |
| 4               | Description                  | The description of the AMR reference |
| 5               | Percent Kmer Identity        | Kmer identity the sample shares with the AMR reference |
| 6               | Coverage Change from Genome  | The coverage change observed for the resistant gene. This is an indication of the number of copies of the gene compared to the chromosome coverage |
| 7               | Originating Species          | The species of origin for the AMR reference |

### Example of Antibiotic Resistance Genes

<kbd>
  <img src="docs/antibiotic_gene_results.png" border="1">
</kbd>

---

# ARUP StrainCompareMer
__ARUP StrainCompareMer__ is the second plugin in the __ARUP StrainTypeMer__ analysis.

## Features
* Includes multiple comparisons
    * Full Genome comparison
        * Similarity matrix showing percent kmer identity shared between each to the samples.
    * Full Genome comparison with select references
        * Same as above but includes selected NCBI references.  References selected based on sample ID.  
    * Core genome comparison for organism in table below
        * Compares the core kmer subset between samples (most conserved kmers). Comparison is made if all samples are
        the same organism and the core reference set is available.      
    * Non-core comparison for organism in table below
        * Compares the non-core kmer subset between sample (accessory genome). Comparison is made if all samples are
            the same organism and the core reference set is available.
    * Rescue
        * Modified Full genome comparison that only compares kmers in the smallest of the to samples processed.
        * This is helpful:
            * To compare a sample with low coverage to a sample with higher coverage
            * To observe changes caused a large acquisition of DNA
    * Relationships
        * Interpretation based on cutoffs:
            * Indistinguishable \>99.9% kmer identity
            * Closely Related 98.7-99.9% kmer identity
            * Possibly Related 95.0-98.7% kmer identity
            * Unrelated \<95.0% kmer identity
    * Comparision Table
        * Table shows details about each comparison
    * Strain Summary
        * Similar to the Summary page form __ARUP StrainConstructMer__.  MLST profiles included in output.

# Guide
The plugin requires a CSV file as input. The CSV file indicates the strains to be processed. The strains must exist in
the database and the plugin must be able to locate the files associated with the strain.  The database location is set
on the global configuration page. It must be identical to the location set for __ARUP StrainConstructMer__ (default
location `/results/plugins/scratch/`).

#### Core Genomes

Core genomes were constructed using completed RefSeq genomes from NCBI. A maximum of 25 genomes were used if available.
By default we use an 80% threshold: a kmer needs to be observed in 80% of the genomes to be considered a "core" kmer.
This is modified for some species in order to target a core kmer size near 40%-70% of the overall genome size.

| Organism                     | Average Genome Size   | Genomes Analyzed | Core Kmer Size (Pct of Genome Size) | Percent Cutoff   |
|:-----------------------------|:---------------------:|:----------------:|:-----------------------------------:|:----------------:|
|_Acinetobacter baumannii_     | 3,945,908             | 25               | 2,040,177 (52%)                     | 80               |
|_Enterococcus faecalis_       | 2,909,703             | 20               | 1,812,881 (62%)                     | 80               |
|_Enterococcus faecium_        | 2,827,968             | 25               | 1,833,417 (65%)                     | 95               |
|_Escherichia coli_            | 4,969,315             | 25               | 1,914,103 (39%)                     | 80               |
|_Klebsiella pneumoniae_       | 5,315,713             | 25               | 3,825,265 (71%)                     | 90               |
|_Pseudomonas aeruginosa_      | 6,581,730             | 25               | 4,416,660 (67%)                     | 80               |
|_Serratia marcescens_         | 5,234,322             | 25               | 1,910,225 (36%)                     | 60               |
|_Staphylococcus aureus_       | 2,852,092             | 25               | 1,726,487 (61%)                     | 80               |
|_Staphylococcus epidermidis_  | 2,544,188             | 15               | 1,665,072 (65%)                     | 80               |
|_Stenotrophomonas maltophilia_| 4,726,726             | 20               | 1,240,438 (26%)                     | 50               |


<em> Note: We constructed the core kmer sets using utility script `create_core_reference_set.py`. However users may
construct core reference set as they see fit.  The file must be placed into
`/results/plugins/StrainCompareMer/strain_comparison/resources/core_reference_sets` and must be prefixed with a
organism name. </em>

### Example of Full Genome Similarity Matrix

The full genome similarity matrix shows a pairwise comparisons on each of the input samples.  The samples are clustered and a dendrogram is drawn to show the
relationships.

<kbd>
<img src="docs/full_matrix.png" border="1">
</kbd>


### Example of Relationship Table

The relationship table has a row for each sample or group (Indistinguishable isolates).  The columns display how these groups are related to other samples using interpretive criteria (e.g. Closely Related, Possibly Related, Unrelated).

<kbd>
<img src="docs/relationship_table.png" border="1">
</kbd>

### Example of Comparison Table

The comparison table gives detailed information including error correction for each comparison made.

<kbd>
<img src="docs/comparison_table.png" border="1">
</kbd>

### Example of Reference Matrix

The reference similarity matrix supplements the full genome matrix with selected references from NCBI genomes. A reference genome matrix will only be created if all samples in the comparison are one of the following: _Enterococcus faecalis_, _Enterococcus faecium_, _Staphylococcus aureus_, or _Acinetobacter baumannii_.

<kbd>
<img src="docs/reference_matrix.png" border="1">
</kbd>

### Example of Core Matrix

The core genome similarity matrix will be created if all samples are identified as one of the organisms listed in the Core Genomes table.  The comparison is made only with kmers in the core kmer set.  Details about the comparison can be viewed in the Comparison Table.

<kbd>
<img src="docs/core_matrix.png" border="1">
</kbd>

### Example of Non-Core Matrix

The non-core similarity matrix is the antithesis of the core matrix. Only non-core kmers are compared between samples.  

<kbd>
<img src="docs/non-core_martix.png" border="1">
</kbd>

### Example of Strain Summary

The Strain summary is similar to __ARUP StrainConstructMer's__ output. It also includes a column showing the MLST type identified for the Sample.

<kbd>
<img src="docs/strain_summary.png" border="1">
</kbd>

___

## Installing Plugins

1. Download the `ARUPStrainConstructMer.zip` and `ARUPStrainCompareMER.zip`
2. Install the zip file through the Torrent Server Plugin interface

___

# Limitations and Notes

## Sample Database and Backups
`ARUPStrainCompareMer` relies on samples to be run through `ARUPStrainConstructMer`. `ARUPStrainConstructMer` writes information
to an SQLite database which by default is placed into `/results/plugins/scratch/`.  The database only holds sparse information
about the sample. If the plugin results for a sample are deleted then `ARUPStrainCompareMer` may fail. This can be easily be rectified
by rerunning `ARUPStrainConstructMer` on the project that contains the sample of interest. Alternatively, you can specify a
backup location. The backup location will hold the required data need to perform comparison and will not be affected if
a plugin is run multiple times or deleted.  If desired backup data will need to be deleted manually using the commandline interface.

Both the backup directory and SQLite database can be configured using the global configuration page on the plugins on the Torrent Server page.

By default the backup directory is not set. If you are interested in setting up a backup directory and do not need to have
multiple instruments accessing the same data `/results/plugins/scratch/` is a good option.

If you require multiple instruments accessing to the same `ARUPStrainConstructMer` results, this can be achieved by mounting a
network drive. Each instrument should be configured to point to the identical location on this network drive.
___

## Comparing Strains
* We have found that at least 25X coverage is needed to preform a high quality comparison. Above this coverage it is
easier to remove errors from the data set.

* Determining if estimated genome size is accurate is based on the classifier results. The top hit from the classifier
  is compared to a list of known genome sizes for the species. The genome sizes are based on completed NCBI genomes.

* Access to the raw data files can be achieved by clicking the status icon

* QC parameters are hard coded
    * \>25X coverage
    * Genome size \> 1,000,000bp
    * Genome size within 10% of min and max for species. If the species is not in NCBI genomes, size of genus will be used.
    if the genus is not present the genome size is ignored as a QC parameter.

---


_When using this resource in publications, please cite the following:
Keith Simmon PhD, ARUP Laboratories at the University of Utah, ARUP StrainTypeMer software_

&copy; 2019, ARUP Laboratories

___
