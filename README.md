# StrainContructMer
Ion torrent plugin to perform strain typing. 
This plugin process a sample and transform the data so that comparisons can be performed with the __StrainCompareMer__ plugin

___
# Features
__StrainConstructMer__ and __StrainCompareMer__ together comprise a package called __StrainTypeMer__, which performs strain comparison using
kmers.  Strain Comparsion is performed reference-free creating a universal strain comparison method that will work across many species and genera
of organisms.

### StrainTypeMer
* Calculates coverage
* Determines the genome size of the organism
    * Determines if genome size is expected value
* Identifies the bacteria based on NCBI [16S rRNA RefSeq] (PRJNA33175)
* Determines MLST type using [PUBMLST](https://pubmlst.org/])
* Identify antibiotic resistant genes present (PRJNA313047)
* Includes mulitple comparisons
    * Full Genome comparison
    * Filter Comparison to remove errors
    * Core genome comparison for *Staphylococcus aureus*, *Acinetobacter baumannii*, *Enterococcus facalis*, *E. faecium*
    * Non-core comparison for *Staphylococcus aureus*, *Acinetobacter baumannii*, *Enterococcus facalis*, *E. faecium*
    * Rescue comparison for contaminated samples

___

## Setup Applications and Resources
### jellyfish [required]

jellyfish is a kmer counter need by `StrainConstructMer`
Get the compiled tar ball from github and install.

```bash
mkdir /results/plugins/scratch/apps
cd /results/plugins/scratch/apps
wget https://github.com/gmarcais/Jellyfish/releases/download/v2.2.7/jellyfish-2.2.7.tar.gz
tar xvf jellyfish-2.2.7.tar.gz
cd jellyfish-2.2.7
./configure
make -j 4

# may not be required
sudo make install

# test command
$: /results/plugins/scratch/apps/jellyfish-2.2.7/bin/jellyfish
# should return
...Too few arguments
...Usage: jellyfish <cmd> [options] arg...
...Where <cmd> is one of: count, bc, info, stats, histo, dump, merge, query, cite, mem, jf.
...Options:
...  --version        Display version
...  --help           Display this message
```
____
## Installing

1. Click releases on the github page.
2. Download the `Source Code (zip)`
3. Change the name of the download file to `StrainConstructMerBeta`
    - may need to recompress if the file decompress automatically during download
4. Install the zip file throught the Torrent Server Plugin interface

___
## Usage

__StrainConstructMer__ will process all the sample on a project.  No input are configuration is required.


[]



## Limitations and notes

#### Sample database and backups
`StrainCompareMer` relies on samples to be run through `StrainConstructMer`. `StrainConstructMer` therefore writes information
to an SQLite database which by default is placed into `/results/plugins/scratch/`.  The database only holds sparse information
about the sample. If the plugin results for a sample are deleted then `StrainComparMer` may fail. This can be easlily rectify
by reruning `StrainConstructMer` on the project that contains the sample of interest. Alternatively, you can specify a
backup location. The backup location will hold the required data need to perform comparsion and will not be affected if
a plugin is run mulitple times or deleted.  Backup data will need to be deleted manually using the commandline interface.

Both the backup directory and SQLite database can be configured using the global configuration page on the plugin settings page.

By default the backup directory is not set. If you are interested in setting up a backup directory and do not need to have
mulitple instruments have access to the same data `/results/plugins/scratch/` is a good option.

If you desire to have mulitple instruments have access to the same StrainConstructMer results, this can be achevied by mounting a
network drive and specifiing the the database and backup location point to this network drive.
___
#### Comparing Strains
* We have found that at least 25X coverage is need to preform a high quaility comparison. Above this coverage it is
easier to remove errors from the dataset.

* Determining if Estimated genome size is accurated is based on the classifier results. The top hit from the clasifier
  is compared to a list of known genome size for the species. Currently the genome size list only contains *Acinetobacter*,
  *Staphylococcus*, and *Enterococcus*.
 
* Access to the raw data files, which include antibiotic resistant genes, can be accessed by clicking the status icon

* QC parameters are hard coded
    * \>25X coverage
    * Genome size \> 1,000,000bp
    * Genome size within 10% of min and max for species or genus (only used for *Acinetobacter*, *Staphylococcus*, and *Enterococcus*)
    

## Preparing a Release

Software is released in the Dev, Cert, Prod enviroment.  To release a verion do the following

1. Create a Branch
2. Rename the `StrainConstructMer` python file and `Class` to include the release name
    1. Example. `StrainConstructMer_DEV`
3. Create a release on `github`
4. Download `Zip` File and install through plugin TS interface
5. Run test on command line
___
