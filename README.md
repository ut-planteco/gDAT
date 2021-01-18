# gDAT (Graphical Downstream Analysis Tool)
gDAT is an open-sourced pipeline, which offers an intuitive graphical interface allowing users to conduct reliable HTS read analysis that yields datasets ready for statistical analysis. The main focus of the pipeline is to maximise the ease-of-use of different tools by having most optimal predefined input parameters and allowing selection of files through open box dialogs, thus avoiding input errors, which can arise when working from the command line. The pipeline is a generalisation of the procedures used in several studies of arbuscular mycorrhizal (AM) fungal communities and takes into account best practises while dealing with HTS data. As such, the pipeline parameter defaults are optimised for analysis of AM fungal DNA sequences of small subunit (SSU) ribosomal RNA gene or internal transcribed spacer (ITS) regions. Nonetheless, the workflow is generally applicable to analysis of a broad range of metabarcoding data generated with different sequencing platforms and multitude of organism groups. Main focus group for the pipeline are ecologists who have minimal experience with command line and bioinformatics tools. Pipeline is optimized to be used with commodity hardware.

## Quick start

Binaries for Windows 10, Linux and macOS can be [downloaded here](https://github.com/ut-planteco/gDAT/releases). 

## The interface and main functions of the gDAT.
![main.png](https://github.com/ut-planteco/gDAT/raw/master/manual/assets/pipeline_main.png)

## Prerequisite

If quick start option is not used, user needs to install following programs and packages: [vsearch 2.11.0](https://github.com/torognes/vsearch/releases), [NCBI-BLAST+ 2.8.1+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download), [FLASh 1.2.11](https://ccb.jhu.edu/software/FLASH/), [Python 2.7+ or 3+](https://www.python.org/downloads/) and python-tk. All the programs are distributed for Windows, Linux and macOS. Debian Linux users can use following commands to download all the necessary packages:

```
sudo apt install vsearch
sudo apt install ncbi-blast+
sudo apt install python python-tk
sudo apt install python3 python3-tk
```

## Download the pipeline

Pipeline is bundled to binaries and easiest is to download [the binary set](https://github.com/ut-planteco/gDAT/releases). Users can manually download GitHub repository, release version or from command line with wget or using GitHub clone. 

## Start the pipeline

Pipeline can be executed using `gdat.py` file with python or using double click from file manager (last option is dependant of operating system and configuration).

`python gdat.py`

## Documentation

The gDAT documentation is available under the manual/ folder as HTML files that can be opened with browser or through GitHub when browsing to manual/ folder.

## Database

Pipeline includes MaarjAM database consisting of arbuscular mycorrhizal sequences from SSU gene region. Database uses virtual taxa (VT) concept, which are phylogenetically defined sequence groups roughly corresponding to species-level taxa. Pipeline allows usage of custom databases, copy FASTA files to db/ folder and files will be automatically shown in the GUI. For BLAST+ search, use "Build a BLAST+ database" function in the pipeline to convert the FASTA files into correct format.  

## License

The gDAT is distributed under the GNU General Public License version 3. [CC-BY](https://creativecommons.org/licenses/by/3.0/)

## Examples

Example datasets for 454 and Illumina data in various formats is provided in the example/ folder.

## Releases

* v1.0: Initial release (29.09.2020).

* v1.1: Added support to use vsearch for pairing paired-end reads and vsearch global pairwise alignment search to identify sequences against reference database (18.01.2021).

## Citation

Pending

Please note that citing any of the third party programs, e.g. vsearch, NCBI-BLAST+ and FLASh, is strongly suggested.