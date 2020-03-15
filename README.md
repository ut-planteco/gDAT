# gDAT (Graphical Downstream Analysis Tool)
gDAT is open-sourced pipeline, which offers an intuitive graphical interface allowing users to conduct reliable HTS read analysis that yields datasets ready for statistical analysis. The main focus of the pipeline is to maximise the ease-of-use of different tools by having most optimal predefined input parameters and allowing selection of files through open dialog, thus avoiding input errors, which can arise when working from the command line. The pipeline is a generalisation of the procedures used in several studies of arbuscular mycorrhizal (AM) fungal communities and takes into account best practises while dealing with HTS data. As such, the pipeline parameter defaults are optimised for analysis of AM fungal DNA sequences of small subunit (SSU) ribosomal RNA gene or internal transcribed spacer (ITS) regions. Nonetheless, the workflow is generally applicable to analysis of a broad range of metabarcoding data generated with different sequencing platforms and multitude of organism groups. Main focus group of the pipeline are newcomers with minimal experiments with command line and bioinformatic tools. Pipeline is optimized to be used with commodity hardware. 

The interface and main functions of the gDAT.

![main.png](https://github.com/utplanteco/gDAT/blob/master/manual/assets/pipeline_main.png?raw=true)

## Quick start

Binaries for Windows 10, Linux and macOS can be [downloaded here](https://github.com/utplanteco/gDAT/releases). 

## Prerequisite

If quick start option is not used, user needs to install following programs and packages: [vsearch 2.11.0](https://github.com/torognes/vsearch/releases), [NCBI-BLAST+ 2.8.1+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download), [FLASh 1.2.11](https://ccb.jhu.edu/software/FLASH/), [Python 2.7+ or 3+](https://www.python.org/downloads/) and python-tk. All the programs are distributed for Windows, Linux and macOS. Debian Linux users can use following commands to download all the necessary packages:

```
sudo apt install vsearch
sudo apt install ncbi-blast+
sudo apt install python2 python2-tk
sudo apt install python3 python3-tk
```

## Download pipeline

Pipeline is bundled to binaries and easiest is to download [the binary set](https://github.com/utplanteco/gDAT/releases). Users can manually download GitHub repository, release version or from command line with wget or using GitHub clone. 

## Start pipeline

Pipeline can be executed using `gdat.py` file with python or using double click from file manager (last option dependant of operating system and configuration).

`python gdat.py`

## Documentation

The gDAT documentation is available under the manual/ folder as HTML files that can be opened with browser or through GitHub when browsing to manual/ folder.

## License

The gDAT is distributed under the GNU General Public License version 3.

## Examples

Example datasets for 454 and Illumina data in various formats is provided in the example/ folder.

## Citation

Pending

Please note that citing any of the third party programs, e.g. vsearch, NCBI-BLAST+ and FLASh, is strongly suggested.
