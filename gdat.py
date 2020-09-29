#!/usr/bin/env python

import os
import sys
import signal
import subprocess
import time
import glob
import webbrowser
import signal
from datetime import datetime
#import psutil
import multiprocessing
import gzip
from threading import Thread
try:
    import Tkinter as tk
    from tkFileDialog import askopenfilenames, askopenfilename, askdirectory, asksaveasfilename
    from tkFont import Font as font
    import tkMessageBox
    from Queue import Queue, Empty
    import ConfigParser
except ImportError:
    try:
        import tkinter as tk
        from tkinter.filedialog import askopenfilenames, askopenfilename, askdirectory, asksaveasfilename
        from tkinter import font
        from tkinter import messagebox as tkMessageBox
        from queue import Queue, Empty
        import configparser as ConfigParser
    except ImportError:
        sys.stderr.write("TKinter package is missing. ")
        if os.name == "nt":
            sys.stderr.write("Please reinstall Python and enable tkinter in the installation menu.\n")
            input("Press Enter to continue...")
        # expect that Mac has it by default, only Linux is missing it
        else:
            if sys.version_info[0] >= 3:
                sys.stderr.write("For Debian/Ubuntu/Mint, use following\ncommand to install: sudo apt install python3-tk\n")
                input("Press Enter to continue...")
            else:
                sys.stderr.write("For Debian/Ubuntu/Mint, use following\ncommand to install: sudo apt install python-tk\n")
                raw_input("Press Enter to continue...")
        sys.exit(1)

class CreateToolTip(object):
    """
    create a tooltip for a given widget (label)
    """
    def __init__(self, widget, text):
        self.waittime = 500     #miliseconds
        self.wraplength = 400 # 180   #pixels
        self.widget = widget
        self.text = text
        self.widget.bind("<Enter>", self.enter)
        self.widget.bind("<Leave>", self.leave)
        self.widget.bind("<ButtonPress>", self.leave)
        self.id = None
        self.tw = None

    def enter(self, event=None):
        self.schedule()

    def leave(self, event=None):
        self.unschedule()
        self.hidetip()

    def schedule(self):
        self.unschedule()
        self.id = self.widget.after(self.waittime, self.showtip)

    def unschedule(self):
        id = self.id
        self.id = None
        if id:
            self.widget.after_cancel(id)

    def showtip(self, event=None):
        x = y = 0
        x, y, cx, cy = self.widget.bbox("insert")
        x += self.widget.winfo_rootx() + 25
        y += self.widget.winfo_rooty() + 20
        # creates a toplevel window
        self.tw = tk.Toplevel(self.widget)
        # Leaves only the label and removes the app window
        self.tw.wm_overrideredirect(True)
        self.tw.wm_geometry("+%d+%d" % (x, y))
        label = tk.Label(self.tw, text=self.text, justify='left',
                       background="#ffffff", relief='solid', borderwidth=1,
                       wraplength = self.wraplength)
        label.pack()

    def hidetip(self):
        tw = self.tw
        self.tw= None
        if tw:
            tw.destroy()


class BaseFrame(tk.Frame):
    """An abstract base class for the frames that sit inside PythonGUI.

    Args:
      master (tk.Frame): The parent widget.
      controller (PythonGUI): The controlling Tk object.

    Attributes:
      controller (PythonGUI): The controlling Tk object.

    """
    _form_row = 0
    _form_col = 0
    _max_row = 0
    padx = 10
    pady = 10
    left_button = []
    right_button = []
    params = {}
    configuration = {}
    tooltips = {
        "fasta" : """Select raw reads file for cleaning and quality filtering (.fasta, .fna, .fastq or .fq).
If you select a FASTA file, be sure also to select a quality file in the next field.""",
 
        "input" : "Select an input FASTA file (.fna, .fasta), files packed with gz extension are supported.",
 
        "fqinput" : "Select input FASTQ file (.fq, .fastq), files packed with gz extension are supported.",
 
        "qual" : "Define a 454 quality file (.qual), where a numerical quality value is provided for each base, files packed with gz extension are supported.",
 
        "ffastq" : "Select a forward FASTQ raw reads file for cleaning and quality filtering (R1 .fastq or .fq file), files packed with gz extension are supported, barcodes located in the header are not supported.",
 
        "rfastq" : "Select a reverse FASTQ raw reads file for cleaning and quality filtering (R2 or R4 .fastq or .fq file), files packed with gz extension are supported, barcodes located in the header are not supported.",
 
        "foligo" : "Select a forward FASTQ oligo reads file. Sometimes barcode sequences are provided in a separate file (R2 .fastq or .fq file), files packed with gz extension are supported.",
 
        "roligo" : "Select a reverse FASTQ oligo reads file. Sometimes barcode sequences are provided as separate file (R3 .fastq or .fq file), files packed with gz extension are supported.",
 
        "barcode" : "Select a sample sheet file for demultiplexing sequencing into samples (.tsv or .txt file), files packed with gz extension are supported.",

        "forward" : """Select a forward FASTQ reads file (.fq, .fastq), which will be combined with the reverse reads.
Interleaved reads means that forward and reverse reads are in the same file and come in pairs.
When using interleaved reads, make sure that first part of the header name matches for read pairs; 
otherwise an error may be thrown.""",
 
        "reverse" : """Select reverse FASTQ reads file (.fq, .fastq), which will be  combined with the forward reads.
Leave empty if an interleaved file is used.""",
 
        "interleaved" : "Specify if the forward and reverse reads are interleaved or in separate files.",
 
        "min_length" : "Minimum allowed overlap when combining sequences; sequences with shorter overlap are discarded.",
 
        "max_length" : "Maximum allowed overlap when combining sequences; sequences with larger overlap are discarded.",
 
        "identity" : "Identity in percentage (without percentage mark) for overlap, clustering or BLAST+ alignment. Default values: overlap = 75%, clustering and BLAST+ match = 97%.",
 
        "alignment" : "Minimum allowed alignment length in percentage of the BLAST+ result compared with the shorter of the query or reference sequence length.",
        
        "alignmentbp" : "Minimum allowed alignment length in base pairs of the BLAST+ result allowing removale of short alignments.", 

        "threads" : "Number of threads used to run the script. Some programs will not run if the number of threads exceeds the computer's thread count. It is automatically set to the maximum threads available on the computer minus 1 thread to allow other programs to work in the background.",
 
        "output" : "Define file name where the output is written. If no name is given, the system will automatically generate a name.",
 
        "folder" : "Define a folder where demultiplexed reads are located. Supports Illumina demultiplexed reads. The pipeline searches for .fq and .fastq files, which have R1 and R2 in the name. Sample names are generated from the first part of the file. The pipeline supports packed files with gz extension", 
 
        "sample" : "Sample conversion table, where for each file or pair of files a sample name is defined. Select the conversion option to generate the conversion table for your defined file automatically based on the selected folder. You can then use a spreadsheet program to change the sample names accordingly",
 
        "sample_col" : "Define the column number where sample names are located in the sample sheet file; the left-most column is counted as 1.",
 
        "fbarcode_col" : "Define the column number where forward barcode sequences are located in the sample sheet file.",
 
        "fprimer_col" : "Define the column number where forward primer sequences are located in the sample sheet file.",
 
        "rbarcode_col" : "Define the column where reverse barcode sequences are located in the sample sheet file. If forward and reverse barcodes are fused together ina single column, use the same number as used for the forward barcode.",

        "rprimer_col" : "Define the column where reverse primer sequences are located in the sample sheet file.",
 
        "fprimer" : "Users can override the forward primer sequence from the sample sheet file. Use commas to separate multiple primers.",
 
        "rprimer" : "Users can override the reverse primer sequence from the sample sheet file; Use commas to separate multiple primers.",
 
        "f_len" : "Define the trimming length (number of nucleotides) to remove primers, adapters or low quality regions  in the forward reads.",
 
        "r_len" : "Define the trimming length (number of nucleotides) to remove primers, adapters or low quality regions in the reverse reads.",
 
        "quality_avg" : "Allowed average read quality between 0 to 41. Average quality is calculated after barcode+primer removal and trimming.",
 
        "fmin_len" : "Minimum allowed sequence length after trimming and barcode+primer removal.",
 
        "fmax_len" : "Trim longer sequences to a defined length; useful for removing primers and low quality regions.", 
 
        "trimming_window" : "Trimming window size (nucleotides) used to calculate average quality. If the average quality drops below the threshold, the rest of the sequence is removed from the beginning of the sliding window position.",

        "trimming_avg" : "Trimming window average quality threshold between 0 to 41. Use a lower value than the average quality of the sequence.",
 
        "mismatch_barcode" : "Allow barcode and/or primer mismatches. Only 1 nucleotide mismatches are allowed.",
 
        "mismatch" : "Allow barcode and/or primer mismatches. Only 1 nucleotide mismatches are allowed.",
 
        "ignore_bases" : "Ignore first bases for quality and primer checking. This may be useful when the sequencer has calibration issues and allows including sequences having mismatches at the front of the read.",
 
        "mismatch_indel" : "Allow insertion and/or deletion errors in the primer.",
 
        "zeros" : "Fill empty cells in the output data matrix with 0 or NA .",
 
        "lookup" : "Convert BLAST+ hits to values defined in the lookup. If no match is found in the lookup, the description in the BLAST+ output is used.", 
 
        "qvariable_start" : "Specify the start position of a variable region in the query sequence. This is useful when using tagmentation based Illumina and users want to filter out hits in the invariable region. These features need to have all database sequences trimmed to start from the same position. If users' sequences start after the NS31 primer, nucleotides 70 - 300 represent a variable region.",

        "qvariable_end" : "Specify the end position of a variable region in the query sequence.",
        
        "variable_start" : "Specify the start position of a variable region in the reference database. This is useful when using tagmentation based Illumina and users want to filter out hits in the invariable region. These features needs to have all database sequences trimmed to start from the same position. If users' sequences start after the NS31 primer, nucleotides 70 - 300 represent a variable region.",
 
        "variable_end" : "Specify the end position of a variable regionin the reference database.",
 
        "taxonomynode" : "Specify an NCBI taxonomy node file to place BLAST+ hits onto the taxonomy (nodes.dmp); supports packed .gz files.",
 
        "taxonomynames" : "Specify an NCBI taxonomy name file (names.dmp); supports packed .gz files.",
 
        "taxonomynucl" : "Specify an NCBI taxonomy lookup file (nucl_gb.accession2taxid); supports packed .gz files.",
 
        "mode" : "Specify the BLAST+ database used. For MaarjAM and UNITE, the last part of the hit description is used; GB/NCBI uses the second column GB accession code; and full description uses the third column from the BLAST+ output.",
 
        "clmode" : "Specify which clustering method to use with vsearch program. cluster_fast will first sort reads based on sequence length, while cluster_size will sort based on sequence abundance in decreasing order.",
 
        "cmode" : "Specify which method to use for chimera checking. With reference mode, users need to define the database FASTA file location. In de novo mode, sequences are first clustered with 97 percentage and cluster centroids are used as the reference for chimera checking. Using both methods in parallel is recommended.",
 
        "database" : "Select a database in FASTA format (.fna, .fasta) to be used for chimera checking.",

        "blastdb" : "Select a database in BLAST+ format.",
 
        "blastdbradio" : "Select a prebuilt database in BLAST+ format from the db/ folder.",
 
        "db" : "Select a database in FASTA format for chimera checking",
 
        "dbradio" : "Select a prebuilt database in FASTA format from the db/ folder for chimera checking",
 
        "blastinput" : "Specify a tab delimited BLAST+ output file, where each query sequence alignment hit against the reference database is shown.",
 
        "blast" : "Specify a tab-delimited BLAST+ output file to identify alignment strands and reverse complement sequences with incorrect orientation.",
 
        "lookup_blast": "Specify a tab-delimited BLAST+ output file when using NCBI BLAST+. The second column will be used to build the taxonomy used to construct pivot table.",
        
        "uc" : "If using a clustered file for BLAST+, define the cluster output file *.uc to count sequence occurences in each cluster for presentation in the final data matrix.",
 
        "ncbi" : "If a partitioned run is used for the BLAST+, it will automatically run all the partitions sequentially. This is useful when using the GB/NCBI database, which is partitioned into multiple smaller databases.",

        "muco_search" : "Specify organism names, separated by commas, to identify sequences based on the BLAST+ results.",
 
        "cluster_first" : "Pick the first (centroid) sequences from each cluster.",
 
        "cluster_group" : "The output will group clusters; otherwise clusters will be written as they appear in the input fasta.",
 
        "cluster_sort" : "Sort clusters by their size (number of sequences in the cluster) in descending order. This can helpfor further sorting and filtering of clusters by hand.",
 
        "cluster_pick" : "Select a method for picking sequences from clusters. Sequentially picks sequences in the order in which they were added to clusters; divided equally uses a constant step to pick sequences from highest to lowest identity within the cluster and randomly selects sequences from clusters at random. It is recommended is to use equally since it represents the full diversity of a cluster.",
 
        "multiblast" : "Select multiple BLAST+ outputs to combine into a single file. To select multiple files, hold down SHIFT or CTRL.",
 
        "fetchblast" : "Select GB/NCBI BLAST+ output, where the second column contains the GB accession code that will be used to build a lookup table and fetch sequences from the database.",
 
        "cluster_uc" : "Specify cluster file (*.uc) generated by vsearch containing information about how clusters were generated and information about cluster elements identities.",

        "primers_mixed" : "Specify whether the primers are mixed between the forward and the reverse reads.",

        "homopolymer" : "Truncate sequences with repeating bases by reducing homopolymer length to user specific value.",
	
        "qiime" : "Select an input FASTA file to convert between QIIME v1 and pipeline compatible formats. The pipeline will automatically detect the format of the file based on the header: | (pipe) is used for gDAP and _ (underscore) is used for QIIME v1 to separate sample information from the sequence identifier.",
 
        "phred" : "Define the PHRED quality base score used by the FASTQ file. The default is 33, but older Illumina formats use 64.",
    
        "output_qual" : "Specify whether to generate the QUAL file based on FASTQ qualities.",
 
        "evalue" : "BLAST+ based cut-off expectation value to show the significance of the alignment; use lower numbers to discard shorter alignments.",
 
        "count" : "Specify a number of sequences to pick from each cluster.",
 
        "cluster_size" : "Specify the minimum allowed cluster size (number of sequences in the cluster).",
 
        "cluster_length" : "Specify the minimum allowed sequence length to be included.",
 
        "cluster_blast" : "Specify a BLAST+ file to include the BLAST hit description with identity, alignment and e-value for each sequence header in the FASTA file.",

        "muco_fasta" : "Specify a FASTA file to be used for identifying sequences representing a defined taxon.",
 
        "muco_blast" : "Specify a BLAST+ output file to identify a defined taxon.",
 
        "muco_lookup" : "Specify a LOOKUP file to convert BLAST+ values. Necessary for GB/NCBI BLAST+ to convert accession codes to taxonomy names.",
 
        "subfolders" : "Also search for files in subfolders.",
 
        "outputfasta" : "Write out in FASTA format (default FASTQ).",
 
        "forwardreads" : "Use only forward reads from the set and ignore reverse reads; useful when reads do not overlap, and a variable region is located near the forward primer.",
 
        "reversereads" : "Use only reverse reads from the set and ignore forward reads; useful when reads do not overlap, and a variable region is located near the reverse primer.",
    
        "min_allowed_base" : "Specify the minimum allowed quality for a nucleotide within the trimmed section for a sequence to be retained.",
 
        "min_base_trimmed" : "Specify the minimum allowed quality for a nucleotide before a sequence is truncated at that position.",

        "outies" : "Allows overlap at the front of the sequences; the default is to have overlap at the end of the sequences. This occurs if sequences are not correctly trimmed and the DNA insert size is shorter than the generated sequences.",
        
        "hits" : "Specify the number of multiple best hits for a single query using BLAST+. At least 10 best hits is recommended when using partitioned GB/NCBI BLAST+ to build a common taxonomy.",
 
        "fastainput" : "Specify a clustered FASTA to produce a nohits file based on those reads that did not meet the defined thresholds.",
 
        "cluster_sequences" : "For each out, the centroid sequence is added to the output table.",
 
        "pick_fasta" : "Select the FASTA file that was used for clustering or for the BLAST+ identification",
 
        "pick_blast" : "Select a BLAST+ file.",
 
        "pick_uc" : "Specify a clustered (*.uc) file to map reads from clusters back to individual sequences.",
 
        "pick_lookup" : "Specify a LOOKUP file to convert GB/NCBI BLAST+ accession codes to taxonomy information.",
 
        "pick_search" : "Specify keywords to be searched in the BLAST+ output.",
 
        "pick_count" : "Specify a number of sequences to be picked; use 0 to pick all.",
 
        "pick_description" : "Writes the BLAST+ description for each picked sequence header.",
        
        "cluster_fasta" : "Specify a FASTA file before clustering to map reads from clusters back to individual sequences",
 
        "cluster_len" : "Specify a minimum allowed length for picking sequences.",
    
        "rep_fasta" : "Select the cleaned FASTA file that was used for clustering and/or BLAST+ identification",
 
        "rep_len" : "Specify a minimum allowed length for picking sequences.",
 
        "rep_group" : "Groups hits together for easier filtering.",
 
        "rep_sort" : "Sorts hits by abundance in descending order; common hits are at the beginning and rare hits are at the end.",
 
        "rep_pick" : "Specify how sequences are picked. 1. Sequentially, hits are picked by highest identity. 2. Equally covers all the sequences in the hit. 3. Randomly.",
 
        "simplify" : "Discards uncommon ranks that are not be available for the majority of reference sequences. The following ranks are included: kingdom, phylum, class, order, family, genus, species",
 
        "ranks" : "Output ranks for each classification.",
 
        "adapter" : "Illumina overhang adapter sequences can occur at the 3' ends of reads when read length exceeds the DNA insert size. This option will remove these overhangs.",
 
        "persample" : "Default picks sequences per taxon, but this allows sequences to be picked per sample per taxon.",

        "multifastq" : "Specify multiple FASTQ files to analyse the quality and primer distribution.",

        "cquality" : "Calculate quality distribution (slow process).",

        "cskip" : "To speed up calculations and to get estimates quickly pipeline uses every 10th sequence to analyse.",

        "cout" : "Summarize for each file statistics separately.",

        "kmerlen" : "Specify kmer size to be used to calculate occurrences.",

        "kmerpos" : "Specify at which position kmers should be counted. This is useful when having barcode at the front of the sequence and want to ignore barcode information and check primer information.",

        "kmershow" : "Specify how many kmers to display in the output, use 0 to display all of them.",

        "kmerfront" : "Only calculate kmer occurrences at the front of the sequence, where barcode and/or primers are located.",

        "resolution" : "Specify minimum taxonomic level to be included when constructing common taxonomy. This helps to reduce uncultured hits to improve output of the BLAST."

    }  

    # define file types for open and save dialogs
    file_type_bio = ("Bioinformatics files",("*.fasta","*.fna","*.qual","*.fq",".*fastq","*.cf","*.txt","*.csv","*.tsv","*.blast","*.uc","*.dmp","*.nhr","*.taxa"))

    file_type_all = ("all files","*.*")

    def __init__(self, master, controller):
        tk.Frame.__init__(self, master)
        config_exists = False
        if os.path.isfile("pipeline.ini"):
            config = ConfigParser.RawConfigParser()
            config.read('pipeline.ini')
            try:
                self.configuration['python_exe'] = config.get('EXECUTABLES', 'python_exe')
                self.configuration['vsearch_exe'] = config.get('EXECUTABLES', 'vsearch_exe')
                self.configuration['flash_exe'] = config.get('EXECUTABLES', 'flash_exe')
                self.configuration['blastn_exe'] = config.get('EXECUTABLES', 'blastn_exe')
                self.configuration['makeblastdb_exe'] = config.get('EXECUTABLES', 'makeblastdb_exe')
                self.configuration['blastdbcmd_exe'] = config.get('EXECUTABLES', 'blastdbcmd_exe')
                self.configuration['database_folder'] = config.get('FOLDERS', 'database_folder')
                self.configuration['gb_folder'] = config.get('FOLDERS', 'gb_folder')
                self.configuration['taxonomy_folder'] = config.get('FOLDERS', 'taxonomy_folder')
                config_exists = True
            except:
                config_exists = False
        if not config_exists:
            self.configuration['python_exe'] = "python"
            self.configuration['vsearch_exe'] = "vsearch"
            self.configuration['flash_exe'] = "flash"
            self.configuration['blastn_exe'] = "blastn"
            self.configuration['makeblastdb_exe'] = "makeblastdb"
            self.configuration['blastdbcmd_exe'] = "blastdbcmd"
            self.configuration['database_folder'] = "db/"
            self.configuration['gb_folder'] = "gb/"
            self.configuration['taxonomy_folder'] = "taxonomy/"
        self.controller = controller
        self.grid(padx = 15, pady = 15)
        self.create_widgets()

    def create_widgets(self):
        """Create the widgets for the frame."""
        raise NotImplementedError
        
    def run_command(self, command, **kwargs):
        #sub = tk.Tk()
        #termf = tk.Frame(sub, height=400, width=800)
        #termf.pack(fill=tk.BOTH, expand=tk.YES)
        #wid = termf.winfo_id()
        #os.system('xterm -into %d -geometry 120x20 -sb &' % wid)
        # change predefined command with user defined command or command location
        split = command.split(" ")
        if split[0] == "python":
            split[0] = self.configuration['python_exe']
        if split[0] == "blastn":
            split[0] = self.configuration['blastn_exe']
        if split[0] == "vsearch":
            split[0] = self.configuration['vsearch_exe']
        if split[0] == "makeblastdb":
            split[0] = self.configuration['makeblastdb_exe']
        if split[0] == "blastdbcmd":
            split[0] = self.configuration['blastdbcmd_exe']
        # do not change scripts that are python or located in another folder using relative or absolute path
        if not '/' in split[0] and not '\\' in split[0]:
            if split[0] != "python":
                # if file exists in same folder and is linux/macos, add ./ execution at the front of the command
                if os.path.isfile(split[0]) and os.name != "nt":
                    split[0] = "./" + split[0]
        command = " ".join(split)
        console = ConsoleGUI(command)
        output = kwargs.get('output', "")
        if len(output) > 0:
            console.add_output(output)
        console.mainloop()
        #os.system('xterm -geometry 120x20 -sb -e %s' % (command))

    def create_label(self, label, row, column, **kwargs):
        font = kwargs.get('font', 8)
        tk.Label(self, text = label, font = (None, font)).grid(row = row, column = column, padx = self.padx, pady = self.pady)

    def create_button(self, label, command, row, column, **kwargs):
        width = kwargs.get('width', 35)
        columnspan = kwargs.get('columnspan', 1)
        tk.Button(self, text = label, width = width, command = command).grid(row = row, column = column, columnspan = columnspan, padx = self.padx, pady = self.pady)

    def check_errors(self, error):
        if len(error) > 0:
            tkMessageBox.showwarning("Warning", "Following error(s) occured: \n%s" % "\n".join(error))
            return True
        else:
            return False

    def add_column(self):
        self._max_row = self._form_row
        self._form_row = 1
        self._form_col += 3

    def run_buttons(self):
        self._form_row = max(self._form_row, self._max_row) + 1
        tk.Button(self, text = "Back", width=15, command = lambda: self.controller.switch_frame(HomeFrame)).grid(row = self._form_row, column = 1 + self._form_col, padx = 5, pady = 5)
        tk.Button(self, text = "Run command", width=15, command = lambda: self.check_fields()).grid(row = self._form_row, column = 2 + self._form_col, padx = 5, pady = 5) 

    def form_title(self, label, **kwargs):
        _columnspan = kwargs.get('columnspan', 3)
        heading = kwargs.get('heading', False)
        if heading:
            tk.Label(self, text = label, font=(None, 16)).grid(row = self._form_row, column = 1 + self._form_col, columnspan = _columnspan)
        else:
            tk.Label(self, text = label, font=(None, 12)).grid(row = self._form_row, column = 1 + self._form_col, padx = 5, pady = 5, columnspan = _columnspan, sticky = "w")
        self._form_row += 1

    def form_label(self, label, **kwargs):
        _columnspan = kwargs.get('columnspan', 3)
        tk.Label(self, text = label, font=(None, 8)).grid(row = self._form_row, column = 1 + self._form_col, padx = 5, pady = 5, columnspan = _columnspan, sticky = "w")
        self._form_row += 1

    def form_input_openfolder(self, variable):
        _var = askdirectory(title = "Select a folder")
        if len(_var) > 0:
            variable.set(_var)

    def allowed_files_list(self, allowed_files):
        if allowed_files == "fasta":
            return (("FASTA files", ("*.fasta", "*.fna", "*.fasta.gz", "*.fna.gz")),("all files","*.*"))
        if allowed_files == "fasta+fastq":
            return (("FASTA/FASTQ files", ("*.fasta", "*.fna", "*.fasta.gz", "*.fna.gz", "*.fastq", "*.fq", "*.fastq.gz", "*.fq", "*.fq.gz")),("all files","*.*"))
        elif allowed_files == "qual":
            return (("QUALITY files", ("*.qual", "*.qual.gz")),("all files","*.*"))
        elif allowed_files == "fastq":
            return (("FASTQ files", ("*.fastq", "*.fq", "*.fastq.gz", "*.fq", "*.fq.gz")),("all files","*.*"))
        elif allowed_files == "uc":
            return (("Cluster output", "*.uc"),("all files","*.*"))
        elif allowed_files == "database":
            return (("BLAST database", "*.nhr"),("all files","*.*"))
        elif allowed_files == "taxa":
            return (("Taxonomy files", "*.taxa"),("all files","*.*"))
        elif allowed_files == "txt":
            return (("Text files", ("*.txt", "*.csv", "*.tsv")),("all files","*.*"))
        elif allowed_files == "blast":
            return (("BLAST output", "*.blast"),("all files","*.*"))
        elif allowed_files == "dmp":
            return (("Taxonomy files", ("*.dmp", "*.dmp.gz")),("all files","*.*"))
        else:
            return (("Bioinformatics files",("*.fasta","*.fna","*.qual","*.fq",".*fastq","*.cf","*.txt","*.csv","*.tsv","*.blast","*.uc","*.dmp","*.nhr","*.taxa")),("FASTA",("*.fasta","*.fna")),("FASTQ",("*.fastq","*.fq")),("sample sheet",("*.csv","*.tsv","*.txt")),("all files","*.*"))

    def form_input_openfiles(self, variable, allowed_files):
        _var = askopenfilenames(title = "Select a file", filetypes = self.allowed_files_list(allowed_files))
        if len(_var) > 0:
           variable.set(_var)

    def form_input_savefile(self, variable, allowed_files):
        _var = asksaveasfilename(title = "Save output into file", filetypes = self.allowed_files_list(allowed_files))
        if len(_var) > 0:
           variable.set(_var)

    def form_input_openfile(self, variable, allowed_files):
        _var = askopenfilename(title = "Select a file", filetypes = self.allowed_files_list(allowed_files))
        if len(_var) > 0:
           variable.set(_var)

    def form_get(self, variable_name):
        if variable_name in self.params:
            return self.params[variable_name].get()
        else:
            return ""

    def form_set(self, variable_name, variable_value):
        return self.params[variable_name].set(variable_value)

    def form_hasvalue(self, variable_name):
        if "%s_on" % variable_name in self.params:
            if self.params['%s_on' % variable_name].get() == 0:
                return False
        if variable_name in self.params:
            if self.params[variable_name].get():
                return True
            else:
                return False
        else:
            return False

    def get_command(self, program, params):
        command = [program]
        for i in range(0, len(params), 2):
            if params[i + 1][0] == "#":
                command.append('%s "%s"' % (params[i], params[i + 1][1:]))
            elif self.form_hasvalue(params[i + 1]):
                command.append('%s "%s"' % (params[i], self.form_get(params[i + 1])))
        return command

    def form_select(self, e):
        e.widget.select_range(0, 'end')     

    def form_delete(self, e): 
        e.widget.delete(0, 'end')

    def form_display_dbs(self):
        dbs = glob.glob("%s*.nhr" % self.configuration['database_folder']) + glob.glob("%s*.00.nhr" % self.configuration['gb_folder'])
        first_value = None
        tmp = 0
        for db in dbs:
            title = "Custom database"
            if "maarjam" in db.lower() or "pioneer" in db.lower():
                title = "MaarjAM database"
            elif "unite" in db.lower():
                title = "UNITE database"
            elif "nt" in db.lower():
                title = "GB/NCBI database"
                # only show first database of NCBI
                if "00" not in db:
                    continue
            value = db.replace(".nhr", "")
            shorten = value
            if len(value) > 40:
                shorten = value[:17] + "..." + value[-17:]
            if tmp == 0:
                first_value = value
                _input = self.form_input("Select BLAST+ database file (*.nhr)", "radiostr", "blastdbradio", box_label = "%s (%s)" % (shorten, title), value = value)
            else:
                _input = self.form_input("", "radiostr", "blastdbradio", box_label = "%s (%s)" % (shorten, title), value = value)
            _input.config(command = lambda: self.form_set('blastdb', ''))
            tmp += 1
        if tmp > 0:
            self.form_input("or select custom BLAST+ database file (*.nhr)", "file", 'blastdb', allowed_files = 'database')
        else:
            self.form_title("No installed database found, read manual to install database")
            self.form_input("Select custom BLAST+ database file (*.nhr)", "file", "blastdb", allowed_files = 'database')
        if first_value:
            self.form_set("blastdbradio", first_value)
        return tmp

    def form_display_dbs_fasta(self):
        dbs = glob.glob("%s/*.fasta" % self.configuration['database_folder']) + glob.glob("%s*.fna" % self.configuration['database_folder'])
        first_value = None
        tmp = 0
        for db in dbs:
            title = "Custom database"
            if "maarjam" in db.lower() or "pioneer" in db.lower():
                title = "MaarjAM database"
            elif "unite" in db.lower():
                title = "UNITE database"
            elif "nt" in db.lower():
                title = "GB/NCBI database"
                # only show first database of NCBI
                if "00" not in db:
                    continue
            value = db
            shorten = value
            if len(value) > 40:
                shorten = value[:17] + "..." + value[-17:]
            if tmp == 0:
                first_value = value
                _input = self.form_input("Select database file (*.fasta)", "radiostr", "dbradio", box_label = "%s (%s)" % (shorten, title), value = value)
            else:
                _input = self.form_input("", "radiostr", "dbradio", box_label = "%s (%s)" % (shorten, title), value = value)
            _input.config(command = lambda: self.form_set('db', ''))
            tmp += 1
        if tmp > 0:
            self.form_input("or select custom database file (*.fasta)", "file", 'db', allowed_files = 'fasta')
        else:
            self.form_title("No installed database found, read manual to install database")
            self.form_input("Select custom database file (*.fasta)", "file", "db", allowed_files = 'fasta')
        if first_value:
            self.form_set("dbradio", first_value)
        return tmp

    def form_input(self, label, input_type, variable_name, **kwargs):
        #if self.form_hasvalue(variable_name):
        #    default = self.form_get(variable_name)
        if input_type == "int" or input_type == "radio" or input_type == "check":
            default = kwargs.get('default', 0)
            if variable_name not in self.params:
                self.params[variable_name] = tk.IntVar()
            self.params[variable_name].set(default)
        else:
            default = kwargs.get('default', "")
            if variable_name not in self.params:
                self.params[variable_name] = tk.StringVar()
            self.params[variable_name].set(default)
        variable = self.params[variable_name]
        input_form = None
        tooltip = kwargs.get('tooltip', False)
        if variable_name in self.tooltips:
            tooltip = self.tooltips[variable_name]
        if tooltip and len(label) > 0:
            label = "%s ?" % label
        input_label = tk.Label(self, text = label)
        input_label.grid(row = self._form_row, column = 1 + self._form_col, padx = 5, pady = 5, sticky = "e")

        if input_type == "entry":
            input_form = tk.Entry(self, textvariable = variable, width = 25)
            input_active = kwargs.get('active', -1)
            if input_active >= 0:
                self.params["%s_on" % variable_name] = tk.IntVar()
                self.params["%s_on" % variable_name].set(input_active)
                checkbox = tk.Checkbutton(self, text = "", variable = self.params["%s_on" % variable_name], onvalue = 1, offvalue = 0)
                checkbox.grid(row = self._form_row, column = 3 + self._form_col, padx = 5, pady = 5, sticky = "w")
                CreateToolTip(checkbox, "If checked, this parameter will be used when running the command")
        elif input_type == "file":
            input_form = tk.Entry(self, textvariable = variable, width = 25)
            allowed_files = kwargs.get('allowed_files', 'all')
            tk.Button(self, text = "OPEN", command = lambda: self.form_input_openfile(variable, allowed_files)).grid(row = self._form_row, column = 3 + self._form_col)
        elif input_type == "files":
            input_form = tk.Entry(self, textvariable = variable, width = 25)
            allowed_files = kwargs.get('allowed_files', 'all')
            tk.Button(self, text = "OPEN", command = lambda: self.form_input_openfiles(variable, allowed_files)).grid(row = self._form_row, column = 3 + self._form_col)
        elif input_type == "savefile":
            input_form = tk.Entry(self, textvariable = variable, width = 25)
            allowed_files = kwargs.get('allowed_files', 'all')
            tk.Button(self, text = "SAVE", command = lambda: self.form_input_savefile(variable, allowed_files)).grid(row = self._form_row, column = 3 + self._form_col)
        elif input_type == "folder":
            input_form = tk.Entry(self, textvariable = variable, width = 25)
            tk.Button(self, text = "OPEN", command = lambda: self.form_input_openfolder(variable)).grid(row = self._form_row, column = 3 + self._form_col)
        elif input_type == "radio" or input_type == "radiostr":
            _label = kwargs.get('box_label', '')
            _value = kwargs.get('value', '')
            input_form = tk.Radiobutton(self, text = _label, variable = variable, value = _value)
        elif input_type == "check":
            _label = kwargs.get('box_label', '')
            _onvalue = kwargs.get('onvalue', 1)
            _offvalue = kwargs.get('offvalue', 0)
            input_form = tk.Checkbutton(self, text = _label, variable = variable, onvalue = _onvalue, offvalue = _offvalue)
        elif input_type == "dropdown":
            _label = kwargs.get('box_label', '')
            _values = kwargs.get('values', '')
            input_form = tk.OptionMenu(self, variable, *_values)
        elif input_type == "int":
            input_from = kwargs.get('label_from', 0)
            input_to = kwargs.get('label_to', 10)
            input_active = kwargs.get('active', 0)
            checkbox = kwargs.get('checkbox', "")
            variable.set(default)
            if variable_name == "threads":
                #tmp = psutil.cpu_count(logical = False) - 1
                tmp = multiprocessing.cpu_count() - 1
                if tmp > 1:
                    variable.set(tmp)
            input_form = tk.Spinbox(self, from_ = input_from, to = input_to, textvariable = variable, width = 25)
            if input_active >= 0:
                self.params["%s_on" % variable_name] = tk.IntVar()
                self.params["%s_on" % variable_name].set(input_active)
                checkbox = tk.Checkbutton(self, text = "", variable = self.params["%s_on" % variable_name], onvalue = 1, offvalue = 0)
                checkbox.grid(row = self._form_row, column = 3 + self._form_col, padx = 5, pady = 5, sticky = "w")
                CreateToolTip(checkbox, "If checked, this parameter will be used when running the command")
            #if len(checkbox) > 0:
            #    self.params[checkbox] = tk.BooleanVar
                #self.params[checkbox].set(True)
            #    tk.Checkbutton(self, text = "", variable = self.params[checkbox], onvalue = True, offvalue = False).grid(row = self._form_row, column = 3, padx = 5, pady = 5)
        if input_form is not None:
            if tooltip and len(label) > 0:
                CreateToolTip(input_label, tooltip)
            if input_type != "check" and input_type != "radio":
                input_form.grid(row = self._form_row, column = 2 + self._form_col, padx = 5, pady = 5, sticky = "we")
            else:
                input_form.grid(row = self._form_row, column = 2 + self._form_col, padx = 5, pady = 5, sticky = "w")
        self._form_row += 1
        """
        if os.name != "nt":
            input_form.bind('<Control-a>', self.form_select)
            input_form.bind('<Command-a>', self.form_select)
            input_form.bind('<Control-c>', self.form_select)
            input_form.bind('<Command-c>', self.form_select)
            input_form.bind('<Delete>', self.form_delete)
        """
        return input_form


class HomeFrame(BaseFrame):
    def create_widgets(self):
        """Create the base widgets for the frame."""
        tk.Label(self, text="Pipeline steps", font=(None, 16)).grid(row = 0, column = 0, padx = self.padx, pady = self.pady)
        tk.Label(self, text="Additional tools", font=(None, 16)).grid(row = 0, column = 1, padx = self.padx, pady = self.pady)
        left_btns = [
            ['1. Clean and quality filter reads', DemultiplexFrame],
            ['2. Combine paired-end reads', CombineFrame],
            ['Cluster reads (optional)', ClusterFrame],
            ['3. Remove chimeric reads', ChimeraFrame],
            ['4. Identify reads with BLAST+', BlastFrame],
            ['5. Generate pivot table from BLAST+', PivotFrame],
            ['Generate Pivot table from cluster (optional)', PivotClusterFrame],
            ['6. Pick sequences from BLAST+', PickBlastFrame],
            ['7. Pick sequences from cluster', PickClusterFrame],
            ['8. Pick representative sequences', PickRepresentativeFrame]
        ]
        right_btns = [
            ['Analyse FASTQ files', AnalyseFrame],
            ['Merge NCBI BLAST+ results', MergeFrame],
            ['Fetch taxonomy from NCBI BLAST+', LookupFrame],
            ["Build a BLAST+ database", DatabaseFrame],
            ["Convert FASTA->FASTQ", ConvertFASTAFrame],
            ["Convert FASTQ->FASTA", ConvertFASTQFrame],
            ["Autocorrect read strands", StrandFrame],
            ['Fetch references from NCBI BLAST+', FetchFrame],
            ['Convert between QIIME format', QiimeFrame],
            ['Convert NCBI BLAST+ to common taxonomy', TaxonomyFrame]
        ]

        tk.Frame(self, height=1, width = 600, bg="grey").grid(padx = 5, pady = 5, columnspan = 2)
        
        key = 0
        for btn in left_btns:
            self.create_button(btn[0], lambda btn = btn: self.controller.switch_frame(btn[1]), 2 + key, 0)
            key += 1
        key = 0
        for btn in right_btns:
            self.create_button(btn[0], lambda btn = btn: self.controller.switch_frame(btn[1]), 2 + key, 1)
            key += 1
        
        key = max(len(left_btns), len(right_btns)) + 2

        tk.Frame(self, height=1, width = 600, bg="grey").grid(row = key, padx = 5, pady = 5, columnspan = 2)
        
        self.create_button("Close", lambda: self.controller.close(), key + 3, 0)
        # try to open manual with browser, different methods for different OS
        if os.name == "nt":
            self.create_button("Open manual", lambda: os.system("start %s/manual/manual.html" % (os.path.realpath("."))), key + 1, 1)
        elif sys.platform == "darwin":
            self.create_button("Open manual", lambda: subprocess.Popen(["xdg-open", "%s/manual/manual.html" % (os.path.realpath("."))]), key + 1, 1)
        else:
            self.create_button("Open manual", lambda: webbrowser.open("file:///%s/manual/manual.html" % (os.path.realpath("."))), key + 1, 1)
        self.create_button("Configuration", lambda: self.controller.switch_frame(ConfigurationFrame), key + 3, 1)


class PivotFrame(BaseFrame):

    def check_fields(self):
        errors = []
        _output = os.path.splitext(self.form_get('fasta'))[0]
        #if os.path.isfile(self.form_get('uc')) is False and os.path.isfile(self.form_get('fasta')) is False:
        #    errors.append("FASTA or UC file is not selected or does not exist")
        if os.path.isfile(self.form_get('blastinput')) is False:
            errors.append("BLAST+ result file is not selected or does not exist")

        if self.check_errors(errors) is False:
            params = [
                "-f", "fasta",
                "-b", "blastinput",
                "-uc", "uc",
                "-i", "identity",
                "-l", "alignment",
                "-lbp", "alignmentbp",
                "-t", "#%s" % self.form_get('mode'),
                "-vqs", "qvariable_start",
                "-vqe", "qvariable_end",
                "-vs", "variable_start",
                "-ve", "variable_end",
                "-zeros", "zeros",
                "-lookup", "lookup",
                "-nohits", "nohits",
                "-proportions", "proportions",
                "-commontaxa", "commontaxa"
            ]

            self.run_command(" ".join(self.get_command("python py/pipeline_pivot_table.py", params)), output = "%s.i%s.a%s.tsv" % (self.form_get('blastinput'), self.form_get('identity'), self.form_get('alignment')))

    def create_widgets(self):
        self.form_title("Generate pivot table based on BLAST+ output", heading = True, columnspan = 6)

        self.form_label("Provides a results table with samples in rows and BLAST+ hits in columns (and a transpoed version)", columnspan = 6)
        self.form_label("If you used clustering, use the nonclustered FASTA, but provide the UC file, which allows to be mapped back to sequences when calculating sequence counts", columnspan = 6)

        self.form_title("Input parameters")
        self.form_input("Select nonclustered cleaned FASTA (*.fasta, optional)", "file", 'fasta', allowed_files = 'fasta')
        self.form_input("Select BLAST+ result file (*.blast)", "file", 'blastinput', allowed_files = 'blast')
        self.form_input("Select clustered UC file (*.uc, optional)", "file", 'uc', allowed_files = 'uc')
        self.form_input("Select hit lookup file (*.taxa, optional)", "file", 'lookup', allowed_files = 'taxa')
        
        self.form_title("Options")
        self.form_input("Alignment identity (%)", "int", 'identity', label_from = 0, label_to = 100, default = 97, active = -1)
        self.form_input("Alignment length (%)", "int", 'alignment', label_from = 0, label_to = 100, default = 95, active = -1)
        self.form_input("Alignment length (bp)", "int", 'alignmentbp', label_from = 0, label_to = 1000, default = 450, active = 0)
        self.form_input("Alignment query variable region start (bp)", "int", 'qvariable_start', label_from = 0, label_to = 1000, default = 70, active = 0)
        self.form_input("Alignment query variable region end (bp)", "int", 'qvariable_end', label_from = 0, label_to = 1000, default = 300, active = 0)
        self.form_input("Alignment reference variable region start (bp)", "int", 'variable_start', label_from = 0, label_to = 1000, default = 70, active = 0)
        self.form_input("Alignment reference variable region end (bp)", "int", 'variable_end', label_from = 0, label_to = 1000, default = 300, active = 0)
        self.add_column()
        self._form_row = 3
                
        self.form_title("Options")
        self.form_input("BLAST+ mode", "radio", 'mode', box_label = 'MaarjAM database (VTX)', value = 0)
        self.form_input("", "radio", 'mode', box_label = 'Full description', value = 1)
        self.form_input("", "radio", 'mode', box_label = 'GB/NCBI database', value = 2)
        self.form_input("", "radio", 'mode', box_label = 'UNITE database (SH)', value = 3)
        self.form_input("Write out zeros", "check", 'zeros', box_label = '', default = 1)
        self.form_input("Write out nohit file", "check", 'nohits', box_label = '', default = 1)
        self.form_title("Options for multiple best hits")
        self.form_input("Calculate hits based on proportions", "check", 'proportions', box_label = '', default = 0)
        self.form_input("Construct common taxonomy", "check", 'commontaxa', box_label = 'GB/NCBI BLAST+ results', default = 0)

        #self.form_title("GB/NCBI database options")
        #self.form_input("Taxonomy node file (nodes.dmp)", "file", 'taxonomynode', default = "taxonomy/nodes.dmp" if os.path.isfile("taxonomy/nodes.dmp") else "")
        #self.form_input("Taxonomy description file (names.dmp)", "file", 'taxonomynames', default = "taxonomy/names.dmp" if os.path.isfile("taxonomy/names.dmp") else "")
        #self.form_input("Taxonomy tree file (nucl_gb.accession2taxid)", "file", 'taxonomynucl', default = "taxonomy/nucl_gb.accession2taxid" if os.path.isfile("taxonomy/nucl_gb.accession2taxid") else "")

        self.run_buttons()


class ConfigurationFrame(BaseFrame):

    def save_fields(self):
        config = ConfigParser.RawConfigParser()
        tmp = ['python_exe', 'vsearch_exe', 'flash_exe', 'blastn_exe', 'makeblastdb_exe', 'blastdbcmd_exe', 'database_folder', 'gb_folder', 'taxonomy_folder']
        for i in tmp:
            self.configuration[i] = self.form_get(i)
        config.add_section('EXECUTABLES')
        config.set('EXECUTABLES', 'python_exe', self.form_get('python_exe'))
        config.set('EXECUTABLES', 'vsearch_exe', self.form_get('vsearch_exe'))
        config.set('EXECUTABLES', 'flash_exe', self.form_get('flash_exe'))
        config.set('EXECUTABLES', 'blastn_exe', self.form_get('blastn_exe'))
        config.set('EXECUTABLES', 'makeblastdb_exe', self.form_get('makeblastdb_exe'))
        config.set('EXECUTABLES', 'blastdbcmd_exe', self.form_get('blastdbcmd_exe'))
        config.add_section('FOLDERS')
        config.set('FOLDERS', 'database_folder', self.form_get('database_folder'))
        config.set('FOLDERS', 'gb_folder', self.form_get('gb_folder'))
        config.set('FOLDERS', 'taxonomy_folder', self.form_get('taxonomy_folder'))
        with open('pipeline.ini', 'w') as configfile:
            config.write(configfile)
        self.controller.switch_frame(HomeFrame)
        tkMessageBox.showinfo("Success", "Values successfully updated and saved")

    def reset_fields(self):
        config = ConfigParser.RawConfigParser()
        tmp = {
            'python_exe': 'python', 
            'vsearch_exe': 'vsearch', 
            'flash_exe': 'flash', 
            'blastn_exe': 'blastn', 
            'makeblastdb_exe': 'makeblastdb', 
            'blastdbcmd_exe': 'blastdbcmd', 
            'database_folder': 'db/', 
            'gb_folder': 'gb/', 
            'taxonomy_folder': 'taxonomy/'}
        for key in tmp:
            self.configuration[key] = tmp[key]
        config.add_section('EXECUTABLES')
        config.set('EXECUTABLES', 'python_exe', "python")
        config.set('EXECUTABLES', 'vsearch_exe', "vsearch")
        config.set('EXECUTABLES', 'flash_exe', "flash")
        config.set('EXECUTABLES', 'blastn_exe', "blastn")
        config.set('EXECUTABLES', 'makeblastdb_exe', "makeblastdb")
        config.set('EXECUTABLES', 'blastdbcmd_exe', "blastdbcmd")
        config.add_section('FOLDERS')
        config.set('FOLDERS', 'database_folder', "db/")
        config.set('FOLDERS', 'gb_folder', "gb/")
        config.set('FOLDERS', 'taxonomy_folder', "taxonomy/")
        with open('pipeline.ini', 'w') as configfile:
            config.write(configfile)
        self.controller.switch_frame(HomeFrame)
        tkMessageBox.showinfo("Success", "Values reset to default and saved")

    def create_widgets(self):
        self.form_title("Configure pipeline execution files and folders", heading = True)
        self.form_input("Define python executable", "file", "python_exe", default = self.configuration["python_exe"])
        self.form_input("Define vsearch executable", "file", "vsearch_exe", default = self.configuration["vsearch_exe"])
        self.form_input("Define FLASh executable", "file", "flash_exe", default = self.configuration["flash_exe"])
        self.form_input("Define BLAST+ blastn", "file", "blastn_exe", default = self.configuration["blastn_exe"])
        self.form_input("Define BLAST+ makeblastdb", "file", "makeblastdb_exe", default = self.configuration["makeblastdb_exe"])
        self.form_input("Define BLAST+ blastdbcmd", "file", "blastdbcmd_exe", default = self.configuration["blastdbcmd_exe"])
        self.form_input("Define database folder", "folder", "database_folder", default = self.configuration["database_folder"])
        self.form_input("Define gb folder", "folder", "gb_folder", default = self.configuration["gb_folder"])
        self.form_input("Define taxonomy folder", "folder", "taxonomy_folder", default = self.configuration["taxonomy_folder"])
        self._form_row = max(self._form_row, self._max_row) + 1
        tk.Button(self, text = "Back", width=15, command = lambda: self.controller.switch_frame(HomeFrame)).grid(row = self._form_row, column = 1 + self._form_col, padx = 5, pady = 5)
        tk.Button(self, text = "Reset fields", width=15, command = lambda: self.reset_fields()).grid(row = self._form_row, column = 2 + self._form_col, padx = 5, pady = 5) 
        tk.Button(self, text = "Save fields", width=15, command = lambda: self.save_fields()).grid(row = self._form_row, column = 3 + self._form_col, padx = 5, pady = 5) 


class MergeFrame(BaseFrame):

    def check_fields(self):
        errors = []
        if self.form_get("folder"):
            if not os.path.isdir(self.form_get("folder")):
                errors.append("Folder is not selected or is incorrect")

            if self.check_errors(errors) is False:
                _output = "%s/all.blast" % self.form_get("folder")
                self.run_command('python py/pipeline_merge_blasts.py -f "%s" %s > "%s"' % 
                    (
                    self.form_get("folder"),
                    "" if self.form_get('besthit') else "-allhits 1",
                    _output
                    ), output = _output)
        else:
            _files = self.form_get('multiblast').split("', '")
            files = []
            for _f in _files:
                # some errors with Python 2 and 3 reading multiple files different, only allow in Python 3?
                #files.append(_f.replace("('", "'").replace("',)", "'").replace("')", "'"))#.replace("',", "").replace("'", ""))
                files.append(_f.replace("('", "").replace("',)", "").replace("')", ""))#.replace("',", "").replace("'", ""))

            _output = "%s.all.blast" % files[0].replace(".blast", "").replace(".0", "")
            if os.path.isfile(files[0]) is False and os.path.isfile(files[-1]) is False :
                errors.append("Input files are not selected")

            if self.check_errors(errors) is False:
                self.run_command('python py/pipeline_merge_blasts.py -i "%s" %s > "%s"' %
                    (
                    ",".join(files),
                    "" if self.form_get('besthit') else "-allhits 1",
                    _output
                    ), output = _output)

    def create_widgets(self):
        self.form_title("Merge different BLAST+ results into one file", heading = True)

        self.form_label("When using partitioned databases, convert multiple BLAST+ results into one file for easier data manipulation")

        self.form_title("Input parameters")
        self.form_input("Select folder containing BLAST+ outputs", "folder", "folder")
        self.form_input("or select BLAST+ outputs (hold shift to select mulitple, *.blast)", "files", 'multiblast', allowed_files = 'blast')
        self.form_input("Select only best hit", "check", "besthit", default = 0)
        
        self.run_buttons()


class FetchFrame(BaseFrame):

    def check_fields(self):
        errors = []
        if len(self.form_get('blastdb')) == 0:
            self.form_set('blastdb', self.form_get('blastdbradio'))
        _database = self.form_get('blastdb').replace(".nhr", "").replace(".00", "")
        _files = self.tk.splitlist(self.form_get('fetchblast'))
        files = []
        for _f in _files:
            files.append(_f.replace("('", "").replace("',)", "").replace("')", "").replace("',", "").replace("'", ""))
        _output = "%s.reference.fasta" % files[0].replace(".blast", "").replace(".0", "")
        if os.path.isfile(files[0]) is False and os.path.isfile(files[-1]) is False :
            errors.append("Input files are not selected")
        if os.path.isfile("%s.nhr" % _database) is False:
            errors.append("Database is not selected, use *.nhr or any extension or write without extension")

        if self.check_errors(errors) is False:
            _cmd = self.configuration['blastdbcmd_exe']
            if not '/' in _cmd and not '\\' in _cmd and os.path.isfile(_cmd) and os.name != "nt":
                _cmd = "./" + _cmd
            self.run_command('python py/pipeline_fetch_blast_ids.py -i "%s" > "%s.ids" && %s -entry_batch "%s.ids" -db "%s" -dbtype nucl -out "%s"' %
                (
                ",".join(files),
                files[0],
                _cmd,
                files[0],
                _database,
                _output
                ), output = _output)

    def create_widgets(self):
        self.form_title("Fetch reference sequences from NCBI BLAST+ result", heading = True)
        
        self.form_label("Returns the reference sequences used to identify query sequences using BLAST+")

        self.form_title("Input parameters")
        self.form_input("Select NCBI BLAST+ outputs (hold shift to select mulitple, *.blast)", "files", 'fetchblast', allowed_files = 'blast')
        self.form_display_dbs()
        self.run_buttons()


class QiimeFrame(BaseFrame):

    def check_fields(self):
        errors = []
        if os.path.isfile(self.form_get('qiime')) is False:
            errors.append("Input file is not selected or does not exist")
        _output = "%s.converted.fasta" % self.form_get('qiime').replace(".fasta", "").replace(".fna", "")

        if self.check_errors(errors) is False:
            self.run_command('python py/pipeline_convert_qiime.py -i "%s" > "%s.ids"' %
                (
                self.form_get('qiime'),
                _output
                ), output = _output)

    def create_widgets(self):
        self.form_title("Convert FASTA file to be compatible with QIIME v1 or with the pipeline", heading = True)

        self.form_label("Users can convert cleaned FASTA file between current pipeline and QIIME v1 pipeline,")
        self.form_label("allowing the converted file to be used in either of the pipelines for further analyses.")

        self.form_title("Input parameters")
        self.form_input("Select cleaned FASTA file (auto-detects format, *.fasta)", "file", 'qiime', allowed_files = 'fasta')
        self.run_buttons()


class TaxonomyFrame(BaseFrame):

    def check_fields(self):
        errors = []
        if os.path.isfile(self.form_get('blast')) is False:
            errors.append("BLAST+ file is not selected or does not exist")
        if os.path.isfile(self.form_get('pick_lookup')) is False:
            errors.append("LOOKUP file is not selected or does not exist")
        _output = "%s.converted.blast" % self.form_get('blast').replace(".blast", "")

        if self.check_errors(errors) is False:
            self.run_command('python py/pipeline_convert_common_taxa.py -b "%s" -lookup "%s" -i "%s" -l "%s" -resolution "%s" > "%s"' %
                (
                self.form_get('blast'),
                self.form_get('pick_lookup'),
                self.form_get('identity'),
                self.form_get('alignment'),
                self.form_get('resolution'),
                _output
                ), output = _output)

    def create_widgets(self):
        self.form_title("Convert GB/NCBI BLAST+ to common taxonomy", heading = True)

        self.form_label("Users can convert multiple best hits to include only single hit as common taxonomy.")

        self.form_title("Input parameters")
        self.form_input("Select BLAST+ file (*.blast)", "file", 'blast', allowed_files = 'blast')
        self.form_input("Select hit lookup file (*.taxa, optional)", "file", 'pick_lookup', allowed_files = 'taxa')
        self.form_input("Alignment identity (%)", "int", 'identity', label_from = 0, label_to = 100, default = 95, active = -1)
        self.form_input("Alignment length (%)", "int", 'alignment', label_from = 0, label_to = 100, default = 90, active = -1)
        self.form_input("Minimum taxonomy resolution", "int", 'resolution', label_from = 0, label_to = 100, default = 3, active = -1)
        self.run_buttons()

class DatabaseFrame(BaseFrame):

    def check_fields(self):
        errors = []        
        if len(self.form_get('db')) == 0:
            self.form_set('db', self.form_get('dbradio'))
        _database = self.form_get('db')
        _output = os.path.splitext(_database)[0]
        if os.path.isfile(_database) is False:
            errors.append("Database input file is not selected")

        if self.check_errors(errors) is False:
            self.run_command('makeblastdb -in "%s" -dbtype nucl -out "%s"' %
                (
                _database,
                _output
                ), output = _output)

    def create_widgets(self):
        self.form_title("Convert FASTA file to BLAST+ database", heading = True)

        self.form_label("Converts FASTA into BLAST+ database, which can be used to identify sequences")
        
        self.form_title("Input parameters")
        self.form_display_dbs_fasta()
        
        self.run_buttons()


class StrandFrame(BaseFrame):

    def check_fields(self):
        errors = []
        if len(self.form_get('blastdb')) == 0:
            self.form_set('blastdb', self.form_get('blastdbradio'))
        _database = self.form_get('blastdb').replace(".nhr", "")
        _output = os.path.splitext(self.form_get('input'))[0]
        if os.path.isfile(self.form_get('input')) is False:
            errors.append("FASTA file is not selected")
        if os.path.isfile(self.form_get('blast')) is False and os.path.isfile("%s.nhr" % _database) is False:
            errors.append("BLAST+ database or output is not selected, one of these has to be selected")

        if self.check_errors(errors) is False:
            if os.path.isfile(self.form_get('blast')):
                self.run_command('python py/pipeline_correct_direction.py -f "%s" -b "%s"' %
                    (
                    self.form_get('input'),
                    self.form_get('blast') 
                    ), output = "%s.cs.fasta" % self.form_get('input').replace("fasta", "").replace("fna", ""))
            else:
                self.run_command('blastn -query "%s" -dust no -evalue %s -max_target_seqs 1 -num_threads %d -db "%s" -outfmt "6 qseqid sseqid stitle evalue pident nident length frames qstart qend sstart send qlen slen score" > tmp.blast && python py/pipeline_correct_direction.py -f "%s" -b tmp.blast' %
                    (self.form_get('input'),
                    self.form_get('evalue'),
                    self.form_get('threads'),
                    _database,
                    _output
                    ), output = "%s.cs.fasta" % self.form_get('input').replace("fasta", "").replace("fna", ""))

    def create_widgets(self):
        self.form_title("Correct FASTA file strand orientation", heading = True)

        self.form_label("Uses reference database or BLAST+ results to correct strands using reverse complement")
        self.form_label("Clustering and chimera checking needs to have sequences in correct strands")
        
        self.form_title("Input parameters")
        self.form_input("Select FASTA file (*.fasta)", "file", 'input', allowed_files = 'fasta')
        self.form_title("Select database or BLAST+ result file")
        self.form_display_dbs()
        self.form_input("E-value", "entry", 'evalue', default = "1e-50")
        self.form_input("Number of threads", "int", 'threads', label_from = 0, label_to = 64, default = 1, active = -1)
        self.form_input("Select BLAST+ result file (*.blast)", "file", 'blast', allowed_files = 'blast')
        
        self.run_buttons()


class ConvertFASTAFrame(BaseFrame):

    def check_fields(self):
        errors = []
        _output = os.path.splitext(self.form_get('input'))[0]
        if os.path.isfile(self.form_get('input')) is False:
            errors.append("FASTA file is not selected")

        if self.check_errors(errors) is False:
            command = ['python', 'py/pipeline_fasta_fastq.py',
                        '-f "%s"' % self.form_get('input'),
                        '-phred %d' % self.form_get('phred')]
            if len(self.form_get('qual')) > 0:
                command.append('-qf "%s"' % self.form_get('qual'))

            self.run_command(" ".join(command), output = self.form_get('input').replace(".fasta", ".fastq"))

    def create_widgets(self):
        self.form_title("Convert FASTA to FASTQ", heading = True)

        self.form_label("Converts FASTA or FASTA and QUAL file to FASTQ")
        
        self.form_title("Input parameters")
        self.form_input("Select FASTA file (*.fasta)", "file", 'input', allowed_files = 'fasta')
        self.form_input("Select QUAL file (*.qual)", "file", 'qual', allowed_files = 'qual')
        self.form_input("Phred score", "int", 'phred', label_from = 33, label_to = 64, default = 33, active = -1)
        # TODO add options for input file manipulations: ADAPTER REMOVAL, PRIMER REMOVAL, QUALITY CHECKING
        
        self.run_buttons()


class ConvertFASTQFrame(BaseFrame):

    def check_fields(self):
        errors = []
        _output = os.path.splitext(self.form_get('fqinput'))[0]
        if os.path.isfile(self.form_get('fqinput')) is False:
            errors.append("FASTQ file is not selected")

        if self.check_errors(errors) is False:
            command = ['python', 'py/pipeline_fastq_fasta.py',
                        '-f "%s"' % self.form_get('fqinput'), 
                        '-phred %d' % self.form_get('phred')]
            if self.params['output_qual'].get() == 1:
                command.append('-output_qual 1')

            self.run_command(" ".join(command), output = self.form_get('fqinput').replace("fastq", "fasta"))

    def create_widgets(self):
        self.form_title("Convert FASTQ to FASTA", heading = True)

        self.form_label("Converts FASTQ to FASTA or FASTA and QUAL file, which can be used for clustering, chimera checking and identification")
                
        self.form_title("Input parameters")
        self.form_input("Select fastq read (*.fastq)", "file", 'fqinput', allowed_files = 'fastq')
        self.form_input("Phred score", "int", 'phred', label_from = 33, label_to = 64, default = 33, active = -1)
        self.form_input("Output QUAL file", "check", 'output_qual', box_label = "")
        # TODO add options for input file manipulations: ADAPTER REMOVAL, PRIMER REMOVAL, QUALITY CHECKING
        
        self.run_buttons()


class DereplicateFrame(BaseFrame):

    def check_fields(self):
        errors = []
        _output = os.path.splitext(self.form_get('input'))[0]
        if os.path.isfile(self.form_get('input')) is False:
            errors.append("FASTA file is not selected")

        if self.check_errors(errors) is False:
            _mode = ['derep_fulllength', 'derep_prefix']
            self.run_command('vsearch --%s "%s" --output "%s.dereplicated.fasta" --threads %s' %
                (_mode[self.form_get('demode')],
                self.form_get('input'),
                _output,
                self.form_get('threads')
                ), output = "%s.dereplicated.fasta" % (_output))

    def create_widgets(self):
        self.form_title("Dereplicate reads", heading = True)
        
        self.form_title("Input parameters")
        self.form_input("Select FASTA file (*.fasta)", "file", 'input', allowed_files = 'fasta')
        
        self.form_title("Options")
        self.form_input("Number of threads", "int", 'threads', label_from = 0, label_to = 64, default = 1, active = -1)
        self.form_input("Dereplicate mode", "radio", 'demode', box_label = 'full length', value = 0)
        self.form_input("", "radio", 'demode', box_label = 'prefix', value = 1)           
        self.run_buttons()


class ClusterFrame(BaseFrame):

    def check_fields(self):
        errors = []
        _output = os.path.splitext(self.form_get('input'))[0]
        if os.path.isfile(self.form_get('input')) is False:
            errors.append("FASTA file is not selected")

        if self.check_errors(errors) is False:
            _mode = ['cluster_fast', 'cluster_size']
            _identity = self.form_get('identity') / 100.0
            self.run_command('vsearch --%s "%s" --id %f --centroids "%s.cluster%s.fasta" --threads %s --qmask none --uc "%s.cluster%s.uc" && python py/pipeline_convert_uc.py -uc "%s.cluster%s.uc" -sort t > "%s.cluster%s.clusters"' %
                (_mode[self.form_get('clmode')],
                self.form_get('input'),
                _identity,
                _output,
                self.form_get('identity'),
                self.form_get('threads'),
                _output,
                self.form_get('identity'),
                _output,
                self.form_get('identity'),
                _output,
                self.form_get('identity')
                ), output = "%s.cluster%s.fasta" % (_output, self.form_get('identity')))

    def create_widgets(self):
        self.form_title("Cluster reads", heading = True)

        self.form_label("Clustering reads is a way of reducing the computational complexity and time taken to run subsequent steps")
        self.form_label("FASTA format required (use FASTQ->FASTA conversion if needed)")
         
        self.form_title("Input parameters")
        self.form_input("Select FASTA file (*.fasta)", "file", 'input', allowed_files = 'fasta')
        
        self.form_title("Options")
        self.form_input("Identity (%)", "int", 'identity', label_from = 0, label_to = 100, default = 97, active = -1)
        self.form_input("Number of threads", "int", 'threads', label_from = 0, label_to = 64, default = 1, active = -1)
        self.form_input("Cluster mode", "radio", 'clmode', box_label = 'cluster_fast', value = 0)
        self.form_input("", "radio", 'clmode', box_label = 'cluster_size', value = 1) 
        self.run_buttons()


class PickBlastFrame(BaseFrame):

    def check_fields(self):
        errors = []
        _output = self.form_get('pick_fasta').replace(".fasta", "")
        if os.path.isfile(self.form_get('pick_fasta')) is False:
            errors.append("FASTA file is not selected")
        if os.path.isfile(self.form_get('pick_blast')) is False:
            errors.append("BLAST file is not selected")

        if self.check_errors(errors) is False:

            command = ['python', 'py/pipeline_pick_sequences_from_blast.py',
                        '-f "%s"' % self.form_get('pick_fasta'),
                        '-b "%s"' % self.form_get('pick_blast'),
                        '-i %d' % self.form_get('identity'),
                        '-a %d' % self.form_get('alignment'),
                        '-c %d' % self.form_get('pick_count'),
                        '-s "%s"' % self.form_get('pick_search')
                        ]
            if self.form_get('pick_lookup'):
                command.append('-l "%s"' % self.form_get('pick_lookup'))
            if self.form_get('pick_uc'):
                command.append('-uc "%s"' % self.form_get('pick_uc'))

            if self.form_get('pick_description'):
                command.append('-d 1')

            command.append('> "%s.filtered.fasta"' % _output)
            self.run_command(" ".join(command), output = "%s.filtered.fasta" % _output)

    def create_widgets(self):
        self.form_title("Pick sequences from BLAST+", heading = True)

        self.form_label("Pick sequences with a given identity and taxa on the basis of BLAST+ results")
                
        self.form_title("Input parameters")
        self.form_input("Select nonclustered FASTA file (*.fasta)", "file", 'pick_fasta', allowed_files = 'fasta')
        self.form_input("Select BLAST+ file (*.blast)", "file", 'pick_blast', allowed_files = 'blast')
        self.form_input("Select clustered UC file (*.uc, optional)", "file", 'pick_uc', allowed_files = 'uc')
        # if lookup is not defined, hit names are looked from BLAST results
        self.form_input("Select hit lookup file (*.taxa, optional)", "file", 'pick_lookup', allowed_files = 'taxa')
        
        self.form_title("Options")
        self.form_input("Return the following taxa, separate by comma", "entry", 'pick_search', default = 'Glomeromycetes')
        self.form_input("Identity (%)", "int", 'identity', label_from = 0, label_to = 100, default = 90, active = -1)
        self.form_input("Alignment (%)", "int", 'alignment', label_from = 0, label_to = 100, default = 90, active = -1)
        self.form_input("Maximum number of occurances", "int", 'pick_count', label_from = 0, label_to = 100, default = 0, active = -1)
        self.form_input("Write out BLAST description", "check", 'pick_description', box_label = '', default = 1)
        
        self.run_buttons()


class PickClusterFrame(BaseFrame):

    def check_fields(self):
        errors = []
        _output = os.path.splitext(self.form_get('cluster_uc'))[0]
        if os.path.isfile(self.form_get('cluster_fasta')) is False:
            errors.append("FASTA file is not selected")
        if os.path.isfile(self.form_get('cluster_uc')) is False:
            errors.append("Cluster file is not selected")

        if self.check_errors(errors) is False:

            command = ['python', 'py/pipeline_pick_sequences_from_cluster.py',
                        '-f "%s"' % self.form_get('cluster_fasta'),
                        '-c "%s"' % self.form_get('cluster_uc'),
                        '-s %d' % self.form_get('count'),
                        '-m %d' % self.form_get('cluster_size'),
                        '-ml %d' % self.form_get('cluster_len'),
                        '-mi %d' % self.form_get('identity'),
                        '-pick %d' % self.form_get('cluster_pick')
                        ]
            if self.form_get('cluster_first'):
                command.append('-centroid 1')
            else:
                command.append('-centroid ""')
            if self.form_get('cluster_group'):
                command.append('-group 1')
            else:
                command.append('-group ""')
            if self.form_get('cluster_sort'):
                command.append('-sort 1')
            else:
                command.append('-sort ""')
            if self.form_get('cluster_blast') and os.path.isfile(self.form_get('cluster_blast')):
                command.append('-blast "%s"' % self.form_get('cluster_blast'))

            command.append('> "%s.cpicked.fasta"' % _output)
            self.run_command(" ".join(command), output = "%s.cpicked.fasta" % _output)

    def create_widgets(self):
        self.form_title("Pick sequences from cluster", heading = True)
        
        self.form_label("Pick sequences from a clustered FASTA file for use in further analysis (e.g. phylogenetic analysis)")
        
        self.form_title("Input parameters")
        self.form_input("Select nonclustered FASTA file (*.fasta)", "file", 'cluster_fasta', allowed_files = 'fasta')
        self.form_input("Select clustered UC file (*.uc)", "file", 'cluster_uc', allowed_files = 'uc')
        self.form_input("Select BLAST (*.blast, optional)", "file", 'cluster_blast', allowed_files = 'blast')
        
        self.form_title("Options")
        self.form_input("BLAST+ identity to show hit description (%)", "int", 'identity', label_from = 0, label_to = 100, default = 97, active = -1)
        self.form_input("Number of sequences to pick", "int", 'count', label_from = 0, label_to = 10000000, default = 4, active = -1)
        
        self.form_input("Minimum cluster size", "int", 'cluster_size', label_from = 0, label_to = 10000000, default = 100, active = -1)
        self.form_input("Minimum sequence length", "int", 'cluster_len', label_from = 0, label_to = 500, default = 100, active = -1)
        
        self.form_input("Pick centroid", "check", "cluster_first", box_label = "pick first item from cluster", default = 1)
        self.form_input("Group", "check", "cluster_group", box_label = "group clusters together for output", default = 1)
        self.form_input("Sort", "check", "cluster_sort", box_label = "sort clusters by size for output", default = 1)

        self.form_input("Pick sequences", "radio", "cluster_pick", box_label = "sequentially", value = 1)
        self.form_input("", "radio", "cluster_pick", box_label = "divided equally", value = 2)
        self.form_input("", "radio", "cluster_pick", box_label = "randomly", value = 3, default = 2)

        self.run_buttons()


class PickRepresentativeFrame(BaseFrame):

    def check_fields(self):
        errors = []
        _output = os.path.splitext(self.form_get('rep_fasta'))[0]
        if os.path.isfile(self.form_get('rep_fasta')) is False:
            errors.append("FASTA file is not selected")
        if os.path.isfile(self.form_get('rep_blast')) is False:
            errors.append("BLAST file is not selected")

        if self.check_errors(errors) is False:

            command = ['python', 'py/pipeline_pick_representative_sequences.py',
                        '-f "%s"' % self.form_get('rep_fasta'),
                        '-b "%s"' % self.form_get('rep_blast'),
                        '-s %d' % self.form_get('count'),
                        '-ml %d' % self.form_get('rep_len'),
                        '-i %d' % self.form_get('identity'),
                        '-a %d' % self.form_get('alignment'),
                        '-t %d' % self.form_get('mode'),
                        '-pick %d' % self.form_get('rep_pick')
                        ]
            if self.form_get('uc'):
                command.append('-uc "%s"' % self.form_get('uc'))
            if self.form_get('rep_group'):
                command.append('-group 1')
            else:
                command.append('-group ""')
            if self.form_get('rep_sort'):
                command.append('-sort 1')
            else:
                command.append('-sort ""')
            if self.form_get('persample'):
                command.append('-persample 1')
            else:
                command.append('-persample ""')
            if self.form_get('outputfasta'):
                command.append('-outputfasta 1')
                command.append('> "%s.cpicked.fasta"' % _output)
                self.run_command(" ".join(command), output = "%s.cpicked.fasta" % _output)
            else:
                command.append('-outputfasta ""')
                command.append('> "%s.cpicked.tsv"' % _output)
                self.run_command(" ".join(command), output = "%s.cpicked.tsv" % _output)


    def create_widgets(self):
        self.form_title("Pick representative sequences from BLAST+ results", heading = True, colspan = 6)
        
        self.form_label("Generates representative sequence set for publication or sequence archiving", colspan = 6)
        
        self.form_title("Input parameters")
        self.form_input("Select nonclustered FASTA file (*.fasta)", "file", 'rep_fasta', allowed_files = 'fasta')
        self.form_input("Select BLAST+ (*.blast)", "file", 'rep_blast', allowed_files = 'blast')
        self.form_input("Select clustered UC file (*.uc, optional)", "file", 'uc', allowed_files = 'uc')
        
        self.form_title("Options")
        self.form_input("Identity (%)", "int", 'identity', label_from = 0, label_to = 100, default = 97, active = -1)
        self.form_input("Alignment (%)", "int", 'alignment', label_from = 0, label_to = 100, default = 95, active = -1)
        self.form_input("Number of sequences to pick", "int", 'count', label_from = 0, label_to = 10000000, default = 2, active = -1)
        
        self.form_input("Minimum sequence length", "int", 'rep_len', label_from = 0, label_to = 500, default = 100, active = -1)
        
        self.add_column()
        self._form_row += 1

        self.form_title("Options")
        self.form_input("Group", "check", "rep_group", box_label = "group hits together", default = 1)
        self.form_input("Sort", "check", "rep_sort", box_label = "sort hits by abundance", default = 1)
        self.form_input("Pick", "check", "persample", box_label = "pick sequences per sample per taxon", default = 0)
        self.form_input("Output", "check", "outputfasta", box_label = "output as FASTA file, otherwise tabulated text", default = 0)
        
        
        self.form_input("BLAST+ mode", "radio", 'mode', box_label = 'MaarjAM database (VTX)', value = 0)
        self.form_input("", "radio", 'mode', box_label = 'Full description', value = 1)
        self.form_input("", "radio", 'mode', box_label = 'UNITE database (SH)', value = 2)

        # sequences are ordered and picked based on BLAST+ score
        self.form_input("Pick sequences", "radio", "rep_pick", box_label = "sequentially", value = 1)
        self.form_input("", "radio", "rep_pick", box_label = "divided equally", value = 2)
        self.form_input("", "radio", "rep_pick", box_label = "randomly", value = 3, default = 1)

        self.run_buttons()


class PivotClusterFrame(BaseFrame):

    def check_fields(self):
        errors = []
        _output = self.form_get('cluster_uc')
        if os.path.isfile(self.form_get('cluster_uc')) is False:
            errors.append("Cluster file not selected")
        if os.path.isfile(self.form_get('blastinput')) is False:
            errors.append("BLAST file is not selected")

        if self.check_errors(errors) is False:
            command = ['python', 'py/pipeline_cluster_pivot_table.py',
                        '-uc "%s"' % self.form_get('cluster_uc'),
                        '-b "%s"' % self.form_get('blastinput'),
                        '-i %s' % self.form_get('identity'),
                        '-l %s' % self.form_get('alignment')
                        ]
            if self.form_get('fastainput'):
                command.append('-f "%s"' % self.form_get('fastainput'))
            if self.form_get('cluster_sort'):
                command.append('-sort t')
            if self.form_get('lookup'):
                command.append('-lookup "%s"' % self.form_get('lookup'))
            if self.form_get('cluster_sequences'):
                command.append('-s 1')
            if self.form_get('commontaxa'):
                command.append('-commontaxa 1')

            command.append('> %s.tsv' % _output)
            self.run_command(" ".join(command), output = "%s.tsv" % _output)

    def create_widgets(self):
        self.form_title("Generate pivot table from cluster and BLAST+ result", heading = True)
        
        self.form_label("Provides a results table with samples in rows and clusters (OTUs) with BLAST+ hits in columns")
        self.form_label("If you want to combine clusters into hits, use option '5. Generate pivot table from BLAST+'")
        
        self.form_title("Input parameters")
        self.form_input("Select clustered UC file (*.uc)", "file", 'cluster_uc', allowed_files = 'uc')
        self.form_input("Select BLAST+ file (*.blast)", "file", 'blastinput', allowed_files = 'blast')
        self.form_input("Select clustered FASTA file (*.fasta, optional for nohits)", "file", 'fastainput', allowed_files = 'fasta')
        self.form_input("Select lookup file (*.taxa, optional)", "file", 'lookup', allowed_files = 'taxa')
        
        self.form_title("Options")
        self.form_input("Sort", "check", "cluster_sort", box_label = "sort clusters by size for output", default = 1)
        self.form_input("Alignment identity (%)", "int", 'identity', label_from = 0, label_to = 100, default = 97, active = -1)
        self.form_input("Alignment length (%)", "int", 'alignment', label_from = 0, label_to = 100, default = 95, active = -1)
        self.form_input("Construct common taxonomy based on multiple best hits", "check", 'commontaxa', box_label = 'GB/NCBI BLAST+ results', default = 0)
        self.form_input("Write out centroids for each OTUs at the end of the table", "check", "cluster_sequences", box_label = "", default = 0)

        self.run_buttons()


class LookupFrame(BaseFrame):

    def check_fields(self):
        errors = []
        _output = os.path.splitext(self.form_get('lookup_blast'))[0]
        if os.path.isfile(self.form_get('lookup_blast')) is False:
            errors.append("BLAST file is not selected")

        if self.check_errors(errors) is False:

            command = ['python', 'py/pipeline_build_taxa.py',
                        '-b "%s"' % self.form_get('lookup_blast'),
                        '-ti "%s"' % self.form_get('taxonomynucl'),
                        '-tt "%s"' % self.form_get('taxonomynames'),
                        '-tn "%s"' % self.form_get('taxonomynode'),
                        "-simplify 1" if self.form_get('simplify') else "",
                        "-ranks 1" if self.form_get('ranks') else "",
                        ]

            command.append('> "%s.taxa"' % _output)
            self.run_command(" ".join(command), output = self.form_get('lookup_blast').replace(".blast", ".taxa"))

    def create_widgets(self):
        self.form_title("Fetch taxa from GenBank/NCBI BLAST+", heading = True)
        
        self.form_label("Use the script in taxonomy/ folder to download correct files. Converts GB/NCBI identities into taxonomy tree")

        self.form_title("Input parameters")
        self.form_input("Select BLAST+ output (*.blast)", "file", 'lookup_blast', allowed_files = 'blast')

        self.form_title("Select taxonomy files")
        
        self.form_input("Taxonomy node file (nodes.dmp)", "file", 'taxonomynode', default = "%snodes.dmp" % self.configuration['taxonomy_folder'] if os.path.isfile("%snodes.dmp" % self.configuration['taxonomy_folder']) else "", allowed_files = 'dmp')
        self.form_input("Taxonomy description file (names.dmp)", "file", 'taxonomynames', default = "%snames.dmp" % self.configuration['taxonomy_folder'] if os.path.isfile("%snames.dmp" % self.configuration['taxonomy_folder']) else "", allowed_files = 'dmp')
        if os.path.isfile("%snucl_gb.accession2taxid.gz" % self.configuration['taxonomy_folder']):
            tmp_file = "%snucl_gb.accession2taxid.gz" % self.configuration['taxonomy_folder']
        elif os.path.isfile("%snucl_gb.accession2taxid" % self.configuration['taxonomy_folder']):
            tmp_file = "%snucl_gb.accession2taxid" % self.configuration['taxonomy_folder']
        else:
            tmp_file = ""
        self.form_input("Taxonomy tree file (nucl_gb.accession2taxid)", "file", "taxonomynucl", default = tmp_file, allowed_files = 'dmp')

        self.form_title("Output options")
        self.form_input("Use simplified output", "check", "simplify", box_label = "", default = 1)
        self.form_input("Display ranks for categories", "check", "ranks", box_label = "", default = 0)
        self.run_buttons()


class BlastFrame(BaseFrame):

    def check_fields(self):
        errors = []
        if len(self.form_get('blastdb')) == 0:
            self.form_set('blastdb', self.form_get('blastdbradio'))
        _database = self.form_get('blastdb').replace(".nhr", "")
        _output = self.form_get('input')
        if os.path.isfile(self.form_get('input')) is False:
            errors.append("Input file is not selected")
        if os.path.isfile("%s.nhr" % _database) is False:
            errors.append("Database is not selected, use *.nhr or any extension or write without extension")

        if self.check_errors(errors) is False:
            if self.form_get('ncbi'):
                """
                # run BLASTs in partitions / not compatible with Windows as limited parameter length
                # uses additional python script to run one by one partitions, can also continue
                file_name = os.path.basename(self.form_get('blastdb'))
                directory_name = os.path.dirname(self.form_get('blastdb'))
                tmp = file_name.replace(".nhr", "").split(".")
                tmp.pop()   # remove index
                partitions = glob.glob("%s/%s.*.nhr" % (directory_name, ".".join(tmp)))
                ids = {}
                ids_sorted = []
                for partition in partitions:
                    partition_i = int(partition.split(".")[-2])
                    ids[partition_i] = partition.replace(".nhr", "")
                    ids_sorted.append(partition_i)
                ids_sorted.sort()
                predefined_id = int(self.form_get('blastdb').replace(".nhr", "").split(".")[-1])
                commands = []
                for i in ids:
                    if predefined_id <= i:
                        # prebuild BLAST commands
                        _cmd = self.configuration['blastn_exe']
                        if not '/' in _cmd and not '\\' in _cmd and os.path.isfile(_cmd) and os.name != "nt":
                            _cmd = "./" + _cmd
                        commands.append('%s -query "%s" -dust no -evalue %s -max_target_seqs %s -num_threads %d -db "%s" -outfmt "6 qseqid sseqid stitle evalue pident nident length frames qstart qend sstart send qlen slen score" > "%s.%d.blast"' %
                        (_cmd,
                        self.form_get('input'),
                        self.form_get('evalue'),
                        self.form_get('hits'),
                        self.form_get('threads'),
                        ids[i],
                        _output,
                        i))
                self.run_command(" && ".join(commands), output = "%s.0.blast" % _output)
                """
                _cmd = self.configuration['blastn_exe']
                if not '/' in _cmd and not '\\' in _cmd and os.path.isfile(_cmd) and os.name != "nt":
                    _cmd = "./" + _cmd
                params = [
                    "-f", "input",
                    "-db", "blastdb",
                    "-hits", "hits",
                    "-evalue", "evalue",
                    "-t", "threads",
                    "-i", "identity",
                    "-py", "#%s" % self.configuration['python_exe'],
                    "-program", "#%s" % _cmd
                ]
                self.run_command(" ".join(self.get_command("python py/pipeline_runblast.py", params)), output = "%s.#.blast" % _output)

            else:
                self.run_command('blastn -query "%s" -dust no -evalue %s -max_target_seqs %s -num_threads %d -db "%s" -outfmt "6 qseqid sseqid stitle evalue pident nident length frames qstart qend sstart send qlen slen score" | %s py/pipeline_filter.py -f "%s" -i "%s" > "%s.blast"' %
                    (self.form_get('input'),
                    self.form_get('evalue'),
                    self.form_get('hits'),
                    self.form_get('threads'),
                    _database,
                    self.configuration['python_exe'],
                    self.form_get('input'),
                    self.form_get('identity'),
                    _output), output = "%s.blast" % _output)

    def create_widgets(self):
        self.form_title("Identify reads with BLAST+", heading = True)
        
        self.form_label("Use cleaned, combined or clustered reads for identification. FASTA format required (use FASTQ->FASTA conversion if needed)")
        
        
        self.form_title("Input parameters")
        self.form_input("Select FASTA file (*.fasta)", "file", 'input', allowed_files = 'fasta')
        
        self.form_title("Options")
        self.form_display_dbs()
        self.form_input("E-value", "entry", 'evalue', default = "1e-10")
        self.form_input("Identity threshold", "int", 'identity', label_from = 0, label_to = 100, default = 95, active = -1)
        self.form_input("Best hits", "int", 'hits', label_from = 1, label_to = 100, default = 1, active = -1)
        self.form_input("NCBI BLAST+ or partitioned database", "check", "ncbi", box_label = 'automatically run multiple partitions')
        self.form_input("Number of threads", "int", 'threads', label_from = 0, label_to = 64, default = 1, active = -1)
        
        self.run_buttons()


class CombineFrame(BaseFrame):

    def check_fields(self):
        errors = []
        if self.form_get('interleaved') == 1:
            if os.path.isfile(self.form_get('forward')) is False:
                errors.append("Interleaved file is not selected")
        else:
            if os.path.isfile(self.form_get('forward')) is False:
                errors.append("Forward read is not selected")
            if os.path.isfile(self.form_get('reverse')) is False:
                errors.append("Reverse read is not selected")

        _mismatch = 1.0 - self.form_get('identity') / 100.0
        _output = os.path.splitext(self.form_get('forward'))[0]

        if self.check_errors(errors) is False:
            _convert = ""
            if self.params['outputfasta'].get():
                _convert = " && python py/pipeline_fastq_fasta.py -fq %s.combined.fastq" % _output
            if self.params['interleaved'].get() == 1:
                self.run_command('flash %s -c -m %d -M %d -x %f -t %s --interleaved-input "%s" > "%s.combined.fastq"%s' % 
                    (
                    "-O" if self.form_get('outies') else "",
                    self.form_get('min_length'),
                    self.form_get('max_length'),
                    _mismatch,
                    self.form_get('threads'),
                    self.form_get('forward'),
                    _output,
                    _convert
                    ), output = "%s.combined.fastq" % _output)
            else:
                self.run_command('flash %s -c -m %d -M %d -x %f -t %s "%s" "%s" > "%s.combined.fastq"%s' % 
                    (
                    "-O" if self.form_get('outies') else "",
                    self.form_get('min_length'),
                    self.form_get('max_length'),
                    _mismatch,
                    self.form_get('threads'),
                    self.form_get('forward'),
                    self.form_get('reverse'),
                    _output,
                    _convert
                    ), output = "%s.combined.fastq" % _output)

    def create_widgets(self):
        self.form_title("Combine paired-end reads", heading = True)

        self.form_label("Select cleaned Illumina reads in interleaved format (i.e. *.demultiplexed.cleaned.fastq)")
        
        self.form_title("Input parameters")
        self.form_input("Select forward/interleaved reads (*.fastq)", "file", 'forward', allowed_files = 'fastq')
        self.form_input("Select reverse reads (*.fastq)", "file", 'reverse', allowed_files = 'fastq')
        
        self.form_title("Options")
        self.form_input("Input file type", "radio", 'interleaved', box_label = 'interleaved', value = 1)
        self.form_input("", "radio", 'interleaved', box_label = 'separate files', value = 0, default = 1)
        self.form_input("Minimum overlap (bp)", "int", 'min_length', label_from = 0, label_to = 1000, default = 10, active = -1)
        self.form_input("Maximum overlap (bp)", "int", 'max_length', label_from = 0, label_to = 1000, default = 300, active = -1)
        self.form_input("Overlap identity (%)", "int", 'identity', label_from = 0, label_to = 100, default = 75, active = -1)
        self.form_input("Number of threads", "int", 'threads', label_from = 0, label_to = 64, default = 1, active = -1)
        self.form_input("Allow outies alignments", "check", 'outies', box_label = '', default = 0)
        self.form_input("Output FASTA file", "check", 'outputfasta', box_label = '', default = 1)
        
        self.run_buttons()


class ChimeraFrame(BaseFrame):
    
    def check_fields(self):
        errors = []
        if len(self.form_get('db')) == 0:
            self.form_set('db', self.form_get('dbradio'))
        _database = self.form_get('db')
        _output = "%s.cf.fasta" % self.form_get('input').replace(".fna", "").replace(".fasta", "")

        if os.path.isfile(self.form_get('input')) is False:
            errors.append("FASTA file is not selected")
        if self.form_get('cmode') == 0 and os.path.isfile(_database) is False:
            errors.append("Reference database is not selected")

        if self.check_errors(errors) is False:
            _extracommand = ""
            if self.form_get('uc') and os.path.isfile(self.form_get('uc')):
                _extracommand = (" && python py/pipeline_uc_chimeras.py -f \"%s\" -uc \"%s\"" % 
                    (
                    _output,
                    self.form_get('uc')
                    ))
            if self.form_get('cmode') == 1:
                self.run_command('vsearch --uchime_denovo "%s" --threads %s --nonchimeras "%s" %s' % 
                    (
                    self.form_get('input'), 
                    self.form_get('threads'), 
                    _output,
                    "" if not _extracommand else _extracommand
                    ), output = _output)
            else:
                self.run_command('vsearch --uchime_ref "%s" --db "%s" --threads %s --nonchimeras "%s" %s' % 
                    (
                    self.form_get('input'), 
                    _database, 
                    self.form_get('threads'), 
                    _output,
                    "" if not _extracommand else _extracommand
                    ), output = _output)
    
    def create_widgets(self):
        self.form_title("Remove chimeric reads", heading = True)
        
        self.form_label("Select cleaned or combined reads (i.e. *.demultiplexed.cleaned.fastq or *.demulitplexed.cleaned.combined.fastq)")
        
        self.form_title("Input parameters")
        self.form_input("Select file (*.fasta)", "file", 'input', allowed_files = 'fasta')
        self.form_input("Select clustered UC file (*.uc, optional)", "file", 'uc', allowed_files = 'uc')
        
        self.form_title("De novo or reference")
        self.form_input("Chimera checking mode", "radio", 'cmode', box_label = 'de novo', value = 1)
        self.form_input("", "radio", 'cmode', box_label = 'reference', value = 0)
        self.form_display_dbs_fasta()
        self.form_input("Number of threads", "int", 'threads', label_from = 0, label_to = 64, default = 1, active = -1)
        
        self.run_buttons()


class AnalyseFrame(BaseFrame):

    def check_fields(self):
        errors = []

        files = None
        _output = None
        if self.form_get("folder"):
            if os.path.isdir(self.params['folder'].get()) is False:
                errors.append("Folder does not exist")
            else:
                _output = "%s/all.kmer%s.stats.txt" % (self.params['folder'].get(), self.params['kmerlen'].get())
        else:
            if not self.params['multifastq'].get():
                errors.append("Files not selected")
            else:
                _files = self.form_get('multifastq').split("', '")
                files = []
                for _f in _files:
                    if _output is None:
                        tmp = _f.replace("('", "").replace("',)", "").replace("')", "")
                        dir = os.path.dirname(tmp)
                        _output = "%s/all.kmer%s.stats.txt" % (dir, self.params['kmerlen'].get())
                    # some errors with Python 2 and 3 reading multiple files different, only allow in Python 3?
                    files.append(_f.replace("('", "").replace("',)", "").replace("')", ""))

        if self.check_errors(errors) is False:
            params = [
                    "-f", "folder",
                    "-skip", "cskip",
                    "-separatestats", "cout",
                    "-quality", "cquality",
                    "-kmerfront", "kmerfront",
                    "-kmerlen", "kmerlen",
                    "-kmerpos", "kmerpos",
                    "-kmershow", "kmershow"
                ]
            if files:
                params.append("-i")
                params.append("#%s" % ",".join(files))
            self.run_command("%s > %s" % (" ".join(self.get_command("python py/pipeline_statistics.py", params)), _output), output = _output)
        
    def create_widgets(self):
        self.form_title("Calculate the quality and barcodes/primers/kmers distribution of FASTQ files", heading = True)
        
        self.form_title("Select input folder or files")

        self.form_input("Select folder containing FASTQ files", "folder", "folder")
        self.form_input("or select FASTQ files separately (hold shift to select multiple, *.fastq)", "files", 'multifastq', allowed_files = 'fastq')
  
        self.form_input("Calculate quality", "check", 'cquality', box_label = '', default = 1)
        self.form_input("Check every 10th sequence", "check", 'cskip', box_label = '', default = 1)
        self.form_input("Write separate statistics for each file", "check", 'cout', box_label = '', default = 1)
        self.form_input("Check kmers at the front", "check", 'kmerfront', box_label = '', default = 1)
        self.form_input("Specify kmer position", "int", 'kmerpos', label_from = 0, label_to = 1000, default = 1, active = -1)
        self.form_input("Specify kmers length", "int", 'kmerlen', label_from = 0, label_to = 100, default = 20, active = -1)
        self.form_input("Specify number of kmers to show", "int", 'kmershow', label_from = 0, label_to = 1000, default = 50, active = -1)
        
        self.run_buttons()     

class DemultiplexFrame(BaseFrame):

    mode = 0

    def check_fields(self):
        errors = []
        if self.form_get('fprimer'):
            self.params['fprimer'].set(self.params['fprimer'].get().split("[")[0])
        if self.form_get('rprimer'):
            self.params['rprimer'].set(self.params['rprimer'].get().split("[")[0])
        if not self.form_get('adapter_on'):
            self.form_set('adapter', '')

        if self.mode == 1:
            # 454/IonTorrent sequences cleaning block
            if os.path.isfile(self.params['fasta'].get()) is False:
                errors.append("Input FASTA/FASTQ not found")
            if os.path.isfile(self.params['barcode'].get()) is False:
                errors.append("Sample sheet not found")
            if self.form_hasvalue("qual") and os.path.isfile(self.params['qual'].get()) is False:
                errors.append("Quality file not found")

            if self.check_errors(errors) is False:
                mismatch = 0
                if self.form_get("mismatch_barcode") == 1:
                    mismatch += 1
                if self.form_get("mismatch_primer") == 1:
                    mismatch += 2
                params = [
                        "-f", "fasta",
                        "-b", "barcode",
                        "-bs", "sample_col",
                        "-bb", "fbarcode_col",
                        "-bp", "fprimer_col",
                        "-q", "quality_avg",
                        "-ml", "fmin_len",
                        "-tl", "fmax_len",
                        "-trimq", "trimming_avg",
                        "-trimw", "trimming_window",
                        "-qf", "qual",
                        "-homopolymer", "homopolymer",
                        "-mismatch", "#%s" % mismatch,
                        "-allow_indel", "mismatch_indel",
                        "-adapter", "adapter"
                    ]
                if len(self.params['fprimer'].get()) > 0:
                    params.append("-primer")
                    params.append("#%s" % self.params['fprimer'].get())
                self.run_command(" ".join(self.get_command("python py/pipeline_clean_454.py", params)), output = "%s.cleaned.fasta" % self.form_get('fasta').replace(".fasta", "").replace(".fna", ""))
        elif self.mode == 2:
            # demultiplexed cleaning block
            if os.path.isdir(self.params['folder'].get()) is False:
                errors.append("Input folder not defined")
            if self.check_errors(errors) is False:
                params = [
                    "-folder", "folder",
                    "-q", "quality_avg",
                    "-forward_primer", "fprimer",
                    "-reverse_primer", "rprimer",
                    "-forward_trim", "f_len",
                    "-reverse_trim", "r_len",
                    "-ml", "fmin_len",
                    "-min_allowed_base", "min_allowed_base",
                    "-min_base_trimmed", "min_base_trimmed",
                    "-trimq", "trimming_avg",
                    "-trimw", "trimming_window",
                    "-homopolymer", "homopolymer",
                    "-mismatch", "mismatch",
                    "-allow_indel", "mismatch_indel",
                    "-ignore_bases", "ignore_bases",
                    "-adapter", "adapter"
                ]
                if self.form_hasvalue("subfolders"):
                    params.append("-subfolders")
                    params.append("#1")
                if self.form_hasvalue("primers_mixed"):
                    params.append("-primers_mixed")
                    params.append("#1")
                if self.form_hasvalue("outputfasta"):
                    params.append("-fasta")
                    params.append("#1")
                if self.form_hasvalue("forwardreads"):
                    params.append("-force_forward")
                    params.append("#1")
                elif self.form_hasvalue("reversereads"):
                    params.append("-force_reverse")
                    params.append("#1")
                self.run_command(" ".join(self.get_command("python py/pipeline_clean_demultiplexed.py", params)), output = "demultiplexed.cleaned.fastq")
        elif self.mode >= 3:
            # Illumina/PacBio sequences cleaning block
            if os.path.isfile(self.params['ffastq'].get()) is False:
                errors.append("Forward reads not selected")
            if os.path.isfile(self.params['barcode'].get()) is False:
                errors.append("Sample sheet not defined")
            if self.check_errors(errors) is False:
                params = [
                        "-fr", "ffastq",
                        "-rr", "rfastq",
                        "-b", "barcode",
                        "-bs", "sample_col",
                        "-bfb", "fbarcode_col",
                        "-q", "quality_avg",
                        "-fo", "foligo",
                        "-ro", "roligo",
                        "-forward_primer", "fprimer",
                        "-reverse_primer", "rprimer",
                        "-bfp", "fprimer_col",
                        "-brb", "rbarcode_col",
                        "-brp", "rprimer_col",
                        "-forward_trim", "f_len",
                        "-reverse_trim", "r_len",
                        "-ml", "fmin_len",
                        "-min_allowed_base", "min_allowed_base",
                        "-min_base_trimmed", "min_base_trimmed",
                        "-trimq", "trimming_avg",
                        "-trimw", "trimming_window",
                        "-homopolymer", "homopolymer",
                        "-mismatch", "mismatch",
                        "-allow_indel", "mismatch_indel",
                        "-ignore_bases", "ignore_bases",
                        "-adapter", "adapter"
                    ]
                if self.form_hasvalue("primers_mixed"):
                    params.append("-primers_mixed")
                    params.append("#1")
                self.run_command(" ".join(self.get_command("python py/pipeline_clean.py", params)), output = "%s.cleaned.fastq" % self.form_get('ffastq').replace(".fastq", "").replace(".fq", ""))

    def create_widgets(self):
        self.change_frame(0)

    def change_frame(self, mode):
        for widget in self.winfo_children():
            widget.destroy()

        self._form_row = 0
        self.form_title("Demultiplex and quality filter reads", heading = True)
        
        self.mode = mode
        if self.mode == 0:
            self.create_widgets_home()
        else:
            self.create_widgets_input()

        if self.mode == 0:
            self.create_button("Back", lambda: self.controller.switch_frame(HomeFrame), self._form_row, 1, width = 35)

    def create_widgets_home(self):
        self.form_title("Select input data format")
        self.create_label("Select demultiplexed, if you have samples divided into multiple .fastq or .fastq.gz files.", 2, 1)
        self.create_button("Demultiplexed FASTA/FASTQ (454/Illumina/PacBio) reads", lambda: self.change_frame(2), 3, 1, width = 45)
        self.create_label("Select single/paired-end if you have one or two .fastq or .fastq.gz raw files containing all the samples", 4, 1)
        self.create_button("Single/Paired-end FASTQ (Illumina/PacBio) reads", lambda: self.change_frame(4), 5, 1, width = 45)
        self.create_label("For 454 and IonTorrent use FASTA + QUAL or FASTQ", 6, 1)
        self.create_button("FASTA+QUAL/FASTQ (454/IonTorrent) reads", lambda: self.change_frame(1), 7, 1, width = 45)
        self._form_row = 8

    def create_widgets_inputblock(self):
        self.form_title("Input parameters")
        if self.mode == 1:
            self.form_input("Select FASTA/FASTQ file (*.fasta/*.fastq)", "file", 'fasta', allowed_files = 'fasta+fastq')
            self.form_input("Select QUAL file (*.qual)", "file", 'qual', allowed_files = 'qual')
        if self.mode == 2:
            self.form_input("Select folder with FASTQ files", "folder", 'folder')
        if self.mode >= 3:
            self.form_input("Select forward reads file (*.fastq)", "file", 'ffastq', allowed_files = 'fastq')
            self.form_input("Select forward oligos file (*.fastq)", "file", 'foligo', allowed_files = 'fastq')
        if self.mode == 4:
            self.form_input("Select reverse reads file (*.fastq)", "file", 'rfastq', allowed_files = 'fastq')
            self.form_input("Select reverse oligos file (*.fastq)", "file", 'roligo', allowed_files = 'fastq')            
        if self.mode != 2:
            self.form_input("Select sample sheet file (*.txt, *.tsv, *.csv)", "file", 'barcode', allowed_files = 'txt')
            self.create_button("Check sample sheet file", lambda: PreviewGUI(self.form_get("barcode")), self._form_row, 2, width = 18)
            self._form_row += 1

    def create_widgets_barcodeblock(self):

        # predefine some primers, [] indicates primer name and region, but is removed while executing command
        fprimers = [
            "CAGCCGCGGTAATTCCAGCT[WANDA-SSU]",
            "TTGGAGGGCAAGTCTGGTGCC[NS31-SSU]",
            "GTGARTCATCGAATCTTTG[fITS7-ITS]"
        ]
        rprimers = [
            "GAACCCAAACACTTTGGTTTCC[AML2-SSU]",
            "GTTTCCCGTAAGGCGCCGAA[AM1-SSU]",
            "TCCTCCGCTTATTGATATGC[ITS4-ITS]"
        ]
        self.form_title("Sample sheet and primer parameters")
        if self.mode != 2:
            self.form_input("Define sample column", "int", 'sample_col', label_from = 1, label_to = 20, default = 1, checkbox = 'sample_col_chk', checkbox_default = True, active = 1)
            self.form_input("Define forward barcode column", "int", 'fbarcode_col', label_from = 1, label_to = 20, default = 2, active = 1)
            self.form_input("Define forward primer column", "int", 'fprimer_col', label_from = 1, label_to = 20, default = 3, active = 1)
        if self.mode == 4:
            self.form_input("Define reverse barcode column", "int", 'rbarcode_col', label_from = 1, label_to = 20, default = 4)
            self.form_input("Define reverse primer column", "int", 'rprimer_col', label_from = 1, label_to = 20, default = 5)
        self.form_input("Select forward primer", "dropdown", 'fprimer', values = fprimers)
        self.form_input("or provide forward primer sequence", "entry", 'fprimer')
        if self.mode == 2 or self.mode == 4:
            self.form_input("Select reverse primer", "dropdown", 'rprimer', values = rprimers)
            self.form_input("or provide reverse primer sequence", "entry", 'rprimer')

    def create_widgets_filterblock(self):

        if self.mode != 2:
            self.add_column()
        self.form_title("Filtering parameters")
        if self.mode == 1:
            self.form_input("Average quality (q)", "int", 'quality_avg', label_from = 0, label_to = 41, default = 25, active = -1)
        else:
            self.form_input("Average quality (q)", "int", 'quality_avg', label_from = 0, label_to = 41, default = 30, active = -1)
        self.form_input("Minimum allowed quality for a base (q)", "int", 'min_allowed_base', label_from = 0, label_to = 41, default = 5)
        self.form_input("Minimum quality for a base before trimming the end (q)", "int", 'min_base_trimmed', label_from = 0, label_to = 41, default = 10)
        
        self.form_input("Average trimming window quality (q)", "int", 'trimming_avg', label_from = 0, label_to = 41, default = 20)
        self.form_input("Trimming window length (bp)", "int", 'trimming_window', label_from = 0, label_to = 1000, default = 50)
        if self.mode == 1:
            self.form_input("Minimum allowed read length", "int", 'fmin_len', label_from = 0, label_to = 1000, default = 170, active = 1)
            self.form_input("Maximum trimming length", "int", 'fmax_len', label_from = 0, label_to = 1000, default = 520, active = 1)
        else:
            self.form_input("Minimum allowed read length", "int", 'fmin_len', label_from = 0, label_to = 1000, default = 100)
        if self.mode > 1:
            self.form_input("Forward read trim length", "int", 'f_len', label_from = 0, label_to = 1000, default = 220)
        if self.mode == 2 or self.mode == 4:
            self.form_input("Reverse read trim length", "int", 'r_len', label_from = 0, label_to = 1000, default = 180)
        self.form_input("Truncate homopolymers to this length (bp)", "int", 'homopolymer', label_from = 4, label_to = 100, default = 8)
        
        if self.mode != 1:
            self.form_input("Ignore first bases", "int", 'ignore_bases', label_from = 0, label_to = 100, default = 3)
        
        self.form_input("Overhang adapter sequence", "entry", 'adapter', default = 'CTGTCTCTTAT', active = 0)

        if self.mode == 1:
            self.form_input("Allow mismatch", "check", 'mismatch_barcode', box_label = "barcode")
            self.form_input("", "check", 'mismatch_primer', box_label = "primer")
            #self.form_input("", "check", 'mismatch_indel', box_label = "indel")
        elif self.mode == 2:
            self.add_column()
            self.form_title("Filtering parameters")
            self.form_input("Allow mismatch", "check", 'mismatch', box_label = "primer")
            #self.form_input("", "check", 'mismatch_indel', box_label = "indel")
            self.form_input("Primers are mixed between forward and reverse", "check", 'primers_mixed', box_label = "")
            self.form_input("Check subfolders", "check", 'subfolders', box_label = "", default = 1)
            self.form_input("Output FASTA", "check", "outputfasta", box_label = "")
            self.form_input("Use only forward reads", "check", "forwardreads", box_label = "")
            self.form_input("Use only reverse reads", "check", "reversereads", box_label = "")
        else:
            self.form_input("Allow mismatch", "check", 'mismatch', box_label = "barcode and primer")
            #self.form_input("", "check", 'mismatch_indel', box_label = "indel")
            self.form_input("Primers are mixed between forward and reverse", "check", 'primers_mixed', box_label = "")

    def create_widgets_input(self):
        if self.mode == 2:
            self.create_widgets_inputblock()
            self.create_widgets_filterblock()
            self.create_widgets_barcodeblock()
        else:
            self.create_widgets_inputblock()
            self.create_widgets_barcodeblock()
            self.create_widgets_filterblock()

        self.run_buttons()


class PythonGUI(tk.Tk):
    """The main window of the GUI.

    Attributes:
      container (tk.Frame): The frame container for the sub-frames.
      frames (dict of tk.Frame): The available sub-frames.

    """

    def __init__(self):
        tk.Tk.__init__(self)
        self.title("gDAT pipeline")
        self.geometry("+0+5")
        try:
            self.iconbitmap("gdat.ico")
        except:
            sys.stderr.write("Cannot open gdat.ico\n")
        self._frame = None
        #default_font = font.nametofont("TkDefaultFont")
        #default_font.configure(family='Helvetica', size=12)
        #self.option_add("*Font", default_font)
        self.option_add("*Font", "lucida 10")
        self.switch_frame(HomeFrame)
        self.create_widgets()
        self.resizable(0, 0)

    def close(self):
        self.destroy()

    def create_widgets(self):
        """Create the widgets for the frame."""             
        #   Frame Container
        self.container = tk.Frame(self)
        #self.container.grid(row=0, column=0, sticky=tk.W+tk.E)

    def switch_frame(self, frame_class):
       """ Destroys current frame and replace it with new frame """
       new_frame = frame_class(self, self)
       if self._frame is not None:
          self._frame.destroy()
       self._frame = new_frame
       #self._frame.pack() 


class PreviewGUI(tk.Tk):

    filename = None
    text = None

    def __init__(self, filename):
        tk.Tk.__init__(self)
        self.title("File preview")
        self.filename = filename
        self.create_widgets()
        self.resizable(0, 0)

    def create_widgets(self):
        """Create the widgets for the frame."""     
        self.container = tk.Frame(self)
        self.text = tk.Text(self, width = 80, height = 20)
        self.text.pack()
        tk.Button(self, text="Close window", width = 70, command = lambda: self.kill_process()).pack()
        self.text.config(background="black", foreground="white")
        if os.path.isfile(self.filename):
            if self.filename.split(".")[-1].lower == "gz":
                fh = gzip.open(self.filename, 'r')
            else:
                fh = open(self.filename, 'r')
            with fh as f:
                i = 0
                for r in f:
                    i += 1
                    if i > 1000:
                        break
                    self.text.insert(tk.END, "%s" % r)
        else:
            self.text.insert(tk.END, "Error, file '%s' does not exist\n" % self.filename)

    def kill_process(self):
        self.destroy()


class ConsoleGUI(tk.Tk):

    process = None
    text = None
    queue = None
    command = None
    start_time = None
    output = None
    terminate = False
    fh = None

    def __init__(self, command):
        tk.Tk.__init__(self)
        self.start_time = time.time()
        self.fh = open("gdat.log", "a+")
        if self.fh:
            if sys.version_info[0] >= 3:
                command = command.encode()
            self.fh.write("[%s] Executing command '%s'\n" % (datetime.now().strftime("%d/%m/%Y %H:%M:%S"), command))
        print(command)
        self.title("Command output")
        self.command = command
        self._frame = None
        self.create_widgets()
        self.resizable(0, 0)

    def create_widgets(self):
        """Create the widgets for the frame."""             
        #   Frame Container
        self.container = tk.Frame(self)
        self.text = tk.Text(self, width = 80, height = 20)
        self.text.pack()
        tk.Button(self, text="Terminate process and close window", width = 70, command = lambda: self.kill_process()).pack()
        #self.protocol("WM_DELETE_WINDOW", self.kill_process())
        self.text.config(background="black", foreground="white")
        self.append_line("Executing command:\n%s" % self.command)
        self.append_line("Please wait for the 'process finished' text to close the window")
        if os.name == "nt" or sys.version_info[0] >= 3:
            if sys.version_info[0] >= 3:
                self.command = self.command.decode()
            self.process = subprocess.Popen(self.command, 
                stdout=subprocess.PIPE, 
                stderr=subprocess.PIPE,
                #encoding="utf8",
                shell=True
                #start_new_session=True
            )
        else:
            self.process = subprocess.Popen(self.command, 
                stdout=subprocess.PIPE, 
                stderr=subprocess.PIPE,
                shell=True,
                preexec_fn=os.setsid
            )

        self.queue = Queue()
        self.thread = Thread(target=self.readlines, args=(self.process, self.queue, 0))
        self.thread.start()

        self.queue = Queue()
        self.thread = Thread(target=self.readlines, args=(self.process, self.queue, 1))
        self.thread.start()

        self.after(100, self.updateLines)

    def kill_process(self):
        self.terminate = True
        if self.fh:
            self.fh.close()
        if self.process is not None and self.process.poll() is None:
            try:
                if os.name == "nt":
                    os.kill(self.process.pid, signal.CTRL_C_EVENT)
                else:
                    os.killpg(os.getpgid(self.process.pid), signal.SIGTERM)
                #self.process.kill()
            except OSError:
                pass
        else:
            self.destroy()

    def add_output(self, output):
        self.output = output

    def append_line(self, line):
        elapsed_time = time.time() - self.start_time
        m, s = divmod(elapsed_time, 60)
        h, m = divmod(m, 60)
        try:
            line = line.decode("utf-8")
        except:
            pass
        # disabling does not allow to copy text by the user
        #self.text.config(state=tk.NORMAL)
        self.text.insert(tk.END, "[%02d:%02d:%02d] %s\n" % (h, m, s, line.replace("\r", "").replace("\n", "")))
        #self.text.config(state=tk.DISABLED)

    def updateLines(self):
        if self.terminate:
            self.destroy()
            return

        while True:
            try:
                line = self.queue.get(False)
                if len(line) > 0:
                    if self.fh:
                        if sys.version_info[0] >= 3:
                            line = line.decode()
                        self.fh.write("[%s] %s\n" % (datetime.now().strftime("%d/%m/%Y %H:%M:%S"), line.strip()))
                    self.append_line(line)
            except Empty:
                break

        if self.process.poll() is None:
            self.after(100, self.updateLines)
        else:
            self.append_line("Process finished, you can close the window\n")
            if self.output != None:
                self.append_line("Output written into file: %s\n" % self.output)
                if self.fh:
                    self.fh.write("[%s] Output written into file '%s'\n" % (datetime.now().strftime("%d/%m/%Y %H:%M:%S"), self.output))
                    self.fh.close()
        self.text.see("end")

    def readlines(self, process, queue, type):
        while process.poll() is None:
            if type == 0:
                queue.put(process.stdout.readline())
            elif type == 1:
                queue.put(process.stderr.readline())


if __name__ == "__main__":
    app = PythonGUI()
    app.mainloop()
    exit()
