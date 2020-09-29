#!/bin/bash

wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz 
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
tar xzvf taxdump.tar.gz
# no need to unpack accession2taxid, script can read packed files
#tar xzvf nucl_gb.accession2taxid.gz
