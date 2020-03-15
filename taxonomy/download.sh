#!/bin/bash

wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz 
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/gi_taxid_nucl.dmp.gz
tar xzvf taxdump.tar.gz 
tar xzvf gi_taxid_nucl.dmp.gz
