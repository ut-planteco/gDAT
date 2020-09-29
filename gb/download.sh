#!/bin/bash

for i in {00..27}; do echo "Downloading NCBI nt.$i"; wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.$i.tar.gz; tar xzvf nt.$i.tar.gz; rm nt.$i.tar.gz; done
