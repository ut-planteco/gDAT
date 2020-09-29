try:
    from urllib.request import urlretrieve
except:
    from urllib import urlretrieve

files = ["ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz", "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz"]

for f in files:
	print("Downloading file %s" % f)
	urlretrieve(f, f.split("/")[-1])
	print("File %s downloaded" % f)
