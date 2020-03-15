try:
    from urllib.request import urlretrieve
except:
    from urllib import urlretrieve

files = ["ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz", "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/gi_taxid_nucl.dmp.gz"]

for f in files:
	print("Downloading file %s" % f)
	urlretrieve(f, f.split("/")[-1])
	print("File %s downloaded" % f)
