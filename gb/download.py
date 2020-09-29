try:
    from urllib.request import urlretrieve
except:
    from urllib import urlretrieve

for i in range(0, 27):
	nr = "0" + str(i) if i < 10 else i
	print("Downloading file nt.%s.tar.gz" % nr)
	urlretrieve("ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.%s.tar.gz" % nr, 'nt.%s.tar.gz' % nr)
	print("File nt.%s.tar.gz downloaded" % nr)
