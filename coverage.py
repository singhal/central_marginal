import gzip
import re
import sys
import glob

quals = 0
cov = 0
num_sites_q = 0
num_sites_c = 0

files = glob.glob("/Users/singhal/Dropbox (Personal)/publications/Center_Marginal/data/variants/keep/*gz")

for file in files:
	print(file)
	v = gzip.open(file)
	for l in v:
		l = l.decode("utf-8")

		if not re.search('#', l) and not re.search('INDEL', l):
			d = re.split("\t", l.rstrip())

			if len(re.split(',', d[4])) == 1:
				num_sites_q += 1
				quals += float(d[5])
				dp = int(re.search("DP=(\d+)", d[7]).group(1))


				num_inds = 0
				genos = d[9:]
				genos = [re.search('^(\S\/\S)', x).group(1) for x in genos]
				for geno in genos:
					if geno in ['0/0', '0/1', '1/1']:
						num_inds += 1

				if num_inds > 0:
					dp = dp / float(num_inds)
					cov += dp
					num_sites_c += 1

print(quals)
print(cov)
print(num_sites_q)
print(num_sites_c)