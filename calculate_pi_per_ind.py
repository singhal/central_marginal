import re
import subprocess
import pandas as pd
import os
import numpy as np
import argparse
import gzip

parser = argparse.ArgumentParser(description="Calculate het per individual")
parser.add_argument('--vcf', help="VCF file for which to run.")
parser.add_argument('--outdir', help="directory to dump results")
args = parser.parse_args()
vcf = args.vcf

def get_diversity(vcffile, outdir):

	f = gzip.open(vcffile, 'r')

	for l in f:
		l = l.decode('utf-8')
		if re.search('^#CHROM', l):
			d = re.split('\t', l.rstrip())
			inds = d[9:]

			# initialize counter
			hets = dict([(ind, {'diff': 0, 'denom': 0}) for ind in inds])

		else:
			if not re.search('#', l) and not re.search('INDEL', l):
				d = re.split('\s+', l.rstrip())
				# don't mess with multiallelics
				if len(re.split(',', d[4])) == 1:
					genos = d[9:]
					genos = [re.search('^(\S\/\S)', x).group(1) for x in genos]

					for ind, geno in zip(inds, genos):
						if geno in ['0/0', '0/1', '1/1']:
							hets[ind]['denom'] += 1
							if geno == '0/1':
								hets[ind]['diff'] += 1
	f.close()

	stem = re.sub('^.*/', '', vcffile)
	stem = re.sub('.vcf.gz', '.diversity.csv')
	out_file = os.path.join(outdir, stem)

	o = open(out_file, 'w')
	o.write('cluster,individual,denominator,pi\n')
	for ind in hets:
		if hets[ind]['denom'] > 0:
			pi = hets[ind]['diff'] / float(hets[ind]['denom'])
		else:
			pi = np.nan
		o.write('%s,%s,%s,%.6f\n' % (cl, ind, hets[ind]['denom'], pi))
	o.close()

get_diversity(vcf, args.outdir)