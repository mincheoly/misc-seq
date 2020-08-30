"""
	Author: Min Cheol Kim
	
	Simple set of tools for analyzing CHIP-seq data from ENCODE.
"""


from pybedtools import BedTool
import os
import pandas as pd


class Encode():
	
	def __init__(self, gene_locations_path):

		self.gene_locations = pd.read_csv(gene_locations_path, sep='\t')
		self.gene_locations.columns = ['chrom', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'symbol', 'name', 'strand']
		
		
	def get_tss_window(self, window_size):
		
		gene_sites = self.gene_locations.copy()
		gene_sites['tss_site'] = gene_sites['txStart']*(gene_sites['strand'] == '+') + gene_sites['txEnd']*(gene_sites['strand'] == '-') 
		gene_sites['tss_window_start'] = gene_sites['tss_site'] - window_size
		gene_sites['tss_window_end'] = gene_sites['tss_site'] + window_size
		gene_sites['tss_window_start'] = gene_sites['tss_window_start']*(gene_sites['tss_window_start'] > 0) 
		gene_sites = gene_sites[['chrom', 'tss_window_start', 'tss_window_end', 'symbol','name']].sort_values(['chrom', 'tss_window_start'])
		
		return gene_sites
	
	
	def get_encode_peaks(self, encode_link):
		
		fname = encode_link.split('/')[-1]
		dl_status = os.system('wget {}'.format(encode_link))

		assert dl_status == 0
		peaks = BedTool(fname).sort()

		rm_status = os.system('rm {}'.format(fname))
		assert rm_status == 0
		
		return peaks
	
	
	def get_encode_peaks_union(self, encode_links):
		
		beds = [self.get_encode_peaks(link) for link in encode_links]
		
		concat  = None
		for bed in beds:
			
			if concat is None:
				concat = bed
			else:
				concat = concat.cat(bed, postmerge=False)
		
		return concat.sort()
				

		
	def get_peak_genes_encode(self, encode_link, window_size):
		
		tss_sites = self.get_tss_window(window_size)
		tss_sites_bed = BedTool.from_dataframe(tss_sites).sort()
		
		peaks = self.get_encode_peaks(encode_link)
		
		closest = tss_sites_bed.closest(peaks, d=True)
		closest_df = closest.to_dataframe(header=None).iloc[:, [3, -1]].copy()
		closest_df.columns = ['gene', 'distance']
		closest_df = closest_df.sort_values('distance').drop_duplicates(['gene'], keep='first')

		return closest_df

	
	def get_peak_genes_bed(self, peaks, window_size):
		
		tss_sites = self.get_tss_window(window_size)
		tss_sites_bed = BedTool.from_dataframe(tss_sites).sort()
		
		closest = tss_sites_bed.closest(peaks.sort(), d=True)
		closest_df = closest.to_dataframe(header=None).iloc[:, [3, -1]].copy()
		closest_df.columns = ['gene', 'distance']
		closest_df = closest_df.sort_values('distance').drop_duplicates(['gene'], keep='first')

		return closest_df