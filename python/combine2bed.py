#usr/bin/python3


import portion as po
import pandas as pd
from collections import defaultdict


def read_bed(file):
	bed = pd.read_csv(file, sep='\t', names=["Interval_name", "Chr", "Target_Interval_Start", "Target_Interval_End", \
										     "Target_Interval_Fraction", "gene", "Exon", "Exon_Start", "Exon_End"])
	return bed

def level_bed(bed):
    bed_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(dict))))
    for index, row in bed.iterrows():
        interval_name = row.iloc[0]
        chrom = row.iloc[1]
        gene_name = row.iloc[5]
        exon_name = row.iloc[6]
        info = row[1:].tolist()
        bed_dict[chrom][gene_name][exon_name].update({interval_name: info})
    return bed_dict

def merge_all_interval(interval_list):
	interval = interval_list[0]
	for i in range(len(interval_list)):
		interval = interval | interval_list[i]
	return interval
		

def combine2bed(bed_dict1, bed_dict2):
	combined_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(dict))))
	chromosomes1 = set(bed_dict1.keys())
	chromosomes2 = set(bed_dict2.keys())
	
	uncovered_chromosomes1 = chromosomes1 - chromosomes2
	uncovered_chromosomes2 = chromosomes2 - chromosomes1
	
    # 处理未覆盖的染色体
	for chrom in uncovered_chromosomes1:
		combined_dict[chrom] = bed_dict1[chrom]
	for chrom in uncovered_chromosomes2:
		combined_dict[chrom] = bed_dict2[chrom]

    # 合并共同的染色体
	for chrom in (chromosomes1 & chromosomes2):
		genes_bed1 = bed_dict1[chrom]
		genes_bed2 = bed_dict2[chrom]
		genes1 = set(genes_bed1.keys())
		genes2 = set(genes_bed2.keys())
        
		for gene in genes1 - genes2:
			combined_dict[chrom][gene] = genes_bed1[gene]

		for gene in genes2 - genes1:
			combined_dict[chrom][gene] = genes_bed2[gene]
			
		for gene in (genes1 & genes2):
			exons_bed1 = genes_bed1[gene]
			exons_bed2 = genes_bed2[gene]
			exons1 = set(exons_bed1.keys())
			exons2 = set(exons_bed2.keys())
            
			for exon in exons1 - exons2:
				combined_dict[chrom][gene][exon] = exons_bed1[exon]
				
                
			for exon in exons2 - exons1:
				combined_dict[chrom][gene][exon] = exons_bed2[exon]
		
			for exon in (exons1 & exons2):
				interval_list = []
				exon_info = next(iter(exons_bed1[exon].values()))[6:8]
				
				for interval_name, info_list in exons_bed1[exon].items():
					interval_list.append(po.closed(info_list[1], info_list[2]))
				for interval_name, info_list in exons_bed2[exon].items():
					interval_list.append(po.closed(info_list[1], info_list[2]))
				merged_interval = merge_all_interval(interval_list)
				
				for i in range(len(merged_interval)):
					name = "interval" + str(i)
					combined_dict[chrom][gene][exon].update({name: [chrom, merged_interval[i].lower, merged_interval[i].upper, 
                                                                     (merged_interval[i].upper - merged_interval[i].lower + 1) / (exon_info[1] - exon_info[0] + 1),
                                                                     gene, exon, exon_info[0], exon_info[1]]})
	return combined_dict

def main():
	file1 = open(r"C:\Users\geneplus\Desktop\lung_domestic_final.bed")
	file2 = open(r"C:\Users\geneplus\Desktop\lung_oversea_final.bed")
	
	bed1 = read_bed(file1)
	bed2 = read_bed(file2)
	
	file1.close()
	file2.close()
	
	bed_dict1 = level_bed(bed1)
	bed_dict2 = level_bed(bed2)
	
	combined_dict = combine2bed(bed_dict1, bed_dict2)
	with open (r"D:\BaiduNetdiskDownload\test.tsv", 'w') as a:
		for chrom in combined_dict:
			for gene in combined_dict[chrom]:
				for exon in combined_dict[chrom][gene]:
					for interval_name, interval_info in combined_dict[chrom][gene][exon].items():
						line = ""
						for i in range(len(interval_info)):
							line += str(interval_info[i])
							line += '\t'
						line += '\n'
						a.writelines(line)
	a.close()


if __name__ == '__main__':
	main()

					
				
						
		
		