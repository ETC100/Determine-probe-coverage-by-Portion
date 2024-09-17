import re
import os
import itertools
import portion
import numpy as np
import math
import copy
from collections import defaultdict
import time
import multiprocessing
from multiprocessing import Pool
import pandas as pd


def portion_open_to_close(pIntervals):
	Transfered_pIntervals = portion.empty()
	for interval in pIntervals:
		if interval.left == portion.OPEN:
			interval = interval.replace(portion.CLOSED,interval.lower+1,interval.upper,interval.right)
		if not interval:
			continue
		if interval.right == portion.OPEN:
			interval = interval.replace(interval.left,interval.lower,interval.upper-1,portion.CLOSED)
		if not interval:
			continue
		Transfered_pIntervals = Transfered_pIntervals | interval
	return Transfered_pIntervals


def portion_interval_merge(Interval_List):
	Merged_pIntervals = portion.empty()
	for list_interval in Interval_List:
		Merged_pIntervals = Merged_pIntervals | portion.closed(list_interval[0],list_interval[1])
	Merged_Interval_List = []
	for interval in Merged_pIntervals:
		Merged_Interval_List.append([interval.lower,interval.upper])
	return Merged_Interval_List


def portion_interval_intersect(Interval_List1,Interval_List2):
	Intersected_pIntervals = portion.empty()
	for list_interval1 in Interval_List1:
		for list_interval2 in Interval_List2:
			if portion.closed(list_interval1[0], list_interval1[1]) & portion.closed(list_interval2[0], list_interval2[1]):
				Intersected_pIntervals = Intersected_pIntervals | ( portion.closed(list_interval1[0], list_interval1[1]) & portion.closed(list_interval2[0], list_interval2[1]) )
	Intersected_Interval_List = []
	for interval in Intersected_pIntervals:
		Intersected_Interval_List.append([interval.lower,interval.upper])
	return Intersected_Interval_List


def portion_interval_differ(Interval_List1,Interval_List2):
	Differ_pIntervals1 = portion.empty()
	for list_interval1 in Interval_List1:
		interval1 = portion.closed(list_interval1[0], list_interval1[1])
		for list_interval2 in Interval_List2:
			interval2 = portion.closed(list_interval2[0], list_interval2[1])
			interval1 = interval1 - interval2
		if interval1:
	
			Differ_pIntervals1 = Differ_pIntervals1 | interval1
	Differ_pIntervals1 = portion_open_to_close(Differ_pIntervals1)
	Differed_Interval_List1 = []	
	for interval in Differ_pIntervals1:	
		Differed_Interval_List1.append([interval.lower,interval.upper])
	Differ_pIntervals2 = portion.empty()
	for list_interval2 in Interval_List2:
		interval2 = portion.closed(list_interval2[0], list_interval2[1])
		for list_interval1 in Interval_List1:
			interval1 = portion.closed(list_interval1[0], list_interval1[1])
			interval2 = interval2 - interval1
		if interval2:
			Differ_pIntervals2 = Differ_pIntervals2 | interval2
	Differ_pIntervals2 = portion_open_to_close(Differ_pIntervals2)
	Differed_Interval_List2 = []
	for interval in Differ_pIntervals2:
		Differed_Interval_List2.append([interval.lower,interval.upper])
	return Differed_Interval_List1,Differed_Interval_List2


def portion_interval_extend(Interval_List1,Interval_List2):
	pIntervals1 = portion.empty()
	for list_interval1 in Interval_List1:
		pIntervals1 = pIntervals1 | portion.closed(list_interval1[0], list_interval1[1])
	pIntervals2 = portion.empty()
	for list_interval2 in Interval_List2:
		pIntervals2 = pIntervals2 | portion.closed(list_interval2[0], list_interval2[1])
	if pIntervals1.upper == pIntervals2.upper:
		pIntervals1 = pIntervals1.replace(portion.CLOSED,pIntervals1.lower,pIntervals1.upper+3,portion.CLOSED)
	if pIntervals1.lower == pIntervals2.lower:
		pIntervals1 = pIntervals1.replace(portion.CLOSED,pIntervals1.lower-3,pIntervals1.upper,portion.CLOSED)
	Extended_Interval_List1 = []
	for interval in pIntervals1:
		Extended_Interval_List1.append([interval.lower,interval.upper])
	return Extended_Interval_List1


def have_neighbour_interval(Interval_List1,Interval_List2):
	Have_Neighbour = False
	start1, end1 = Interval_List1[0]
	for start2, end2 in Interval_List2:
		if not ( start1 - end2 -1 > 80 or start2 - end1 -1 > 80):
			Have_Neighbour = True
	return Have_Neighbour
			

def calculate_neighbour_intervals_length(Interval_List):
	Length = 0
	original_pIntervals = portion.empty()
	for list_interval in Interval_List:
		original_pIntervals = original_pIntervals | portion.closed(list_interval[0],list_interval[1])
	range_pIntervals = portion.closed(original_pIntervals.lower,original_pIntervals.upper)
	insert_pIntervals = range_pIntervals - original_pIntervals
	insert_pIntervals = portion_open_to_close(insert_pIntervals)
	for interval in insert_pIntervals:
		if interval.upper - interval.lower + 1 <= 80:
			original_pIntervals = original_pIntervals | interval
	for interval in original_pIntervals:
		Length = Length + ( interval.upper - interval.lower + 1 )
	return Length


def merge_neighbour_intervals(Interval_List):
	original_pIntervals = portion.empty()
	for list_interval in Interval_List:
		original_pIntervals = original_pIntervals | portion.closed(list_interval[0],list_interval[1])
	range_pIntervals = portion.closed(original_pIntervals.lower,original_pIntervals.upper)
	insert_pIntervals = range_pIntervals - original_pIntervals
	insert_pIntervals = portion_open_to_close(insert_pIntervals)
	for interval in insert_pIntervals:
		if interval.upper - interval.lower + 1 <= 80:
			original_pIntervals = original_pIntervals | portion.closed(interval.lower-1,interval.upper+1)
	Merged_Interval_List = []
	Length = 0
	for interval in original_pIntervals:
		Length = Length + ( interval.upper - interval.lower + 1 )
		Merged_Interval_List.append([interval.lower,interval.upper])

	return Merged_Interval_List,Length

def get_probe_target_include_and_real_coverage():
	pattern = re.compile("EX\\d+E?")
	All_Gene_Exon = {}
	ncbi_anno_file_in = open("../../1021/config/ncbi_anno_rel104_db_b37_TXS_IVS.bed","r")
	for line in ncbi_anno_file_in:
		line = line.strip()
		if not line or line.startswith("#"):
			continue
		gene_annotation = line.split("\t")[3]
		if pattern.search(gene_annotation.split("|")[3]) and gene_annotation.split("|")[6] == "CDS":
			gene_info = gene_annotation.split("|")[1]+":"+gene_annotation.split("|")[4]
			exon_info = gene_annotation.split("|")[1]+":"+gene_annotation.split("|")[4]+":"+gene_annotation.split("|")[3]
			exon_start = int(line.split("\t")[1]) + 1
			exon_end = int(line.split("\t")[2])
			exon_chr = line.split("\t")[0]
			if gene_info not in All_Gene_Exon:
				All_Gene_Exon[gene_info] = {}
			All_Gene_Exon[gene_info][exon_info] = [exon_start,exon_end,exon_chr]
	ncbi_anno_file_in.close()

	Double_Probe_Include_Gene = {}
	Probe_Gene_RawCoveranges_Dict = {}

	for probe in ["cd3","cd4a"]:
		probe_include_file_in = open("../../1021/config/ncbi_anno_rel104_db_b37_%s_chip_cover.bed"%(probe),"r")
		for line in probe_include_file_in:
			line = line.strip()
			if not line or line.startswith("#"):
				continue
			include_gene_annotation = line.split("\t")[3]
			if pattern.search(include_gene_annotation.split("|")[3]) and include_gene_annotation.split("|")[6] == "CDS":
				include_gene_info = include_gene_annotation.split("|")[1]+":"+include_gene_annotation.split("|")[4]
				if not include_gene_info in Probe_Gene_RawCoveranges_Dict:
					Probe_Gene_RawCoveranges_Dict[include_gene_info] = {}
				include_exon_info = include_gene_annotation.split("|")[1]+":"+include_gene_annotation.split("|")[4]+":"+include_gene_annotation.split("|")[3]
				if not include_exon_info in Probe_Gene_RawCoveranges_Dict[include_gene_info]:
					Probe_Gene_RawCoveranges_Dict[include_gene_info][include_exon_info] = {}
				include_exon_ann_start = int(line.split("\t")[1]) + 1
				include_exon_ann_end = int(line.split("\t")[2])
				probe_rowcover_start = int(line.split("\t")[5]) + 1
				probe_rowcover_end = int(line.split("\t")[6])
				temp_overlap_list = [include_exon_ann_start,include_exon_ann_end,probe_rowcover_start,probe_rowcover_end]
				temp_overlap_list.sort()
				Probe_Gene_RawCoveranges_Dict[include_gene_info][include_exon_info].setdefault(probe,[]).append([temp_overlap_list[1],temp_overlap_list[2]])
				temp_raw_coveranges = Probe_Gene_RawCoveranges_Dict[include_gene_info][include_exon_info][probe]
				Probe_Gene_RawCoveranges_Dict[include_gene_info][include_exon_info][probe] = portion_interval_merge(temp_raw_coveranges)
				if include_gene_info not in Double_Probe_Include_Gene:
					Double_Probe_Include_Gene[include_gene_info] = {}
				Double_Probe_Include_Gene[include_gene_info][include_exon_info] = [[include_exon_ann_start,include_exon_ann_end]]
		probe_include_file_in.close()

	Probe_Gene_Coveranges_Dict = {}
	for probe in ["cd3","cd4a"]:
		probe_flank50_cover_file_in = open("../../1021/config/ncbi_anno_rel104_db_b37_%s_chip_flank_cover.bed"%(probe),"r")
		for line in probe_flank50_cover_file_in:
			line = line.strip()
			if not line or line.startswith("#"):
				continue
			cover_gene_annotation = line.split("\t")[3]
			if pattern.search(cover_gene_annotation.split("|")[3]) and cover_gene_annotation.split("|")[6] == "CDS":
				cover_gene_info = cover_gene_annotation.split("|")[1]+":"+cover_gene_annotation.split("|")[4]
				cover_exon_info = cover_gene_annotation.split("|")[1]+":"+cover_gene_annotation.split("|")[4]+":"+cover_gene_annotation.split("|")[3]
				if cover_gene_info in Double_Probe_Include_Gene and cover_exon_info in Double_Probe_Include_Gene[cover_gene_info]:
					cover_exon_start = int(line.split("\t")[1]) + 1
					cover_exon_end = int(line.split("\t")[2])
					if cover_gene_info not in Probe_Gene_Coveranges_Dict:
						Probe_Gene_Coveranges_Dict[cover_gene_info] = {}
					if cover_exon_info not in Probe_Gene_Coveranges_Dict[cover_gene_info]:
						Probe_Gene_Coveranges_Dict[cover_gene_info][cover_exon_info] = {}
					#提取探针真实覆盖区间
					temp_overlap_list = [int(line.split("\t")[1]) + 1,int(line.split("\t")[2]),int(line.split("\t")[5])+1,int(line.split("\t")[6])]
					temp_overlap_list.sort()
					#将同一外显子内所有的探针真实覆盖区间，进行合并
					Probe_Gene_Coveranges_Dict[cover_gene_info][cover_exon_info].setdefault(probe,[]).append([temp_overlap_list[1],temp_overlap_list[2]])
					temp_coveranges = Probe_Gene_Coveranges_Dict[cover_gene_info][cover_exon_info][probe]
					###Probe_Gene_Coveranges_Dict赋值行
					Probe_Gene_Coveranges_Dict[cover_gene_info][cover_exon_info][probe] = portion_interval_merge(temp_coveranges)
		probe_flank50_cover_file_in.close()

	return Double_Probe_Include_Gene,Probe_Gene_RawCoveranges_Dict,Probe_Gene_Coveranges_Dict,All_Gene_Exon


def probe_different_coverange():
	Double_Probe_Include_Gene,Probe_Gene_RawCoveranges_Dict,Probe_Gene_Coveranges_Dict,All_Gene_Exon = get_probe_target_include_and_real_coverage()
	Probe_Gene_Differange_Dict = {}
	Probe_Gene_Differange_Extend_Dict = {}
	for gene_info in Probe_Gene_Coveranges_Dict.keys():
		Probe_Gene_Differange_Dict[gene_info] = {}
		Probe_Gene_Differange_Extend_Dict[gene_info] = {}
		for exon_info in Probe_Gene_Coveranges_Dict[gene_info].keys():
			Probe_Gene_Differange_Dict[gene_info][exon_info] = {}
			Probe_Gene_Differange_Extend_Dict[gene_info][exon_info] = {}
			if "cd3" in Probe_Gene_Coveranges_Dict[gene_info][exon_info] and not "cd4a" in Probe_Gene_Coveranges_Dict[gene_info][exon_info]:
				Probe_Gene_Differange_Dict[gene_info][exon_info]["cd3_only"] = Probe_Gene_Coveranges_Dict[gene_info][exon_info]["cd3"]
				Probe_Gene_Differange_Extend_Dict[gene_info][exon_info]["cd3_only"] = portion_interval_extend(  Probe_Gene_Differange_Dict[gene_info][exon_info]["cd3_only"],Double_Probe_Include_Gene[gene_info][exon_info]  )
			elif "cd4a" in Probe_Gene_Coveranges_Dict[gene_info][exon_info] and "cd3" not in Probe_Gene_Coveranges_Dict[gene_info][exon_info]:
				Probe_Gene_Differange_Dict[gene_info][exon_info]["cd4a_only"] = Probe_Gene_Coveranges_Dict[gene_info][exon_info]["cd4a"]
				Probe_Gene_Differange_Extend_Dict[gene_info][exon_info]["cd4a_only"] = portion_interval_extend( Probe_Gene_Differange_Dict[gene_info][exon_info]["cd4a_only"],Double_Probe_Include_Gene[gene_info][exon_info] )
			else:
				Probe_Gene_Differange_Dict[gene_info][exon_info]["cd3-cd4a"] = portion_interval_intersect( Probe_Gene_Coveranges_Dict[gene_info][exon_info]["cd3"], Probe_Gene_Coveranges_Dict[gene_info][exon_info]["cd4a"] )
				Probe_Gene_Differange_Extend_Dict[gene_info][exon_info]["cd3-cd4a"] = portion_interval_extend( Probe_Gene_Differange_Dict[gene_info][exon_info]["cd3-cd4a"], Double_Probe_Include_Gene[gene_info][exon_info] )
				cd3_only,cd4a_only = portion_interval_differ( Probe_Gene_Coveranges_Dict[gene_info][exon_info]["cd3"], Probe_Gene_Coveranges_Dict[gene_info][exon_info]["cd4a"] )
				if cd3_only:
					Probe_Gene_Differange_Dict[gene_info][exon_info]["cd3_only"] = cd3_only
					Probe_Gene_Differange_Extend_Dict[gene_info][exon_info]["cd3_only"] = portion_interval_extend( cd3_only , Double_Probe_Include_Gene[gene_info][exon_info] )
				if cd4a_only:
					Probe_Gene_Differange_Dict[gene_info][exon_info]["cd4a_only"] = cd4a_only
					Probe_Gene_Differange_Extend_Dict[gene_info][exon_info]["cd4a_only"] = portion_interval_extend( cd4a_only , Double_Probe_Include_Gene[gene_info][exon_info] )

	return Probe_Gene_Differange_Dict,Probe_Gene_Differange_Extend_Dict


def read_Panel_data_file():
	Panel_Dataset = {}
	Panel_data_file_in = open("../../1021/output/retained_mutation1024.tsv","r")
	for line in Panel_data_file_in:
		line = line.strip()
		if not line or line.startswith("#") or line.startswith("patientDNA_ID"):
			continue
		gene_info = ":".join(line.split("\t")[1:3])
		exon_info = ":".join(line.split("\t")[1:4])
		if gene_info not in Panel_Dataset:
			Panel_Dataset[gene_info] = {}
		if exon_info not in Panel_Dataset[gene_info]:
			Panel_Dataset[gene_info][exon_info] = []
		Panel_Dataset[gene_info][exon_info].append(":".join(line.split("\t")[:4])+":"+":".join(line.split("\t")[5:12]) + ":" + line.split("\t")[-2])
		
	Panel_data_file_in.close()
	
	return Panel_Dataset


def panel_mutation_remapping(OUTPUT_DIR):
	Pan_Probe_Differange_CoverMut = {}
	if os.path.exists("../output/BINGXING/%s/Pan_Probe_Differange_CoverMut.tsv"%(OUTPUT_DIR)):
		Pan_Probe_Differange_CoverMut_file_in = open("../output/BINGXING/%s/Pan_Probe_Differange_CoverMut.tsv"%(OUTPUT_DIR),"r")
		for line in Pan_Probe_Differange_CoverMut_file_in:
			line = line.strip()
			if not line or line.startswith("#"):
				continue
			gene_info = line.split("\t")[0]
			exon_info = line.split("\t")[1]
			probe = line.split("\t")[2]
			mutations = line.split("\t")[3:]
			if gene_info not in Pan_Probe_Differange_CoverMut:
				Pan_Probe_Differange_CoverMut[gene_info] = {}
			if exon_info not in Pan_Probe_Differange_CoverMut[gene_info]:
				Pan_Probe_Differange_CoverMut[gene_info][exon_info] = {}
			Pan_Probe_Differange_CoverMut[gene_info][exon_info][probe] = mutations
		Pan_Probe_Differange_CoverMut_file_in.close()
		
		return Pan_Probe_Differange_CoverMut,Probe_Gene_Differange_Extend_Dict
	
	Probe_Gene_Differange_Dict,Probe_Gene_Differange_Extend_Dict = probe_different_coverange()
	Panel_Dataset = read_Panel_data_file()
	remapping_failed_mutation = open("../output/BINGXING/%s/panel_remapping_failed_mutation.tsv"%(OUTPUT_DIR),"w")
	for gene_info in Probe_Gene_Differange_Extend_Dict.keys():
		Pan_Probe_Differange_CoverMut[gene_info] = {}
		for exon_info in Probe_Gene_Differange_Extend_Dict[gene_info].keys():
			Pan_Probe_Differange_CoverMut[gene_info][exon_info] = {}
			for probe in Probe_Gene_Differange_Extend_Dict[gene_info][exon_info].keys():
				Pan_Probe_Differange_CoverMut[gene_info][exon_info][probe] = []


	for gene_info in Panel_Dataset.keys():
		if gene_info in Probe_Gene_Differange_Extend_Dict:
			for exon_info in Panel_Dataset[gene_info]:
				if exon_info in Probe_Gene_Differange_Extend_Dict[gene_info]:
					for mutation_info in Panel_Dataset[gene_info][exon_info]:
						tag = True
						for probe,differanges in Probe_Gene_Differange_Extend_Dict[gene_info][exon_info].items():
							for probe_range in differanges:
								if not ( int(mutation_info.split(":")[6]) <= probe_range[0]-1 or int(mutation_info.split(":")[5]) >= probe_range[1] ):
									Pan_Probe_Differange_CoverMut[gene_info][exon_info].setdefault(probe,[]).append(mutation_info)
									tag = False
						if tag:
							remapping_failed_mutation.write(mutation_info+"\tmapping"+"\n")
				else:
					for mutation_info in Panel_Dataset[gene_info][exon_info]:
						remapping_failed_mutation.write(mutation_info+"\texon_info"+"\n")
		else:
			for exon_info in Panel_Dataset[gene_info]:
				for mutation_info in Panel_Dataset[gene_info][exon_info]:
					remapping_failed_mutation.write(mutation_info+"\texon_info"+"\n")
	
	Pan_Probe_Differange_CoverMut_file_out = open("../output/BINGXING/%s/Pan_Probe_Differange_CoverMut.tsv"%(OUTPUT_DIR),"w")
	for gene_info in Pan_Probe_Differange_CoverMut.keys():
		for exon_info in Pan_Probe_Differange_CoverMut[gene_info]:
			for probe in Pan_Probe_Differange_CoverMut[gene_info][exon_info]:
				Pan_Probe_Differange_CoverMut_file_out.write(gene_info+"\t"+exon_info+"\t"+probe+"\t"+"\t".join(Pan_Probe_Differange_CoverMut[gene_info][exon_info][probe]) + "\n")	
	Pan_Probe_Differange_CoverMut_file_out.close()

	return Pan_Probe_Differange_CoverMut,Probe_Gene_Differange_Extend_Dict


def read_ncbi_anno_file():
	pattern = re.compile("EX\\d+E?")
	wes_transcript_list = []
	tr_nc_file_in = open("/root/probe_update/1021/tr_nc.lst","r")
	for line in tr_nc_file_in:
		line = line.strip()
		if not line or line.startswith("#"):
			continue
		wes_transcript_list.append(line)
	tr_nc_file_in.close()
	WES_All_Gene_Exon = {}
	ncbi_anno_file_in = open("../../1021/config/ncbi_anno_rel104_db_b37_TXS_IVS.bed","r")
	for line in ncbi_anno_file_in:
		line = line.strip()
		if not line or line.startswith("#"):
			continue
		gene_annotation = line.split("\t")[3]
		if gene_annotation.split("|")[4] in wes_transcript_list and pattern.search(gene_annotation.split("|")[3]) and gene_annotation.split("|")[6] == "CDS":
			gene_info = gene_annotation.split("|")[1]+":"+gene_annotation.split("|")[4]
			exon_info = gene_annotation.split("|")[1]+":"+gene_annotation.split("|")[4]+":"+gene_annotation.split("|")[3]
			exon_start = int(line.split("\t")[1]) + 1
			exon_end = int(line.split("\t")[2])
			exon_chr = line.split("\t")[0]
			if gene_info not in WES_All_Gene_Exon:
				WES_All_Gene_Exon[gene_info] = {}
			WES_All_Gene_Exon[gene_info][exon_info] = [exon_start,exon_end,exon_chr]
	ncbi_anno_file_in.close()

	return WES_All_Gene_Exon


def read_WES_file(TARGET_TUMOR,TARGET_TUMOR_WES_INPUT_LIST):
	pattern = re.compile("EX\\d+E?")
	WES_Dataset_By_Gene = {}
	WES_Dataset_By_Patient = {}
	for INPUT in TARGET_TUMOR_WES_INPUT_LIST:
		WES_file_in = open(INPUT,"r")
		for line in WES_file_in:
			line = line.strip()
			if not line or line.startswith("#") or line.startswith("Sample_ID") :
				continue
			gene_info = line.split("\t")[0]
			exon_info = line.split("\t")[1]
			patient_id = line.split("\t")[2].split(":")[0]
			tumor_type = line.split("\t")[2].split(":")[13]
			if gene_info not in WES_Dataset_By_Gene:		
				WES_Dataset_By_Gene[gene_info] = {}
			if exon_info not in WES_Dataset_By_Gene[gene_info]:
				WES_Dataset_By_Gene[gene_info][exon_info] = []
			WES_Dataset_By_Gene[gene_info][exon_info].append( ":".join(line.split("\t")[2].split(":")[:14]) )
		
			if tumor_type not in WES_Dataset_By_Patient:
				WES_Dataset_By_Patient[tumor_type] = {}
			if patient_id not in WES_Dataset_By_Patient[tumor_type]:
				WES_Dataset_By_Patient[tumor_type][patient_id] = []
			WES_Dataset_By_Patient[tumor_type][patient_id].append( ":".join(line.split("\t")[2].split(":")[:14]) )
		WES_file_in.close()

	return WES_Dataset_By_Gene,WES_Dataset_By_Patient


def read_determined_gene_exon(WES_All_Gene_Exon,OLD_TUMOR,OUTPUT_DIR):
	pattern = re.compile("EX\\d+E?")
	WES_Determined_Gene_Exon = {}
	#read_old
	if os.path.exists("../output/BINGXING/%s/%s_WES_Determined_Gene_Exon.tsv"%(OUTPUT_DIR,OLD_TUMOR)):
		WES_Determined_Gene_Exon_file_in = open("../output/BINGXING/%s/%s_WES_Determined_Gene_Exon.tsv"%(OUTPUT_DIR,OLD_TUMOR))
		for line in WES_Determined_Gene_Exon_file_in:
			line = line.strip()
			if not line or line.startswith("#"):
				continue
			gene_info = line.split("\t")[0]
			exon_info = line.split("\t")[1]
			interval_start = int(line.split("\t")[5])
			interval_end = int(line.split("\t")[6])
			if gene_info not in WES_Determined_Gene_Exon:
				WES_Determined_Gene_Exon[gene_info] = {}
			WES_Determined_Gene_Exon[gene_info].setdefault(exon_info,[]).append([interval_start,interval_end])
		WES_Determined_Gene_Exon_file_in.close()
		return WES_Determined_Gene_Exon

	#de_novo
	determined_interval_file_in = open("../input/lung_determined_Y12.bed","r")
	for line in determined_interval_file_in:
		line = line.strip()
		if not line or line.startswith("#"):
			continue
		cover_gene_annotation = line.split("\t")[3]
		if pattern.search(cover_gene_annotation.split("|")[3]) and cover_gene_annotation.split("|")[6] == "CDS":
			cover_gene_info = cover_gene_annotation.split("|")[1]+":"+cover_gene_annotation.split("|")[4]
			cover_exon_info = cover_gene_annotation.split("|")[1]+":"+cover_gene_annotation.split("|")[4]+":"+cover_gene_annotation.split("|")[3]
			temp_overlap_list = [int(line.split("\t")[1]) + 1,int(line.split("\t")[2]),int(line.split("\t")[5])+1,int(line.split("\t")[6])]
			temp_overlap_list.sort()
			if not cover_gene_info in WES_Determined_Gene_Exon:
				WES_Determined_Gene_Exon[cover_gene_info] = {}
			WES_Determined_Gene_Exon[cover_gene_info].setdefault(cover_exon_info,[]).append([temp_overlap_list[1],temp_overlap_list[2]])
			temp_coveranges = WES_Determined_Gene_Exon[cover_gene_info][cover_exon_info]
			WES_Determined_Gene_Exon[cover_gene_info][cover_exon_info] = portion_interval_merge(temp_coveranges)
	determined_interval_file_in.close()
	WES_Determined_Gene_Exon_file_out = open("../output/BINGXING/%s/Part1_2_WES_Determined_Gene_Exon.tsv"%(OUTPUT_DIR),"w")
	for gene_info in WES_Determined_Gene_Exon.keys():
		for exon_info in WES_Determined_Gene_Exon[gene_info]:
			for interval in WES_Determined_Gene_Exon[gene_info][exon_info]:
				WES_Determined_Gene_Exon_file_out.write(gene_info+"\t"+exon_info+"\t"+str(WES_All_Gene_Exon[gene_info][exon_info][2])+"\t"+str(WES_All_Gene_Exon[gene_info][exon_info][0])+"\t"+str(WES_All_Gene_Exon[gene_info][exon_info][1])+"\t"+str(interval[0])+"\t"+str(interval[1])+"\n")
	WES_Determined_Gene_Exon_file_out.close()
	return WES_Determined_Gene_Exon


def candidate_gene_exon_generator(WES_Determined_Gene_Exon,WES_All_Gene_Exon,TARGET_TUMOR,OUTPUT_DIR):
	WES_Candidate_Gene_Exon = {}
	#read_old
	if os.path.exists("../output/BINGXING/%s/%s_WES_Candidate_Gene_Exon.tsv"%(OUTPUT_DIR,TARGET_TUMOR)):
		WES_Candidate_Gene_Exon_file_in = open("../output/BINGXING/%s/%s_WES_Candidate_Gene_Exon.tsv"%(OUTPUT_DIR,TARGET_TUMOR),"r")
		for line in WES_Candidate_Gene_Exon_file_in:
			line = line.strip()
			if not line or line.startswith("#"):
				continue
			gene_info = line.split("\t")[0]
			exon_info = line.split("\t")[1]
			interval_start = int(line.split("\t")[5])
			interval_end = int(line.split("\t")[6])
			if gene_info not in WES_Candidate_Gene_Exon:
				WES_Candidate_Gene_Exon[gene_info] = {}
			WES_Candidate_Gene_Exon[gene_info].setdefault(exon_info,[]).append([interval_start,interval_end])
		WES_Candidate_Gene_Exon_file_in.close()
		
		return WES_Candidate_Gene_Exon

	for gene_info in WES_All_Gene_Exon.keys():
		WES_Candidate_Gene_Exon[gene_info] = {}
		save_gene_info = False
		for exon_info in WES_All_Gene_Exon[gene_info]:
			if gene_info in WES_Determined_Gene_Exon and exon_info in WES_Determined_Gene_Exon[gene_info]:
				differed_interval_list1, differed_interval_list2 = portion_interval_differ( [ [WES_All_Gene_Exon[gene_info][exon_info][0],WES_All_Gene_Exon[gene_info][exon_info][1]] ] , WES_Determined_Gene_Exon[gene_info][exon_info] )
				if [ x for x in differed_interval_list1 if x ]:
					WES_Candidate_Gene_Exon[gene_info][exon_info] = portion_interval_extend( differed_interval_list1 , [ [WES_All_Gene_Exon[gene_info][exon_info][0],WES_All_Gene_Exon[gene_info][exon_info][1]] ] )
					save_gene_info = True
			else:
				WES_Candidate_Gene_Exon[gene_info][exon_info] = [ [WES_All_Gene_Exon[gene_info][exon_info][0]-3,WES_All_Gene_Exon[gene_info][exon_info][1]+3] ]
				save_gene_info = True
		if not save_gene_info:
			del WES_Candidate_Gene_Exon[gene_info]
	
	WES_Candidate_Gene_Exon_file_out = open("../output/BINGXING/%s/%s_WES_Candidate_Gene_Exon.tsv"%(OUTPUT_DIR,TARGET_TUMOR),"w")
	for gene_info in WES_Candidate_Gene_Exon.keys():
		for exon_info in WES_Candidate_Gene_Exon[gene_info]:
			for interval in WES_Candidate_Gene_Exon[gene_info][exon_info]:
				WES_Candidate_Gene_Exon_file_out.write(gene_info+"\t"+exon_info+"\t"+str(WES_All_Gene_Exon[gene_info][exon_info][2])+"\t"+str(WES_All_Gene_Exon[gene_info][exon_info][0])+"\t"+str(WES_All_Gene_Exon[gene_info][exon_info][1])+"\t"+str(interval[0])+"\t"+str(interval[1])+"\n")
	WES_Candidate_Gene_Exon_file_out.close()
	return WES_Candidate_Gene_Exon


def wes_mutation_remapping(WES_Dataset_By_Gene,WES_Tag_Gene_Exon,Tag,TARGET_TUMOR,OUTPUT_DIR):
	Tag_Gene_Exon_CoverMut = {}
	#read_old
	if os.path.exists("../output/BINGXING/%s/%s_%s_remapping_successed_mutation.tsv"%(OUTPUT_DIR,TARGET_TUMOR,Tag)):
		Tag_Gene_Exon_CoverMut_file_in = open("../output/BINGXING/%s/%s_%s_remapping_successed_mutation.tsv"%(OUTPUT_DIR,TARGET_TUMOR,Tag),"r")
		for line in Tag_Gene_Exon_CoverMut_file_in:
			line = line.strip()
			if not line or line.startswith("#"):
				continue
			gene_info = line.split("\t")[0]
			exon_info = line.split("\t")[1]
			exon_covermut = line.split("\t")[2:]
			if gene_info not in Tag_Gene_Exon_CoverMut:
				Tag_Gene_Exon_CoverMut[gene_info] = {}
			Tag_Gene_Exon_CoverMut[gene_info][exon_info] = exon_covermut
		Tag_Gene_Exon_CoverMut_file_in.close()

		return Tag_Gene_Exon_CoverMut

	remapping_failed_mutation = open("../output/BINGXING/%s/%s_%s_remapping_failed_mutation.tsv"%(OUTPUT_DIR,TARGET_TUMOR,Tag),"w")
	Tag_Gene_Exon_CoverMut = {}
	for gene_info in WES_Tag_Gene_Exon.keys():
		Tag_Gene_Exon_CoverMut[gene_info] = {}
		for exon_info in WES_Tag_Gene_Exon[gene_info].keys():
			Tag_Gene_Exon_CoverMut[gene_info][exon_info] = []
	for gene_info in WES_Dataset_By_Gene.keys():
		if gene_info in WES_Tag_Gene_Exon:
			for exon_info in WES_Dataset_By_Gene[gene_info]:
				if exon_info in WES_Tag_Gene_Exon[gene_info]:
					for mutation_info in WES_Dataset_By_Gene[gene_info][exon_info]:
						for start,end in WES_Tag_Gene_Exon[gene_info][exon_info]:
							if not ( int(mutation_info.split(":")[6]) < start or int(mutation_info.split(":")[5]) > end ):
								Tag_Gene_Exon_CoverMut[gene_info][exon_info].append(mutation_info)
							else:
								remapping_failed_mutation.write(mutation_info+"\tmapping"+"\n")
				else:
					for mutation_info in WES_Dataset_By_Gene[gene_info][exon_info]:
						remapping_failed_mutation.write(mutation_info+"\texon_info"+"\n")
		else:
			for exon_info in WES_Dataset_By_Gene[gene_info]:
				for mutation_info in WES_Dataset_By_Gene[gene_info][exon_info]:
					remapping_failed_mutation.write(mutation_info+"\texon_info"+"\n")
	remapping_failed_mutation.close()
	
	remapping_successed_mutation = open("/mnt/GenePlus002/genecloud/Org_terminal/org_52/terminal/wangym_16601250368/probe_update/internal_wes/output/BINGXING/%s/%s_%s_remapping_successed_mutation.tsv"%(OUTPUT_DIR,TARGET_TUMOR,Tag),"w")
	for gene_info in Tag_Gene_Exon_CoverMut.keys():
		for exon_info in Tag_Gene_Exon_CoverMut[gene_info].keys():
			remapping_successed_mutation.write(gene_info+"\t"+exon_info+"\t"+"\t".join(Tag_Gene_Exon_CoverMut[gene_info][exon_info])+"\n")
	remapping_successed_mutation.close()

	return	Tag_Gene_Exon_CoverMut


def interval_generator_bystep(Candidate_Gene_Exon_CoverMut,WES_Candidate_Gene_Exon,TARGET_TUMOR,OUTPUT_DIR):
	Candidate_Interval_CoverMut = {}
	Ori_Candidate_Interval_CoverMut = {}
	#read_old
	if os.path.exists("../output/BINGXING/%s/%s_Candidate_Interval_CoverMut.tsv"%(OUTPUT_DIR,TARGET_TUMOR)):
		black_list = []
		
		black_list_file_in = open("../input/intervals_blacklist/lung_LINCHUANG_blacklist.bed","r")
		for line in black_list_file_in:
			line = line.strip()
			if not line or line.startswith("#"):
				continue
			black_list.append( ( int(line.split("\t")[1]),int(line.split("\t")[2]) ) )
		black_list_file_in.close()

		black_list_file_in = open("../input/intervals_blacklist/lung_paper_domestic_blacklist.bed","r")
		for line in black_list_file_in:
			line = line.strip()
			if not line or line.startswith("#"):
				continue
			black_list.append( ( int(line.split("\t")[1]),int(line.split("\t")[2]) ) )
		black_list_file_in.close()

		black_list_file_in = open("../input/intervals_blacklist/lung_LINCHUANG_blacklist.bed","r")
		for line in black_list_file_in:
			line = line.strip()
			if not line or line.startswith("#"):
				continue
			black_list.append( ( int(line.split("\t")[1]),int(line.split("\t")[2]) ) )
		black_list_file_in.close()

		black_list_file_in = open("../input/intervals_blacklist/lung_paper_domestic_blacklist.bed","r")
		for line in black_list_file_in:
			line = line.strip()
			if not line or line.startswith("#"):
				continue
			black_list.append( ( int(line.split("\t")[1]),int(line.split("\t")[2]) ) )
		black_list_file_in.close()

		selected_interval_list = []
		Selected_Interval_File_in = open("../output/BINGXING/%s/%s_WES_interval_improve_coverage_selected_interval.tsv"%(OUTPUT_DIR,TARGET_TUMOR),"r")
		for line in Selected_Interval_File_in:
			line = line.strip()
			if not line or line.startswith("#") or line.startswith("Interval_ID"):
				continue
			interval_id = line.split("\t")[0]
			selected_interval_list.append(interval_id)
		Selected_Interval_File_in.close()

		Candidate_Interval_CoverMut_file_in = open("../output/BINGXING/%s/%s_Candidate_Interval_CoverMut.tsv"%(OUTPUT_DIR,TARGET_TUMOR),"r")
		for line in Candidate_Interval_CoverMut_file_in:
			line = line.strip()
			if not line or line.startswith("#"):
				continue
			gene_info = line.split("\t")[0]
			exon_info = line.split("\t")[1]	
			interval_info = line.split("\t")[2]
			interval_covermut = line.split("\t")[3:]
			if (int(interval_info.split(":")[5]), int(interval_info.split(":")[6])) not in black_list and interval_info not in selected_interval_list:
				if gene_info not in Candidate_Interval_CoverMut:
					Candidate_Interval_CoverMut[gene_info] = {}
				if exon_info not in Candidate_Interval_CoverMut[gene_info]:
					Candidate_Interval_CoverMut[gene_info][exon_info] = {}
				Candidate_Interval_CoverMut[gene_info][exon_info][interval_info] = interval_covermut
			if gene_info not in Ori_Candidate_Interval_CoverMut:
				Ori_Candidate_Interval_CoverMut[gene_info] = {}
			if exon_info not in Ori_Candidate_Interval_CoverMut[gene_info]:
				Ori_Candidate_Interval_CoverMut[gene_info][exon_info] = {}
			Ori_Candidate_Interval_CoverMut[gene_info][exon_info][interval_info] = interval_covermut
		Candidate_Interval_CoverMut_file_in.close()
		return Candidate_Interval_CoverMut,Ori_Candidate_Interval_CoverMut
	
	for gene_info in Candidate_Gene_Exon_CoverMut.keys():
		for exon_info,mutation_list in Candidate_Gene_Exon_CoverMut[gene_info].items():
			if mutation_list:
				if gene_info not in Candidate_Interval_CoverMut:
					Candidate_Interval_CoverMut[gene_info] = {}
				Candidate_Interval_CoverMut[gene_info][exon_info] = {}
				long_mutation_list = [ mutation for mutation in mutation_list if ( int(mutation.split(":")[6]) - int(mutation.split(":")[5]) + 1 ) > 40 ]
				if long_mutation_list:
					for long_mutation in long_mutation_list:
						lung_mutation_cluster = []
						for mutation_info in mutation_list:
							if not ( int(mutation_info.split(":")[6]) < int(long_mutation.split(":")[5]) or int(mutation_info.split(":")[5]) > int(long_mutation.split(":")[6]) ):
								lung_mutation_cluster.append(mutation_info)
						long_interval_info =  gene_info + ":" + exon_info + ":" + str(min([ int(x.split(":")[5]) for x in lung_mutation_cluster  ])) + ":" + str(max([ int(x.split(":")[6]) for x in lung_mutation_cluster ]))
						Candidate_Interval_CoverMut[gene_info][exon_info][long_interval_info] = lung_mutation_cluster

				normal_mutation_list = [ mutation for mutation in mutation_list if mutation not in long_mutation_list ]
				if normal_mutation_list:
					sorted_byend_mutation_list = sorted(normal_mutation_list, key=lambda x: int(x.split(":")[6]))
					sorted_bystart_mutation_list = sorted(normal_mutation_list, key=lambda x: int(x.split(":")[5]))
					last_mutation = sorted_byend_mutation_list[-1]				
					a_bed_file = open("../output/BINGXING/%s/temp_a.bed"%(OUTPUT_DIR),"w")
					b_bed_file = open("../output/BINGXING/%s/temp_b.bed"%(OUTPUT_DIR),"w")
					for mutation_info in  sorted_bystart_mutation_list:
						a_bed_file.write(mutation_info.split(":")[4]+"\t"+str(int(mutation_info.split(":")[5])-1)+"\t"+str(mutation_info.split(":")[5])+"\t"+mutation_info+"\n")
					for mutation_info in  sorted_byend_mutation_list:
						b_bed_file.write(mutation_info.split(":")[4]+"\t"+str(int(mutation_info.split(":")[6])-1)+"\t"+str(mutation_info.split(":")[6])+"\t"+mutation_info+"\n")
					a_bed_file.close()
					b_bed_file.close()
					window_command = "bedtools window -a ../output/BINGXING/%s/temp_a.bed -b ../output/BINGXING/%s/temp_b.bed -l 0 -r 39 > ../output/BINGXING/%s/temp_window.bed"%(OUTPUT_DIR,OUTPUT_DIR,OUTPUT_DIR)
					os.system(window_command)
					step_follow_mutation = {}
					window_file_in = open("../output/BINGXING/%s/temp_window.bed"%(OUTPUT_DIR),"r")
					for line in window_file_in: 
						line = line.strip()
						if not line :
							continue
						step_mutation = line.split("\t")[3]
						follow_mutation = line.split("\t")[7]
						step_follow_mutation.setdefault(step_mutation,[]).append(follow_mutation)
						if last_mutation in follow_mutation:
							break
					window_file_in.close()
					for interval_covermut in step_follow_mutation.values():
						interval_info =  gene_info + ":" + exon_info + ":" + str(min([ int(x.split(":")[5]) for x in interval_covermut  ])) + ":" + str(max([int(x.split(":")[6]) for x in interval_covermut]))
						Candidate_Interval_CoverMut[gene_info][exon_info][interval_info] = interval_covermut
	Candidate_Interval_CoverMut_file_out = open("../output/BINGXING/%s/%s_Candidate_Interval_CoverMut.tsv"%(OUTPUT_DIR,TARGET_TUMOR),"w")
	for gene_info in Candidate_Interval_CoverMut.keys():
		for exon_info in Candidate_Interval_CoverMut[gene_info]:
			for interval_info in Candidate_Interval_CoverMut[gene_info][exon_info]:
				Candidate_Interval_CoverMut_file_out.write(gene_info+"\t"+exon_info+"\t"+interval_info+"\t"+"\t".join(Candidate_Interval_CoverMut[gene_info][exon_info][interval_info])+"\n")	
	Candidate_Interval_CoverMut_file_out.close()

	return Candidate_Interval_CoverMut


def interval_panel_coverage(Candidate_Interval_CoverMut,TARGET_TUMOR,OUTPUT_DIR):
	Candidate_Interval_CoverPanelMut = {}
	#read_old
	if os.path.exists("../output/BINGXING/%s/%s_Candidate_Interval_CoverPanelMut.tsv"%(OUTPUT_DIR,TARGET_TUMOR)):
		Candidate_Interval_CoverPanelMut_file_in = open("../output/BINGXING/%s/%s_Candidate_Interval_CoverPanelMut.tsv"%(OUTPUT_DIR,TARGET_TUMOR),"r")
		for line in Candidate_Interval_CoverPanelMut_file_in:
			line = line.strip()
			if not line or line.startswith("#"):
				continue
			gene_info = line.split("\t")[0]
			exon_info = line.split("\t")[1]
			interval_info = line.split("\t")[2]
			probe_type = line.split("\t")[3]
			interval_coverpanelmut = line.split("\t")[4:]
			if gene_info not in Candidate_Interval_CoverPanelMut:
				Candidate_Interval_CoverPanelMut[gene_info] = {}
			if exon_info not in Candidate_Interval_CoverPanelMut[gene_info]:
				Candidate_Interval_CoverPanelMut[gene_info][exon_info] = {}
			if interval_info not in Candidate_Interval_CoverPanelMut[gene_info][exon_info]:
				Candidate_Interval_CoverPanelMut[gene_info][exon_info][interval_info] = {}
			Candidate_Interval_CoverPanelMut[gene_info][exon_info][interval_info][probe_type] = interval_coverpanelmut
		Candidate_Interval_CoverPanelMut_file_in.close()

		return Candidate_Interval_CoverPanelMut

	#de novo
	Pan_Probe_Differange_CoverMut,Probe_Gene_Differange_Extend_Dict = panel_mutation_remapping(OUTPUT_DIR)
	for gene_info in Candidate_Interval_CoverMut.keys():
		Candidate_Interval_CoverPanelMut[gene_info] = {}
		for exon_info in Candidate_Interval_CoverMut[gene_info]:
			Candidate_Interval_CoverPanelMut[gene_info][exon_info] = {}
			for interval_info in Candidate_Interval_CoverMut[gene_info][exon_info]:
				interval_start = int(interval_info.split(":")[-2])
				interval_end = int(interval_info.split(":")[-1])
				Candidate_Interval_CoverPanelMut[gene_info][exon_info][interval_info] = {}
				if gene_info in Probe_Gene_Differange_Extend_Dict and exon_info in Probe_Gene_Differange_Extend_Dict[gene_info]:
					for probe_type in Probe_Gene_Differange_Extend_Dict[gene_info][exon_info].keys():
						differed_interval_list1,differed_interval_list2 = portion_interval_differ( [ [interval_start,interval_end] ],Probe_Gene_Differange_Extend_Dict[gene_info][exon_info][probe_type] )
						if not [ x for x in differed_interval_list1 if x ] :
							Candidate_Interval_CoverPanelMut[gene_info][exon_info][interval_info][probe_type] = []
							if Pan_Probe_Differange_CoverMut[gene_info][exon_info][probe_type]:
								for mutation_info in Pan_Probe_Differange_CoverMut[gene_info][exon_info][probe_type]:
									if not ( int(mutation_info.split(":")[6]) < interval_start or int(mutation_info.split(":")[5]) > interval_end ):
										Candidate_Interval_CoverPanelMut[gene_info][exon_info][interval_info][probe_type].append(mutation_info)
	Candidate_Interval_CoverPanelMut_file_out = open("../output/BINGXING/%s/%s_Candidate_Interval_CoverPanelMut.tsv"%(OUTPUT_DIR,TARGET_TUMOR),"w")
	for gene_info in Candidate_Interval_CoverPanelMut.keys():
		for exon_info in Candidate_Interval_CoverPanelMut[gene_info]:
			for interval_info in Candidate_Interval_CoverPanelMut[gene_info][exon_info]:
				for probe_type in Candidate_Interval_CoverPanelMut[gene_info][exon_info][interval_info]:
					Candidate_Interval_CoverPanelMut_file_out.write(gene_info+"\t"+exon_info+"\t"+interval_info+"\t"+probe_type+"\t"+"\t".join(Candidate_Interval_CoverPanelMut[gene_info][exon_info][interval_info][probe_type])+"\n")
	Candidate_Interval_CoverPanelMut_file_out.close()

	return Candidate_Interval_CoverPanelMut


def calculate_determined_interval_performance(Selected_Interval_Determined_CoverMut,WES_Dataset_By_Patient,Tumor_Type):
	Patient_Covered_Mutation = {}
	for tumor_type in WES_Dataset_By_Patient:
		if tumor_type in Tumor_Type:
			for patient_id in WES_Dataset_By_Patient[tumor_type]:
				Patient_Covered_Mutation[patient_id] = set([])
	for gene_info in Selected_Interval_Determined_CoverMut.keys():
		for exon_info in Selected_Interval_Determined_CoverMut[gene_info]:
			for mutation_info in Selected_Interval_Determined_CoverMut[gene_info][exon_info]:
				if mutation_info.split(":")[-1] in Tumor_Type:
					Patient_Covered_Mutation[mutation_info.split(":")[0]].add(mutation_info)
	return Patient_Covered_Mutation

def process_sub_dict(args):
	(sub_dict, Target_Mutation_Count, Gene_Interval_Improvement_input,
	 gene_coverage_count_input, WES_Determined_Gene_Exon, Patient_Covered_Mutation, WES_Patients_Count) = args
	Candidate_Interval_Improvement = {}
	for gene_info in sub_dict.keys():
		for exon_info in sub_dict[gene_info]:
			for interval_info in sub_dict[gene_info][exon_info]:
				Gene_Interval_Improvement = Gene_Interval_Improvement_input.copy()
				gene_coverage_count = gene_coverage_count_input.copy()
				Candidate_Interval_Improvement[interval_info] = {}
				interval_cover_patients = set([x.split(":")[0] for x in sub_dict[gene_info][exon_info][interval_info]])
				interval_coverage = round(len(interval_cover_patients) / WES_Patients_Count, 6)

				interval_start, interval_end = map(int, interval_info.split(":")[-2:])
				if ((gene_info not in WES_Determined_Gene_Exon) or (exon_info not in WES_Determined_Gene_Exon[gene_info]) or 
				(not have_neighbour_interval([[interval_start, interval_end]], WES_Determined_Gene_Exon[gene_info][exon_info]))):
					add_length = math.ceil((interval_end - interval_start + 1) / 40) * 40 + 80
				else:
					exon_interval_info = WES_Determined_Gene_Exon[gene_info][exon_info]
					add_length = calculate_neighbour_intervals_length(exon_interval_info + [[interval_start, interval_end]]) - calculate_neighbour_intervals_length(exon_interval_info)
					add_length = max(add_length, 1)
					
				addition_mutation_count = 0
				addition_clonal_mutation_count = 0
				addition_nb_mutation_count = 0
				Interval_Patient_Mutation = {}
				
				old_gene_coverage_count = gene_coverage_count.copy()
				added_patient_set = set()

				for mutation in sub_dict[gene_info][exon_info][interval_info]:
					Interval_Patient_Mutation.setdefault(mutation.split(":")[0],[]).append(mutation)
					patient_name, gene_name = mutation.split(":")[0], mutation.split(":")[1] + ":" + mutation.split(":")[2]
					numpy_index = Gene_Interval_Improvement.loc[Gene_Interval_Improvement['Gene_Name'] == gene_name, 'Gene_Index'].values[0]
					current_patient_set = Gene_Interval_Improvement.at[numpy_index, 'Patient_Set']
					if patient_name not in current_patient_set:
						new_patient_set = current_patient_set.union({patient_name})
						Gene_Interval_Improvement.at[numpy_index, 'Patient_Set'] = frozenset(new_patient_set)
						gene_coverage_count[numpy_index] += 1
						added_patient_set.add(patient_name)

				mask = gene_coverage_count_initial <= Target_Coverage_Ratio
				Covered_gene_Ratio_improve = round(np.sum((gene_coverage_count[mask] - gene_coverage_count_initial[mask]) / WES_Patients_Count), 6)
				gene_coverage_improve = np.sum(gene_coverage_count - old_gene_coverage_count)
				new_gene_coverage = np.count_nonzero(gene_coverage_count) - np.count_nonzero(old_gene_coverage_count)
				improve_covered_gene = gene_coverage_improve / np.count_nonzero(gene_coverage_count)
				
				for patient_id in Interval_Patient_Mutation.keys():
					patient_covered_all_mutation = list(Patient_Covered_Mutation[patient_id])
					patient_covered_clonal_mutation = [mutation for mutation in patient_covered_all_mutation if mutation.split(":")[-3] == "major"]
					patient_covered_nb_mutation = [mutation for mutation in patient_covered_all_mutation if mutation.split(":")[-2] == "nonsynonymous"]
					interval_patient_all_mutation = Interval_Patient_Mutation[patient_id]
					interval_patient_clonal_mutation = [mutation for mutation in interval_patient_all_mutation if mutation.split(":")[-3] == "major"]
					interval_patient_nb_mutation = [mutation for mutation in interval_patient_all_mutation if mutation.split(":")[-2] == "nonsynonymous"]
					
					patient_covered_all_mutation_length = len(patient_covered_all_mutation)
					combined_interval_patients_length = len(set(interval_patient_all_mutation + patient_covered_all_mutation))
					
					if patient_covered_all_mutation_length >= Target_Mutation_Count:
						pass
					elif combined_interval_patients_length < Target_Mutation_Count:
						addition_mutation_count += combined_interval_patients_length - patient_covered_all_mutation_length
					elif combined_interval_patients_length >= Target_Mutation_Count:
						addition_mutation_count += Target_Mutation_Count - patient_covered_all_mutation_length
					else:
						pass
					
					patient_covered_clonal_mutation_length = len(patient_covered_clonal_mutation)
					combined_interval_patients_clonal_length = len(set(interval_patient_clonal_mutation + patient_covered_clonal_mutation))
					if patient_covered_clonal_mutation_length >= Target_Mutation_Count:
						pass
					elif combined_interval_patients_clonal_length < Target_Mutation_Count:
						addition_clonal_mutation_count += combined_interval_patients_clonal_length - len(patient_covered_clonal_mutation)
					elif combined_interval_patients_clonal_length >= Target_Mutation_Count:
						addition_clonal_mutation_count += Target_Mutation_Count - patient_covered_clonal_mutation_length
					else:
						pass
					
					patient_covered_nb_mutation_length = len(patient_covered_nb_mutation)
					combined_interval_patients_nb_length = len(set(interval_patient_nb_mutation + patient_covered_nb_mutation))
					if patient_covered_nb_mutation_length >= Target_Mutation_Count:
						pass
					elif combined_interval_patients_nb_length < Target_Mutation_Count:
						addition_nb_mutation_count += combined_interval_patients_nb_length - len(patient_covered_clonal_mutation)
					elif combined_interval_patients_nb_length >= Target_Mutation_Count:
						addition_nb_mutation_count += Target_Mutation_Count - patient_covered_nb_mutation_length
					else:
						pass
					
				mutation_improvement = round(addition_mutation_count / add_length, 6)
				clonal_mutation_improvement = round(addition_clonal_mutation_count / add_length, 6)
				nb_mutation_improvement = round(addition_nb_mutation_count / add_length, 6)
				
				Candidate_Interval_Improvement[interval_info]["Added_covered_patients"] = added_patient_set
				Candidate_Interval_Improvement[interval_info]["improved_patients_coverage"] = gene_coverage_improve
				Candidate_Interval_Improvement[interval_info]["improved_gene_coverage"] = new_gene_coverage
				Candidate_Interval_Improvement[interval_info]["improve_covered_gene"] = improve_covered_gene
				Candidate_Interval_Improvement[interval_info]["interval_coverage"] = interval_coverage
				Candidate_Interval_Improvement[interval_info]["mutation_improvement"] = mutation_improvement
				Candidate_Interval_Improvement[interval_info]["clonal_mutation_improvement"] = clonal_mutation_improvement
				Candidate_Interval_Improvement[interval_info]["nb_mutation_improvement"] = nb_mutation_improvement
				Candidate_Interval_Improvement[interval_info]["add_length"] = add_length
				Candidate_Interval_Improvement[interval_info]["Covered_gene_Ratio_improve"] = Covered_gene_Ratio_improve
	
	def sort_key(x):
		return (-x[1]['mutation_improvement'], -x[1]['interval_coverage'],
		        -x[1]['clonal_mutation_improvement'], -x[1]["nb_mutation_improvement"],
				x[1]["add_length"])
    
	sorted_Candidate_Interval_Improvement = sorted(Candidate_Interval_Improvement.items(), key=sort_key)
    
	return sorted_Candidate_Interval_Improvement[0]
					

def candidate_interval_improvement_evaluation(Target_Mutation_Count, Candidate_Interval_CoverMut, Gene_Interval_Improvement_input, gene_coverage_count_input, WES_Determined_Gene_Exon, Patient_Covered_Mutation, WES_Patients_Count):
	"""主函数"""
	max_workers = multiprocessing.cpu_count()
	num_workers = 15
	sub_dicts = split_dict(Candidate_Interval_CoverMut, num_workers)

	tasks = [
		(sub_dict, Target_Mutation_Count, Gene_Interval_Improvement_input, gene_coverage_count_input,
		 WES_Determined_Gene_Exon, Patient_Covered_Mutation, WES_Patients_Count) for sub_dict in sub_dicts
	]

	with multiprocessing.Pool(processes=num_workers) as pool:
		results = pool.map(process_sub_dict, tasks)

	all_results = []
	for result in results:
		all_results.append(result)

	sorted_results = sorted(all_results, key=lambda x: (-x[1]['mutation_improvement'], -x[1]['interval_coverage'],
														-x[1]['improved_patients_coverage'],-x[1]['clonal_mutation_improvement'],
														-x[1]["nb_mutation_improvement"],x[1]["add_length"]))
	return sorted_results[0]
		        

def add_selected_interval(Determined_Gene_Exon_CoverMut,Ori_Candidate_Interval_CoverMut,TARGET_TUMOR,OUTPUT_DIR):
	Selected_Interval_Determined_CoverMut = Determined_Gene_Exon_CoverMut 
	black_list = []

	black_list_file_in = open("../input/intervals_blacklist/lung_paper_oversea_blacklist.bed","r")
	for line in black_list_file_in:
		line = line.strip()
		if not line or line.startswith("#"):
			continue
		black_list.append( ( int(line.split("\t")[1]),int(line.split("\t")[2]) ) )
	black_list_file_in.close()

	black_list_file_in = open("../input/intervals_blacklist/lung_TCGA_blacklist.bed","r")
	for line in black_list_file_in:
		line = line.strip()
		if not line or line.startswith("#"):
			continue
		black_list.append( ( int(line.split("\t")[1]),int(line.split("\t")[2]) ) )
	black_list_file_in.close()

	black_list_file_in = open("../input/intervals_blacklist/lung_LINCHUANG_blacklist.bed","r")
	for line in black_list_file_in:
		line = line.strip()
		if not line or line.startswith("#"):
			continue
		black_list.append( ( int(line.split("\t")[1]),int(line.split("\t")[2]) ) )
	black_list_file_in.close()

	black_list_file_in = open("../input/intervals_blacklist/lung_paper_domestic_blacklist.bed","r")
	for line in black_list_file_in:
		line = line.strip()
		if not line or line.startswith("#"):
			continue
		black_list.append( ( int(line.split("\t")[1]),int(line.split("\t")[2]) ) )
	black_list_file_in.close()

	Selected_Interval_File_in = open("../output/BINGXING/%s/%s_WES_interval_improve_coverage_selected_interval.tsv"%(OUTPUT_DIR,TARGET_TUMOR),"r")
	for line in Selected_Interval_File_in:
		line = line.strip()
		if not line or line.startswith("#") or line.startswith("Interval_ID"):
			continue
		interval_id = line.split("\t")[0]
		gene_info = interval_id.split(":")[0]+":"+interval_id.split(":")[1] 
		exon_info = interval_id.split(":")[2]+":"+interval_id.split(":")[3]+":"+interval_id.split(":")[4]
		if ( int(interval_id.split(":")[5]) , int(interval_id.split(":")[6]) ) in black_list:
			continue
		if not Ori_Candidate_Interval_CoverMut[gene_info][exon_info][interval_id]:
			print("###################Interval ID ERROR######################")
		if gene_info not in Selected_Interval_Determined_CoverMut:
			Selected_Interval_Determined_CoverMut[gene_info] = {}
		if exon_info not in Selected_Interval_Determined_CoverMut[gene_info]:
			Selected_Interval_Determined_CoverMut[gene_info][exon_info] = []
		for mutation in Ori_Candidate_Interval_CoverMut[gene_info][exon_info][interval_id]:
			if mutation not in Selected_Interval_Determined_CoverMut[gene_info][exon_info]:
				Selected_Interval_Determined_CoverMut[gene_info][exon_info].append(mutation)
	Selected_Interval_File_in.close()

	return Selected_Interval_Determined_CoverMut


def find_shortest_list_key(patients_dict_input, black_list, max_covered_for_patients, target_count):
    shortest_key = None
    shortest_length = float('inf')

    patients_dict = {key: value for key, value in patients_dict_input.items() if key not in black_list}

    for key, value in patients_dict.items():
        current_length = len(value)

        max_covered_value = max_covered_for_patients[max_covered_for_patients["names"] == key].values[0][1]  # 获取单一值
        
        # 检查条件
        if current_length + 2 >= target_count and current_length < target_count and max_covered_value >= target_count:
            return key
        elif max_covered_value >= target_count:
            if current_length < shortest_length:
                shortest_length = current_length
                shortest_key = key

    return shortest_key


def convert_to_regular_dict(dd):
    if isinstance(dd, defaultdict):
        dd = {k: convert_to_regular_dict(v) for k, v in dd.items()}
    return dd


def find_satisfied_candidate_genes(Candidate_Interval_CoverMut, Target_Coverage_Ratio, WES_Patient_Count):
	sub_Candidate_Interval_CoverMut = defaultdict(lambda: defaultdict(lambda: defaultdict()))
	target_count = round(Target_Coverage_Ratio * WES_Patient_Count)
	count_dict = defaultdict(int)
	name_dict = defaultdict(set)
	for gene_id in Candidate_Interval_CoverMut.keys():
		name_dict[gene_id] = set()
		count_dict[gene_id] = 0
		for exon_name in Candidate_Interval_CoverMut[gene_id]:
			for interval_name in Candidate_Interval_CoverMut[gene_id][exon_name]:
				for mutation in Candidate_Interval_CoverMut[gene_id][exon_name][interval_name]:
					patient_name = mutation.split(":")[0]
					if patient_name not in name_dict[gene_id]:
						name_dict[gene_id].add(patient_name)
						count_dict[gene_id] += 1
	
	satisfied_genes = []
	for gene_name, count in count_dict.items():
		if count >= target_count:
			satisfied_genes.append(gene_name)
	
	sub_Candidate_Interval_length = 0
	for gene_id in satisfied_genes:
		for exon_name in Candidate_Interval_CoverMut[gene_id]:
			for interval_name in Candidate_Interval_CoverMut[gene_id][exon_name]:
				sub_Candidate_Interval_length += 1
				sub_Candidate_Interval_CoverMut[gene_id][exon_name][interval_name] = Candidate_Interval_CoverMut[gene_id][exon_name][interval_name]

	return sub_Candidate_Interval_length, sub_Candidate_Interval_CoverMut


def find_best_adding_gene(Gene_Interval_Improvement_initial, gene_covarage_count_initial, Target_Coverage_Ratio, WES_Patients_Count, low_limit_ratio):					
    sub_Candidate_Interval_CoverMut = defaultdict(lambda: defaultdict(lambda: defaultdict()))
    goal_gene_cover_count = round(Target_Coverage_Ratio * WES_Patients_Count)
    Gene_Interval_Improvement_input = Gene_Interval_Improvement_initial.copy()
    gene_covarage_count_input = gene_covarage_count_initial.copy()
    
    mask = (gene_covarage_count_input < goal_gene_cover_count) & (gene_covarage_count_input >= round(goal_gene_cover_count * low_limit_ratio))
    sub_Gene_name_list = Gene_Interval_Improvement_input[mask]["Gene_Name"].tolist()
    sub_Candidate_Interval_length = 0
    
    for gene_id in sub_Gene_name_list:
        try:
            for exon_name in Candidate_Interval_CoverMut[gene_id]:
                for interval_name in Candidate_Interval_CoverMut[gene_id][exon_name]:
                    sub_Candidate_Interval_length += 1
                    sub_Candidate_Interval_CoverMut[gene_id][exon_name][interval_name] = Candidate_Interval_CoverMut[gene_id][exon_name][interval_name]
        except:
            continue
		
    return sub_Candidate_Interval_length, sub_Candidate_Interval_CoverMut
    

def process_sub_uniformity(sub_dict, Gene_Interval_Improvement_input, gene_coverage_count_input, WES_Determined_Gene_Exon,Patient_Covered_Mutation, WES_Patients_Count, evaluation_stop):
	Candidate_Interval_CoverMut_Patient = sub_dict
	Candidate_Interval_Improvement = {}
	for gene_info in Candidate_Interval_CoverMut_Patient.keys():
		for exon_info in Candidate_Interval_CoverMut_Patient[gene_info]:
			for interval_info in Candidate_Interval_CoverMut_Patient[gene_info][exon_info]:
				Gene_Interval_Improvement_cp = Gene_Interval_Improvement_input.copy()
				gene_coverage_count_cp = gene_coverage_count_input.copy()
				interval_cover_patients = set([x.split(":")[0] for x in Candidate_Interval_CoverMut_Patient[gene_info][exon_info][interval_info]])
				Candidate_Interval_Improvement[interval_info] = {}
				interval_coverage = round(len(interval_cover_patients)/WES_Patients_Count,6)

				interval_start, interval_end = map(int, interval_info.split(":")[-2:])
				if (gene_info not in WES_Determined_Gene_Exon) or (exon_info not in WES_Determined_Gene_Exon[gene_info]) or ( not have_neighbour_interval([ [interval_start,interval_end] ], WES_Determined_Gene_Exon[gene_info][exon_info] ) ) :
					add_length = math.ceil((interval_end-interval_start+1)/40)*40+80
				else:
					add_length = calculate_neighbour_intervals_length( WES_Determined_Gene_Exon[gene_info][exon_info] + [ [interval_start,interval_end] ] )  - calculate_neighbour_intervals_length( WES_Determined_Gene_Exon[gene_info][exon_info] )
					add_length = max(add_length, 1)
					
				addition_mutation_count = 0
				addition_clonal_mutation_count = 0
				addition_nb_mutation_count = 0
				Interval_Patient_Mutation = {}
				old_gene_coverage_count = gene_coverage_count_cp.copy()
				added_patient_set = set()

				for mutation in Candidate_Interval_CoverMut_Patient[gene_info][exon_info][interval_info]:
					Interval_Patient_Mutation.setdefault(mutation.split(":")[0],[]).append(mutation)
					patient_name, gene_name = mutation.split(":")[0], mutation.split(":")[1] + ":" + mutation.split(":")[2]
					numpy_index = Gene_Interval_Improvement_cp.loc[Gene_Interval_Improvement_cp['Gene_Name'] == gene_name, 'Gene_Index'].values[0]
					current_patient_set = Gene_Interval_Improvement_cp.at[numpy_index, 'Patient_Set']
					if patient_name not in current_patient_set:
						new_patient_set = current_patient_set.union({patient_name})
						Gene_Interval_Improvement_cp.at[numpy_index, 'Patient_Set'] = frozenset(new_patient_set)
						gene_coverage_count_cp[numpy_index] += 1
						added_patient_set.add(patient_name)

				mask = gene_coverage_count_initial <= Target_Coverage_Ratio
				Covered_gene_Ratio_improve = round(np.sum(gene_coverage_count_cp[mask] / WES_Patients_Count - gene_coverage_count_initial[mask] / WES_Patients_Count), 6)
				gene_coverage_improve = np.sum(gene_coverage_count_cp - old_gene_coverage_count)
				new_gene_coverage = np.count_nonzero(gene_coverage_count_cp) - np.count_nonzero(old_gene_coverage_count)
				improve_covered_gene = gene_coverage_improve / np.count_nonzero(gene_coverage_count_cp)

				for patient_id in Interval_Patient_Mutation.keys():
					patient_covered_all_mutation = list(Patient_Covered_Mutation[patient_id])
					patient_covered_clonal_mutation = [ mutation for mutation in patient_covered_all_mutation if mutation.split(":")[-3] == "major" ]
					patient_covered_nb_mutation = [ mutation for mutation in patient_covered_all_mutation if mutation.split(":")[-2] == "nonsynonymous" ]
					interval_patient_all_mutation = Interval_Patient_Mutation[patient_id]
					interval_patient_clonal_mutation = [mutation for mutation in interval_patient_all_mutation if mutation.split(":")[-3] == "major" ]
					interval_patient_nb_mutation = [mutation for mutation in interval_patient_all_mutation if mutation.split(":")[-2] == "nonsynonymous" ]
										
					addition_mutation_count = addition_mutation_count + (len(set(interval_patient_all_mutation + patient_covered_all_mutation)) - len(patient_covered_all_mutation))/max(len(patient_covered_all_mutation), 1) * 10
					addition_clonal_mutation_count = addition_clonal_mutation_count + (len(set(interval_patient_clonal_mutation + patient_covered_clonal_mutation)) - len(patient_covered_clonal_mutation)) /max(len(patient_covered_clonal_mutation), 1) * 10
					addition_nb_mutation_count = addition_nb_mutation_count + (len(set(interval_patient_nb_mutation + patient_covered_nb_mutation)) -  len(patient_covered_clonal_mutation)) /max(len(patient_covered_nb_mutation), 1) * 10

				mutation_improvement = round(addition_mutation_count/add_length, 6)
				clonal_mutation_improvement = round(addition_clonal_mutation_count/add_length, 6)
				nb_mutation_improvement = round(addition_nb_mutation_count/add_length, 6)
				Candidate_Interval_Improvement[interval_info]["Added_covered_patients"] = added_patient_set
				Candidate_Interval_Improvement[interval_info]["interval_coverage"] = interval_coverage
				Candidate_Interval_Improvement[interval_info]["mutation_improvement"] = mutation_improvement
				Candidate_Interval_Improvement[interval_info]["clonal_mutation_improvement"] = clonal_mutation_improvement
				Candidate_Interval_Improvement[interval_info]["nb_mutation_improvement"] = nb_mutation_improvement
				Candidate_Interval_Improvement[interval_info]["add_length"] = add_length
				Candidate_Interval_Improvement[interval_info]["improved_patients_coverage"] = gene_coverage_improve
				Candidate_Interval_Improvement[interval_info]["improved_gene_coverage"] = new_gene_coverage
				Candidate_Interval_Improvement[interval_info]["improve_covered_gene"] = improve_covered_gene
				Candidate_Interval_Improvement[interval_info]["gene_info"] = gene_info
				Candidate_Interval_Improvement[interval_info]["exon_info"] = exon_info
				Candidate_Interval_Improvement[interval_info]["Covered_gene_Ratio_improve"] = Covered_gene_Ratio_improve

	if evaluation_stop == 1:
		if not Candidate_Interval_Improvement:
			return "", False
		sorted_Candidate_Interval_Improvement = sorted(Candidate_Interval_Improvement.items(), key=lambda x: (-x[1]["Covered_gene_Ratio_improve"], -x[1]["improve_covered_gene"],
																										-x[1]['mutation_improvement'], -x[1]['interval_coverage'],
																										-x[1]["improved_patients_coverage"], -x[1]['clonal_mutation_improvement'],
																										-x[1]["nb_mutation_improvement"], x[1]["add_length"],
																										x[1]["improved_gene_coverage"]))
	else:
		if not Candidate_Interval_Improvement:
			return "", False
		sorted_Candidate_Interval_Improvement = sorted(Candidate_Interval_Improvement.items(), key=lambda x: (-x[1]['mutation_improvement'], -x[1]['interval_coverage'],
																										-x[1]['clonal_mutation_improvement'],-x[1]["nb_mutation_improvement"],
																										x[1]["add_length"], x[1]["improved_gene_coverage"]))

	return sorted_Candidate_Interval_Improvement[0]


def candidate_interval_improvement_uniformity(sub_Candidate_Interval_CoverMut,Gene_Interval_Improvement_input, gene_coverage_count_input, WES_Determined_Gene_Exon,Patient_Covered_Mutation,WES_Patients_Count, evaluation_stop, CoverMut_length_org, CoverMut_length):
	"""主函数"""
	max_workers = multiprocessing.cpu_count()
	#print(max_workers)
	num_workers = 15
	if evaluation_stop == 2:
		num_workers = 0
	elif evaluation_stop == 1 and CoverMut_length <= 2000 and CoverMut_length == CoverMut_length_org * 0.8:
		num_worker = max(num_worker * 0.8, 2)
		CoverMut_length_org = CoverMut_length_org * 0.8
	sub_dicts = split_dict(sub_Candidate_Interval_CoverMut, num_workers)

	tasks = [
        (sub_dict, Gene_Interval_Improvement_input, gene_coverage_count_input, WES_Determined_Gene_Exon,
        Patient_Covered_Mutation, WES_Patients_Count, evaluation_stop) for sub_dict in sub_dicts
    ]

	with multiprocessing.Pool(processes=num_workers) as pool:
		results = pool.starmap(process_sub_uniformity, tasks)

	all_results = []
	for result in results:
		all_results.append(result)

	all_results = [result for result in results if result]
	if not all_results:
		return "", False

	if evaluation_stop == 1:
		sorted_results = sorted(all_results, key=lambda x: (-x[1]["Covered_gene_Ratio_improve"], -x[1]["improve_covered_gene"],
															-x[1]['mutation_improvement'], -x[1]['interval_coverage'],
															-x[1]["improved_patients_coverage"], -x[1]['clonal_mutation_improvement'],
															-x[1]["nb_mutation_improvement"], x[1]["add_length"],
															x[1]["improved_gene_coverage"]))
	else:
		sorted_results = sorted(all_results, key=lambda x: (-x[1]['mutation_improvement'], -x[1]['interval_coverage'],
															-x[1]['clonal_mutation_improvement'],-x[1]["nb_mutation_improvement"],
															x[1]["add_length"], x[1]["improved_gene_coverage"]))

	return sorted_results[0]


def calculate_patient_coverage(patient_covered_mutation, total_patients):
	return round(len([x for x in patient_covered_mutation.keys() if len(patient_covered_mutation[x]) > 0]) / total_patients, 6)


def calculate_patient_mutation_counts(patient_covered_mutation):
	return {
		'1st': np.percentile([len(x) for x in patient_covered_mutation.values()], 1),
		'5th': np.percentile([len(x) for x in patient_covered_mutation.values()], 5),
		'25th': np.percentile([len(x) for x in patient_covered_mutation.values()], 25),
		'50th': np.percentile([len(x) for x in patient_covered_mutation.values()], 50),
		'75th': np.percentile([len(x) for x in patient_covered_mutation.values()], 75),
	}


def split_dict(d, n):
    """将字典 d 切分为 n 个子字典"""
    items = list(d.items())
    chunks = [dict(items[i::n]) for i in range(n)]
    return chunks


def final_determined_output(OLD_TUMOR,TARGET_TUMOR,OUTPUT_DIR):
	WES_All_Gene_Exon = read_ncbi_anno_file()
	Final_Determined_Gene_Exon = read_determined_gene_exon(WES_All_Gene_Exon,OLD_TUMOR,OUTPUT_DIR)
	if os.path.exists("../output/BINGXING/%s/%s_WES_interval_improve_coverage_selected_interval_v4.tsv"%(OUTPUT_DIR,TARGET_TUMOR)):
		selected_interval_file_in = open("../output/BINGXING/%s/%s_WES_interval_improve_coverage_selected_interval_v4.tsv"%(OUTPUT_DIR,TARGET_TUMOR),"r")
		for line in selected_interval_file_in:
			line = line.strip()
			if not line or line.startswith("Interval_ID"):
				continue
			gene_info = line.split("\t")[1]
			exon_info = line.split("\t")[2]
			interval_start = int(line.split("\t")[6])
			interval_end = int(line.split("\t")[7])
			if gene_info not in Final_Determined_Gene_Exon:
				Final_Determined_Gene_Exon[gene_info] = {}
			Final_Determined_Gene_Exon[gene_info].setdefault(exon_info,[]).append([interval_start,interval_end])
		selected_interval_file_in.close()
	for gene_info in Final_Determined_Gene_Exon.keys():
		for exon_info in Final_Determined_Gene_Exon[gene_info]:
			temp_coveranges = Final_Determined_Gene_Exon[gene_info][exon_info]
			Final_Determined_Gene_Exon[gene_info][exon_info] = portion_interval_merge(temp_coveranges)

	final_determined_gene_exon_file_out = open("../output/BINGXING/%s/%s_WES_Determined_Gene_Exon_v4.tsv"%(OUTPUT_DIR,TARGET_TUMOR),"w")
	final_determined_gene_exon_file_out.write("Gene\tExon\tChr\tExon_Start\tExon_End\tIntervals\tMerged_Intervals\tMerged_Length\tMerged_Intervals_Fraction\n")
	final_target_interval_file_out = open("../output/BINGXING/%s/%s_WES_Determined_Gene_Exon_BY_TARGET_INTERVAL_v4.tsv"%(OUTPUT_DIR,TARGET_TUMOR),"w")
	final_target_interval_file_out.write("gene\tExon\tChr\tExon_Start\tExon_End\tTarget_Interval_Start\tTarget_Interval_End\tTarget_Interval_Fraction\n")
	for gene_info in Final_Determined_Gene_Exon.keys():
		for exon_info in Final_Determined_Gene_Exon[gene_info]:
			merged_intervals,merged_intervals_length = merge_neighbour_intervals(Final_Determined_Gene_Exon[gene_info][exon_info])			
			final_determined_gene_exon_file_out.write(gene_info+"\t"+exon_info+"\t"+str(WES_All_Gene_Exon[gene_info][exon_info][2])+"\t"+str(WES_All_Gene_Exon[gene_info][exon_info][0])+"\t"+str(WES_All_Gene_Exon[gene_info][exon_info][1])+"\t"+str(Final_Determined_Gene_Exon[gene_info][exon_info])+"\t"+str(merged_intervals)+"\t"+str(merged_intervals_length)+"\t"+str(round( merged_intervals_length/(WES_All_Gene_Exon[gene_info][exon_info][1]-WES_All_Gene_Exon[gene_info][exon_info][0]+1),6 )) +"\n")
			
			for interval in merged_intervals:
				final_target_interval_file_out.write(gene_info+"\t"+exon_info+"\t"+str(WES_All_Gene_Exon[gene_info][exon_info][2])+"\t"+str(WES_All_Gene_Exon[gene_info][exon_info][0])+"\t"+str(WES_All_Gene_Exon[gene_info][exon_info][1])+"\t"+str(interval[0])+"\t"+str(interval[1])+"\t"+str(round( (interval[1]-interval[0]+1)/(WES_All_Gene_Exon[gene_info][exon_info][1]-WES_All_Gene_Exon[gene_info][exon_info][0]+1),6 ))+"\n")
	final_determined_gene_exon_file_out.close()
	final_target_interval_file_out.close()


def restructure_candidate_interval_cover_mut(Candidate_Interval_CoverMut):
	patient_to_interval = defaultdict(set)
	interval_to_patient = defaultdict(set)

	for gene_info in Candidate_Interval_CoverMut.keys():
		for exon_info in Candidate_Interval_CoverMut[gene_info]:
			for interval_info in Candidate_Interval_CoverMut[gene_info][exon_info]:
				for mutation in Candidate_Interval_CoverMut[gene_info][exon_info][interval_info]:
					temp = mutation.split(":")
					patient_info = temp[0]
					patient_to_interval[patient_info].add(interval_info)
					interval_to_patient[interval_info].add(patient_info)

	return patient_to_interval, interval_to_patient

def calculate_max_mut_for_patient(Candidate_Interval_CoverMut, Patient_Covered_Mutation):
	max_count_dict = defaultdict(set)
	for gene_info in Candidate_Interval_CoverMut.keys():
		for exon_info in Candidate_Interval_CoverMut[gene_info]:
			for interval_info in Candidate_Interval_CoverMut[gene_info][exon_info]:
				for mutation in Candidate_Interval_CoverMut[gene_info][exon_info][interval_info]:
					temp = mutation.split(":")
					patient_info = temp[0]
					max_count_dict[patient_info].add(mutation)

	for patient_name, mutation_set in Patient_Covered_Mutation.items():
		max_count_dict[patient_name].update(mutation_set)

	pandas_data = []
	for patient_name, mutation_set in max_count_dict.items():
		pandas_data.append([patient_name, len(mutation_set)])
	max_cover_patient_df = pd.DataFrame(pandas_data, columns=["names", "max_count"])

	return max_cover_patient_df


def add_interval_by_patient(patient_id):
	sub_Candidate_Interval_CoverMut = defaultdict(lambda: defaultdict(lambda: defaultdict()))
	for interval in patient_to_interval[patient_id]:
		temp = interval.split(":")
		gene_info = temp[0] + ":" + temp[1]
		exon_info = temp[2] + ":" + temp[3] + ":" + temp[4]
		interval_info = interval
		sub_Candidate_Interval_CoverMut[gene_info][exon_info][interval_info] = Candidate_Interval_CoverMut[gene_info][exon_info][interval_info]
	return sub_Candidate_Interval_CoverMut


def delete_interval_for_uniformity(add_interval_info):
	rm_patient_to_interval = {}
	for patient_name in interval_to_patient[add_interval_info]:
		for interval in patient_to_interval[patient_name]:
			if interval == add_interval_info:
				rm_patient_to_interval.update({patient_name:add_interval_info})

	rm_interval_to_patient = set()
	for patient_name, interval_info in rm_patient_to_interval.items():
		patient_to_interval[patient_name].remove(interval_info)
		rm_interval_to_patient.add(patient_name)
		if len(patient_to_interval[patient_name]) == 0:
			del patient_to_interval[patient_name]

	for patient_name in rm_interval_to_patient:
		interval_to_patient[add_interval_info].remove(patient_name)
	del interval_to_patient[add_interval_info]



def main():
	WES_All_Gene_Exon = read_ncbi_anno_file()
	WES_Dataset_By_Gene,WES_Dataset_By_Patient = read_WES_file("lung", ("train/lingchuang_lung_clear_train3.tsv","train/paper_lung34_1207_remapping_successed_mutation_train3.tsv",
																		"train/paper_tracerX_1207_remapping_successed_mutation_train3.tsv","train/keyan_lung_clear_train3.tsv",
																		"train/database_TCGA_lung_1207_remapping_successed_mutation_train3.tsv",))
	TARGET_MUTATION_COUNT = 5
	TARGET_PATIENT_COVERAGE = 0.99
	WES_Patients_Count = 3168
	Target_condition_Covered_Gene = 500 #How many genes to desired coverage?
	Target_Coverage_Ratio = 0.02 #How much coverage you wat?
	gene_coverage_count_initial = np.zeros(10)
	OUT_PUT_DIR = "lung_train3"
	WES_Determined_Gene_Exon = read_determined_gene_exon(WES_All_Gene_Exon,"Part1_2",OUT_PUT_DIR)
	WES_Candidate_Gene_Exon = candidate_gene_exon_generator(WES_Determined_Gene_Exon,WES_All_Gene_Exon,"lung",OUT_PUT_DIR)
	Candidate_Gene_Exon_CoverMut = wes_mutation_remapping(WES_Dataset_By_Gene,WES_Candidate_Gene_Exon,"Candidate","lung",OUT_PUT_DIR)
	Candidate_Gene_Exon_CoverMut = wes_mutation_remapping(WES_Dataset_By_Gene,WES_Candidate_Gene_Exon,"Candidate","lung",OUT_PUT_DIR)
	Determined_Gene_Exon_CoverMut = wes_mutation_remapping(WES_Dataset_By_Gene,WES_Determined_Gene_Exon,"Determined","lung",OUT_PUT_DIR)
	Candidate_Interval_CoverMut,Ori_Candidate_Interval_CoverMut = interval_generator_bystep(Candidate_Gene_Exon_CoverMut,WES_Candidate_Gene_Exon,"lung",OUT_PUT_DIR)
	Candidate_Interval_CoverPanelMut =  interval_panel_coverage(Candidate_Interval_CoverMut,"lung", OUT_PUT_DIR)
	Selected_Interval_Determined_CoverMut = add_selected_interval(Determined_Gene_Exon_CoverMut,Ori_Candidate_Interval_CoverMut,"lung", OUT_PUT_DIR)

	selected_interval_file_out = open("../output/BINGXING/%s/%s_WES_interval_improve_coverage_selected_interval_v4.tsv"%(OUT_PUT_DIR,"lung"),"w")
	selected_interval_file_out.write("Interval_ID\tGene\tExon\tchr\texon_start\texon_end\tinterval_start\tinterval_end\tadd_length\tmutation_improvement\timproved_patients_coverage\timproved_gene_coverage\timprove_covered_gene\tinterval_coverage\tclonal_mutation_improvement\tnb_mutation_improvement\tPatient_Coverage\tDetected_Mutation_Count(1)\tDetected_Mutation_Count(5)\tDetected_Mutation_Count(25)\tDetected_Mutation_Count(50)\tDetected_Mutation_Count(75)\n")
	selected_interval_file_out.flush()

	performance_file_out = open("../output/BINGXING/%s/%s_WES_interval_improve_performance_v4.tsv"%(OUT_PUT_DIR,"lung"),"w")
	performance_file_out.write("Interval_Total_Length\tInterval_Total_Length(Include_determined_interval)" + "\tCovered_Patient_Count\tTotal_Patient_Count\tPatient_Coverage\tDetected_Mutation_Count(1)\tDetected_Mutation_Count(5)\tDetected_Mutation_Count(25)\tDetected_Mutation_Count(50)\tDetected_Mutation_Count(75)"*1 + "\n")
	performance_file_out.flush()

	Determined_Gene_Exon_Length = 0
	for gene_info in WES_Determined_Gene_Exon.keys():
		for exon_info in WES_Determined_Gene_Exon[gene_info]:
			for interval in WES_Determined_Gene_Exon[gene_info][exon_info]:
				Determined_Gene_Exon_Length = Determined_Gene_Exon_Length + (interval[1] - interval[0] + 1)

	performance_file_out.write(str(Determined_Gene_Exon_Length)+"\t"+str(Determined_Gene_Exon_Length)+"\t")

	Patient_Covered_Mutation = calculate_determined_interval_performance(Selected_Interval_Determined_CoverMut,WES_Dataset_By_Patient,"lung")
	max_covered_for_patients = calculate_max_mut_for_patient(Candidate_Interval_CoverMut, Patient_Covered_Mutation)
	Covered_Patients_Count = len([ x for x in Patient_Covered_Mutation.keys() if len(Patient_Covered_Mutation[x]) > 0 ])
	Patient_Coverage = round( len([ x for x in Patient_Covered_Mutation.keys() if len(Patient_Covered_Mutation[x]) > 0 ])/WES_Patients_Count,6 )
	print(Patient_Covered_Mutation.keys())
	Patient_CoverMut_Count_25 = np.percentile([len(x) for x in list(Patient_Covered_Mutation.values())],25)
	Patient_CoverMut_Count_50 = np.percentile([len(x) for x in list(Patient_Covered_Mutation.values())],50)
	Patient_CoverMut_Count_75 = np.percentile([len(x) for x in list(Patient_Covered_Mutation.values())],75)
	Patient_CoverMut_Count_5 = np.percentile([len(x) for x in list(Patient_Covered_Mutation.values())],5)
	Patient_CoverMut_Count_1 = np.percentile([len(x) for x in list(Patient_Covered_Mutation.values())],1)

	performance_file_out.write("\t".join([str(x) for x in [Covered_Patients_Count,WES_Patients_Count,Patient_Coverage,Patient_CoverMut_Count_1,Patient_CoverMut_Count_5,Patient_CoverMut_Count_25,Patient_CoverMut_Count_50,Patient_CoverMut_Count_75] ]) + "\n" )
	performance_file_out.flush()

	#保存一下初始的WES_Determined_Gene_Exon
	Original_WES_Determined_Gene_Exon = copy.deepcopy(WES_Determined_Gene_Exon)

	def default_set():
		return defaultdict(set)

	for TUMOR_TYPE in [tuple(["lung"])]:
		print("ready to add")
		stop_flag, evaluation_stop = 0, 0

		gene_data = []
		Gene_name_index = 0
		for gene_name in WES_Dataset_By_Gene.keys():
			gene_data.append([gene_name, Gene_name_index, frozenset()])
			Gene_name_index += 1
		Gene_Interval_Improvement_initial = pd.DataFrame(gene_data, columns=['Gene_Name', 'Gene_Index', 'Patient_Set'])
		gene_coverage_count_initial = np.zeros(len(Gene_Interval_Improvement_initial))

		for key, mutations in Patient_Covered_Mutation.items():
			all_mutations = list(mutations)
			All_genes = set(map(lambda x: x.split(":")[1] + ":" + x.split(":")[2], all_mutations))
			for gene_info in All_genes:
				numpy_index = Gene_Interval_Improvement_initial.loc[Gene_Interval_Improvement_initial['Gene_Name'] == gene_info, 'Gene_Index'].values[0]
				current_patient_set = Gene_Interval_Improvement_initial.at[numpy_index, 'Patient_Set']
				new_patient_set = current_patient_set.union({key})
				Gene_Interval_Improvement_initial.at[numpy_index, 'Patient_Set'] = frozenset(new_patient_set)
				gene_coverage_count_initial[numpy_index] += 1

		patient_black_list = set()
		target_complete_count = 5
		find_best_adding_gene_flag = 1
		CoverMut_length, sub_Candidate_Interval_CoverMut = 0, {}
		no_add_according_genecov = 0

		while True:
			Covered_gene_Count_pre = np.sum(gene_coverage_count_initial)
			if Patient_Coverage >= TARGET_PATIENT_COVERAGE and Patient_CoverMut_Count_1 >= TARGET_MUTATION_COUNT and stop_flag == 1:
				break
			if evaluation_stop == 0:
				add_interval_info, add_interval_improvement = candidate_interval_improvement_evaluation(TARGET_MUTATION_COUNT,Candidate_Interval_CoverMut,
																							Gene_Interval_Improvement_initial,gene_coverage_count_initial,
																							WES_Determined_Gene_Exon,Patient_Covered_Mutation,WES_Patients_Count)
				if not add_interval_improvement["mutation_improvement"] > 0:
					print("the first step come to an end, let's start the second stage")
					evaluation_stop = 1
					continue

			elif evaluation_stop == 1 and no_add_according_genecov == 0:
				if find_best_adding_gene_flag == 1:
					CoverMut_length_org, sub_Candidate_Interval_CoverMut = find_satisfied_candidate_genes(Candidate_Interval_CoverMut, Target_Coverage_Ratio, WES_Patients_Count)
					convert_Candidate_Interval_CoverMut = convert_to_regular_dict(sub_Candidate_Interval_CoverMut)
					CoverMut_length = CoverMut_length_org
					find_best_adding_gene_flag == 0
				else:
					pass
				add_interval_info, add_interval_improvement = candidate_interval_improvement_uniformity(convert_Candidate_Interval_CoverMut,Gene_Interval_Improvement_initial, gene_coverage_count_initial, WES_Determined_Gene_Exon,
																										Patient_Covered_Mutation, WES_Patients_Count, evaluation_stop, CoverMut_length_org, CoverMut_length)
				CoverMut_length -= 1
				if add_interval_improvement == False:
					no_add_according_genecov = 1
					patient_to_interval, interval_to_patient = restructure_candidate_interval_cover_mut(Candidate_Interval_CoverMut)
					find_best_adding_gene_flag = 2
					continue
				del sub_Candidate_Interval_CoverMut[add_interval_improvement["gene_info"]][add_interval_improvement["exon_info"]][add_interval_info]

			else:
				min_patient = find_shortest_list_key(Patient_Covered_Mutation, patient_black_list, max_covered_for_patients, target_complete_count)
				sub_Candidate_Interval_CoverMut = add_interval_by_patient(min_patient)
				add_interval_info, add_interval_improvement = candidate_interval_improvement_uniformity(sub_Candidate_Interval_CoverMut,Gene_Interval_Improvement_initial, gene_coverage_count_initial,
																										WES_Determined_Gene_Exon,Patient_Covered_Mutation,TUMOR_TYPE,WES_Patients_Count, evaluation_stop)
				if not add_interval_improvement:
					patient_black_list.add(min_patient)
					continue
				if add_interval_improvement["mutation_improvement"] <= 0:
					patient_black_list.add(min_patient)
					if len(patient_black_list) >= WES_Patients_Count / 2 or min_patient == None:
						break
					else:
						continue
				try:
					delete_interval_for_uniformity(add_interval_info)
				except:
					pass


			add_gene_info = add_interval_info.split(":")[0]+":"+add_interval_info.split(":")[1]
			add_exon_info = add_interval_info.split(":")[2]+":"+add_interval_info.split(":")[3]+":"+add_interval_info.split(":")[4]
			print("the previous numpy: ", np.sum(gene_coverage_count_initial))
			for patient_name in add_interval_improvement["Added_covered_patients"]:
				numpy_index = Gene_Interval_Improvement_initial.loc[Gene_Interval_Improvement_initial['Gene_Name'] == add_gene_info, 'Gene_Index'].values[0]
				if patient_name not in Gene_Interval_Improvement_initial.at[numpy_index, 'Patient_Set']:
					current_patient_set = Gene_Interval_Improvement_initial.at[numpy_index, 'Patient_Set']
					new_patient_set = current_patient_set.union({patient_name})
					Gene_Interval_Improvement_initial.at[numpy_index, 'Patient_Set'] = frozenset(new_patient_set)
					gene_coverage_count_initial[numpy_index] += 1
			print("the present numpy: ", np.sum(gene_coverage_count_initial))
			add_interval_start = int(add_interval_info.split(":")[-2])
			add_interval_end = int(add_interval_info.split(":")[-1])
			WES_Determined_Gene_Exon.setdefault(add_gene_info, {}).setdefault(add_exon_info, []).append([add_interval_start, add_interval_end])
			temp_coveranges = WES_Determined_Gene_Exon[add_gene_info][add_exon_info]
			WES_Determined_Gene_Exon[add_gene_info][add_exon_info] = portion_interval_merge(temp_coveranges)
			add_interval_mutation = Candidate_Interval_CoverMut[add_gene_info][add_exon_info][add_interval_info]
			Selected_Interval_Determined_CoverMut.setdefault(add_gene_info, {}).setdefault(add_exon_info, [])
			for add_mutation in add_interval_mutation:
				if add_mutation not in Selected_Interval_Determined_CoverMut[add_gene_info][add_exon_info]:
					Selected_Interval_Determined_CoverMut[add_gene_info][add_exon_info].append(add_mutation)

			selected_interval_file_out.write("\t".join([
				str(x) for x in [
					add_interval_info, add_gene_info, add_exon_info,
					WES_All_Gene_Exon[add_gene_info][add_exon_info][2],
					WES_All_Gene_Exon[add_gene_info][add_exon_info][0],
					WES_All_Gene_Exon[add_gene_info][add_exon_info][1],
					add_interval_start, add_interval_end,
					add_interval_improvement["add_length"],
					add_interval_improvement["mutation_improvement"],
					add_interval_improvement["improved_patients_coverage"],
					add_interval_improvement["improved_gene_coverage"],
					add_interval_improvement["improve_covered_gene"],
					add_interval_improvement["interval_coverage"],
					add_interval_improvement["clonal_mutation_improvement"],
					add_interval_improvement["nb_mutation_improvement"]
				]
			]))

			selected_interval_file_out.flush()

			del Candidate_Interval_CoverMut[add_gene_info][add_exon_info][add_interval_info]
			if len(Candidate_Interval_CoverMut[add_gene_info][add_exon_info]) == 0:
				del Candidate_Interval_CoverMut[add_gene_info][add_exon_info]
			if len(Candidate_Interval_CoverMut[add_gene_info]) == 0:
				del Candidate_Interval_CoverMut[add_gene_info]

			Patient_Covered_Mutation = calculate_determined_interval_performance(Selected_Interval_Determined_CoverMut, WES_Dataset_By_Patient, TUMOR_TYPE)
			Patient_Coverage = calculate_patient_coverage(Patient_Covered_Mutation, WES_Patients_Count)

			Patient_CoverMut_Count_1 = np.percentile([len(x) for x in list(Patient_Covered_Mutation.values())],1)
			Patient_Coverage_5 = round(len([x for x in Patient_Covered_Mutation.keys() if len(Patient_Covered_Mutation[x]) >= 5]) / WES_Patients_Count, 6)
			Patient_Coverage_10 = round(len([x for x in Patient_Covered_Mutation.keys() if len(Patient_Covered_Mutation[x]) >= 10]) / WES_Patients_Count, 6)
			Patient_Coverage_15 = round(len([x for x in Patient_Covered_Mutation.keys() if len(Patient_Covered_Mutation[x]) >= 15]) / WES_Patients_Count, 6)
			Patient_Coverage_20 = round(len([x for x in Patient_Covered_Mutation.keys() if len(Patient_Covered_Mutation[x]) >= 20]) / WES_Patients_Count, 6)
			mutation_counts = calculate_patient_mutation_counts(Patient_Covered_Mutation)
			selected_interval_file_out.write("\t" + str(Patient_Coverage) + "\t" + str(mutation_counts["1st"]) + "\t" + str(mutation_counts["5th"]) + "\t" + str(mutation_counts["25th"]) + "\t"+str(mutation_counts["50th"]) + "\t"+str(mutation_counts["75th"]) + "\n")
			selected_interval_file_out.flush()

			print("Difference for each loop: ", np.sum(gene_coverage_count_initial) - Covered_gene_Count_pre)
			Covered_gene_Ratio = np.count_nonzero((gene_coverage_count_initial / WES_Patients_Count) >= Target_Coverage_Ratio)
			print(Patient_Coverage, "\t", Patient_CoverMut_Count_1, "\t", Patient_Coverage_5, "\t", Patient_Coverage_10, "\t", Patient_Coverage_15, "\t", Patient_Coverage_20, "\t", Covered_gene_Ratio)

			if Patient_Coverage_5 > 0.8:
				target_complete_count = 10
				if Patient_Coverage_10 > 0.5:
					target_complete_count = 15
					if Patient_Coverage_15 > 0.3:
						target_complete_count = 20
						if Patient_Coverage_20 > 0.2:
							stop_flag = 1

		for i in range(21):
			Patient_Coverage = round(len([ x for x in Patient_Covered_Mutation.keys() if len(Patient_Covered_Mutation[x]) >= i])/WES_Patients_Count,6)
			print("Patient_Coverage" + str(i), "\t", Patient_Coverage)
	selected_interval_file_out.close()

	Determined_Gene_Exon_Length = 0
	for gene_info in WES_Determined_Gene_Exon.keys():
		for exon_info in WES_Determined_Gene_Exon[gene_info]:
			for interval in WES_Determined_Gene_Exon[gene_info][exon_info]:
				Determined_Gene_Exon_Length = Determined_Gene_Exon_Length + (interval[1] - interval[0] + 1)

	performance_file_out.write(str(Determined_Gene_Exon_Length)+"\t"+str(Determined_Gene_Exon_Length)+"\t")

	for TUMOR_TYPE in [tuple(["lung"])]:
		Patient_Covered_Mutation = calculate_determined_interval_performance(Selected_Interval_Determined_CoverMut,WES_Dataset_By_Patient,TUMOR_TYPE)
		Covered_Patients_Count = len([ x for x in Patient_Covered_Mutation.keys() if len(Patient_Covered_Mutation[x]) > 0 ])
		Patient_Coverage = round( len([ x for x in Patient_Covered_Mutation.keys() if len(Patient_Covered_Mutation[x]) > 0 ])/WES_Patients_Count,6 )

		Patient_CoverMut_Count_25 = np.percentile([len(x) for x in list(Patient_Covered_Mutation.values())],25)
		Patient_CoverMut_Count_50 = np.percentile([len(x) for x in list(Patient_Covered_Mutation.values())],50)
		Patient_CoverMut_Count_75 = np.percentile([len(x) for x in list(Patient_Covered_Mutation.values())],75)
		Patient_CoverMut_Count_5 = np.percentile([len(x) for x in list(Patient_Covered_Mutation.values())],5)
		Patient_CoverMut_Count_1 = np.percentile([len(x) for x in list(Patient_Covered_Mutation.values())],1)

		performance_file_out.write("\t".join([ str(x) for x in [Covered_Patients_Count,WES_Patients_Count,Patient_Coverage,Patient_CoverMut_Count_1,Patient_CoverMut_Count_5,Patient_CoverMut_Count_25,Patient_CoverMut_Count_50,Patient_CoverMut_Count_75] ]) + "\n" )
		performance_file_out.flush()
	performance_file_out.close()


if __name__ == "__main__":
	final_determined_output("Part1_2", "lung", OUT_PUT_DIR)

