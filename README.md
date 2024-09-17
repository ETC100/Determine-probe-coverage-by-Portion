# Determine-probe-coverage-by-Portion
# What is portion?
Portion is a python package for set operation. If you use it well, you can even achieve a "pseudo bedtools" in your python program.  
The first thing we need to do is to know how to use portion, and then, we will learn how to determine the probe region for disease monitoring.  
Let's install the portion now.
```
pip install portion
```
and import it
```python
import portion as P # import portion
interval = P.closedopen(2, 4) #(CLOSED, 2, 4, OPEN), four args here, two for open/closed, two for number recording.
```
The document for portion is here. I will help you to remember the functions in portion package.  
https://pypi.org/project/portion/  
```python
## "closed" on the left or right means closed interval on the left or right, as for "open", it's the same.
P.open(1, 2) # (1,2)
P.closed(1, 2) # [1,2]
P.openclosed(1, 2) # (1,2]
P.closedopen(1, 2) # {1,2)
P.singleton(1) # only to contain one number
P.empty() # create empty P object
interval = P.open(1, 2)
interval.empty # judge if the interval contain nothing, here return False P.open(0, 0) True
```
The most widly used function are "and"/"or"(intersection/union) operation.
```python
# "|" means or, "&" means and.
(P.open(1, 11) | P.closed(0, 1) | P.closed(20, 21)) # [0, 11) | [20, 21]
P.open(1, 11) & P.closed(0, 1) # () minus P.empty()
P.closed(1, 11) & P.closed(0, 1) # [1], note that I use P.closed for (1, 11)
~P.closed(0, 1) # reverse operation, (-inf,0) | (1,+inf) inf, P.inf
P.closed(1,2) in P.closed(0, 3) # True, the left interval surely in the right
```
A portion object may contain one or several intervals, we can use atomic to judge if there is only one in it.  
This is usefu
```python
interval = (P.open(1, 11) | P.closed(0, 1) | P.closed(20, 21))
interval.atomic # return False, here are two intervals
P.open(0, 1).atomic # return True
```
Acquire the smallest interval covering provided portion object
```python
(P.closed(0, 1) | P.open(2, 3)).enclosure #[0,3)
interval = (P.open(1, 11) | P.closed(0, 1) | P.closed(20, 21))
interval.enclosure # [0 ,21]
```
Modfiy portion object needs replace function, no exterior setter.
```python
i = P.closed(0, 2)
i.replace(P.OPEN, -1, 3, P.CLOSED) # (-1,3]
```
Portion can also work as the dictionary.
```python
## dictionary
d = P.IntervalDict()
d[P.closed(0, 3)] = 'banana'
d[4] = 'apple'
{[0,3]: 'banana', [4]: 'apple'}
d.find('apple') # 4

list(d.keys()) # [[0,2), [2,4]]
list(d.values()) # ['banana', 'orange']
list(d.items()) # [([0,2), 'banana'), ([2,4], 'orange')]
```
Output and Input
```python
s = P.to_string(P.closedopen(0, 1)) # '[0,1)', we can store this string in somewhere
interval = P.from_string(s) # read the string
interval # [0,1), back to portion object
```

Okay, let's start to design the probe for disease monitoring. Here we use pan-cancer as the exmaple.  
The first thing is to find out the regions we should monitor.  
Data we need:  
1. Known conductive regions, from reference or disease consensus (Fixed data)
2. Your in-house regions, from your experience (Candidate data)

# Step1 Read the NCBI
```python
'''
1	11873	12227	+|DDX11L1|R1|EX1|NR_046018.2|Y|UTR
1	12227	12612	+|DDX11L1|IVS|IVS1|NR_046018.2|Y|IVS
1	12612	12721	+|DDX11L1|R2|EX2|NR_046018.2|Y|UTR
1	12721	13220	+|DDX11L1|IVS|IVS2|NR_046018.2|Y|IVS
1	13220	14409	+|DDX11L1|R3E|EX3|NR_046018.2|Y|UTR
'''
ncbi_anno_file_in = open()
for line in ncbi_anno_file_in:
		line = line.strip()
		if not line or line.startswith("#"):
			continue
		gene_annotation = line.split("\t")[3]
		#select CDS
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
```

# Step2, read and initialize your own data, both determined and candidate data
Generally, there are two kinds of mutation data. Mutation on conductive cancer genes and mutation on other genes with unclear correlation with cancer. The coding for this part is in "probe_design_by_portion.py".
1. As for the first class, we only need to join the mutation togther into determined interval, with several tens externsion.
2. As for the second class, read the data and also join them as candidate interval, while cut all determined intervals out.
```python
#data fomat: gene name\texon name\tinfo
#info: patient_name:gene_name:transcript_name:exon_num:chrom:start:end:base:base:base_alt:aa_alt:clonal:is_synonymous:tissue:gene_name
'''
KCTD15:NM_001129994.1	KCTD15:NM_001129994.1:EX3	CRUK0005:KCTD15:NM_001129994.1:EX3:19:34291428:34291428:G:T:NA:NA:minor:nonsynonymous:lung:KCTD15:
KCTD15:NM_001129994.1	KCTD15:NM_001129994.1:EX7	CRUK0269:KCTD15:NM_001129994.1:EX7:19:34303796:34303796:G:T:c.G795T:p.E265D:minor:nonsynonymous:lung:KCTD15:
KCTD15:NM_001129994.1	KCTD15:NM_001129994.1:EX6	CRUK0418:KCTD15:NM_001129994.1:EX6:19:34302440:34302440:C:T:c.C676T:p.R226W:minor:nonsynonymous:lung:KCTD15:
KCTD15:NM_001129994.1	KCTD15:NM_001129994.1:EX6	CRUK0418:KCTD15:NM_001129994.1:EX6:19:34302180:34302180:G:A:c.G416A:p.R139H:minor:nonsynonymous:lung:KCTD15:
KCTD15:NM_001129994.1	KCTD15:NM_001129994.1:EX4	CRUK0533:KCTD15:NM_001129994.1:EX4:19:34292172:34292172:C:A:c.C167A:p.A56E:minor:nonsynonymous:lung:KCTD15:
KCTD15:NM_001129994.1	KCTD15:NM_001129994.1:EX6	CRUK0274:KCTD15:NM_001129994.1:EX6:19:34302319:34302319:C:T:c.C555T:p.G185G:minor:synonymous:lung:KCTD15:
KCTD15:NM_001129994.1	KCTD15:NM_001129994.1:EX7	CRUK0060:KCTD15:NM_001129994.1:EX7:19:34303797:34303797:C:T:c.C796T:p.R266W:minor:nonsynonymous:lung:KCTD15:
'''
```

# Step3, select candidate intervals and add them to your interval set.
1. Calculate how many patients covered by the determined intervals, and start to calculate the contribution of candidate intervals to your patients.
2. Reach the threshold, such as 99% patients with at least 5 mutations detected.

# Step4, extend the interval and obtain the final probes.
1. 3X method
2. Superposition method

# How to extract probes for sole cancer type from pan-cancer probes?
This is a problem given by my colleague， which is quite easy from the aspect of computer science.  
If the patient coverage for pan-cancer probes is 99% (5 mutation for each patient) in test data, we can also extract a subset with as few as probes for sole cancer type with 99% coverage.
1. the first thing is to get the pan-cancer probe-patient coverage. The test data is in the data cache.
```bash
bedtools -a pan-cancer.bed -b patient_mut.bed -wa -wb > probe_mut_patient.bed
```
I recommend to name each probe and only remain the names of probes in the following analysis.  
2. Sort the probes, find the best one to start. All the mutation count for patients are 0 now, so the probe covering the most patients is the best one. This rule works for at least 5 turns.  
3. According to the patient-mutation, sort the probe again, because patients with 5 mutations covered have no priority temporarily.  
4. Find the patients with relatively few mutation coaverage, although the have reached the coverage goals. Use more probes to cover them.   
5. Add some meaningful probes back from pan-caner probes set again.  
```python
import sys
import numpy as np
import pandas as pd
from collections import defaultdict
from multiprocessing import Pool, cpu_count

def read_mutation(file):
    probe_mutation = pd.read_csv(file, names=["probe_name", "mutation"], sep="\t", header=None)
    
    # extract patient names
    probe_mutation["patient_name"] = probe_mutation["mutation"].str.split(":").str[0]
    
    # unique the names
    patients = list(probe_mutation["patient_name"].unique())
    records = len(patients)
    numbers = np.arange(records).tolist()
    
    # create patient_mutation DataFrame
    patient_mutation = pd.DataFrame({"patient_name": patients, "numpy_index": numbers})
    patient_mutation["covered_mut"] = [frozenset() for _ in range(records)]
    patient_mut_count = np.zeros(records, dtype=int)

    # clustered by probe_name
    probe_mutation = probe_mutation.groupby('probe_name').agg({
        'patient_name': lambda x: set(x),
        'mutation': lambda x: set(x)
    }).reset_index()

    return patient_mutation, patient_mut_count, probe_mutation

def add_the_first(probe_mutation, patient_mutation, patient_mut_count):
    probe_mutation['patient_num'] = probe_mutation['patient_name'].apply(len)
    probe_mutation['mut_num'] = probe_mutation['mutation'].apply(len)

    # Greed method, get the top 20000, save the time
    sorted_probe_mutation = (probe_mutation.sort_values(by=['patient_num', 'mut_num'], ascending=[False, False])
                             .head(20000))
    result_dict = defaultdict(set)

    for index, row in sorted_probe_mutation.iterrows():
        for mutation in row['mutation']:
            patient = mutation.split(":")[0]
            result_dict[patient].add(mutation)
    
    for key, value in result_dict.items():
        # obtain numpy_index
        numpy_index = patient_mutation.loc[patient_mutation["patient_name"] == key, "numpy_index"]
        index = numpy_index.iloc[0]
        patient_mutation.at[index, "covered_mut"] = frozenset(value)  
        # refresh the patient count
        patient_mut_count[index] = len(value)

    top_20000_probes = sorted_probe_mutation['probe_name']
    remaining_probe_mutation = probe_mutation[~probe_mutation['probe_name'].isin(top_20000_probes)]

    return remaining_probe_mutation, patient_mutation, patient_mut_count


def add_the_best(sub_df, patient_mutation, patient_mut_count):
    improvement = []
    improvement1 = []

    for index, row in sub_df.iterrows():
        patient_mutation_cp = patient_mutation.copy()
        patient_mut_count_cp = patient_mut_count.copy()
        patient_mut_count_cp1 = patient_mut_count.copy()
        
        for mutation in row['mutation']:
            patient = mutation.split(":")[0]
            present = patient_mutation_cp.loc[patient_mutation_cp['patient_name'] == patient, 'covered_mut'].values
            
            if present.size > 0:
                present = present[0]
                if mutation not in present:
                    patient_num = patient_mutation_cp.loc[patient_mutation_cp['patient_name'] == patient, 'numpy_index'].values[0]
                    if len(present) < 5:
                        patient_mut_count_cp[patient_num] += 1
                    else:
                        pass
                    patient_mut_count_cp1[patient_num] += 1
        
        improvement.append(np.sum(patient_mut_count_cp - patient_mut_count))
        improvement1.append(np.sum(patient_mut_count_cp1 - patient_mut_count))

    sub_df["improvement"] = improvement
    sub_df["improvement1"] = improvement1
    top_probe = sub_df.sort_values(by=["improvement", "improvement1", "patient_num", "mut_num"], ascending=[False, False, False, False]).head(1).iloc[0]
    
    return top_probe

def find_best_probe_multiprocess(probe_mutation, patient_mutation, patient_mut_count, num_processes=None):
    if num_processes is None:
        num_processes = 15

    sub_dfs = np.array_split(probe_mutation, num_processes)
    
    with Pool(processes=num_processes) as pool:
        results = pool.starmap(add_the_best, [(sub_df, patient_mutation, patient_mut_count) for sub_df in sub_dfs])
    
    best_probe = max(results, key=lambda x: (x["improvement"], x["improvement1"], x["patient_num"], x["mut_num"]))
    
    probe_mutation = probe_mutation[probe_mutation["probe_name"] != best_probe["probe_name"]]

    for mutation in best_probe["mutation"]:
        patient = mutation.split(":")[0]
        present = patient_mutation.loc[patient_mutation['patient_name'] == patient, 'covered_mut'].values
        
        if present.size > 0:
            present = present[0]
            if mutation not in present:
                numpy_index = patient_mutation.loc[patient_mutation["patient_name"] == patient, "numpy_index"]
                
                if not numpy_index.empty:
                    index = numpy_index.iloc[0]
                    patient_mut_count[index] += 1
                    updated_present = set(present)
                    updated_present.add(mutation)
                    patient_mutation.at[index, "covered_mut"] = frozenset(updated_present)  

    return probe_mutation, patient_mutation, patient_mut_count



def main():
    patient_mutation, patient_mut_count, probe_mutation = read_mutation("./mutation_probes.tsv")
    probe_mutation, patient_mutation, patient_mut_count = add_the_first(probe_mutation, patient_mutation,
                                                                        patient_mut_count)

    stop_ratio = np.count_nonzero(patient_mut_count >= 5) / len(patient_mut_count)
    stop_ratio1 = np.count_nonzero(patient_mut_count >= 6) / len(patient_mut_count)
    max_iter = 15000
    times = 0
    while stop_ratio <= 0.99 and stop_ratio1 <= 0.9:
        if times >= max_iter:
            break
        probe_mutation, patient_mutation, patient_mut_count = find_best_probe_multiprocess(probe_mutation, patient_mutation,
                                                                           patient_mut_count)
        stop_ratio = np.count_nonzero(patient_mut_count >= 5) / len(patient_mut_count)
        stop_ratio1 = np.count_nonzero(patient_mut_count >= 6) / len(patient_mut_count)
        print(stop_ratio, stop_ratio1, sep="\t")
        times += 1

    probe_mutation.drop(columns=["patient_name", "mutation"]).to_csv("./remained_probed.tsv", index=False, sep="\t")


if __name__ == "__main__":
    main()

```
