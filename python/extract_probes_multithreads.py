import sys
import numpy as np
import pandas as pd
from collections import defaultdict
from multiprocessing import Pool, cpu_count

def read_mutation(file):
    probe_mutation = pd.read_csv(file, names=["probe_name", "mutation"], sep="\t", header=None)
    
    # 从 mutation 列中提取患者名称
    probe_mutation["patient_name"] = probe_mutation["mutation"].str.split(":").str[0]
    
    # 获取唯一的患者名称列表
    patients = list(probe_mutation["patient_name"].unique())
    records = len(patients)
    numbers = np.arange(records).tolist()
    
    # 创建 patient_mutation DataFrame
    patient_mutation = pd.DataFrame({"patient_name": patients, "numpy_index": numbers})
    patient_mutation["covered_mut"] = [frozenset() for _ in range(records)]
    patient_mut_count = np.zeros(records, dtype=int)

    # 按照 probe_name 分组并聚合
    probe_mutation = probe_mutation.groupby('probe_name').agg({
        'patient_name': lambda x: set(x),
        'mutation': lambda x: set(x)
    }).reset_index()

    return patient_mutation, patient_mut_count, probe_mutation

def add_the_first(probe_mutation, patient_mutation, patient_mut_count):
    probe_mutation['patient_num'] = probe_mutation['patient_name'].apply(len)
    probe_mutation['mut_num'] = probe_mutation['mutation'].apply(len)

    # 按患者数量和突变数量排序，并取前20000个
    sorted_probe_mutation = (probe_mutation.sort_values(by=['patient_num', 'mut_num'], ascending=[False, False])
                             .head(20000))
    result_dict = defaultdict(set)

    for index, row in sorted_probe_mutation.iterrows():
        for mutation in row['mutation']:
            patient = mutation.split(":")[0]
            result_dict[patient].add(mutation)
    
    for key, value in result_dict.items():
        # 找到患者的 numpy_index
        numpy_index = patient_mutation.loc[patient_mutation["patient_name"] == key, "numpy_index"]
        index = numpy_index.iloc[0]
        # 赋值 frozenset
        patient_mutation.at[index, "covered_mut"] = frozenset(value)  
        # 更新患者突变计数
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
    # 如果没有指定进程数，使用可用 CPU 核心数
    if num_processes is None:
        num_processes = 15

    # 将 probe_mutation 切割成子集
    sub_dfs = np.array_split(probe_mutation, num_processes)
    
    # 使用进程池执行 add_the_best
    with Pool(processes=num_processes) as pool:
        results = pool.starmap(add_the_best, [(sub_df, patient_mutation, patient_mut_count) for sub_df in sub_dfs])
    
    # 找到全局最优的 probe
    best_probe = max(results, key=lambda x: (x["improvement"], x["improvement1"], x["patient_num"], x["mut_num"]))
    
    # 从 probe_mutation 中删除最优的 probe
    probe_mutation = probe_mutation[probe_mutation["probe_name"] != best_probe["probe_name"]]

    # 更新 patient_mutation 和 patient_mut_count
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
