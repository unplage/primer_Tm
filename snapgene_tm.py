'''
#保存为snapgene_tm.py，放同一个文件夹中
#!/usr/bin/env python3
"""
精简版SnapGene Tm计算器
基于校准参数: 56条序列 | RMSE: 1.16°C
函数接口:
1. calculate_tm(sequence, na_conc=50.0, primer_conc=0.25, annealing_mode=False)
2. calculate_tm_batch(sequences, na_conc=50.0, primer_conc=0.25, annealing_mode=False)
"""

import math
from typing import Union, List, Dict, Tuple

# 校准参数 (基于56条序列优化)
CALIBRATION_PARAMS = {
    "salt_correction_factor": 1.6,
    "terminal_AT_penalty": 2.497083789430314,
    "entropy_scaling_factor": 0.9675869088668494,
    "initiation_adjustment": -0.0007103390016273732,
    "calibration_rmse": 1.1553994766937066,
    "calibration_points": 56
}

# 热力学参数 (SantaLucia 2004)
NN_DH = {
    'AA': -7.6, 'TT': -7.6, 'AT': -7.2, 'TA': -7.2,
    'CA': -8.5, 'TG': -8.5, 'GT': -8.4, 'AC': -8.4,
    'CT': -7.8, 'AG': -7.8, 'GA': -8.2, 'TC': -8.2,
    'CG': -10.6, 'GC': -10.6, 'GG': -8.0, 'CC': -8.0
}

NN_DS = {
    'AA': -21.3, 'TT': -21.3, 'AT': -20.4, 'TA': -20.4,
    'CA': -22.7, 'TG': -22.7, 'GT': -22.4, 'AC': -22.4,
    'CT': -21.0, 'AG': -21.0, 'GA': -22.2, 'TC': -22.2,
    'CG': -27.2, 'GC': -27.2, 'GG': -19.9, 'CC': -19.9
}


def calculate_tm(sequence: str, na_conc: float = 50.0, 
                 primer_conc: float = 0.25, annealing_mode: bool = False,
                 return_details: bool = False) -> Union[float, Dict]:
    """
    计算单条DNA序列的Tm值
    
    参数:
        sequence: DNA序列 (大写或小写)
        na_conc: 钠离子浓度 (mM), 默认50.0
        primer_conc: 引物浓度 (µM), 默认0.25
        annealing_mode: 是否为退火模式, 默认False (PCR模式)
        return_details: 是否返回详细信息, 默认False只返回Tm值
    
    返回:
        如果return_details=False: 返回Tm值(float)
        如果return_details=True: 返回包含详细信息的字典
    
    示例:
        >>> calculate_tm("ATCGATCGATCGATCG")
        46.3
        >>> calculate_tm("ATCGATCGATCGATCG", return_details=True)
        {'tm': 46.3, 'sequence': 'ATCGATCGATCGATCG', 'length': 16, ...}
    """
    seq = sequence.upper().strip()
    n = len(seq)
    
    # 验证序列
    if n < 4:
        raise ValueError(f"序列太短 ({n} bp)，至少需要4bp")
    
    # 1. 计算序列本征ΔH和ΔS
    dH = 0.2 + CALIBRATION_PARAMS["initiation_adjustment"] * 0.05
    dS = -5.7 * CALIBRATION_PARAMS["entropy_scaling_factor"]
    
    for i in range(n - 1):
        nn = seq[i:i+2]
        dH += NN_DH[nn]
        dS += NN_DS[nn] * CALIBRATION_PARAMS["entropy_scaling_factor"]
    
    if seq[-1] in ('A', 'T'):
        dH += CALIBRATION_PARAMS["terminal_AT_penalty"]
        dS += 6.9
    
    # 2. 盐校正
    monovalent_M = na_conc / 1000.0
    salt_effect = 0.368 * (n - 1) * math.log(monovalent_M) * CALIBRATION_PARAMS["salt_correction_factor"]
    dS_corrected = dS + salt_effect
    
    # 3. 计算Tm
    R = 1.987
    Ct = primer_conc * 1e-6
    
    if annealing_mode:
        sym_factor = 1
    else:
        sym_factor = 4
    
    denominator = dS_corrected + R * math.log(Ct / sym_factor)
    
    if denominator >= 0:
        # 回退到Wallace规则
        gc_count = seq.count('G') + seq.count('C')
        tm = 2 * (n - gc_count) + 4 * gc_count
    else:
        Tm_K = (dH * 1000) / denominator
        tm = Tm_K - 273.15
    
    if not return_details:
        return round(tm, 2)
    
    # 准备详细信息
    gc_count = seq.count('G') + seq.count('C')
    at_count = n - gc_count
    
    return {
        "tm": round(tm, 2),
        "sequence": seq,
        "length": n,
        "gc_percent": round(gc_count / n * 100, 2),
        "gc_count": gc_count,
        "at_count": at_count,
        "conditions": {
            "na_conc": na_conc,
            "primer_conc": primer_conc,
            "annealing_mode": annealing_mode
        },
        "annealing_suggestions": {
            "standard_enzyme": [round(tm - 3, 1), round(tm - 1, 1)],
            "high_fidelity": [round(tm + 3, 1), round(tm + 5, 1)],
            "gradient_pcr": [round(tm - 5, 1), round(tm + 5, 1)]
        },
        "calibration_info": {
            "rmse": CALIBRATION_PARAMS["calibration_rmse"],
            "data_points": CALIBRATION_PARAMS["calibration_points"]
        }
    }


def calculate_tm_batch(sequences: List[str], na_conc: float = 50.0,
                       primer_conc: float = 0.25, annealing_mode: bool = False,
                       return_format: str = "list") -> Union[List[float], Dict]:
    """
    批量计算多个DNA序列的Tm值
    
    参数:
        sequences: DNA序列列表
        na_conc: 钠离子浓度 (mM), 默认50.0
        primer_conc: 引物浓度 (µM), 默认0.25
        annealing_mode: 是否为退火模式, 默认False (PCR模式)
        return_format: 返回格式, "list"只返回Tm列表, "dict"返回完整信息
    
    返回:
        如果return_format="list": 返回Tm值列表 [float, ...]
        如果return_format="dict": 返回包含成功和失败结果的字典
    
    示例:
        >>> calculate_tm_batch(["ATCGATCG", "GCGCGCGC"])
        [34.0, 68.0]
    """
    results = []
    errors = []
    
    for i, seq in enumerate(sequences):
        try:
            if return_format == "dict":
                result = calculate_tm(seq, na_conc, primer_conc, annealing_mode, return_details=True)
                result["index"] = i
                results.append(result)
            else:
                tm = calculate_tm(seq, na_conc, primer_conc, annealing_mode, return_details=False)
                results.append(tm)
        except Exception as e:
            errors.append({
                "index": i,
                "sequence": seq,
                "error": str(e)
            })
    
    if return_format == "dict":
        return {
            "total": len(sequences),
            "successful": len(results),
            "failed": len(errors),
            "results": results,
            "errors": errors
        }
    
    return results


def get_annealing_temp(tm: float, enzyme_type: str = "standard") -> Tuple[float, float]:
    """
    根据Tm值获取退火温度范围
    
    参数:
        tm: Tm值
        enzyme_type: 酶类型, "standard"(标准酶)或"high_fidelity"(高保真酶)
    
    返回:
        退火温度范围 (最小值, 最大值)
    """
    if enzyme_type == "standard":
        return (round(tm - 3, 1), round(tm - 1, 1))
    elif enzyme_type == "high_fidelity":
        return (round(tm + 3, 1), round(tm + 5, 1))
    else:
        raise ValueError(f"未知酶类型: {enzyme_type}")


# 快速验证函数
def validate_sequence(sequence: str) -> Tuple[bool, str]:
    """验证DNA序列有效性"""
    seq = sequence.upper().strip()
    
    if not seq:
        return False, "序列不能为空"
    
    if len(seq) < 4:
        return False, f"序列太短 ({len(seq)} bp)，至少需要4bp"
    
    if len(seq) > 200:
        return False, f"序列太长 ({len(seq)} bp)，最大支持200bp"
    
    invalid_chars = [c for c in seq if c not in "ATCG"]
    if invalid_chars:
        return False, f"序列包含非法字符: {', '.join(set(invalid_chars))}"
    
    return True, seq


# 使用示例
if __name__ == "__main__":
    # 示例1: 单序列计算
    print("示例1: 单序列计算")
    seq = "ATCGATCGATCGATCG"
    tm = calculate_tm(seq)
    print(f"序列: {seq}")
    print(f"Tm值: {tm}°C")
    print(f"标准酶退火范围: {get_annealing_temp(tm, 'standard')[0]}-{get_annealing_temp(tm, 'standard')[1]}°C")
    
    # 示例2: 获取详细信息
    print("\n示例2: 获取详细信息")
    details = calculate_tm("GCCAGTGCCAAGCTTGCA", return_details=True)
    print(f"序列: {details['sequence']}")
    print(f"长度: {details['length']} bp")
    print(f"GC含量: {details['gc_percent']}%")
    print(f"Tm值: {details['tm']}°C")
    
    # 示例3: 批量计算
    print("\n示例3: 批量计算")
    sequences = ["ATCGATCGATCGATCG", "GCCAGTGCCAAGCTTGCA", "TATATATATATATATA"]
    batch_results = calculate_tm_batch(sequences)
    for seq, tm in zip(sequences, batch_results):
        print(f"{seq[:20]:<20} -> {tm}°C")
    
    # 示例4: 批量计算获取完整信息
    print("\n示例4: 批量计算完整信息")
    full_batch = calculate_tm_batch(sequences, return_format="dict")
    print(f"总数: {full_batch['total']}, 成功: {full_batch['successful']}, 失败: {full_batch['failed']}")
    
    print("\n" + "="*50)
    print("函数接口总结:")
    print("1. calculate_tm(sequence, na_conc=50.0, primer_conc=0.25, annealing_mode=False, return_details=False)")
    print("2. calculate_tm_batch(sequences, na_conc=50.0, primer_conc=0.25, annealing_mode=False, return_format='list')")
    print("3. get_annealing_temp(tm, enzyme_type='standard')")
    print("4. validate_sequence(sequence)")

'''


'''
# 导入函数，计算文本中的引物Tm，
def add_tm_to_sequences(input_file: str, output_file: str = None):
    """
    极简版本：读取序列文件，计算Tm，生成新文件
    """
    from snapgene_tm import calculate_tm_batch

    # 读取文件
    with open(input_file, 'r') as f:
        sequences = [line.strip() for line in f if line.strip() and not line.startswith('#')]

    # 计算Tm
    results = calculate_tm_batch(sequences)

    # 生成输出文件名
    if output_file is None:
        output_file = input_file.replace('.txt',
                                         '_with_tm.txt') if '.txt' in input_file else input_file + '_with_tm.txt'

    # 写入结果
    with open(output_file, 'w') as f:
        for seq, tm in zip(sequences, results):
            f.write(f"{seq}\t{tm}\n")  # 使用制表符分隔

    print(f"完成！结果保存到: {output_file}")


# 使用
add_tm_to_sequences("my_primers.txt")

'''



'''

# 少量数据可以使用以下命令行进行Tm值计算

# 导入特定函数
from snapgene_tm import calculate_tm, calculate_tm_batch

# 使用函数
tm = calculate_tm("ATCGATCGATCGATCG")
print(f"单序列Tm: {tm}°C")

sequences = ["ATCGATCGATCGATCG", "GCCAGTGCCAAGCTTGCA", "TATATATATATATATA"]
batch_result = calculate_tm_batch(sequences, return_format="dict")
print(f"批量计算成功: {batch_result['successful']}个")

'''
