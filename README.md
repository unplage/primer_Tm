# primer_Tm
primer Tm value calculate,based in python,simulate Snap Gene

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
