#! /usr/bin/python
import sys
import re

# 检查命令行参数
if len(sys.argv) != 3:
    print("Usage: python intron_stats.py genomic.fasta genes.gtf")
    sys.exit(1)

# 从命令行参数中获取输入文件名
fasta_file = sys.argv[1]
gtf_file = sys.argv[2]

# 打开FASTA文件并读取基因组序列
def read_fasta(fasta_file):
    genome = {}
    current_sequence = None
    with open(fasta_file, 'r') as fasta:
        for line in fasta:
            line = line.strip()
            if line.startswith('>'):
                current_sequence = line[1:]
                genome[current_sequence] = ""
            else:
                genome[current_sequence] += line
    return genome

genome = read_fasta(fasta_file)

# 定义函数以提取内含子位置和计算相关统计信息
def extract_intron_info(gtf_file, genome):
    intron_lengths = []
    intron_sequences = []

    with open(gtf_file, 'r') as gtf:
        current_transcript_id = None
        current_sequence = None
        for line in gtf:
            if not line.startswith('#'):
                fields = line.strip().split('\t')
                feature_type = fields[2]
                if feature_type == 'exon':
                    attributes = fields[8]
                    match = re.search(r'transcript_id "([^"]+)"', attributes)
                    if match:
                        transcript_id = match.group(1)
                        start = int(fields[3]) - 1
                        end = int(fields[4])
                        if current_transcript_id != transcript_id:
                            if current_sequence is not None:
                                intron_sequences.append(current_sequence)
                                intron_lengths.append(len(current_sequence))
                            current_transcript_id = transcript_id
                            current_sequence = genome.get(transcript_id, "")[start:end]
                        else:
                            current_sequence += genome.get(transcript_id, "")[start:end]
        
        # 处理最后一个内含子
        if current_sequence is not None:
            intron_sequences.append(current_sequence)
            intron_lengths.append(len(current_sequence))

    return intron_lengths, intron_sequences

# 提取内含子信息并计算相关统计信息
intron_lengths, intron_sequences = extract_intron_info(gtf_file, genome)

# 写入内含子序列到sre.introns.fasta文件
with open('introns.fasta', 'w') as intron_fasta:
    for i, sequence in enumerate(intron_sequences):
        intron_fasta.write(f">Intron_{i+1}\n{sequence}\n")

# 计算内含子总长度、内含子个数、平均长度、中位数和GC含量
total_intron_length = sum(intron_lengths)
num_introns = len(intron_lengths)
average_intron_length = total_intron_length / num_introns if num_introns > 0 else 0
sorted_intron_lengths = sorted(intron_lengths)
median_intron_length = (
    sorted_intron_lengths[len(sorted_intron_lengths) // 2]
    + sorted_intron_lengths[(len(sorted_intron_lengths) - 1) // 2]
) / 2 if num_introns > 0 else 0
total_gc_count = sum(sequence.count('G') + sequence.count('C') for sequence in intron_sequences)
total_base_count = sum(len(sequence) for sequence in intron_sequences)
gc_content = (total_gc_count / total_base_count) * 100 if total_base_count > 0 else 0

# 打印统计信息
print(f"内含子总长度: {total_intron_length} bp")
print(f"内含子个数: {num_introns}")
print(f"内含子平均长度: {average_intron_length:.2f} bp")
print(f"内含子中位数长度: {median_intron_length:.2f} bp")
print(f"内含子GC含量: {gc_content:.2f}%")
