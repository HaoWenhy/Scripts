#! /usr/bin/python
import sys
import re
import concurrent.futures
import argparse
import pandas as pd

# 解析命令行参数
parser = argparse.ArgumentParser(description='提取内含子信息并计算统计信息')
parser.add_argument('fasta_file', type=str, help='FASTA格式的基因组文件')
parser.add_argument('gtf_file', type=str, help='GTF格式的注释文件')
parser.add_argument('--threads', type=int, default=4, help='线程数量 (默认为4)')
parser.add_argument('--output', type=str, default='intron_stats.csv', help='输出文件名 (默认为intron_stats.csv)')
args = parser.parse_args()

# 从命令行参数中获取输入文件名和输出文件名
fasta_file = args.fasta_file
gtf_file = args.gtf_file
num_threads = args.threads
output_file = args.output

# 打开FASTA文件并读取基因组序列
def read_fasta(fasta_file):
    genome = {}
    current_sequence = None
    with open(fasta_file, 'r') as fasta:
        for line in fasta:
            line = line.strip()
            if line.startswith('>'):
                current_sequence = line[1:]
                genome[current_sequence] = []
            else:
                genome[current_sequence].append(line)
    return genome

genome = read_fasta(fasta_file)

# 使用集合来存储已知的转录本ID
known_transcripts = set()

# 定义函数以提取内含子位置和计算相关统计信息
def extract_intron_info(gtf_file, genome, transcript_id_set):
    intron_lengths = []
    intron_sequences = []

    with open(gtf_file, 'r') as gtf:
        current_transcript_id = None
        current_sequence = None
        for line in gtf:
            if not line.startswith('#'):
                fields = line.strip().split('\t')
                if len(fields) >= 9:  # 确保至少有9个字段
                    feature_type = fields[2]
                    if feature_type == 'exon':
                        attributes = fields[8]
                        match = re.search(r'transcript_id "([^"]+)"', attributes)
                        if match:
                            transcript_id = match.group(1)
                            if current_transcript_id != transcript_id:
                                if current_sequence is not None:
                                    intron_sequences.append(''.join(current_sequence))
                                    intron_lengths.append(len(current_sequence))
                                current_transcript_id = transcript_id
                                current_sequence = genome.get(transcript_id, [])
                                transcript_id_set.add(transcript_id)
                            else:
                                current_sequence.extend(genome.get(transcript_id, []))
        
        # 处理最后一个内含子
        if current_sequence is not None:
            intron_sequences.append(''.join(current_sequence))
            intron_lengths.append(len(current_sequence))

    return intron_lengths, intron_sequences

# 提取内含子信息并计算相关统计信息
intron_lengths = []
intron_sequences = []

# 创建线程池
with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as executor:
    future_to_transcript_id_set = {}
    for _ in range(num_threads):
        transcript_id_set = set()
        future = executor.submit(extract_intron_info, gtf_file, genome, transcript_id_set)
        future_to_transcript_id_set[future] = transcript_id_set

    for future in concurrent.futures.as_completed(future_to_transcript_id_set):
        intron_lengths_part, intron_sequences_part = future.result()
        intron_lengths.extend(intron_lengths_part)
        intron_sequences.extend(intron_sequences_part)
        known_transcripts.update(future_to_transcript_id_set[future])

# 写入内含子序列到sre.introns.fasta文件
with open('sre.introns.fasta', 'w') as intron_fasta:
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

# 创建统计结果的DataFrame
data = {
    'Total Intron Length (bp)': [total_intron_length],
    'Number of Introns': [num_introns],
    'Average Intron Length (bp)': [average_intron_length],
    'Median Intron Length (bp)': [median_intron_length],
    'Intron GC Content (%)': [gc_content]
}

# 将统计结果写入
