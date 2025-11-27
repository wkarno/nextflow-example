#!/bin/bash
# Script to generate toy data for testing the FASTQ pipeline

set -e

# Configuration
REF_LENGTH=5000
NUM_READS=100
READ_LENGTH=150
INSERT_SIZE=300
SAMPLE_ID="test_sample"

echo "======================================"
echo "Generating Test Data"
echo "======================================"
echo "Reference length: ${REF_LENGTH} bp"
echo "Number of reads: ${NUM_READS}"
echo "Read length: ${READ_LENGTH} bp"
echo "Insert size: ${INSERT_SIZE} bp"
echo ""

# Generate a simple reference genome
echo "Creating reference genome..."
python3 - <<'EOF'
import random
import sys

# Set seed for reproducibility
random.seed(42)

# Generate random DNA sequence
bases = ['A', 'C', 'G', 'T']
length = 5000

# Create a more realistic reference with some structure
sequence = []
for i in range(length):
    # Add some GC-rich and AT-rich regions
    if i % 1000 < 500:
        # GC-rich region
        sequence.append(random.choice(['G', 'C', 'A', 'T'] if random.random() < 0.6 else ['A', 'T']))
    else:
        # AT-rich region
        sequence.append(random.choice(['A', 'T', 'G', 'C'] if random.random() < 0.6 else ['G', 'C']))

# Write FASTA
print(">test_reference synthetic 5kb genome for testing")
# Print in 80 character lines
seq_str = ''.join(sequence)
for i in range(0, len(seq_str), 80):
    print(seq_str[i:i+80])

EOF
python3 - > reference.fa <<'PYEOF'
import random
import sys

random.seed(42)
bases = ['A', 'C', 'G', 'T']
length = 5000

sequence = []
for i in range(length):
    if i % 1000 < 500:
        sequence.append(random.choice(['G', 'C', 'A', 'T'] if random.random() < 0.6 else ['A', 'T']))
    else:
        sequence.append(random.choice(['A', 'T', 'G', 'C'] if random.random() < 0.6 else ['G', 'C']))

print(">test_reference synthetic 5kb genome for testing")
seq_str = ''.join(sequence)
for i in range(0, len(seq_str), 80):
    print(seq_str[i:i+80])
PYEOF

echo "Reference genome created: reference.fa"

# Generate simulated paired-end reads
echo "Generating simulated paired-end reads..."
python3 - <<'EOF'
import random
import gzip

# Set seed for reproducibility
random.seed(42)

# Read the reference
with open('reference.fa', 'r') as f:
    lines = f.readlines()
    ref_seq = ''.join([line.strip() for line in lines if not line.startswith('>')])

ref_len = len(ref_seq)
num_reads = 100
read_len = 150
insert_size = 300

def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join(complement.get(base, 'N') for base in reversed(seq))

def add_errors(seq, error_rate=0.01):
    """Add random sequencing errors"""
    bases = ['A', 'C', 'G', 'T']
    seq_list = list(seq)
    for i in range(len(seq_list)):
        if random.random() < error_rate:
            seq_list[i] = random.choice(bases)
    return ''.join(seq_list)

def generate_quality_string(length, mean_qual=35):
    """Generate realistic quality scores"""
    quals = []
    for i in range(length):
        # Quality decreases toward end of read
        pos_factor = 1.0 - (i / length) * 0.3
        qual = int(mean_qual * pos_factor + random.randint(-5, 5))
        qual = max(20, min(40, qual))  # Clamp between 20 and 40
        quals.append(chr(qual + 33))  # Phred+33 encoding
    return ''.join(quals)

# Generate reads
reads_r1 = []
reads_r2 = []

for i in range(num_reads):
    # Random position in reference (ensure we can get both reads)
    max_start = ref_len - insert_size
    if max_start < 0:
        max_start = 0

    start_pos = random.randint(0, max_start)

    # Extract forward read
    r1_seq = ref_seq[start_pos:start_pos + read_len]

    # Extract reverse read (from the other end of the insert)
    r2_start = start_pos + insert_size - read_len
    if r2_start < 0:
        r2_start = 0
    if r2_start + read_len > ref_len:
        r2_start = ref_len - read_len

    r2_seq = ref_seq[r2_start:r2_start + read_len]
    r2_seq = reverse_complement(r2_seq)

    # Add sequencing errors
    r1_seq = add_errors(r1_seq)
    r2_seq = add_errors(r2_seq)

    # Generate quality scores
    r1_qual = generate_quality_string(len(r1_seq))
    r2_qual = generate_quality_string(len(r2_seq))

    # Store reads
    reads_r1.append((f"@read_{i+1}/1", r1_seq, "+", r1_qual))
    reads_r2.append((f"@read_{i+1}/2", r2_seq, "+", r2_qual))

# Write R1
with gzip.open('test_sample_R1.fastq.gz', 'wt') as f:
    for header, seq, plus, qual in reads_r1:
        f.write(f"{header}\n{seq}\n{plus}\n{qual}\n")

# Write R2
with gzip.open('test_sample_R2.fastq.gz', 'wt') as f:
    for header, seq, plus, qual in reads_r2:
        f.write(f"{header}\n{seq}\n{plus}\n{qual}\n")

print(f"Generated {num_reads} paired-end reads")
print(f"R1: test_sample_R1.fastq.gz")
print(f"R2: test_sample_R2.fastq.gz")
EOF

echo ""
echo "======================================"
echo "Test Data Generation Complete!"
echo "======================================"
echo ""
echo "Files created:"
ls -lh reference.fa test_sample_R*.fastq.gz 2>/dev/null || echo "Check for errors above"
echo ""
echo "To use this data, create a samplesheet:"
echo "sampleID,fastq1,fastq2"
echo "test_sample,$(pwd)/test_sample_R1.fastq.gz,$(pwd)/test_sample_R2.fastq.gz"
