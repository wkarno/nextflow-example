#!/usr/bin/env python3
"""
FASTQ Samplesheet Creator

This script recursively searches a directory for FASTQ files (paired or unpaired),
extracts sample names, and generates a CSV samplesheet compatible with the
Nextflow FASTQ processing pipeline.

The output samplesheet is named: <directory_name>_SampleSheet.csv

Usage:
    python create_samplesheet.py <input_directory> [options]

Example:
    python create_samplesheet.py /path/to/fastq/dir
    python create_samplesheet.py /path/to/fastq/dir --output custom_output.csv
"""

import argparse
import csv
import logging
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

# Common FASTQ file extensions
FASTQ_EXTENSIONS = {'.fastq', '.fq', '.fastq.gz', '.fq.gz'}

# Patterns for paired-end reads
PAIRED_PATTERNS = [
    # Illumina-style: _R1_, _R2_
    (re.compile(r'(.+?)_R1(_.*)?$'), re.compile(r'(.+?)_R2(_.*)?$')),
    # Simple style: _1, _2
    (re.compile(r'(.+?)_1(_.*)?$'), re.compile(r'(.+?)_2(_.*)?$')),
]


def is_fastq_file(filepath: Path) -> bool:
    """
    Check if a file is a FASTQ file based on extension.

    Args:
        filepath: Path to the file

    Returns:
        True if the file has a FASTQ extension
    """
    # Check for compound extensions like .fastq.gz
    for ext in FASTQ_EXTENSIONS:
        if filepath.name.endswith(ext):
            return True
    return False


def extract_sample_name(filepath: Path) -> str:
    """
    Extract sample name from FASTQ filename by removing standard Illumina suffixes
    and read pair indicators.

    Examples:
        sample1_R1.fastq.gz -> sample1
        sample1_S01_L001_R1_001.fastq.gz -> sample1
        sample1_1.fq -> sample1

    Args:
        filepath: Path to the FASTQ file

    Returns:
        Extracted sample name
    """
    # Remove all known FASTQ extensions
    name = filepath.name
    for ext in sorted(FASTQ_EXTENSIONS, key=len, reverse=True):
        if name.endswith(ext):
            name = name[:-len(ext)]
            break

    # Try to match paired-end patterns and extract base name
    for r1_pattern, r2_pattern in PAIRED_PATTERNS:
        r1_match = r1_pattern.match(name)
        r2_match = r2_pattern.match(name)
        if r1_match:
            return r1_match.group(1)
        if r2_match:
            return r2_match.group(1)

    # If no paired pattern matched, return the name as-is
    # This handles unpaired reads or non-standard naming
    return name


def find_fastq_files(directory: Path) -> List[Path]:
    """
    Recursively find all FASTQ files in a directory.

    Args:
        directory: Root directory to search

    Returns:
        List of Path objects for all FASTQ files found
    """
    fastq_files = []

    for filepath in directory.rglob('*'):
        if filepath.is_file() and is_fastq_file(filepath):
            fastq_files.append(filepath)

    logger.info(f"Found {len(fastq_files)} FASTQ files in {directory}")
    return fastq_files


def pair_fastq_files(fastq_files: List[Path]) -> Dict[str, Dict[str, Optional[Path]]]:
    """
    Group FASTQ files by sample name and detect paired vs unpaired reads.

    Args:
        fastq_files: List of FASTQ file paths

    Returns:
        Dictionary mapping sample names to their R1 and R2 file paths
        Format: {sample_name: {'R1': path, 'R2': path or None}}
    """
    samples = {}

    for filepath in fastq_files:
        sample_name = extract_sample_name(filepath)

        # Initialize sample entry if not exists
        if sample_name not in samples:
            samples[sample_name] = {'R1': None, 'R2': None}

        # Determine if this is R1, R2, or unpaired
        name = filepath.name
        is_r1 = False
        is_r2 = False

        for r1_pattern, r2_pattern in PAIRED_PATTERNS:
            # Remove extension for pattern matching
            name_no_ext = filepath.name
            for ext in sorted(FASTQ_EXTENSIONS, key=len, reverse=True):
                if name_no_ext.endswith(ext):
                    name_no_ext = name_no_ext[:-len(ext)]
                    break

            if r1_pattern.match(name_no_ext):
                is_r1 = True
                break
            elif r2_pattern.match(name_no_ext):
                is_r2 = True
                break

        # Assign to appropriate slot
        if is_r1:
            if samples[sample_name]['R1'] is not None:
                logger.warning(f"Multiple R1 files found for sample '{sample_name}': "
                             f"{samples[sample_name]['R1'].name} and {filepath.name}")
            samples[sample_name]['R1'] = filepath
        elif is_r2:
            if samples[sample_name]['R2'] is not None:
                logger.warning(f"Multiple R2 files found for sample '{sample_name}': "
                             f"{samples[sample_name]['R2'].name} and {filepath.name}")
            samples[sample_name]['R2'] = filepath
        else:
            # Unpaired or couldn't determine pairing - put in R1 slot
            if samples[sample_name]['R1'] is not None:
                logger.warning(f"Multiple files found for sample '{sample_name}': "
                             f"{samples[sample_name]['R1'].name} and {filepath.name}")
            samples[sample_name]['R1'] = filepath

    # Log pairing status
    paired_count = sum(1 for s in samples.values() if s['R1'] and s['R2'])
    unpaired_count = sum(1 for s in samples.values() if s['R1'] and not s['R2'])
    incomplete_count = sum(1 for s in samples.values() if s['R2'] and not s['R1'])

    logger.info(f"Detected {paired_count} paired samples")
    logger.info(f"Detected {unpaired_count} unpaired samples")

    if incomplete_count > 0:
        logger.warning(f"Found {incomplete_count} samples with R2 but no R1 - "
                      "these will be placed in the fastq1 column")
        # Fix incomplete pairs by moving R2 to R1 slot
        for sample_name, files in samples.items():
            if files['R2'] and not files['R1']:
                files['R1'] = files['R2']
                files['R2'] = None

    return samples


def write_samplesheet(
    samples: Dict[str, Dict[str, Optional[Path]]],
    output_path: Path,
    use_absolute_paths: bool = True
) -> None:
    """
    Write samples to a CSV samplesheet.

    Args:
        samples: Dictionary of samples and their FASTQ files
        output_path: Path to write the samplesheet CSV
        use_absolute_paths: If True, use absolute paths in the CSV
    """
    # Sort samples by name for consistent output
    sorted_samples = sorted(samples.items())

    with open(output_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)

        # Write header
        writer.writerow(['sampleID', 'fastq1', 'fastq2'])

        # Write sample rows
        for sample_name, files in sorted_samples:
            if files['R1'] is None:
                logger.warning(f"Skipping sample '{sample_name}' - no FASTQ files found")
                continue

            fastq1 = files['R1'].resolve() if use_absolute_paths else files['R1']
            fastq2 = files['R2'].resolve() if (files['R2'] and use_absolute_paths) else files['R2']

            # Write row (empty string for missing R2)
            writer.writerow([
                sample_name,
                str(fastq1),
                str(fastq2) if fastq2 else ''
            ])

    logger.info(f"Samplesheet written to: {output_path}")
    logger.info(f"Total samples: {len([s for s in sorted_samples if s[1]['R1']])}")


def main():
    """Main entry point for the script."""
    parser = argparse.ArgumentParser(
        description='Create a samplesheet CSV from a directory of FASTQ files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    parser.add_argument(
        'input_dir',
        type=Path,
        help='Directory to search for FASTQ files (searched recursively)'
    )

    parser.add_argument(
        '-o', '--output',
        type=Path,
        help='Output samplesheet path (default: <dir_name>_SampleSheet.csv)'
    )

    parser.add_argument(
        '--relative-paths',
        action='store_true',
        help='Use relative paths instead of absolute paths in the samplesheet'
    )

    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Enable verbose logging'
    )

    args = parser.parse_args()

    # Set logging level
    if args.verbose:
        logger.setLevel(logging.DEBUG)

    # Validate input directory
    if not args.input_dir.exists():
        logger.error(f"Input directory does not exist: {args.input_dir}")
        sys.exit(1)

    if not args.input_dir.is_dir():
        logger.error(f"Input path is not a directory: {args.input_dir}")
        sys.exit(1)

    # Determine output path
    if args.output:
        output_path = args.output
    else:
        dir_name = args.input_dir.resolve().name
        output_path = args.input_dir / f"{dir_name}_SampleSheet.csv"

    logger.info(f"Searching for FASTQ files in: {args.input_dir}")

    # Find and process FASTQ files
    fastq_files = find_fastq_files(args.input_dir)

    if not fastq_files:
        logger.error("No FASTQ files found in the specified directory")
        sys.exit(1)

    # Pair files and extract sample information
    samples = pair_fastq_files(fastq_files)

    if not samples:
        logger.error("No samples could be extracted from FASTQ files")
        sys.exit(1)

    # Write samplesheet
    write_samplesheet(
        samples,
        output_path,
        use_absolute_paths=not args.relative_paths
    )

    logger.info("Samplesheet creation complete!")


if __name__ == '__main__':
    main()
