#!/usr/bin/env python3
"""
ALU Insertion Event analyser
analyses VCF variants and BAM alignments to identify ALU insertion events
based on coverage patterns and polyA/T tail signatures.
"""

import pysam
pysam.set_verbosity(0)

import argparse
import re
import sys
import logging
import statistics
from typing import List, Dict, Tuple
from dataclasses import dataclass
from pathlib import Path
from collections import defaultdict
import pandas as pd

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# -------------------------------
# Data Classes
# -------------------------------
@dataclass
class CoverageTract:
    """ A contiguous region where coverage is either abnormally high or low.
    
    Attributes:
        start: inclusive start coordinate of the region.
        end: Inclusive end coordinatge of the region.
        direction: Indicates whether the coverage is abnormally HIGH or LOW compared to surroundings.
        length: Length of the region in bases.
    """
    
    start: int
    end: int
    direction: str  # 'HIGH' or 'LOW'
    length: int
    
    def __str__(self) -> str:
        return f"{self.start}-{self.end}({self.direction},len={self.length})"


@dataclass
class EventAnalysis:
    """ Analysis information for a single potential ALU.
    
    Attributes:
        has_coverage: Whether there is a nearby significant coverage change from one base to the next. 
        coverage_tracts: List of contiguous regions of abnormal coverage found near the potential ALU.
        polyA_reads: Number of polyA tracts found near the potential ALU.
        high_tracts: List of contiguous regions of abnormally high coverage found near the potential ALU.
        evidence: A string declaring whether the potential ALU has COVERAGE, POLY-A or BOTH nearby.
    """

    has_coverage: bool
    coverage_tracts: List[CoverageTract]
    polyA_reads: int
    high_tracts: List[CoverageTract]
    evidence: str


# -------------------------------
# Cluster VCF Variants
# -------------------------------
def cluster_variants(records: List, max_gap: int) -> List[List]:
    """ Combines multiple adjacent scramble calls into a cluster.
    
    Parameters:
        records: A list of pysam VCF record objects.
        max_gap: The maximum distance in bases between variants for them to be considered adjacent. 

    Returns:
        clusters: A list of pysam VCF record objects making up a cluster of adjacent variants.
    """

    if not records:
        return []
    
    clusters = []
    current_cluster = [records[0]]
    
    for rec in records[1:]:
        prev = current_cluster[-1]
        if rec.chrom == prev.chrom and (rec.pos - prev.pos) <= max_gap:
            current_cluster.append(rec)
        else:
            clusters.append(current_cluster)
            current_cluster = [rec]
    
    if current_cluster:
        clusters.append(current_cluster)

    return clusters


def representative_position(cluster: List) -> int:
    """ Creates a single representative genomic position for a cluster of adjacent scramble calls.

    Parameters:
        cluster: A list of pysam VCF record objects making up a cluster of adjacent variants.

    Returns:
        The median genomic position of the cluster.
    """
    
    positions = sorted(rec.pos for rec in cluster)
    n = len(positions)
    if n % 2 == 0:
        return (positions[n//2 - 1] + positions[n//2]) // 2
    return positions[n//2]


# -------------------------------
# Coverage Analysis
# -------------------------------
def compute_coverage(
    bam: pysam.AlignmentFile,
    chrom: str,
    start: int,
    end: int,
    min_base_quality: int = 0
) -> Dict[int, int]:
    """ Calculates coverage per base in a specified region.

    Parameters:
        bam: Pysam alignment file object for the sample.
        chrom: Chromosome to calculate coverage for.
        start: Inclusive start coordinate to calculate coverage from. 
        end: Inclusive end coordinate to calculated coverage to. of the candidate ALU.
        min_base_quality: The minimum base quality a base must have before being counted towards coverage at a position.

    Returns:
        coverage: A dictionary containing the calculated coverage values for each base position.
    """

    coverage = defaultdict(int)
    
    try:
        for col in bam.pileup(
            chrom,
            start - 1,
            end,
            truncate=True,
            stepper="all",
            min_base_quality=min_base_quality
        ):
            ref_pos = col.reference_pos + 1
            depth = sum(
                1 for pr in col.pileups
                if not pr.is_del and not pr.is_refskip
            )
            coverage[ref_pos] = depth
    except Exception as e:
        logger.warning(f"Error computing coverage for {chrom}:{start}-{end}: {e}")

    return coverage

def detect_coverage_tracts(
    coverage: Dict[int, int],
    pos: int,
    window: int,
    threshold: float,
    min_length: int = 3,
    min_coverage: int = 30
) -> List[CoverageTract]:
    """ Detect regions of significant coverage change where coverage change occurs within a window.
    Tracts are defined by sharp increases (start) and sharp decreases (end).
    
    Parameters:
        coverage: A dictionary containing the calculated coverage values for a series of base positions.
        pos: Genomic coordinate of a candidate ALU.
        window: Distance either side of the candidate ALU coordinate to detect abnormal coverage tracts.
        threshold: Percentage value used as the threshold at which base-to-base coverage changes are considered significant.

    Returns:
        tracts: A list of CoverageTract objects (see CoverageTract class definition above).
    """

    for p in range(pos - window, pos + window + 1):
        coverage.setdefault(p, 0)
    
    positions = sorted(coverage.keys())
    tracts = []
    i = 1
    
    while i < len(positions):
        prev_cov = coverage[positions[i - 1]]
        curr_cov = coverage[positions[i]]
        
        if prev_cov == 0:
            i += 1
            continue
        
        pct_change = ((curr_cov - prev_cov) / prev_cov) * 100
        
        if abs(pct_change) >= threshold:
            direction = "HIGH" if pct_change > 0 else "LOW"
            
            # For HIGH tracts: start at positions[i] (where coverage becomes high)
            # For LOW tracts: start at positions[i] (where coverage becomes low)
            tract_start = positions[i]
            
            tract_end = positions[i]
            last_cov = curr_cov
            j = i + 1
            
            while j < len(positions):
                next_cov = coverage[positions[j]]
                if next_cov == 0:
                    break
                
                step = ((next_cov - last_cov) / last_cov) * 100 if last_cov else 0
                
                # Check for the reverse change
                if (direction == "HIGH" and step <= -threshold) or \
                   (direction == "LOW" and step >= threshold):
                    # For HIGH tracts: end at positions[j-1] (last high position)
                    # For LOW tracts: end at positions[j-1] (last low position)
                    tract_end = positions[j - 1]
                    break
                
                tract_end = positions[j]
                last_cov = next_cov
                j += 1
            
            length = tract_end - tract_start + 1
            max_cov = max(coverage[p] for p in range(tract_start, tract_end + 1))
            
            if length >= min_length and max_cov >= min_coverage:
                tracts.append(
                    CoverageTract(
                        start=tract_start,
                        end=tract_end,
                        direction=direction,
                        length=length
                    )
                )
            
            i = j
        else:
            i += 1
    
    return tracts

def is_true_high_tract(
    tract: CoverageTract,
    coverage: Dict[int, int],
    threshold: float,
    debug: bool = False,
) -> bool:
    """
    Validates a HIGH coverage tract by verifying sharp step changes at boundaries.
    """
    left_cov = coverage.get(tract.start - 1, 0)
    entry_cov = coverage.get(tract.start, 0)
    right_cov = coverage.get(tract.end + 1, 0)
    exit_cov = coverage.get(tract.end, 0)

    if debug:
        print(f"Tract: {tract.start}-{tract.end}")
        print(f"Left boundary: pos={tract.start-1} cov={left_cov}, pos={tract.start} cov={entry_cov}")
        print(f"Right boundary: pos={tract.end} cov={exit_cov}, pos={tract.end+1} cov={right_cov}")

    if left_cov == 0 or right_cov == 0:
        if debug:
            print(f"REJECTED: zero coverage on boundary (left_cov={left_cov}, right_cov={right_cov})")
        return False

    left_step = ((entry_cov - left_cov) / left_cov) * 100
    
    # This gives the percentage DROP from the tract to the right flank
    right_step = ((right_cov - exit_cov) / exit_cov) * 100

    if debug:
        print(f"left_step={left_step:.1f}% (need >= {threshold})")
        print(f"right_step={right_step:.1f}% (need <= -{threshold})")

    result = left_step >= threshold and right_step <= -threshold
    if debug:
        print(f"RESULT: {result}")
    
    return result


# -------------------------------
# PolyA/T Detection
# -------------------------------
def count_polyA_reads(
    bam: pysam.AlignmentFile,
    chrom: str,
    start: int,
    end: int,
    min_polyA_len: int
) -> int:
    """
    Counts reads with soft-clipped polyA/T tails near the candidate position.
    Only soft-clipped sequence is examined, as novel ALU tails appear there
    rather than in the reference-aligned portion of the read.
    """
    poly_pattern = re.compile(f"(?:A{{{min_polyA_len},}}|T{{{min_polyA_len},}})")
    polyA_count = 0

    try:
        for read in bam.fetch(chrom, start, end):
            if read.is_unmapped or read.query_sequence is None or not read.cigartuples:
                continue

            seq = read.query_sequence.upper()
            cigar = read.cigartuples  # list of (op, length) tuples

            # Extract soft-clipped segments (cigar op 4 = SOFT_CLIP)
            # Left clip
            if cigar[0][0] == 4:
                clip_len = cigar[0][1]
                clipped_seq = seq[:clip_len]
                if poly_pattern.search(clipped_seq):
                    polyA_count += 1
                    continue

            # Right clip
            if cigar[-1][0] == 4:
                clip_len = cigar[-1][1]
                clipped_seq = seq[-clip_len:]
                if poly_pattern.search(clipped_seq):
                    polyA_count += 1

    except Exception as e:
        logger.warning(f"Error counting polyA reads for {chrom}:{start}-{end}: {e}")

    return polyA_count


# -------------------------------
# Event Analysis
# -------------------------------
def analyse_event(
    bam: pysam.AlignmentFile,
    chrom: str,
    pos: int,
    coverage_window: int,
    coverage_threshold: float,
    min_polyA_len: int,
    polyA_window: int
) -> EventAnalysis:
    """ Performs nearby coverage and PolyA analysis on a specified genomic coordinate.

    Parameters:
        bam: Pysam alignment file object for the sample.
        chrom: Chromosome to perform analysis on.
        pos: Chromosome coordinate to perform analysis on.
        coverage_window: Window size used for plus/minus coverage calculations around a position.
        coverage_threshold: Percentage value used as the threshold at which base-to-base coverage changes are considered significant.
        min_polyA_len: The minimum length a polyA tract must be to be counted.   
        polyA_window: Window size used for plus/minus PolyA counting around a position.

    Returns:
        An EventAnalysis object containing information as described in the EventAnalysis class definition.
    """
    cov = compute_coverage(
        bam,
        chrom,
        max(pos - coverage_window, 1),
        pos + coverage_window
    )

    coverage_tracts = detect_coverage_tracts(
        cov, pos, coverage_window, coverage_threshold
    )

    high_tracts = [
        t for t in coverage_tracts
        if t.direction == "HIGH" and is_true_high_tract(t, cov, coverage_threshold)
    ]

    polyA_reads = count_polyA_reads(
        bam,
        chrom,
        max(pos - polyA_window, 1),
        pos + polyA_window,
        min_polyA_len
    )

    if high_tracts and polyA_reads > 50: # do we want a pure number of polyA tails, or as a percentage of present reads? what is the percentage for our known ALUs?
        evidence = "BOTH"
    elif high_tracts:
        evidence = "COVERAGE_ONLY"
    elif polyA_reads > 50:
        evidence = "POLYA_ONLY"
    else:
        evidence = "NONE"

    return EventAnalysis(
        has_coverage=bool(high_tracts),
        coverage_tracts=coverage_tracts,
        polyA_reads=polyA_reads,
        high_tracts=high_tracts,
        evidence=evidence
    )

# -------------------------------
# Putting event analysis VCF into dataframe
# -------------------------------
def analyse_vcf_to_dataframe(
    vcf_path: str,
    bam_path: str,
    window: int = 50,
    polyA_window: int = 10,
    threshold: float = 10.0,
    min_polyA_len: int = 10,
    merge_gap: int = 1,
    verbose: bool = False
) -> pd.DataFrame:
    """ Returns analysis results as a pandas DataFrame.

    Parameters:
        vcf_path: Path to VCF file
        bam_path: Path to BAM file
        window: Coverage analysis window size (bp)
        polyA_window: PolyA detection window size (bp)
        threshold: Percentage value used as the threshold at which base-to-base coverage changes are considered significant.
        min_polyA_len: The minimum length a polyA tract must be to be counted.
        merge_gap: Maximum gap for merging nearby variants
        verbose: Enable verbose logging

    Returns:
        pandas DataFrame with analysis results

    """
    if verbose:
        logger.setLevel(logging.DEBUG)

    # Validate inputs
    valid, msg = validate_inputs(vcf_path, bam_path)
    if not valid:
        raise FileNotFoundError(msg)

    # Open files
    vcf = pysam.VariantFile(vcf_path)
    bam = pysam.AlignmentFile(bam_path, "rb")

    # Initiate results list for later collection
    results = []
    # Extract list of variants in VCF file
    records = list(vcf.fetch())
    # If variants are within {merge_gap} bp of each other, combine them
    clusters = cluster_variants(records, merge_gap)

    logger.info(f"Processing {len(clusters)} ALU insertion clusters...")

    for cluster in clusters:
        chrom = cluster[0].chrom
        start = min(r.pos for r in cluster)
        end = max(r.pos for r in cluster)
        rep_pos = representative_position(cluster)

        analysis = analyse_event(
            bam,
            chrom,
            rep_pos,
            window,
            threshold,
            min_polyA_len,
            polyA_window
        )

        # Format coverage tract info for output
        coverage_str = (
            ",".join(str(t) for t in analysis.high_tracts)
            if analysis.high_tracts else "NA"
        )

        # Add analysis info for cluster to results
        results.append({
            'Chrom': chrom,
            'Start': start,
            'End': end,
            'Representative_Pos': rep_pos,
            'Cluster_Size': len(cluster),
            'HighCoverageTracts': coverage_str,
            'PolyAReads': analysis.polyA_reads,
            'Evidence': analysis.evidence,
            'HasCoverage': analysis.has_coverage,
            'TotalCoverageTracts': len(analysis.coverage_tracts)
        })

    # Close files
    bam.close()
    vcf.close()

    # Create DataFrame
    df = pd.DataFrame(results)

    logger.info(f"Analysis complete. Found {len(df)} ALU insertion events.")

    return df

# -------------------------------
# Terminal printing function
# -------------------------------
def print_results_table(df: pd.DataFrame):
    """Print DataFrame as a formatted table to the terminal.

    Parameters:
        df: A pandas DataFrame
    """

    print("Chrom\tStat\tEnd\tHighCoverageTracts\tPolyAReads\tEvidence")
    for _, row in df.iterrows():
        print(f"{row['Chrom']}\t{row['Start']}\t{row['End']}\t{row['HighCoverageTracts']}\t{row['PolyAReads']}\t{row['Evidence']}")


# -------------------------------
# Output Functions
# -------------------------------
# def print_header():
#     """ Prints header of output table. """
#     print("Chrom\tStart\tEnd\tHighCoverageTracts\tPolyAReads\tEvidence")


# def print_result(chrom: str, start: int, end: int, analysis: EventAnalysis):
#      """ Prints results of candidate ALU analysis

#     Parameters:
#         chrom: Chromosome of candidate ALU.
#         start: Start coordinate of candidate ALU.
#         end: End coordinate of candidate ALU.
#         analysis: EventAnalysis object for candidate ALU, containing information on whether nearby coverage and/or polyA tracts were detected. 
#     """
#     coverage_str = (
#         ",".join(str(t) for t in analysis.high_tracts)
#         if analysis.high_tracts else "NA"
#     )
#     print(f"{chrom}\t{start}\t{end}\t{coverage_str}\t{analysis.polyA_reads}\t{analysis.evidence}")


# -------------------------------
# Input Validation
# -------------------------------
def validate_inputs(vcf_path: str, bam_path: str) -> Tuple[bool, str]:
    """ Validates whether a vcf and bam file exist and have indexes.

    Parameters:
        vcf_path: File path of VCF file.
        bam_path: file path of BAM file.

    Returns:
        Returns a Boolean indicating whether the files and indexes were found. 
        If False, includes error message indicating files were not found or indexing needed.
    """

    vcf_file = Path(vcf_path)
    bam_file = Path(bam_path)
    
    if not vcf_file.exists():
        return False, f"VCF file not found: {vcf_path}"
    if not bam_file.exists():
        return False, f"BAM file not found: {bam_path}"
    
    if not (Path(str(bam_file) + ".bai").exists() or
            Path(str(bam_file).replace(".bam", ".bai")).exists()):
        return False, f"BAM index not found. Please run: samtools index {bam_path}"

    return True, ""


# -------------------------------
# Main
# -------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="analyse ALU insertion events from VCF and BAM files"
    )
    """ 
    Runs analysis of candidate ALUs present in the input VCF file.
    Detects whether there are nearby high coverage regions or polyA tracts for any of the candidate ALUs.
    """

    parser.add_argument("--vcf", required=True)
    parser.add_argument("--window", type=int, default=50)
    parser.add_argument("--polyA_window", type=int, default=10)
    parser.add_argument("--threshold", type=float, default=10.0)
    parser.add_argument("--min_polyA_len", type=int, default=10)
    parser.add_argument("--merge_gap", type=int, default=1)
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--id", type=str)

    args = parser.parse_args()
    if args.verbose:
        logger.setLevel(logging.DEBUG)

    valid, msg = validate_inputs(args.vcf, args.bam)
    if not valid:
        logger.error(msg)
        sys.exit(1)

    # Run analysis
    df = analyse_vcf_to_dataframe(
            vcf_path=args.vcf,
            bam_path=args.bam,
            window=args.window,
            polyA_window=args.polyA_window,
            threshold=args.threshold,
            min_polyA_len=args.min_polyA_len,
            merge_gap=args.merge_gap,
            verbose=args.verbose,
    )

    # Print results to terminal
    print_results_table(df)

    # Save results to CSV
    csv_path = f"/app/output/{args.id}_ALU_analysis.csv"
    df.to_csv(csv_path, index=False)

    # Creating a separate CSV containing only ALUs where coverage abnormalities and polyAs were detected
    df_high_confidence = df[df['Evidence'] == 'BOTH']
    csv_hc_path = f"/app/output/{args.id}_ALU_analysis_high_confidence.csv"
    df_high_confidence.to_csv(csv_hc_path, index=False)

    # Print high confidence results to terminal
    print_results_table(df_high_confidence)



    # vcf = pysam.VariantFile(args.vcf)
    # bam = pysam.AlignmentFile(args.bam, "rb")

    # records = list(vcf.fetch())
    # clusters = cluster_variants(records, args.merge_gap)

    # print_header()
    # for cluster in clusters:
    #     chrom = cluster[0].chrom
    #     start = min(r.pos for r in cluster)
    #     end = max(r.pos for r in cluster)
    #     rep_pos = representative_position(cluster)

    #     analysis = analyse_event(
    #         bam,
    #         chrom,
    #         rep_pos,
    #         args.window,
    #         args.threshold,
    #         args.min_polyA_len,
    #         args.polyA_window,
    #         args.proximity
    #     )

    #     print_result(chrom, start, end, analysis)

    # bam.close()
    # vcf.close()
    logger.info("Analysis complete!")


if __name__ == "__main__":
    main()
