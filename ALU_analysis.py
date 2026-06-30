import subprocess
import sys
import argparse
import re
import pysam
import os
from collections import Counter
import threading
import queue
from concurrent.futures import ThreadPoolExecutor

def run_command(command, check=True):
    """ Executes a command and returns the result."""
    try:
        result = subprocess.run(
            command,
            shell=True,
            check=check
        )
        return result
    except subprocess.CalledProcessError as e:
        print(f"Error executing command: {e}", file=sys.stderr)
        print(f"stdout: {e.stdout}", file=sys.stderr)
        print(f"stderr: {e.stderr}", file=sys.stderr)
        raise

def project_scan(dx_project_id):
    # extract run name from project id
    project_describe = subprocess.run(["dx", "describe", dx_project_id], capture_output=True, text=True)
    project_name = re.search(r"(NGS\d+[AB]?)", project_describe.stdout)

    # open the dx project via the CLI
    cmd = ["dx", "select", dx_project_id]
    subprocess.run(cmd)

    result = subprocess.run(["dx", "ls", "output/"], capture_output=True, text=True)
    all_files = result.stdout.splitlines()
    r134 = [f for f in all_files if "R134" in f]
    r134.sort()

    return r134

def process_bam(r134_file):
    # save bam and bai file names and paths to variables.
    bam_name = r134_file
    bai_name = bam_name.replace("bam","bai")

    bam_path = "output/" + bam_name
    bai_path = "output/" + bai_name

    # extract bam file dx id
    bam_describe = subprocess.run(["dx", "describe", bam_path], capture_output=True, text=True)
    bam_describe = bam_describe.stdout
    pattern = "file-[A-Za-z0-9]{24}"
    bam_id = re.findall(pattern,bam_describe)

    # extract bai file dx id
    bai_describe = subprocess.run(["dx", "describe", bai_path], capture_output=True, text=True)
    bai_describe = bai_describe.stdout
    bai_id = re.findall(pattern,bai_describe)

    return bam_name, bai_name, bam_id[0], bai_id[0]

def sequence_search(bam_id, bai_id,sample_id, bam_name,bai_name):
    # download bam and bai files
    #cmd = ["dx", "download", bam_id, "--no-progress"]
    #subprocess.run(cmd)
    #print("bam file successfully downloaded")
    #cmd = ["dx", "download", bai_id, "--no-progress"]
    #subprocess.run(cmd)
    #print("bai file successfully downloaded")

    # search for the ALU right flanking sequence in the sample bam file.
    # if counts at a single position exceed 100, save to output file
    right_seq = "GGCCGGGCGCGGTGGCTCACGCCTGTAATCC"
            
    # longer right flank seq
    long_right_seq = "TGGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG"

    # left flank seq we see in the ng683 samples, and not the ngs625 sample
    ngs683_38_left_flank_seq = "CGACATTGCGCCACTGCAGTCGGCAGTCCGGCCTGGGCGACAGAGCGAGACTCCATCTCAAAAAAAAAATAATAATAA"
            
    #seq_used_by_oxford 
    ox_seq = "AGAACTGGCGGCTTAAGAACATCAACAGCATTGGCCGGGCGCGGTGGCTCACG"

    sequences = {
        "500": (right_seq, 500),
        "long_right": (long_right_seq, 100),
        "ox": (ox_seq, 100),
        "left_flank": (ngs683_38_left_flank_seq, 100),
    }

    positions = {name: [] for name in sequences}

    with pysam.AlignmentFile(bam_name, "rb") as bam:
        for read in bam.fetch():
            if not read.query_sequence:
                continue

            for name, (seq, threshold) in sequences.items():
                if seq not in read.query_sequence:
                    continue

                match_start = read.query_sequence.index(seq)  # 0-based read position

                # get_aligned_pairs maps each read position to its genomic position
                aligned_pairs = read.get_aligned_pairs(matches_only=False)

                last_genome_pos = None
                last_chrom = None
                for read_pos, genome_pos in aligned_pairs:
                    if genome_pos is not None:
                        last_genome_pos = genome_pos
                        last_chrom = read.reference_name
                    if read_pos == match_start:
                        if genome_pos is not None:
                            positions[name].append(f"{read.reference_name}:{genome_pos + 1}")
                        elif last_genome_pos is not None:
                            # seq falls in an insertion — use the last mapped genomic position
                            positions[name].append(f"{last_chrom}:{last_genome_pos + 1}")
                        break

    for name, (seq, threshold) in sequences.items():
        results = [
            (pos, count)
            for pos, count in Counter(positions[name]).most_common()
            if count > threshold
        ]

        if results:
            with open(f"/app/output/{sample_id}_output_{name}.txt", "w") as out:
                for pos, count in results:
                    out.write(f"{count} {pos}\n")

        print(f"--- subanalysis {name} complete for {sample_id}!")
                

def scramble_analysis(bam_name,sample_id,bed,window,polyA_window,threshold,min_polyA_len,merge_gap,work_dir):
    os.makedirs(work_dir, exist_ok=True)
    print(f"Running cluster_identifier...")
    run_command(
        f"cluster_identifier -m 10 -s 3 {bam_name} > {work_dir}/{sample_id}.clusters.txt"
    )

    print(f"Running SCRAMble.R...")
    run_command(
        f"Rscript --vanilla /app/cluster_analysis/bin/SCRAMble.R "
        f"--out-name {work_dir}/{sample_id} "
        f"--cluster-file {work_dir}/{sample_id}.clusters.txt "
        f"--install-dir /app/cluster_analysis/bin "
        f"--mei-refs /app/cluster_analysis/resources/MEI_consensus_seqs.fa "
        f"--ref /app/data/reference.fa "
        f"--eval-meis"
    )

    print(f"Running bcftools filtering...")
    run_command(f"bgzip {work_dir}/{sample_id}.vcf -f")
    run_command(f"bcftools index {work_dir}/{sample_id}.vcf.gz")
    run_command(f"bcftools view -i 'ALT=\"<INS:ME:ALU>\" && QUAL>=100' {work_dir}/{sample_id}.vcf.gz -o {work_dir}/{sample_id}_ALU_ins.vcf")

    if bed:
        run_command(f"bgzip {work_dir}/{sample_id}_ALU_ins.vcf -f")
        run_command(f"bcftools index {work_dir}/{sample_id}_ALU_ins.vcf.gz")
        run_command(f"bcftools view -R {bed} {work_dir}{sample_id}_ALU_ins.vcf.gz -o {work_dir}/{sample_id}_specified_region_ALU_ins.vcf")
        vcf_path = f"{work_dir}/{sample_id}_specified_region_ALU_ins.vcf"
    else:
        vcf_path = f"{work_dir}/{sample_id}_ALU_ins.vcf"

    print(f"Running ALU analysis...")
    python_command = (
        f"python /app/scramble_filtering_vcf_updated_v2.py "
        f"--vcf {vcf_path} "
        f"--bam {bam_name} "
        f"--window {window} "
        f"--polyA_window {polyA_window} "
        f"--threshold {threshold} "
        f"--min_polyA_len {min_polyA_len} "
        f"--merge_gap {merge_gap} "
        f"--id {sample_id}"
    )

    run_command(python_command)

def downloader_thread(bam_files, result_queue):
    """Runs continuously in the background: downloads one sample after another."""
    for r134_file in bam_files:
        bam_name, bai_name, bam_id, bai_id = process_bam(r134_file)
        subprocess.run(["dx", "download", bam_id, "--no-progress"])
        print(f"bam file successfully downloaded: {bam_name}")
        subprocess.run(["dx", "download", bai_id, "--no-progress"])
        print(f"bai file successfully downloaded: {bai_name}")
        result_queue.put((bam_name, bai_name, bam_id, bai_id))
    result_queue.put(None)  # sentinel: signals "no more downloads coming"

import shutil

def process_one_sample(item, bed, window, polyA_window, threshold, min_polyA_len, merge_gap, scramble_failed):
    bam_name, bai_name, bam_id, bai_id = item
    match = re.search(r"(NGS[^_]+_\d+)", bam_name)
    if match:
        sample_id = match.group(1)
        print(sample_id)
    else:
        raise ValueError(f"Could not extract sample ID from BAM file: {bam_name}")

    work_dir = f"/app/work/{sample_id}"

    print(f"starting sequence search for {sample_id}")
    sequence_search(bam_id, bai_id, sample_id, bam_name, bai_name)
    print(f"sequence search done for {sample_id}")
    try:
        scramble_analysis(bam_name, sample_id, bed, window, polyA_window, threshold, min_polyA_len, merge_gap,work_dir)
    except:
        scramble_failed.append(sample_id)
    os.remove(bam_name)
    os.remove(bai_name)

    # remove intermediate files generated by scramble_analysis for this sample
    shutil.rmtree(work_dir, ignore_errors=True)

    print(f"Analysis of {sample_id} completed!")

def alu_analysis(dx_project_id, bed, window, polyA_window, threshold, min_polyA_len, merge_gap, max_concurrent=5):
    r134_list = project_scan(dx_project_id)
    bam_files = [f for f in r134_list if "bam" in f and "refined" not in f]
    scramble_failed = []

    if not bam_files:
        print("No matching BAM files found.")
        return

    result_queue = queue.Queue(maxsize=max_concurrent)
    dl_thread = threading.Thread(target=downloader_thread, args=(bam_files, result_queue))
    dl_thread.start()

    futures = []
    with ThreadPoolExecutor(max_workers=max_concurrent) as executor:
        while True:
            item = result_queue.get()
            if item is None:
                break  # all downloads finished and consumed
            futures.append(
                executor.submit(
                    process_one_sample, item, bed, window, polyA_window, threshold, min_polyA_len, merge_gap, scramble_failed
            ))

        for f in futures:
            f.result()

    dl_thread.join()

    if scramble_failed:
        with open ('/app/output/failed_scramble.txt', 'w') as f:
            for sample in scramble_failed:
                f.write('%s\n' %sample)
    
        f.close()
    


def main():
    """ Runs ALU detection and analysis on LDLR for the provided sample BAM. """

    parser = argparse.ArgumentParser(description='Run Scramble analysis and ALU filtering')
    #parser.add_argument('--bam', help='Input BAM file path')
    #parser.add_argument('--bai', help='Input BAI file path')
    parser.add_argument("--window", type=int, default=50, help="Coverage analysis window size (bp)")
    parser.add_argument("--polyA_window", type=int, default=10, help="PolyA detection window size (bp)")
    parser.add_argument("--threshold", type=float, default=12.0, help="Coverage change threshold (%)")
    parser.add_argument("--min_polyA_len", type=int, default=10, help="Minimum polyA/T stretch length")
    parser.add_argument("--merge_gap", type=int, default=1, help="Maximum gap for merging nearby variants")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging")
    parser.add_argument("--bed", help='Optional input BED file to limit what genome regions are analysed', default=None)
    parser.add_argument("--dx_project_id", help='ID for DNAnexus project we want to run the ALU analysis on', default = None)

    args = parser.parse_args()

    # create output directory
    try:
        os.makedirs("/app/output/", exist_ok=True)
        print(f"/app/output/ exists and is writable")
    except Exception as e:
        print(f"Cannot create /app/output/: {e}")

    # if the above functions are all cool, this might be all that's required...
    alu_analysis(args.dx_project_id,args.bed,args.window,args.polyA_window,args.threshold,args.min_polyA_len,args.merge_gap)


if __name__ == "__main__":
    main()
