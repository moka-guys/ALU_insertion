import subprocess
import sys
import argparse
import re

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

    # capture a list of R134 files present in the output directory of the dx project
    result = subprocess.run(["dx", "ls", "output/*R134*"], capture_output=True, text=True)
    r134 = result.stdout.splitlines()
    r134.sort()

    return r134

def sequence_search(r134_file):
    # filter R134 files to those that are non-refined bam files
    if ("bam" in file and "refined" not in file):
        # save bam and bai file names and paths to variables.
        bam_name = file
        bai_name = bam_name.replace("bam","bai")

        bam_path = "output/" + bam_name
        bai_path = "output/" + bai_name

        print(bam_path)


def main():
    """ Runs ALU detection and analysis on LDLR for the provided sample BAM. """

    # put project selection here, I think. move the below functionality to another python script. call that one here,
    # once per bam file.

    parser = argparse.ArgumentParser(description='Run Scramble analysis and ALU filtering')
    parser.add_argument('--bam', help='Input BAM file path')
    parser.add_argument('--bai', help='Input BAI file path')
    parser.add_argument("--window", type=int, default=50, help="Coverage analysis window size (bp)")
    parser.add_argument("--polyA_window", type=int, default=10, help="PolyA detection window size (bp)")
    parser.add_argument("--threshold", type=float, default=15.0, help="Coverage change threshold (%)")
    parser.add_argument("--min_polyA_len", type=int, default=10, help="Minimum polyA/T stretch length")
    parser.add_argument("--merge_gap", type=int, default=1, help="Maximum gap for merging nearby variants")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging")
    parser.add_argument("--bed", help='Optional input BED file to limit what genome regions are analysed', default=None)
    parser.add_argument("--dx_project_id", help='ID for DNAnexus project we want to run the ALU analysis on', default = None)

    args = parser.parse_args()

    # project scan using project id parameter
    r134_list = project_scan(args.dx_project_id)
    for file in r134_list:
        sequence_search(file)

    # Extract sample ID from BAM file name
    match = re.search(r"(NGS[^_]+_\d+)", args.bam)
    if match:
        sample_id = match.group(1)
        print(sample_id)
    else:
        raise ValueError(f"Could not extract sample ID from BAM file: {args.bam}")

    print(f"Running cluster_identifier...")
    run_command(
        f"cluster_identifier -m 10 -s 3 {args.bam} > /app/output/{sample_id}.clusters.txt"
    )

    print(f"Running SCRAMble.R...")
    run_command(
        f"Rscript --vanilla /app/cluster_analysis/bin/SCRAMble.R "
        f"--out-name /app/output/{sample_id} "
        f"--cluster-file /app/output/{sample_id}.clusters.txt "
        f"--install-dir /app/cluster_analysis/bin "
        f"--mei-refs /app/cluster_analysis/resources/MEI_consensus_seqs.fa "
        f"--ref /app/data/reference.fa "
        f"--eval-meis"
    )

    print(f"Running bcftools filtering...")
    run_command(f"bgzip /app/output/{sample_id}.vcf -f")
    run_command(f"bcftools index /app/output/{sample_id}.vcf.gz")
    run_command(f"bcftools view -i 'ALT=\"<INS:ME:ALU>\"' /app/output/{sample_id}.vcf.gz -o /app/output/{sample_id}_ALU_ins.vcf")

    if args.bed:
        run_command(f"bgzip /app/output/{sample_id}_ALU_ins.vcf -f")
        run_command(f"bcftools index /app/output/{sample_id}_ALU_ins.vcf.gz")
        run_command(f"bcftools view -R {args.bed} /app/output/{sample_id}_ALU_ins.vcf.gz -o /app/output/{sample_id}_specified_region_ALU_ins.vcf")
        vcf_path = f"/app/output/{sample_id}_specified_region_ALU_ins.vcf"
    else:
        vcf_path = f"/app/output/{sample_id}_ALU_ins.vcf"

    print(f"Running ALU analysis...")
    python_command = (
        f"python /app/scramble_filtering_vcf_updated_v2.py "
        f"--vcf {vcf_path} "
        f"--bam {args.bam} "
        f"--window {args.window} "
        f"--polyA_window {args.polyA_window} "
        f"--threshold {args.threshold} "
        f"--min_polyA_len {args.min_polyA_len} "
        f"--merge_gap {args.merge_gap} "
        f"--id {sample_id}"
    )

    if args.verbose:
        python_command += " --verbose"

    run_command(python_command)


if __name__ == "__main__":
    main()
