import subprocess

class ReadProcessor:
    def __init__(self, genome_index):
        self.genome_index = genome_index
    
    def align_reads(self, fastq_file, output_prefix):
        bam_file = f"{output_prefix}.bam"
        # Handles compressed gz files, encode_downloader.py outputs replicates all in fastq.gz
        cmd = f"zcat {fastq_file} | bwa mem -t 8 {self.genome_index} - | samtools sort -@ 8 -o {bam_file} -"
        subprocess.run(cmd, shell=True, check=True)
        subprocess.run(f"samtools index {bam_file}", shell=True, check=True)
        return bam_file
