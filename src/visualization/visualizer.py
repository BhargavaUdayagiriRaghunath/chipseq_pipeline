import subprocess

class ChIPSeqVisualizer:
    def __init__(self, output_dir):
        self.output_dir = output_dir
    
    def create_coverage_track(self, bam_file):
        output = f"{self.output_dir}/coverage.bw"
        cmd = f"bamCoverage -b {bam_file} -o {output} --normalizeUsing RPKM --binSize 10"
        subprocess.run(cmd, shell=True, check=True)
        return output
