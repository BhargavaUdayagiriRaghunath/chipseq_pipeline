import subprocess
import os

class PeakCaller:
    def __init__(self, output_dir):
        self.output_dir = output_dir
    
    def call_peaks(self, treatment_bam, control_bam):
        output_prefix = os.path.join(self.output_dir, "peaks")
        cmd = f"macs2 callpeak -t {treatment_bam} -c {control_bam} -f BAM -g hs --outdir {self.output_dir} -n {output_prefix}"
        subprocess.run(cmd, shell=True, check=True)
        return f"{output_prefix}_peaks.narrowPeak"
