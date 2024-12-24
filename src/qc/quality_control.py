import subprocess

class QualityControl:
    def __init__(self, workdir):
        self.workdir = workdir
    
    def run_fastqc(self, fastq_files):
        cmd = f"fastqc --outdir {self.workdir}/qc_results --threads 8 {' '.join(fastq_files)}"
        subprocess.run(cmd, shell=True, check=True)
    
    def run_multiqc(self):
        cmd = f"multiqc {self.workdir}/qc_results -o {self.workdir}/multiqc_report"
        subprocess.run(cmd, shell=True, check=True)
