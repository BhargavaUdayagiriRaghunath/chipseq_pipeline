from src.download.encode_downloader import parse_arguments, download_encode_files
from src.qc.quality_control import QualityControl
from src.alignment.read_processor import ReadProcessor
from src.peaks.peak_caller import PeakCaller
from src.visualization.visualizer import ChIPSeqVisualizer
import yaml

def load_config(config_file):
    with open(config_file, 'r') as f:
        return yaml.safe_load(f)

def run_pipeline(args, config):
    qc = QualityControl(args.dir)
    processor = ReadProcessor(config['genome_index'])
    peak_caller = PeakCaller(f"{args.dir}/peaks")
    visualizer = ChIPSeqVisualizer(f"{args.dir}/visualization")
    
    fastq_files = download_encode_files(args)
    qc.run_fastqc(fastq_files)
    qc.run_multiqc()
    bam_files = [processor.align_reads(f, f"{args.dir}/aligned_{i}") for i, f in enumerate(fastq_files)]
    peaks = peak_caller.call_peaks(bam_files[0], bam_files[1])
    visualizer.create_coverage_track(bam_files[0])

if __name__ == "__main__":
    args = parse_arguments()
    config = load_config('config/pipeline_config.yaml')
    run_pipeline(args, config)
