#!/usr/bin/env python
import os
import time
import json
import requests
import argparse
from collections import OrderedDict

ENCODE_BASE_URL = 'https://www.encodeproject.org'

def parse_arguments():
    parser = argparse.ArgumentParser(prog='ENCODE downloader', description='Downloads genome data files from the ENCODE portal \
                                        and save them under [WORK_DIR]/[ACCESSION_ID]/. \
                                        If authentication information (--encode-access-key-id and --encode-secret-key) is given, \
                                        unpublished files only visible to submitters with valid authentication \
                                        can be downloaded.')
    parser.add_argument('url_or_file', metavar='url-or-file', nargs='+', type=str, help='List of URLs/files/accession_ids.')
    parser.add_argument('--dir', default='.', type=str, help='Root directory for all downloaded genome data.')
    parser.add_argument('--file-types', nargs='+', default=['fastq'], type=str,
                        help='List of file types to be downloaded. Default: fastq.')
    parser.add_argument('--assemblies', nargs='+', default=['all'], type=str,
                        help='Assemblies allowed. e.g. --assemblies GRCh38 hg19.')
    parser.add_argument('--encode-access-key-id', type=str, help='ENCODE access key ID.')
    parser.add_argument('--encode-secret-key', type=str, help='ENCODE secret key.')
    parser.add_argument('--ignored-accession-ids-file', type=str, help='Text file with ignored accession IDs.')
    parser.add_argument('--pooled-rep-only', action="store_true", help='Download genome data from pooled replicates only.')
    parser.add_argument('--dry-run', action="store_true", help='Dry-run: downloads nothing.')
    parser.add_argument('--dry-run-list-accession-ids', action="store_true", help='Dry-run: lists accession IDs.')
    parser.add_argument('--max-download', type=int, default=8, help='Maximum number of concurrent downloads.')
    parser.add_argument('--assembly-map', nargs='+', default=['Mus+musculus:mm10', 'Homo+sapiens:GRCh38'], type=str,
                        help='List of strings to infer ENCODE assembly from species name.')
    group_ignore_status = parser.add_mutually_exclusive_group()
    group_ignore_status.add_argument('--ignore-released', action='store_true', help='Ignore released data.')
    group_ignore_status.add_argument('--ignore-unpublished', action='store_true', help='Ignore unpublished data.')
    args = parser.parse_args()

    if (args.encode_access_key_id and not args.encode_secret_key) or (not args.encode_access_key_id and args.encode_secret_key):
        raise ValueError("Both --encode-access-key-id and --encode-secret-key must be specified together.")

    args.dir = os.path.abspath(args.dir)
    args.file_types = [file_type.lower() for file_type in args.file_types]
    args.assemblies = ['GRCh38' if assembly == 'hg38' else assembly for assembly in args.assemblies]
    return args

def is_encode_search_query_url(url):
    return url.startswith(f'{ENCODE_BASE_URL}/search/?')


def is_encode_exp_url(url):
    return url.startswith(f'{ENCODE_BASE_URL}/experiments/ENCSR')


def get_accession_id_from_encode_exp_url(url):
    for part in url.split('/'):
        if part.startswith('ENC'):
            return part
    return None


def download_file(url, destination):
    response = requests.get(url, stream=True)
    response.raise_for_status()
    with open(destination, 'wb') as file:
        for chunk in response.iter_content(chunk_size=8192):
            file.write(chunk)


def get_accession_ids_from_file(filepath):
    with open(filepath, 'r') as f:
        return [line.strip() for line in f if line.strip() and not line.startswith("#")]

def main():
    args = parse_arguments()
    if args.encode_access_key_id:
        auth = (args.encode_access_key_id, args.encode_secret_key)
    else:
        auth = None

    headers = {'accept': 'application/json'}
    accession_ids = []

    for item in args.url_or_file:
        if is_encode_search_query_url(item):
            if 'limit=all' not in item:
                item += '&limit=all'
            if 'format=json' not in item:
                item += '&format=json'
            response = requests.get(item, headers=headers, auth=auth)
            response.raise_for_status()
            json_data = response.json()
            accession_ids.extend([exp['accession'] for exp in json_data['@graph']])
        elif is_encode_exp_url(item):
            accession_id = get_accession_id_from_encode_exp_url(item)
            accession_ids.append(accession_id)
        elif os.path.isfile(item):
            accession_ids.extend(get_accession_ids_from_file(item))
        elif item.startswith('ENCSR'):
            accession_ids.append(item)
        else:
            raise ValueError(f"Invalid input: {item}")

    if args.dry_run_list_accession_ids:
        print(accession_ids)
        return

    os.makedirs(args.dir, exist_ok=True)
    ignored_accession_ids = get_accession_ids_from_file(args.ignored_accession_ids_file) if args.ignored_accession_ids_file else []

    all_file_metadata = OrderedDict()
    
    for accession_id in accession_ids:
        if accession_id in ignored_accession_ids:
            print(f"Ignored accession: {accession_id}")
            continue

        exp_url = f'{ENCODE_BASE_URL}/experiments/{accession_id}?format=json'
        response = requests.get(exp_url, headers=headers, auth=auth)
        response.raise_for_status()
        exp_data = response.json()

        if exp_data.get('status') == 'error':
            print(f"Error accessing accession {accession_id}")
            continue

        assembly = None
        for mapping in args.assembly_map:
            species, asm = mapping.replace('+', ' ').split(':')
            if species in str(exp_data):
                assembly = asm
                break

        if not assembly:
            print(f"Could not infer assembly for {accession_id}")
            continue

        metadata = {k: v for k, v in exp_data.items() if not isinstance(v, (dict, list))}
        metadata['files'] = {}
        downloaded = False

        for file_info in exp_data.get('original_files', []):
            file_response = requests.get(f"{ENCODE_BASE_URL}{file_info}?format=json", headers=headers, auth=auth)
            file_response.raise_for_status()
            file_data = file_response.json()
            status = file_data.get('status', '').lower().replace(' ', '_')
            if status == 'error':
                continue
            file_type = file_data.get('file_type', '').lower()
            file_format = file_data.get('file_format', '').lower()
            output_type = file_data.get('output_type', '').lower()
            file_assembly = file_data.get('assembly', '')

            if not any(ft in f"{file_type}:{file_format}:{output_type}" for ft in args.file_types):
                continue
            if args.ignore_released and status == 'released':
                continue
            if args.ignore_unpublished and status != 'released':
                continue
            if file_type == 'fastq' and 'paired_end' in file_data:
                pair = int(file_data['paired_end'])
            else:
                pair = -1

            bio_rep_id = file_data.get('biological_replicates', [])
            if args.pooled_rep_only and len(bio_rep_id) < 2:
                continue

            dir_suffix = f"{accession_id}/{status}/{file_assembly}/{output_type.replace(' ', '_')}/{file_type.replace(' ', '_')}"
            if file_type != file_format:
                dir_suffix += f"/{file_format}"
            if bio_rep_id:
                dir_suffix += f"/rep{'_rep'.join(map(str, bio_rep_id))}"
            if pair > 0:
                dir_suffix += f"/pair{pair}"

            download_dir = os.path.join(args.dir, dir_suffix)
            os.makedirs(download_dir, exist_ok=True)

            url_file = f"{ENCODE_BASE_URL}{file_data['href']}"
            filename = os.path.join(download_dir, os.path.basename(url_file))

            if os.path.exists(filename):
                print(f"File exists: {filename}")
            elif args.dry_run:
                print(f"Dry-run: {filename}")
            else:
                print(f"Downloading: {url_file}")
                download_file(url_file, filename)

            rel_file = os.path.relpath(filename, args.dir)
            metadata['files'][file_data['accession']] = {
                'file_type': file_type,
                'file_format': file_format,
                'output_type': output_type,
                'status': status,
                'bio_rep_id': bio_rep_id,
                'pair': pair,
                'rel_file': rel_file
            }
            downloaded = True

        if downloaded:
            exp_dir = os.path.join(args.dir, accession_id)
            os.makedirs(exp_dir, exist_ok=True)
            with open(os.path.join(exp_dir, 'metadata.org.json'), 'w') as fp:
                json.dump(exp_data, fp, indent=4)
            with open(os.path.join(exp_dir, 'metadata.json'), 'w') as fp:
                json.dump(metadata, fp, indent=4)
            all_file_metadata[accession_id] = metadata['files']

    if not args.dry_run and all_file_metadata:
        max_files = max(len(files) for files in all_file_metadata.values())
        header = 'accession\tdescription\t' + '\t'.join(f"file{i+1}" for i in range(max_files)) + '\n'
        rows = []

        for accession_id, files in all_file_metadata.items():
            desc = ','.join(f"{acc_id}:{meta['status']}:{meta['file_type']}:{meta['file_format']}:{meta['output_type']}:{meta['bio_rep_id']}:{meta['pair']}"
                            for acc_id, meta in files.items())
            file_list = '\t'.join(meta['rel_file'] for meta in files.values())
            rows.append(f"{accession_id}\t{desc}\t{file_list}")

        with open(os.path.join(args.dir, 'all_files.tsv'), 'w') as fp:
            fp.write(header)
            fp.write('\n'.join(rows))

if __name__ == '__main__':
    main()
