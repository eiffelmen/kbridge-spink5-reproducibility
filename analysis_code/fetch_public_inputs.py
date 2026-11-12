#!/usr/bin/env python3
"""Fetch third-party public inputs from the archives that host them."""
import argparse
import hashlib
from pathlib import Path
import shutil
from tempfile import NamedTemporaryFile
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parent
GEO = 'https://ftp.ncbi.nlm.nih.gov/geo'
XENA = 'https://tcga.xenahubs.net/download/TCGA.HNSC.sampleMap'
MCP = 'https://raw.githubusercontent.com/ebecht/MCPcounter/6b85a671cfa8385798cea6147b6613a172e47ff5/Signatures/genes.txt'
SOURCE_DATA = {
    'kbridge/data/GSE214265_expression.txt.gz': f'{GEO}/series/GSE214nnn/GSE214265/suppl/GSE214265_R1_expression.matrix.txt.gz',
    'kbridge/data/GSE214265_metadata.csv.gz': f'{GEO}/series/GSE214nnn/GSE214265/suppl/GSE214265_R1_metadata.csv.gz',
}
for accession, number in [('GSM4074307', 1), ('GSM4074308', 2)]:
    for suffix in ['endogenes_UMI.csv.gz', 'metadata.xlsx']:
        filename = f'{accession}_Tonsil{number}Bcells_filtered_{suffix}'
        SOURCE_DATA[f'kbridge/external/{filename}'] = f'{GEO}/samples/GSM4074nnn/{accession}/suppl/{filename}'
SOURCE_DATA.update({
    'external_hnsc/GSE65858_series_matrix.txt.gz': f'{GEO}/series/GSE65nnn/GSE65858/matrix/GSE65858_series_matrix.txt.gz',
    'external_hnsc/GPL10558_HumanHT-12_V4_0_R2_15002873_B.txt.gz': f'{GEO}/platforms/GPL10nnn/GPL10558/suppl/GPL10558_HumanHT-12_V4_0_R2_15002873_B.txt.gz',
    'public_hnsc_sensitivity/inputs/mcp_counter_genes.tsv': MCP,
    'public_hnsc_sensitivity/source_inputs/HiSeqV2.gz': f'{XENA}/HiSeqV2.gz',
    'public_hnsc_sensitivity/source_inputs/GSE103322_HNSCC_all_data.txt.gz': f'{GEO}/series/GSE103nnn/GSE103322/suppl/GSE103322_HNSCC_all_data.txt.gz',
})


def digest(path):
    sha = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            sha.update(chunk)
    return sha.hexdigest()


def fetch(relative, url, expected):
    destination = ROOT / relative
    if destination.is_file():
        if expected and digest(destination) != expected:
            raise ValueError(f'Checksum mismatch for existing input {relative}')
        print(f'Already present: {relative}')
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    with NamedTemporaryFile(dir=destination.parent, delete=False) as temporary:
        temp_path = Path(temporary.name)
        try:
            with urlopen(url, timeout=90) as response:
                shutil.copyfileobj(response, temporary)
            temporary.flush()
            if expected and digest(temp_path) != expected:
                raise ValueError(f'Checksum mismatch for {relative} from {url}')
            temp_path.replace(destination)
            print(f'Downloaded: {relative}')
        finally:
            temp_path.unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('module', choices=['kbridge', 'external_hnsc', 'public_hnsc_sensitivity'])
    args = parser.parse_args()
    checksums = {}
    for line in (ROOT / 'source_input_checksums.sha256').read_text().splitlines():
        checksum, relative = line.split(maxsplit=1)
        checksums[relative] = checksum
    for relative, url in SOURCE_DATA.items():
        if relative.startswith(args.module + '/'):
            fetch(relative, url, checksums.get(relative))


if __name__ == '__main__':
    main()
