#!/usr/bin/env python3
"""Fetch only clinical metadata, expression header and one expression row.
Uses public UCSC Xena hub; no keys or private data are transmitted.
"""
from pathlib import Path
from urllib.request import Request,urlopen
import json, gzip, hashlib, datetime
P=Path(__file__).resolve().parent/'data'; P.mkdir(exist_ok=True)
HUB='https://tcga.xenahubs.net'
urls={
'HNSC_clinicalMatrix.tsv':HUB+'/download/TCGA.HNSC.sampleMap/HNSC_clinicalMatrix',
'HiSeqV2.json':HUB+'/download/TCGA.HNSC.sampleMap/HiSeqV2.json',
'HNSC_clinicalMatrix.json':HUB+'/download/TCGA.HNSC.sampleMap/HNSC_clinicalMatrix.json'}
for name,url in urls.items():
    if not (P/name).exists():
        (P/name).write_bytes(urlopen(url,timeout=45).read())
# Gzip stream is closed immediately after the header. Entire expression matrix is not downloaded.
if not (P/'expression_sample_ids.json').exists():
    with urlopen(HUB+'/download/TCGA.HNSC.sampleMap/HiSeqV2.gz',timeout=45) as response:
        with gzip.GzipFile(fileobj=response) as stream:
            samples=stream.readline().decode().strip().split('\t')[1:]
    (P/'expression_sample_ids.json').write_text(json.dumps(samples))
samples=json.loads((P/'expression_sample_ids.json').read_text())
query='(fetch [{:table "TCGA.HNSC.sampleMap/HiSeqV2" :columns ["SPINK5"] :samples '+json.dumps(samples)+'}])'
(P/'spink5_all_samples_query.xq').write_text(query)
if not (P/'SPINK5_all_samples.json').exists():
    req=Request(HUB+'/data/',data=query.encode(),headers={'Content-Type':'text/plain'})
    (P/'SPINK5_all_samples.json').write_bytes(urlopen(req,timeout=45).read())
print('Inputs ready:',len(samples),'samples')
