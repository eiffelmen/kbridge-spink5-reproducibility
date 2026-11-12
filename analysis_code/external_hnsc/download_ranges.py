from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import urllib.request, time, json
base=Path(__file__).resolve().parent
url='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE65nnn/GSE65858/matrix/GSE65858_series_matrix.txt.gz'
N=21329581
op=urllib.request.build_opener(urllib.request.ProxyHandler({}))
def one(i):
 start=i*N//8;stop=(i+1)*N//8-1;p=base/f'matrix.part{i}'
 for attempt in range(4):
  have=p.stat().st_size if p.exists() else 0
  if have==stop-start+1:return str(p)
  req=urllib.request.Request(url,headers={'Range':f'bytes={start+have}-{stop}'})
  try:
   with op.open(req,timeout=20) as r:
    if r.status!=206 or not r.headers.get('Content-Range','').startswith(f'bytes {start+have}-'):raise ValueError('Invalid range response')
    with p.open('ab') as f:
     while True:
      b=r.read(65536)
      if not b:break
      f.write(b)
      if f.tell()>stop-start+1:raise ValueError('Range exceeded')
   print('part',i,'downloaded',p.stat().st_size,flush=True)
  except Exception as e:print('part',i,'attempt',attempt,'error',str(e)[:100],flush=True)
 if p.stat().st_size!=stop-start+1:raise ValueError('Incomplete part')
 return str(p)
with ThreadPoolExecutor(max_workers=8) as ex:
 for f in as_completed([ex.submit(one,i) for i in range(8)]):f.result()
with (base/'GSE65858_series_matrix.txt.gz').open('wb') as out:
 for i in range(8):out.write((base/f'matrix.part{i}').read_bytes())
print('joined',N,flush=True)
