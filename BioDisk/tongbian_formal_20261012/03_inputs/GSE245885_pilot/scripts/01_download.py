#!/usr/bin/env python3
"""Download GEO public processed matrices only; never FASTQ.
Future reruns explicitly disable environment proxies and will fail closed if direct DNS/access is unavailable.
Historical executed download used platform environment defaults before the user no-proxy constraint.
"""
import urllib.request, pathlib, re, hashlib, json, datetime, concurrent.futures, csv
ROOT=pathlib.Path(__file__).resolve().parents[1]
raw=ROOT/'raw'; raw.mkdir(exist_ok=True)
records=[]
def fetch(pair):
    url,name=pair
    path=raw/name
    t=datetime.datetime.now(datetime.timezone.utc).isoformat()
    for attempt in range(3):
        try:
            with urllib.request.build_opener(urllib.request.ProxyHandler({})).open(url,timeout=90) as r:
                data=r.read()
            path.write_bytes(data)
            break
        except Exception:
            if attempt==2: raise
    entry={'url':url,'local_file':str(path.relative_to(ROOT)),'downloaded_utc':t,'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()}
    print(json.dumps(entry),flush=True)
    return entry
series_url='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE245885&targ=self&form=text&view=full'
family_url='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE245nnn/GSE245885/soft/GSE245885_family.soft.gz'
import gzip
records.append(fetch((series_url,'GSE245885_series.soft.txt')))
records.append(fetch((family_url,'GSE245885_family.soft.gz')))
soft=gzip.decompress((raw/'GSE245885_family.soft.gz').read_bytes()).decode()
(raw/'GSE245885_family.soft').write_text(soft)
samples=[]; urls=[]
for block in soft.split('^SAMPLE = ')[1:]:
    gsm=block.splitlines()[0].strip(); tags={}
    for line in block.splitlines()[1:]:
        if line.startswith('^'): break
        if line.startswith('!') and ' = ' in line:
            k,v=line.split(' = ',1);tags.setdefault(k,[]).append(v)
    title=tags['!Sample_title'][0]
    characteristics=tags.get('!Sample_characteristics_ch1',[])
    treatment=[c.split(': ',1)[1] for c in characteristics if c.startswith('treatment:')][0]
    group={'control':'Control','STC':'STC'}[treatment]
    files=[v for k,vs in tags.items() if k.startswith('!Sample_supplementary_file') for v in vs if '.genes.FPKM.' in v]
    assert len(files)==1
    url=files[0].replace('ftp://','https://');filename=url.rsplit('/',1)[-1]
    samples.append({'gsm':gsm,'sample':title,'condition':group,'tissue':tags['!Sample_source_name_ch1'][0],'characteristics':'; '.join(characteristics),'sex':'not_reported','age':'not_reported','colon_segment':'not_reported','batch':'not_reported','supplement_file':filename,'supplement_url':url})
    urls.append((url,filename))
assert len(samples)==10 and len(set(s['gsm'] for s in samples))==10
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex: records.extend(ex.map(fetch,urls))
(ROOT/'metadata'/'download_manifest.json').write_text(json.dumps(records,indent=2))
with (ROOT/'metadata'/'samples.tsv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=list(samples[0]),delimiter='\t');w.writeheader();w.writerows(samples)
print('Downloaded',len(samples),'gene matrices. Expected 5 Control + 5 STC.',flush=True)
