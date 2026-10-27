"""Every network request is direct, TLS-verified, logged, and never retried via a proxy."""
import os
for _key in list(os.environ):
    if _key.lower() in {'http_proxy','https_proxy','all_proxy'}:
        os.environ.pop(_key, None)
os.environ['NO_PROXY']='*'
os.environ['no_proxy']='*'
import requests, hashlib, json, time, threading
from pathlib import Path
from datetime import datetime, timezone
ROOT=Path(__file__).resolve().parents[1]
LOCK=threading.Lock()
for _d in ['raw','logs','tables','figures','results']:
    (ROOT/_d).mkdir(exist_ok=True)

def new_session():
    session=requests.Session()
    session.trust_env=False
    session.proxies={}
    session.headers['User-Agent']='ResearchPreflight/1.0 (noncommercial reproducible research; no proxy)'
    return session

def get(url, label, timeout=(12,25), params=None, session=None, method='GET',json_body=None,data_body=None):
    owned_session=session is None
    session=session or new_session()
    started=time.monotonic()
    record={'time_utc':datetime.now(timezone.utc).isoformat(),'label':label,'url':url,'params':params,'method':method,'json_body':json_body,'proxy_mode':'disabled; trust_env=False; NO_PROXY=*','tls_verify':True}
    if data_body is not None:record['form_body']=data_body
    response=None
    try:
        response=session.request(method,url,params=params,json=json_body,data=data_body,timeout=timeout,verify=True,proxies={})
        data=response.content
        path=ROOT/'raw'/f'{label}.bin'
        path.write_bytes(data)
        record.update(status=response.status_code,final_url=response.url,content_type=response.headers.get('Content-Type',''),bytes=len(data),sha256=hashlib.sha256(data).hexdigest(),raw_file=str(path.relative_to(ROOT)))
        record['source_headers']={k:response.headers[k] for k in ['Date','Last-Modified','X-UniProt-Release','X-UniProt-Release-Date','X-Throttling-Control'] if k in response.headers}
    except requests.RequestException as exc:
        record.update(status=None,error=str(exc))
    finally:
        if owned_session:session.close()
    record['elapsed_seconds']=round(time.monotonic()-started,3)
    with LOCK:
        with (ROOT/'logs'/'http_requests.jsonl').open('a',encoding='utf-8') as fh:
            fh.write(json.dumps(record,ensure_ascii=False)+'\n')
    return record,response
