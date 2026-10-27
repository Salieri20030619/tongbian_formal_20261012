from pathlib import Path
import subprocess,json,urllib.parse,datetime,hashlib,xml.etree.ElementTree as ET
P=Path('/mnt/BioDisk/wangmingyuan/tongbian_formal_20261012');O=P/'03_inputs/references/Ensembl_human_mouse_one2one_20261013';O.mkdir(exist_ok=True)
xml='''<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1" count="" datasetConfigVersion="0.6"><Dataset name="hsapiens_gene_ensembl" interface="default"><Attribute name="ensembl_gene_id"/><Attribute name="external_gene_name"/><Attribute name="mmusculus_homolog_ensembl_gene"/><Attribute name="mmusculus_homolog_associated_gene_name"/><Attribute name="mmusculus_homolog_orthology_type"/><Attribute name="mmusculus_homolog_orthology_confidence"/></Dataset></Query>'''
(O/'query.xml').write_text(xml)
receipts=[]
for host in ['https://www.ensembl.org','https://useast.ensembl.org']:
 url=host+'/biomart/martservice?'+urllib.parse.urlencode({'query':xml});p=O/'human_mouse_all_homology.tsv'
 if p.exists():break
 r=subprocess.run(['curl','--noproxy','*','-fLsS','--connect-timeout','15','--max-time','120',url],capture_output=True)
 valid=r.returncode==0 and r.stdout.startswith(b'Gene stable ID') and b'ortholog_one2one' in r.stdout
 receipts.append({'url':url,'exit_code':r.returncode,'bytes':len(r.stdout),'valid':valid,'error':r.stderr.decode(errors='replace')[:1000]})
 if valid:p.write_bytes(r.stdout)
 else:(O/('failed_'+urllib.parse.urlsplit(host).hostname+'.txt')).write_bytes(r.stdout[:2000])
r=subprocess.run(['curl','--noproxy','*','-fLsS','--max-time','30','https://rest.ensembl.org/info/data?content-type=application/json'],capture_output=True)
if r.returncode==0:(O/'ensembl_release.json').write_bytes(r.stdout)
(O/'retrieval_receipt.json').write_text(json.dumps(receipts,indent=2));print(json.dumps(receipts,indent=2))
if (O/'human_mouse_all_homology.tsv').exists():print((O/'human_mouse_all_homology.tsv').read_text()[:600])
