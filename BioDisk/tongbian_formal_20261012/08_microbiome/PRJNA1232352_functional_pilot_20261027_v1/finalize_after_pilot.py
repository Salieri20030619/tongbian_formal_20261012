import pathlib,json,hashlib,csv,re,subprocess,os,select,datetime
N=pathlib.Path(__file__).parent;L=pathlib.Path(json.loads((N/'plan_frozen.json').read_text())['large_output'])
if not (N/'exit_code.txt').exists():
 pid=int((N/'pid.txt').read_text())
 try:
  fd=os.pidfd_open(pid);select.select([fd],[],[]);os.close(fd)
 except ProcessLookupError:pass
assert (N/'exit_code.txt').read_text().strip()=='0'
p=json.loads((N/'PILOT_COMPLETE.json').read_text());assert p['samples']==3 and p['status']=='TECHNICAL_PILOT_PASS'
def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
checks=0;support={};resources=[]
for x in p['results']:
 o=L/x['run'];assert json.loads((o/'VALIDATED.json').read_text())==x
 for f,h in x['outputs'].items():assert (o/f).stat().st_size==h['bytes'] and sha(o/f)==h['sha256'];checks+=1
 rows=list(csv.DictReader((o/'dbcan/overview.tsv').open(),delimiter='\t'));assert len(rows)==x['dbcan_overview_rows'];support[x['run']]={'any_method':len(rows),'two_methods':sum(int(y['#ofTools'])>=2 for y in rows)}
 for label in ['megahit','prodigal','dbcan']:
  assert (o/(label+'_exit_code.txt')).read_text().strip()=='0';log=(o/(label+'.log')).read_text();rss=re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)',log);wall=re.search(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*(\S+)',log);resources.append([x['run'],label,rss.group(1) if rss else 'NA',wall.group(1) if wall else 'NA'])
with (N/'pilot_resource_summary.tsv').open('w',newline='') as f:
 w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['run','stage','max_RSS_kbytes_GNU_time','elapsed_wall']);w.writerows(resources)
(N/'PILOT_OUTPUT_QA.json').write_text(json.dumps({'status':'PASS','utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'files_rehashed':checks,'stage_exit_codes_checked':9,'method_support':support,'full_sample_peak_validated':False},indent=2)+'\n')
script=N/'make_update_report.py';text=script.read_text();text=text.replace("rows=[['组别 / SRR','读对数','contigs','预测蛋白','注释记录']]","support=json.loads((N/'PILOT_OUTPUT_QA.json').read_text())['method_support']\nrows=[['组别 / SRR','读对数','contigs','预测蛋白','注释记录','双方法']]")
text=text.replace("str(x['dbcan_overview_rows'])])","str(x['dbcan_overview_rows']),str(support[x['run']]['two_methods'])])").replace('tab(rows,[150,80,75,85,75])','tab(rows,[135,70,65,70,65,60])');script.write_text(text,encoding='utf-8')
py='/mnt/BioData/bio/wangmingyuan/tongbian_formal_20261012/software/envs/biobakery_20261020_v1/bin/python'
with (N/'report_build.log').open('w') as f:subprocess.run([py,str(script)],stdout=f,stderr=subprocess.STDOUT,check=True)
q=N/'render_qa';q.mkdir(exist_ok=True);report=N/'tongbian_covariate_functional_update_20261027.pdf';subprocess.run(['pdftoppm','-scale-to','1400','-png',str(report),str(q/'page')],check=True)
from PIL import Image
files=sorted(q.glob('page-*.png')); sheet=Image.new('RGB',(1800,850*((len(files)+2)//3)),'#cccccc')
for i,f in enumerate(files):
 im=Image.open(f);im.thumbnail((580,820));sheet.paste(im,((i%3)*600+10,(i//3)*850+15))
sheet.save(q/'contact.png');print(json.dumps({'status':'READY_FOR_VISUAL_REVIEW','pages':len(files),'support':support,'report_sha256':sha(report)}),flush=True)
