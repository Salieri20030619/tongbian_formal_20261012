from pathlib import Path
import os,signal,json,csv,subprocess,re,collections,hashlib,datetime
P=Path('/mnt/BioDisk/wangmingyuan/tongbian_formal_20261012');O=P/'01_inventory/resume_20261013';I=P/'03_inputs/MTBLS13907_audit_v1'
# Stop only this task's malformed shell probe, whose wildcard was expanded by bash.
killed=[]
for p in Path('/proc').iterdir():
 if not p.name.isdigit():continue
 try:
  if p.stat().st_uid!=os.getuid():continue
  args=(p/'cmdline').read_bytes().split(b'\0');s=b' '.join(args).decode(errors='replace')
  target='https://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/MTBLS13907/FILES/'
  if target in s and '--noproxy' in s and (args[0].endswith(b'bash') or args[0].endswith(b'curl')):
   os.kill(int(p.name),signal.SIGTERM);killed.append(int(p.name))
 except (OSError,IndexError):pass
(O/'malformed_probe_cleanup.json').write_text(json.dumps({'terminated_task_probe_pids':killed,'cause':'PowerShell/ssh quoting removed wildcard quoting; no probe result accepted; use subprocess argument arrays henceforth'},indent=2))
base='https://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/MTBLS13907/FILES/'
for folder in ['DERIVED_FILES','SUPPLEMENTARY_FILES']:
 p=I/(folder+'_index.html');r=subprocess.run(['curl','--noproxy','*','-fLsS','--max-time','45',base+folder+'/'],capture_output=True)
 if not r.returncode:p.write_bytes(r.stdout);print(folder,r.stdout.decode()[:8000])
 else: print(folder,'FAILED',r.stderr.decode())
def rd(n):return list(csv.DictReader((I/n).open(),delimiter='\t'))
s=rd('s_MTBLS13907.txt');a=rd('a_MTBLS13907_LC-MS_alternating_reverse-phase_metabolite_profiling.txt');sm={r['Sample Name']:r for r in s};audit=[]
for row in a:
 sn=row['Sample Name'];an=row['MS Assay Name'];sg=sm.get(sn,{}).get('Factor Value[Group]','UNMATCHED');ag=an.split('_')[0];audit.append({'sample_name':sn,'sample_table_group':sg,'extract_name':row['Extract Name'],'assay_name':an,'assay_prefix_group':ag,'normalization_name':row['Normalization Name'],'transformation_name':row['Data Transformation Name'],'conflict':sg!=ag,'resolved_group':'UNRESOLVED' if sg!=ag else sg})
with (O/'MTBLS13907_sample_assay_conflicts.tsv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(audit[0]),delimiter='\t');w.writeheader();w.writerows(audit)
print('MT_CONFLICTS',json.dumps(audit));summary={'sample_rows':len(s),'assay_rows':len(a),'sample_groups':dict(collections.Counter(r['Factor Value[Group]'] for r in s)),'conflicting_rows':sum(x['conflict'] for x in audit),'metabolite_group_gate':False,'crossomics_pairing_gate':False,'preprocessing':'repository protocol: KNN imputation, total-signal normalization to QC1, QC CV filter already performed; do not repeat'}
(O/'MTBLS13907_group_audit.json').write_text(json.dumps(summary,indent=2));print('MT_SUMMARY',summary)
core=list(csv.DictReader((O/'PRJNA1232352_runs.tsv').open(),delimiter='\t'));print('CORE_GROUPS',collections.Counter(re.sub(r'\d+$','',x['sample_alias']) for x in core))
for n in ['PRJNA1048321','PRJNA1031121','PRJNA1348869']:
 rows=list(csv.DictReader((O/(n+'_runs.tsv')).open(),delimiter='\t'));print(n,'ALIASES',[r['sample_alias'] for r in rows])
