from pathlib import Path
import json,hashlib,zipfile,datetime,csv
r=Path('/mnt/BioDisk/wangmingyuan/tongbian_formal_20261012')
out=r/'12_reports/first_batch_20261012'
assert json.loads((out/'first_batch_completion.json').read_text())['status']=='COMPLETED'
readme='''GSE245885 first formal batch: reproduction and audit entry points

This bundle contains all first-batch model objects/results, full frozen reference mappings, scripts, logs, and module sessions. It does not duplicate the original ten raw input tables or legacy Model A. Their exact server paths, sizes and SHA256 are in 01_inventory/staged_input_manifest.json, with prior validation in server_staged_input_verification.json.

Inputs: /mnt/BioDisk/wangmingyuan/tongbian_formal_20261012/03_inputs/GSE245885_pilot/raw
Runtime: /mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript --vanilla
Execution order: GSE245885_R_DESeq2_bridge.R --run; audit_DE_bridge.R; run_host_functional_scores.R (sources prepare_functional_inputs.R); run_fgsea_pathways.R. Sample QC and technical validation read saved outputs. No new software is required.

All original jobs are complete. Do not resubmit them. The launch.lock, spec.json and status.json in each logs/jobs directory prevent duplicate launch. Bridge refuses to overwrite its output directory. For an intentional independent reproduction, prepare a new isolated project root, point the scripts' explicit root constant to it, stage the same SHA-verified inputs and references, and retain the frozen design/parameters/seed and package versions. Do not run copied scripts against the archived output root.

Unchanged grouping: NP1..NP5 cancer-surgery normal-colon controls; STC1..STC5 STC. Observational host evidence only. No formula treatment contrast, sample exclusion or added FC gate.

The final figure directory under 12_reports contains label/layout corrections only; original figures remain beside their numeric outputs. Tables named ALL are full inferential universes; Top20 tables/figures are display subsets.

The payload manifest excludes the two manifest files themselves. The ZIP receipt supplies the whole-bundle checksum. Windows transfer verification is supplied separately in windows_sync_receipt.json; it is not a statistical recomputation.
'''
(out/'REPRODUCIBILITY_README.txt').write_text(readme)
include=[]
dirs=['05_human_RNA','12_reports/first_batch_20261012','12_reports/first_batch_figures','03_inputs/references/msigdbr_26.1.0_20261012']
for d in dirs:include.extend(p for p in (r/d).rglob('*') if p.is_file() and p.name not in ['deliverable_manifest.json','deliverable_manifest.tsv'])
extra=['03_inputs/references/first_batch_reference_manifest.json','03_inputs/references/dorothea_hs_ABC_v1.22.0.tsv','03_inputs/GSE245885_pilot/metadata/samples.tsv','03_inputs/GSE245885_pilot/genesets/host_modules_frozen.gmt','03_inputs/GSE245885_pilot/genesets/host_modules_frozen.json','03_inputs/GSE245885_pilot/genesets/README.md','00_protocol/functional_plan.json','00_protocol/first_batch_script_freeze.json','00_protocol/project_config.json','00_protocol/execution_state.json','00_protocol/module_status.tsv','00_protocol/blockers.tsv','environment/formal_function_interfaces.txt']
include.extend(r/p for p in extra if (r/p).is_file())
include.extend(r/p for p in ['01_inventory/staged_input_manifest.json','01_inventory/server_staged_input_verification.json'] if (r/p).is_file())
scripts=['GSE245885_R_DESeq2_bridge.R','audit_DE_bridge.R','freeze_msigdb_existing_package.R','prepare_functional_inputs.R','run_host_functional_scores.R','run_fgsea_pathways.R','sample_quality_audit.R','verify_first_batch.R','inspect_score_internals.R','final_figures.R','durable_R_job_runner.py','launch_formal_de_bridge.py','start_remaining_batch.py','start_fgsea.py','start_validation.py','poll_first_batch.py','inspect_completed_outputs.py','audit_output_integrity.py','finish_first_batch.py','package_first_batch_complete.py']
include.extend(r/'environment'/s for s in scripts if (r/'environment'/s).is_file())
jobs=['GSE245885_R_DESeq2_bridge_v1','MSigDB_package_resource_freeze','GSE245885_host_functional_scores_v1','GSE245885_sample_QC_v1','GSE245885_fgsea_v1','GSE245885_technical_validation_v1']
for n in jobs:
 include.extend(p for p in (r/'logs/jobs'/n).iterdir() if p.is_file())
include=sorted(set(include))
manifest=[{'relative_path':str(p.relative_to(r)),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in include]
mp=out/'deliverable_manifest.json';mp.write_text(json.dumps({'generated_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'root':str(r),'files':manifest,'scope':'complete first-batch outputs, model objects, reference maps, scripts, logs; original ten raw inputs remain under frozen project input root and are not duplicated'},indent=2));include.append(mp)
with (out/'deliverable_manifest.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=['relative_path','bytes','sha256'],delimiter='\t');w.writeheader();w.writerows(manifest)
include.append(out/'deliverable_manifest.tsv')
z=r/'12_reports/GSE245885_first_batch_complete.zip'
with zipfile.ZipFile(z,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as f:
 for p in include:f.write(p,str(p.relative_to(r)))
receipt={'path':str(z),'bytes':z.stat().st_size,'sha256':hashlib.sha256(z.read_bytes()).hexdigest(),'payload_files':len(include),'payload_bytes':sum(p.stat().st_size for p in include),'status':'READY_FOR_WINDOWS_CHECKSUM_VERIFIED_TRANSFER'}
(r/'12_reports/GSE245885_first_batch_bundle_receipt.json').write_text(json.dumps(receipt,indent=2))
print(json.dumps(receipt,indent=2))
