from pathlib import Path
import pandas as pd,sys,json,re
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'vendor'))
from rdkit import Chem,RDLogger
from rdkit.Chem import rdMolDescriptors
RDLogger.DisableLog('rdApp.warning')
src=pd.read_csv(ROOT/'tables'/'17_herb_ingredient_details.tsv',sep='\t',dtype=str).fillna('')
cache={p.stem:json.loads(p.read_text(encoding='utf-8')) for p in (ROOT/'raw'/'pubchem_properties').glob('*.json')}
rows=[]
for r in src.to_dict('records'):
 o=dict(r);o.update(structure_status='unresolved',network_structure_eligible=False,structure_source='',canonical_isomeric_smiles='',computed_inchikey='')
 smi=r.get('Molecule smile','');smol=None
 if smi and smi.lower() not in ['na','nan','none']:
  try:smol=Chem.MolFromSmiles(smi)
  except Exception:pass
 skey=Chem.MolToInchiKey(smol) if smol is not None else ''
 o['source_smiles_inchikey']=skey
 assigned=sum(t!='?' for _,t in Chem.FindMolChiralCenters(smol,includeUnassigned=True)) if smol is not None else 0
 ids=re.findall(r'\b\d+\b',r.get('PubChem id',''));candidates=[cache[c] for c in ids if c in cache]
 distinct={c['InChIKey'] for c in candidates};mol=smol
 if len(distinct)==1:
  pc=candidates[0];pkey=pc['InChIKey'];o['referenced_pubchem_cids']=';'.join(ids)
  o['pubchem_inchikey']=pkey;o['pubchem_source_raw']=pc['retrieval_raw']
  if not skey or skey==pkey or (skey[:14]==pkey[:14] and assigned==0):
   mol=Chem.MolFromSmiles(pc['SMILES']);o['structure_source']='explicit_HERB_PubChem_CID';o['structure_status']='cid_structure_consistent' if skey==pkey else ('cid_resolves_unspecified_stereochemistry' if skey else 'cid_only_structure')
  else:o['structure_status']='source_smiles_CID_identity_conflict'
 elif len(distinct)>1:o['structure_status']='multiple_CID_structure_conflict'
 elif smol is not None:o['structure_source']='HERB_source_SMILES';o['structure_status']='source_SMILES_only'
 elif ids:o['structure_status']='CID_not_yet_resolved'
 if mol is not None:
  o.update(canonical_isomeric_smiles=Chem.MolToSmiles(mol,isomericSmiles=True),computed_inchikey=Chem.MolToInchiKey(mol),computed_formula=rdMolDescriptors.CalcMolFormula(mol),fragment_count=len(Chem.GetMolFrags(mol)),has_wildcard=any(a.GetAtomicNum()==0 for a in mol.GetAtoms()),unassigned_stereocenters=sum(t=='?' for _,t in Chem.FindMolChiralCenters(mol,includeUnassigned=True)))
  o['source_formula_agrees']=not r.get('Molecule formula') or r.get('Molecule formula')==o['computed_formula']
  o['network_structure_eligible']=o['structure_status'] in ['cid_structure_consistent','cid_resolves_unspecified_stereochemistry','cid_only_structure','source_SMILES_only'] and o['fragment_count']==1 and not o['has_wildcard'] and bool(o['computed_inchikey']) and o['source_formula_agrees']
 rows.append(o)
df=pd.DataFrame(rows);df.to_csv(ROOT/'tables'/'27_herb_structure_qc.tsv',sep='\t',index=False,encoding='utf-8-sig')
summary={'downloaded_records':len(df),'status_counts':df.structure_status.value_counts().to_dict(),'eligible_single_structure_records':int(df.network_structure_eligible.sum()),'eligible_unique_full_inchikey':df.loc[df.network_structure_eligible,'computed_inchikey'].nunique()}
(ROOT/'results'/'herb_structure_qc_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(summary,ensure_ascii=False))
