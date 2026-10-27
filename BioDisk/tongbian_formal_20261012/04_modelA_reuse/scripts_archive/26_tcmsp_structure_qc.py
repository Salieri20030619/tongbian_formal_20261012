from pathlib import Path
import sys,json,re,pandas as pd
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'vendor'))
from rdkit import Chem,RDLogger
from rdkit.Chem import rdMolDescriptors
RDLogger.DisableLog('rdApp.warning')
src=pd.read_csv(ROOT/'tables'/'11_tcmsp_structure_identifiers.tsv',sep='\t',dtype=str).fillna('')
props=pd.read_csv(ROOT/'tables'/'22_pubchem_structure_properties.tsv',sep='\t',dtype=str).fillna('')
keymap={k:g.to_dict('records') for k,g in props.groupby('InChIKey')};rows=[]
for s in src.to_dict('records'):
 out=dict(s);key=s.get('InChIKey','');candidates=keymap.get(key,[])
 out.update(structure_status='unresolved',structure_source='',canonical_isomeric_smiles='',computed_inchikey='',identity_rule='exact_full_InChIKey agreement required',network_structure_eligible=False)
 mol=None
 if candidates:
  rec=sorted(candidates,key=lambda x:int(x['CID']))[0];mol=Chem.MolFromSmiles(rec['SMILES'])
  out.update(structure_source='PubChem',structure_source_raw=rec['retrieval_raw'],verified_cids=';'.join(sorted({x['CID'] for x in candidates},key=int)),source_cid_exact_key_agreement=bool(set(s.get('pubchem_cids','').split(';'))&{x['CID'] for x in candidates}))
 else:
  p=ROOT/'raw'/f"mol2_{s['MOL_ID']}.bin"
  if p.exists():
   try:mol=Chem.MolFromMol2Block(p.read_text(encoding='utf-8'),sanitize=True,removeHs=True)
   except Exception:mol=None
   out.update(structure_source='TCMSP_MOL2',structure_source_raw=str(p.relative_to(ROOT)))
 if mol is not None:
  computed=Chem.MolToInchiKey(mol);frags=len(Chem.GetMolFrags(mol));wild=any(a.GetAtomicNum()==0 for a in mol.GetAtoms())
  out.update(canonical_isomeric_smiles=Chem.MolToSmiles(mol,isomericSmiles=True),computed_inchikey=computed,computed_formula=rdMolDescriptors.CalcMolFormula(mol),fragment_count=frags,unassigned_stereocenters=sum(x[1]=='?' for x in Chem.FindMolChiralCenters(mol,includeUnassigned=True)),has_wildcard=wild)
  out['structure_status']='exact_key_confirmed' if computed==key else 'structure_key_conflict'
  out['network_structure_eligible']=computed==key and frags==1 and not wild
 elif out['structure_source']:out['structure_status']='structure_parse_failed'
 rows.append(out)
df=pd.DataFrame(rows);df.to_csv(ROOT/'tables'/'26_tcmsp_structure_qc.tsv',sep='\t',index=False,encoding='utf-8-sig')
summary={'source_mol_ids':len(df),'status_counts':df.structure_status.value_counts().to_dict(),'single_structure_eligible_records':int(df.network_structure_eligible.sum()),'unique_eligible_full_inchikey':df.loc[df.network_structure_eligible,'computed_inchikey'].nunique(),'unassigned_stereochemistry_retained_as_annotation':int((pd.to_numeric(df.unassigned_stereocenters,errors='coerce')>0).sum())}
(ROOT/'results'/'tcmsp_structure_qc_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(summary,ensure_ascii=False))
