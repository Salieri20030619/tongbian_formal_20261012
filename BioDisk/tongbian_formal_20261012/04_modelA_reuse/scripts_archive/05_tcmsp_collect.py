from direct_http import *
import re,csv,pandas as pd
from bs4 import BeautifulSoup
catalog=json.loads(re.search(r'data:\s*(\[.*?\])\s*,',(ROOT/'raw'/'tcmsp_browse.bin').read_text(encoding='utf-8'),re.S).group(1))
pd.DataFrame(catalog).to_csv(ROOT/'tables'/'03_tcmsp_catalog.tsv',sep='\t',index=False,encoding='utf-8-sig')
formula=pd.read_csv(ROOT/'tables'/'02_formula_standardization.tsv',sep='\t')
coverage=[];all_arrays={};ingredients=[];targets=[]
session=new_session()
r,resp=get('https://old.tcmsp-e.com/tcmsp.php','tcmsp_collect_session',session=session)
token=BeautifulSoup(resp.content,'html.parser').select_one('input[name=token]')['value']
for herb in formula.to_dict('records'):
 matches=[x for x in catalog if x['herb_cn_name']==herb['search_drug_name']]
 cr={'herb_id':herb['herb_id'],'prescription_label':herb['prescription_label'],'query_name':herb['search_drug_name'],'catalog_records':len(catalog),'exact_name_matches':len(matches),'source':'TCMSP old v2.3','collection_status':'not_in_exact_catalog' if not matches else 'pending','ingredient_rows':None,'target_rows':None,'scope_status':'药材名匹配不等于炮制/部位/基原全部核验'}
 if not matches:
    coverage.append(cr);continue
 for match in matches:
    label='tcmsp_'+herb['herb_id']
    r,resp=get('https://old.tcmsp-e.com/tcmspsearch.php',label,params={'qr':match['herb_en_name'],'qsr':'herb_en_name','token':token},session=session)
    cr.update(source_latin_name=match['herb_en_name'],source_url=r.get('final_url'),raw_file=r.get('raw_file'))
    if resp is None or resp.status_code!=200:
        cr['collection_status']='http_failure';continue
    arrays=[]
    # Each public grid includes the entire JSON array; pagination is client-side.
    for m in re.finditer(r'\bdata:\s*(?=\[)',resp.text):
        try:arrays.append(json.JSONDecoder().raw_decode(resp.text[m.end():])[0])
        except json.JSONDecodeError:pass
    all_arrays[herb['herb_id']]=arrays
    print(herb['herb_id'],herb['prescription_label'],[(len(a),list(a[0]) if a else []) for a in arrays],flush=True)
    cr['collection_status']='public_grid_downloaded'
    for a in arrays:
      if a and 'molecule_name' in a[0] and 'ob' in a[0]:
        cr['ingredient_rows']=len(a)
        ingredients.extend([dict(x,herb_id=herb['herb_id'],prescription_label=herb['prescription_label'],source_url=r['final_url'],source_raw=r['raw_file']) for x in a])
      elif a and 'target_name' in a[0] and 'molecule_ID' in a[0]:
        cr['target_rows']=len(a)
        targets.extend([dict(x,herb_id=herb['herb_id'],source_url=r['final_url'],source_raw=r['raw_file']) for x in a])
 coverage.append(cr)
 time.sleep(.5)
(ROOT/'results'/'tcmsp_grid_arrays.json').write_text(json.dumps(all_arrays,ensure_ascii=False,indent=2),encoding='utf-8')
for name,data in [('04_tcmsp_herb_coverage',coverage),('05_tcmsp_ingredient_rows',ingredients),('06_tcmsp_target_rows_UNCLASSIFIED',targets)]:
 pd.DataFrame(data).to_csv(ROOT/'tables'/f'{name}.tsv',sep='\t',index=False,encoding='utf-8-sig')
print('TOTAL',len(ingredients),len(targets),flush=True)
