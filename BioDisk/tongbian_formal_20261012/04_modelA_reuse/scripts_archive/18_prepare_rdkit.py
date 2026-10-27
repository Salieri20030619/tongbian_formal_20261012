from direct_http import *
import shutil,subprocess,sys
rec,resp=get('https://pypi.org/pypi/rdkit/json','pypi_rdkit_metadata')
if resp is None or resp.status_code!=200:raise SystemExit('Official PyPI metadata unavailable')
data=resp.json();wheels=[x for x in data['urls'] if 'cp312' in x['filename'] and 'win_amd64' in x['filename'] and x['filename'].endswith('.whl')]
if len(wheels)!=1:raise SystemExit('No unique compatible wheel')
wheel=wheels[0];rec,resp=get(wheel['url'],'rdkit_official_wheel',timeout=(12,60))
if resp is None or resp.status_code!=200:raise SystemExit('Wheel download failed')
if hashlib.sha256(resp.content).hexdigest()!=wheel['digests']['sha256']:raise SystemExit('PyPI hash mismatch')
wh=ROOT/'wheelhouse';wh.mkdir(exist_ok=True);path=wh/wheel['filename'];path.write_bytes(resp.content)
(ROOT/'logs'/'rdkit_install_source.json').write_text(json.dumps({'name':'rdkit','version':data['info']['version'],'filename':wheel['filename'],'url':wheel['url'],'sha256':wheel['digests']['sha256'],'install_scope':'workspace/vendor; no system installation'},indent=2),encoding='utf-8')
subprocess.run([sys.executable,'-m','pip','install','--no-index','--no-deps','--no-cache-dir','--target',str(ROOT/'vendor'),str(path)],check=True)
sys.path.insert(0,str(ROOT/'vendor'))
from rdkit import Chem,rdBase
print('rdkit',rdBase.rdkitVersion,Chem.MolToInchiKey(Chem.MolFromSmiles('CCO')),flush=True)
