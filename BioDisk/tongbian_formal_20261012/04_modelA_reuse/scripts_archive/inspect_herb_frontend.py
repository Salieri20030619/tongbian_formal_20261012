from pathlib import Path
import re
h=(Path(__file__).resolve().parents[1]/'raw/herb_js.bin').read_text(encoding='utf-8')
for pattern in ['cd=e=>','cd=','search_api','browse_api','init_api','search_type','Object(g["a"])','base_url','ld=']:
 print('\nPATTERN',pattern)
 for m in list(re.finditer(re.escape(pattern),h))[-10:]:print(h[max(0,m.start()-200):m.start()+1000])
