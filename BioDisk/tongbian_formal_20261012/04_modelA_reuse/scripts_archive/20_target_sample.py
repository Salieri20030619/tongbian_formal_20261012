from direct_http import *
for label,body in [('herb_target_browse_sample',{'func_name':'browse_api','label':'Target','page':1,'page_size':50}),('herb_target_detail_sample',{'func_name':'detail_api','label':'Target','v':'HBTAR003295','key_id':'HBTAR003295'})]:
 r,resp=get('http://herb.ac.cn/chedi/api/',label,method='POST',json_body=body,timeout=(12,35))
 if resp is not None:
  d=resp.json();print(label,json.dumps({k:(v[:3] if isinstance(v,list) else v) for k,v in d.items() if k in ['summary','browse_data','total_num']},ensure_ascii=False),flush=True)
