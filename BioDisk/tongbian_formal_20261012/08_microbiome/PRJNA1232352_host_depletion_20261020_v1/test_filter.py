import importlib.util,io,tempfile,pathlib,gzip,json
p=pathlib.Path(__file__).parent
spec=importlib.util.spec_from_file_location('pairfilter',p/'filter_host_sam.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
def row(name,flag,seq=b'ACGTT',qual=b'ABCDE'):
 return b'\t'.join([name,str(flag).encode(),b'*',b'0',b'0',b'*',b'*',b'0',b'0',seq,qual])+b'\n'
with tempfile.TemporaryDirectory(dir=p) as td:
 d=pathlib.Path(td);a=d/'r1.gz';b=d/'r2.gz'
 sam=row(b'keep',77)+row(b'keep',141)+row(b'drop_one',73)+row(b'drop_one',133)+row(b'drop_both',99)+row(b'drop_both',147)+row(b'reverse',93)+row(b'reverse',173)
 c=m.filter_stream(io.BytesIO(sam),a,b)
 assert c['total_pairs']==4 and c['both_unmapped_pairs']==2 and c['one_mapped_pairs']==1 and c['both_mapped_pairs']==1
 x=gzip.open(a,'rb').read();y=gzip.open(b,'rb').read()
 assert b'@keep/1\nACGTT\n+\nABCDE' in x and b'@reverse/1\nAACGT\n+\nEDCBA' in x
 assert b'drop' not in x and b'drop' not in y
 for i,bad in enumerate([row(b'x',77),row(b'x',77)+row(b'y',141),row(b'x',77)+row(b'x',77)]):
  try:m.filter_stream(io.BytesIO(bad),d/(str(i)+'a.gz'),d/(str(i)+'b.gz'))
  except AssertionError:pass
  else:raise AssertionError('Malformed input not rejected')
(p/'FILTER_TESTS_PASS.json').write_text(json.dumps({'status':'PASS','cases':['both unmapped retained','one mapped pair removed','both mapped removed','reverse complement restored','truncated pair rejected','mismatched IDs rejected','duplicate mate rejected']}),encoding='utf-8')
print('FILTER_TESTS_PASS')
