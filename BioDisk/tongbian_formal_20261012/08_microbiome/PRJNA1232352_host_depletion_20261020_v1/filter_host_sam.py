import sys,gzip,json,hashlib
from pathlib import Path
TRANS=bytes.maketrans(b'ACGTNacgtn',b'TGCANtgcan')
def record(line):
 f=line.rstrip(b'\r\n').split(b'\t')
 assert len(f)>=11,'SAM record has fewer than 11 fields'
 flag=int(f[1]);assert flag&1,'Unpaired SAM record'
 assert not flag&0x900,'Unexpected secondary/supplementary alignment'
 assert bool(flag&64)!=bool(flag&128),'Invalid mate flag'
 return f,flag
def filter_stream(stream,r1path,r2path):
 counts={'total_pairs':0,'both_unmapped_pairs':0,'one_mapped_pairs':0,'both_mapped_pairs':0}
 ids=hashlib.sha256();pending=None
 with gzip.open(r1path,'wb',compresslevel=3) as o1,gzip.open(r2path,'wb',compresslevel=3) as o2:
  for line in stream:
   if line.startswith(b'@'):continue
   f,flag=record(line)
   if pending is None:pending=(f,flag);continue
   g,gflag=pending;pending=None
   assert f[0]==g[0],'Mate names not adjacent/equal'
   assert bool(flag&64)!=bool(gflag&64),'Duplicate mate'
   a,af,b,bf=(f,flag,g,gflag) if flag&64 else (g,gflag,f,flag)
   assert bool(af&8)==bool(bf&4) and bool(bf&8)==bool(af&4),'Inconsistent mate-mapping flags'
   counts['total_pairs']+=1
   nmap=int(not af&4)+int(not bf&4)
   counts[['both_unmapped_pairs','one_mapped_pairs','both_mapped_pairs'][nmap]]+=1
   if nmap==0:
    for rec,flg,out,suffix in [(a,af,o1,b'/1'),(b,bf,o2,b'/2')]:
     seq,qual=rec[9],rec[10]
     assert seq!=b'*' and qual!=b'*' and len(seq)==len(qual),'Missing sequence/quality'
     if flg&16:seq=seq.translate(TRANS)[::-1];qual=qual[::-1]
     out.write(b'@'+rec[0]+suffix+b'\n'+seq+b'\n+\n'+qual+b'\n')
    ids.update(a[0]+b'\n')
 assert pending is None,'Truncated pair at SAM EOF'
 assert counts['total_pairs']==counts['both_unmapped_pairs']+counts['one_mapped_pairs']+counts['both_mapped_pairs']
 counts['retained_pair_id_sha256']=ids.hexdigest()
 return counts
if __name__=='__main__':
 counts=filter_stream(sys.stdin.buffer,sys.argv[1],sys.argv[2])
 Path(sys.argv[3]).write_text(json.dumps(counts,indent=2),encoding='utf-8')
