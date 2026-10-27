import subprocess,pathlib,os
O=pathlib.Path(__file__).parent
with (O/'statistics.log').open('w') as f:
 env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
 p=subprocess.run(['/mnt/BioData/bio/wangmingyuan/tongbian_formal_20261012/software/envs/maaslin2_20261020_v1/bin/Rscript','--vanilla',str(O/'fit_associations.R'),str(O)],stdout=f,stderr=subprocess.STDOUT,env=env)
(O/'statistics_exit_code.txt').write_text(str(p.returncode)+'\n')
