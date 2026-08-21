import sys,time,numpy as np; sys.path.insert(0,'/workspace/code')
import sasp_real as R, sasp_kernels as K, run_figure2a as F
t=time.time(); s=R.load_sample('7250_liver_sham_Male_26-U1'); print('load %.1fs'%(time.time()-t))
sm=s.sender['tierA_p95']; y=s.modules['tnfa_nfkb_proximal']
t=time.time(); d=R.distance_to_set(s.coords,sm); print('dist %.1fs'%(time.time()-t))
keep=(~sm)&np.isfinite(y)&np.isfinite(d); idx=np.flatnonzero(keep)
c=s.coords[idx]; yy=y[idx].astype(float); dd=d[idx].astype(float); X0=np.ones((idx.size,1))
t=time.time(); bid=F.block_ids(c); print('blocks %.1fs n=%d'%(time.time()-t,idx.size))
for f in K.FAMILIES:
    t=time.time(); r=K.fit_family(dd,yy,X0,f,bid,100,n_boot=50,rng=np.random.default_rng(1),p_grid=F.P_GRID_REAL)
    print('%-12s %6.1fs d_half=%6.1f aic=%.0f dAIC=%.0f'%(f,time.time()-t,r['d_half'],r['aic'],r['delta_aic_vs_null']))
