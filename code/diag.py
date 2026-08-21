import time, numpy as np
from sasp_sim import TissueConfig, simulate_tissue, ripley_ratio
from sasp_estimators import analyze
def run(tag,cfg,seeds=range(1,7),nb=300):
    t=time.time(); L=[analyze(simulate_tissue(cfg,[7,s]),[8,s],n_boot=nb) for s in seeds]
    m=lambda k: np.nanmean([r[k] for r in L]); sd=lambda k: np.nanstd([r[k] for r in L])
    print(f"{tag}")
    print(f"   lam  naive {m('lam_naive'):6.1f}+-{sd('lam_naive'):4.1f} | bin {m('lam_naive_bin'):6.1f} | nuis {m('lam_nuis'):6.1f} | decoy {m('lam_decoy'):6.1f}+-{sd('lam_decoy'):4.1f} | dec+nuis {m('lam_decoy_nuis'):6.1f}")
    print(f"   beta naive {m('beta_naive'):6.3f} | nuis {m('beta_nuis'):6.3f} | decoy {m('beta_decoy'):6.3f} | dec+nuis {m('beta_decoy_nuis'):6.3f} | b_t-b_d {m('beta_true_minus_decoy'):6.3f}")
    print(f"   cover naive_iid {m('cover_lam_naive_iid'):.2f} naive_blk {m('cover_lam_naive_blk'):.2f} nuis_blk {m('cover_lam_nuis_blk'):.2f} decoy_blk {m('cover_lam_decoy_blk'):.2f} | SEratio {m('se_ratio_naive'):.1f} | SMD {m('max_smd_before'):.2f}->{m('max_smd_after'):.3f} | {(time.time()-t)/6:.1f}s/rep")
run("PURE  (no nuisance at all)", TissueConfig(clustering=0,conf_strength=0,autocorr_len_um=7.5,base_grf_amp=0.0,gamma_density=0.0,gamma_counts=0.0,density_mod=0.0))
run("EASY  (random senders, no latent confounder, ell=7.5)", TissueConfig(clustering=0,conf_strength=0,autocorr_len_um=7.5))
run("MID   (kappa=3, conf=1, ell=30 = lam_true)", TissueConfig(clustering=3.0,conf_strength=1.0,autocorr_len_um=30.0))
run("HARD  (kappa=3, conf=1, ell=120 = 4*lam)", TissueConfig(clustering=3.0,conf_strength=1.0,autocorr_len_um=120.0))
run("HARD2 (kappa=3, conf=2, ell=120)", TissueConfig(clustering=3.0,conf_strength=2.0,autocorr_len_um=120.0))
