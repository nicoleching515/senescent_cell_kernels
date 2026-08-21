import time, numpy as np
from sasp_sim import TissueConfig, simulate_tissue, ripley_ratio
from sasp_estimators import analyze

def show(tag, cfg, seeds=(1,2,3,4,5)):
    L=[]
    for s in seeds:
        t=time.time(); sim=simulate_tissue(cfg,[7,s]); res=analyze(sim,[8,s],n_boot=60); L.append(res)
    def m(k): 
        v=np.array([r[k] for r in L],float); return np.nanmean(v)
    print(f"--- {tag}  (n={L[0]['n_cells']}, senders={L[0]['n_senders']}, "
          f"ripley50={ripley_ratio(sim['coords'],sim['sender_mask']):.2f}, {time.time()-t:.1f}s/rep)")
    print(f"  lam_true={cfg.lambda_true_um}  lam_naive={m('lam_naive'):.1f}  lam_binned={m('lam_naive_binned'):.1f} "
          f"lam_nuis={m('lam_nuis'):.1f}  lam_decoy={m('lam_decoy_corr'):.1f}")
    print(f"  beta_true={cfg.beta_true}  beta_naive={m('beta_naive_binned'):.3f}  beta_decoyonly={m('beta_decoy_only'):.3f} "
          f" beta_decoycorr={m('beta_decoy_corr'):.3f}")
    print(f"  maxSMD before={m('max_smd_before'):.3f} after={m('max_smd_after'):.3f}  match_rate={m('match_rate'):.2f}")
    print(f"  cover lam naive_iid={m('cover_lam_naive_iid'):.2f} naive_block={m('cover_lam_naive_block'):.2f} decoy_block={m('cover_lam_decoy_block'):.2f}")

# EASY: random senders, no confounder, short baseline autocorrelation
show("EASY", TissueConfig(clustering=0.0, conf_strength=0.0, autocorr_len_um=7.5, base_grf_amp=0.55))
# HARD: clustered senders, strong confounder, long autocorrelation
show("HARD", TissueConfig(clustering=3.0, conf_strength=1.5, autocorr_len_um=120.0))
