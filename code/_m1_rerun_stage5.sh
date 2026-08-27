#!/bin/bash
# Phase 8 / 8.7 -- A7, summaries and the one-pass figure regeneration.
set -x
cd /workspace/code
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
L=/workspace/logs
# A7 depends on the sender call, which changed -> re-run before figure 2h
# A7 was launched earlier in the run, as soon as main_fits.csv landed.
python3 -u summarize_a7.py                      > $L/m1_a7_summary.log 2>&1
# summaries
python3 -u summarize_phase3.py                  > $L/m1_summary_p3.log 2>&1
python3 -u summarize_phase3_c1.py               > $L/m1_summary_p3c1.log 2>&1
python3 -u summarize_phase5.py                  > $L/m1_summary_p5.log 2>&1
python3 -u summarize_super_callers.py           > $L/m1_summary_super.log 2>&1
python3 -u summarize_caller_coverage.py         > $L/m1_caller_gate.log 2>&1
python3 -u caller_disagree_d3.py                 > $L/m1_caller_d3.log 2>&1
echo '--- figure guard, BEFORE the one-pass regeneration ---'
python3 /workspace/code/check_figures_guard.py > $L/m1_guard_before.log 2>&1
echo "guard_before_exit=$?"
# figures -- ONE pass, from the frozen configuration.
# figure2a caches its binned curves; the sender set changed, so the cache must go.
rm -f /workspace/figures/figure2a_stratified_curves.csv
python3 -u make_phase5_figs.py --which 2a,3     > $L/m1_fig2a3.log 2>&1
python3 -u make_figure2bc.py                    > $L/m1_fig2bc.log 2>&1
python3 -u make_figure4.py                      > $L/m1_fig4.log 2>&1
python3 -u make_figure_phase8_callers.py        > $L/m1_fig_callers.log 2>&1
python3 -u make_figure_phase8_d3.py             > $L/m1_fig_d3.log 2>&1
python3 -u make_phase3_figs.py                  > $L/m1_fig_p3.log 2>&1
echo '--- figure guard, AFTER the one-pass regeneration ---'
python3 /workspace/code/check_figures_guard.py > $L/m1_guard_after.log 2>&1
echo "guard_after_exit=$?"
echo M1_STAGE5_DONE
