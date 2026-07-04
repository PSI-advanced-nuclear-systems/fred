#!/bin/bash
cd /home/chan_y/Documents/fred_platform/examples/fred_m_na/timpano_SAS_benchmark/test_200d
source ~/miniforge3/etc/profile.d/conda.sh
conda activate fred-dev

nohup timeout 1200 python3 tuning_run.py loose_tol   1e-3 1e-3 0 200 > tuning_logs/loose_tol.log   2>&1 &
disown
nohup timeout 1200 python3 tuning_run.py maxord1     1e-5 1e-5 1 200 > tuning_logs/maxord1.log     2>&1 &
disown
nohup timeout 1200 python3 tuning_run.py combo       1e-3 1e-3 1 200 > tuning_logs/combo.log       2>&1 &
disown
nohup timeout 1200 python3 tuning_run.py combo_loose 1e-2 1e-2 1 200 > tuning_logs/combo_loose.log 2>&1 &
disown

sleep 2
ps aux | grep tuning_run | grep -v grep
