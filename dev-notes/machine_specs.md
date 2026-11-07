# Machine specs (2026-11-02 22:36:13)

Name                      : 12th Gen Intel(R) Core(TM) i7-12700H
NumberOfCores             : 14
NumberOfLogicalProcessors : 20
MaxClockSpeed             : 2300

RAM         : 16 GB LPDDR5, 4800 MT/s
Manufacturer : ASUSTeK COMPUTER INC.
Model        : Zenbook UP6502ZD_Q539ZD

## R / tools
R version 4.6.1 (2026-06-24 ucrt) 
mirai 2.7.3 
mori 0.2.2 
nanonext 1.10.3 
rstan 2.32.7 
rstanarm 2.32.2 
brms 2.23.0 
posterior 1.7.1 
foil NA 

## WSL
Linux Laptop 6.18.33.2-microsoft-standard-WSL2 #1 SMP PREEMPT_DYNAMIC Thu Jun 18 21:54:43 UTC 2026 x86_64 x86_64 x86_64 GNU/Linux
20
               total        used        free      shared  buff/cache   available
Mem:            8934         739        7067           3        1284        8195
Rscript (R) version 4.5.3 (2026-03-11)
model name	: 12th Gen Intel(R) Core(TM) i7-12700H
[wsl2]
memory=12GB

## Cluster (NRP TIDE), recorded 2026-11-07
Node: Intel Xeon Gold 6430, 2 sockets x 32 cores, 64 CPUs, kernel 5.15
Image: Ubuntu 24.04.5, R 4.5.3 (rocker/r-ver:4.5.3), 6 GiB memory-backed /dev/shm
Each worker count ran as its own pod on its own node (anti-affinity), requests equal to limits. Time only: the cgroup is read-only, so foil's Linux fallback does not measure memory.
