#!/usr/bin/env bash


# Figure 1:
# - compare with/without CMFs
# - compare with/without epidermal pressure

echo "Simulations for figure 1..."

stomasimulator feb -s -w vicia-example        MR    > create-fig-f11.log 2>&1 &
stomasimulator feb -s -w vicia-example        tiMR  > create-fig-f12.log 2>&1 &
stomasimulator feb -s    vicia-example-epi-p  MR    > create-fig-f13.log 2>&1 &
stomasimulator feb -s    vicia-example-epi-p  tiMR  > create-fig-f14.log 2>&1 &

wait

# Figure 2:
# - optimised tiMR parameters

echo "Simulation for figure 2 (and S4)..."

stomasimulator feb -s -w vicia-optimised tiMR  > create-fig-f2.log 2>&1 &

wait

# Figure 3:
# - optimised tiVW parameters

echo "Simulation for figure 3..."

stomasimulator feb -s -w vicia-optimised tiVW  > create-fig-f3.log 2>&1 &

wait

# Table 4:
# - optimised tiVW parameters

echo "Simulations for table 4..."

stomasimulator feb -s ath-cw1  tiVW  > create-fig-t41.log 2>&1 &
stomasimulator feb -s ath-cw16 tiVW  > create-fig-t42.log 2>&1 &
stomasimulator feb -s ath-cw17 tiVW  > create-fig-t43.log 2>&1 &
stomasimulator feb -s ath-cw18 tiVW  > create-fig-t44.log 2>&1 &

wait

# Figure S1:
# - use "vicia-example tiMR" for the figures


# Figure S2:
# - baseline model is "vicia-example MR/tiMR"
# - Compare thicker walls vs. baseline model

echo "Simulations for figure S2..."

stomasimulator feb -s -w vicia-example-thicker          MR    > create-fig-fs21.log 2>&1 &
stomasimulator feb -s    vicia-example-thicker          tiMR  > create-fig-fs22.log 2>&1 &
stomasimulator feb -s -w vicia-example-thicker-ventral  MR    > create-fig-fs23.log 2>&1 &
stomasimulator feb -s    vicia-example-thicker-ventral  tiMR  > create-fig-fs24.log 2>&1 &

wait

# - compare stiffer matrix vs. baseline model
stomasimulator feb -s vicia-example --matrix-factor=2  MR    > create-fig-fs25.log 2>&1 &
stomasimulator feb -s vicia-example --matrix-factor=2  tiMR  > create-fig-fs26.log 2>&1 &
stomasimulator feb -s vicia-example --matrix-factor=4  MR    > create-fig-fs27.log 2>&1 &
stomasimulator feb -s vicia-example --matrix-factor=4  tiMR  > create-fig-fs28.log 2>&1 &

wait

# Figure S3:
# - Compare flattened cross-sections vs. baseline model
# - use whole stoma (for some) to show the initial cross-section

echo "Simulations for figure S3..."

stomasimulator feb -s -w vicia-example-ar-gt-1 MR    > create-fig-fs31.log 2>&1 &
stomasimulator feb -s    vicia-example-ar-gt-1 tiMR  > create-fig-fs32.log 2>&1 &
stomasimulator feb -s -w vicia-example-ar-lt-1 MR    > create-fig-fs33.log 2>&1 &
stomasimulator feb -s    vicia-example-ar-lt-1 tiMR  > create-fig-fs34.log 2>&1 &

wait

# Figure S5
# - Use the inferred tiVW cell wall matrix parameters (no fibres)

echo "Simulations for figure S5..."

stomasimulator feb -s -w vicia-optimised VW  > create-fig-fs5.log 2>&1 &

wait

# Figure S7
# - adjust the cell wall matrix and fibre parameters

echo "Simulations for figure S6..."

stomasimulator feb -s --matrix-factor=0.25  vicia-optimised tiVW  > create-fig-fs61.log 2>&1 &
stomasimulator feb -s --matrix-factor=0.5   vicia-optimised tiVW  > create-fig-fs62.log 2>&1 &
stomasimulator feb -s --matrix-factor=2.0   vicia-optimised tiVW  > create-fig-fs63.log 2>&1 &
stomasimulator feb -s --matrix-factor=4.0   vicia-optimised tiVW  > create-fig-fs64.log 2>&1 &

wait

stomasimulator feb -s --fibre-factor=0.0  vicia-optimised tiVW  > create-fig-fs65.log 2>&1 &
stomasimulator feb -s --fibre-factor=0.01 vicia-optimised tiVW  > create-fig-fs66.log 2>&1 &
stomasimulator feb -s --fibre-factor=0.1  vicia-optimised tiVW  > create-fig-fs67.log 2>&1 &
stomasimulator feb -s --fibre-factor=0.5  vicia-optimised tiVW  > create-fig-fs68.log 2>&1 &

wait

echo "Simulations for Lm > 1..."

# vary LM (strain stiffening fibres)
stomasimulator feb -s vicia-lm-1.01 tiVW  > create-fig-misc1.log 2>&1 &
stomasimulator feb -s vicia-lm-1.05 tiVW  > create-fig-misc2.log 2>&1 &
stomasimulator feb -s vicia-lm-1.1  tiVW  > create-fig-misc3.log 2>&1 &
stomasimulator feb -s vicia-lm-1.2  tiVW  > create-fig-misc4.log 2>&1 &

wait

# Table S1
# - sensitivity analysis

echo "Simulations for table S1..."

stomasimulator feb -s --matrix-factor=0.9  ath-cw1  tiVW  > create-fig-ts111.log 2>&1 &
stomasimulator feb -s --matrix-factor=1.1  ath-cw1  tiVW  > create-fig-ts112.log 2>&1 &
stomasimulator feb -s --matrix-factor=0.9  ath-cw16 tiVW  > create-fig-ts121.log 2>&1 &
stomasimulator feb -s --matrix-factor=1.1  ath-cw16 tiVW  > create-fig-ts122.log 2>&1 &
stomasimulator feb -s --matrix-factor=0.9  ath-cw17 tiVW  > create-fig-ts131.log 2>&1 &
stomasimulator feb -s --matrix-factor=1.1  ath-cw17 tiVW  > create-fig-ts132.log 2>&1 &
stomasimulator feb -s --matrix-factor=0.9  ath-cw18 tiVW  > create-fig-ts141.log 2>&1 &
stomasimulator feb -s --matrix-factor=1.1  ath-cw18 tiVW  > create-fig-ts142.log 2>&1 &

wait

stomasimulator feb -s --fibre-factor=0.9  ath-cw1  tiVW  > create-fig-ts151.log 2>&1 &
stomasimulator feb -s --fibre-factor=1.1  ath-cw1  tiVW  > create-fig-ts152.log 2>&1 &
stomasimulator feb -s --fibre-factor=0.9  ath-cw16 tiVW  > create-fig-ts161.log 2>&1 &
stomasimulator feb -s --fibre-factor=1.1  ath-cw16 tiVW  > create-fig-ts162.log 2>&1 &
stomasimulator feb -s --fibre-factor=0.9  ath-cw17 tiVW  > create-fig-ts171.log 2>&1 &
stomasimulator feb -s --fibre-factor=1.1  ath-cw17 tiVW  > create-fig-ts172.log 2>&1 &
stomasimulator feb -s --fibre-factor=0.9  ath-cw18 tiVW  > create-fig-ts181.log 2>&1 &
stomasimulator feb -s --fibre-factor=1.1  ath-cw18 tiVW  > create-fig-ts182.log 2>&1 &

wait

echo "Finished - phew!"
