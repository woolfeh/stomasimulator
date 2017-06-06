########################################################################################################################
## Figure S3 - Aspect ratio and aperture
########################################################################################################################

reset

########################################################################################################################
## Terminal id and plot placement variables

term_id = 0

sz_w  =  960; sz_h  = 540
scr_x = 1200; scr_y = 100
dx    =   20; dy    = 20

########################################################################################################################
# Name the files

# main iso
f_iso_1 = '../simulations/vicia-example,2gc,d_1.0,v_1.0,MR,c1_9.0000,c2_9.0000.stats.txt'

# iso with AR != 1
f_iso_0p8  = '../simulations/vicia-example-ar-lt-1,2gc,d_1.0,v_1.0,ar_0.8,MR,c1_9.0000,c2_9.0000.stats.txt'
f_iso_1p25 = '../simulations/vicia-example-ar-gt-1,2gc,d_1.0,v_1.0,ar_1.25,MR,c1_9.0000,c2_9.0000.stats.txt'

# main aniso
f_aniso_1 = '../simulations/vicia-example,2gc,d_1.0,v_1.0,tiMR,c1_9.0000,c2_9.0000,c5_1000.00,lm_1.00.stats.txt'

# aniso with AR != 1
f_aniso_0p8  = '../simulations/vicia-example-ar-lt-1,1gc,d_1.0,v_1.0,ar_0.8,tiMR,c1_9.0000,c2_9.0000,c5_1000.00,lm_1.00.stats.txt'
f_aniso_1p25 = '../simulations/vicia-example-ar-gt-1,1gc,d_1.0,v_1.0,ar_1.25,tiMR,c1_9.0000,c2_9.0000,c5_1000.00,lm_1.00.stats.txt'

########################################################################################################################
# Pressure vs. relative change in aperture (AR changes)

set term wxt term_id size sz_w,sz_h position scr_x + dx*term_id, scr_y + dy*term_id
term_id = term_id + 1

set xtics out nomirror
set ytics out
set border 1 + 2 + 8 # Left (2), bottom (1) and right (8)

set key top left Left reverse maxrows 3 width 1 height 0.5 samplen 3

set xlabel 'Turgor pressure (MPa)'
set xrange [0:5]

set ylabel 'Aperture ({/Symbol m}m)'
set yrange [0:8]
#set ytics format "%0.1f"

p \
  f_aniso_0p8   i 0 u 1:'pore-width' w l lc 2      lw 3 t'With CMFs: AR=0.8', \
  f_aniso_1     i 0 u 1:'pore-width' w l lc 1      lw 3 t'With CMFs: AR=1', \
  f_aniso_1p25  i 0 u 1:'pore-width' w l lc 3      lw 3 t'With CMFs: AR=1.25', \
  f_iso_0p8     i 0 u 1:'pore-width' w l lc 2 dt 2 lw 3 t'Isotropic: AR=0.8', \
  f_iso_1       i 0 u 1:'pore-width' w l lc 1 dt 2 lw 3 t'Isotropic: AR=1', \
  f_iso_1p25    i 0 u 1:'pore-width' w l lc 3 dt 2 lw 3 t'Isotropic: AR=1.25'

reset

########################################################################################################################
# Pressure vs. aspect ratio change for the stomata

set term wxt term_id size sz_w,sz_h position scr_x + dx*term_id, scr_y + dy*term_id
term_id = term_id + 1

set xtics out nomirror
set ytics out
set border 1 + 2 + 8 # Left (2), bottom (1) and right (8)

set lmargin screen 0.3

#set key top left # this makes it set the key's coordinates to the top left corner
#set key at 0.32,1.24 Left reverse maxrows 3 width 1 height 0.5 samplen 3

set key top cent Left reverse maxrows 3 width 1 height 0.5 samplen 3

set xlabel 'Turgor pressure (MPa)'
set ylabel 'Aspect ratio (at GC mid-point)'

set xrange [0:5]
set yrange [0.8:1.25]

set ytics format "%0.2f"

p \
  f_aniso_0p8  i 0 u 1:'gc0-ar' w l lc 2      lw 3 t'With CMFs: AR=0.8', \
  f_aniso_1    i 0 u 1:'gc0-ar' w l lc 1      lw 3 t'With CMFs: AR=1', \
  f_aniso_1p25 i 0 u 1:'gc0-ar' w l lc 3      lw 3 t'With CMFs: AR=1.25', \
  f_iso_0p8    i 0 u 1:'gc0-ar' w l lc 2 dt 2 lw 3 t'Isotropic: AR=0.8', \
  f_iso_1      i 0 u 1:'gc0-ar' w l lc 1 dt 2 lw 3 t'Isotropic: AR=1', \
  f_iso_1p25   i 0 u 1:'gc0-ar' w l lc 3 dt 2 lw 3 t'Isotropic: AR=1.25'

reset

########################################################################################################################
