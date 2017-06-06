########################################################################################################################
## Figure 2 - Optimisation results for the tiMR cell wall model for the Vicia faba data of Spence et al. (1986)
########################################################################################################################

reset

########################################################################################################################
## Terminal id and plot placement variables

term_id = 0

sz_w  =  960; sz_h  = 540
scr_x = 1200; scr_y = 100
dx    =   20; dy    = 20

########################################################################################################################
# Name the file

f_aniso_mr = '../simulations/vicia-optimised,2gc,d_1.0,v_1.0,tiMR,c1_4.4546,c2_4.4779,c5_1169.65,lm_1.00.stats.txt'

########################################################################################################################
## Change in aperture and volume

set term wxt term_id size sz_w,sz_h position scr_x + dx*term_id, scr_y + dy*term_id
term_id = term_id + 1

set xtics  out nomirror
set ytics  out nomirror 0,1
set y2tics out nomirror
set mytics 2
set my2tics 2

set border 1 + 2 + 8 # Left (2), bottom (1) and right (8)

set key top left Left reverse

set xlabel 'Turgor pressure (MPa)'

set xrange  [0:5]
set x2range [0:5]

set ylabel 'Relative to initial value'

set ylabel  'Aperture ({/Symbol m}m)'
set y2label 'Guard cell volume ({/Symbol m}m^3)'

#set ytics 1,0.5 format "%.1f"
#set yrange [1:8]

p \
	f_aniso_mr  i 0 u 1:'pore-width' w l lw 3 lc 1 t'Aperture'           axes x1y1, \
    ''          i 0 u 1:'gc0-volume' w l lw 3 lc 7 t'Guard cell volume'  axes x2y2

reset

########################################################################################################################
## Change in stoma length and surface area

set term wxt term_id size sz_w,sz_h position scr_x + dx*term_id, scr_y + dy*term_id
term_id = term_id + 1

set xtics  out nomirror
set ytics  out nomirror 0,1
set y2tics out nomirror
set mytics 2
set my2tics 2

set border 1 + 2 + 8 # Left (2), bottom (1) and right (8)

set key top left Left reverse

set xlabel 'Turgor pressure (MPa)'

set xrange  [0:5]
set x2range [0:5]

set ylabel  'Stoma length ({/Symbol m}m)'
set y2label 'Guard cell surface area ({/Symbol m}m^2)'

p \
	f_aniso_mr  i 0 u 1:'stoma-length'        w l lw 3 lc 2 t 'Stoma length'             axes x1y1, \
	''          i 0 u 1:'gc0-surface-area'    w l lw 3 lc 4 t 'Guard cell surface area'  axes x2y2

reset

########################################################################################################################
