reinitialize

load datapdb2f4k_cut.ent, pdb2f4k
alter pdb2f4k, resi=str(int(resi)-41)
sort

hide everything, all
show cartoon, pdb2f4k
color gray80, pdb2f4k

bg_color white
set cartoon_fancy_helices, 1
set dash_width, 3
set dash_gap, 0.25
set dash_radius, 0.08
set label_size, 18

# Cluster colors
set_color cluster1_col, [1.00, 0.15, 0.15]
set_color cluster2_col, [0.10, 0.35, 1.00]
set_color cluster3_col, [0.10, 0.75, 0.20]
set_color cluster4_col, [1.00, 0.55, 0.00]
set_color cluster5_col, [0.65, 0.20, 0.90]
set_color cluster6_col, [0.00, 0.75, 0.75]
set_color cluster7_col, [0.95, 0.85, 0.10]

# Show contacted residues
#select contacted_residues, pdb2f4k and resi 3+5+6+7+9+10+11+12+13+14+16+17+18+20+24+25+28+29+30+32+33+34+35
#show sticks, contacted_residues
#color white, contacted_residues

# Contact lines use CA atoms. Change /CA to /CB if desired.
distance c1_3_14, resi 3 and name CA, resi 14 and name CA
distance c1_3_13, resi 3 and name CA, resi 13 and name CA
distance c1_6_14, resi 6 and name CA, resi 14 and name CA
distance c1_5_14, resi 5 and name CA, resi 14 and name CA
color cluster1_col, c1_*
disable c1_*

distance c2_7_12, resi 7 and name CA, resi 12 and name CA
distance c2_7_13, resi 7 and name CA, resi 13 and name CA
distance c2_6_12, resi 6 and name CA, resi 12 and name CA
distance c2_7_11, resi 7 and name CA, resi 11 and name CA
distance c2_6_11, resi 6 and name CA, resi 11 and name CA
distance c2_6_17, resi 6 and name CA, resi 17 and name CA
color cluster2_col, c2_*
disable c2_*

distance c3_12_17, resi 12 and name CA, resi 17 and name CA
distance c3_12_16, resi 12 and name CA, resi 16 and name CA
distance c3_12_20, resi 12 and name CA, resi 20 and name CA
distance c3_13_17, resi 13 and name CA, resi 17 and name CA
color cluster3_col, c3_*
disable c3_*

distance c4_18_25, resi 18 and name CA, resi 25 and name CA
distance c4_17_25, resi 17 and name CA, resi 25 and name CA
distance c4_20_25, resi 20 and name CA, resi 25 and name CA
color cluster4_col, c4_*
disable c4_*

distance c5_24_28, resi 24 and name CA, resi 28 and name CA
distance c5_20_28, resi 20 and name CA, resi 28 and name CA
distance c5_25_29, resi 25 and name CA, resi 29 and name CA
color cluster5_col, c5_*
disable c5_*

distance c6_29_35, resi 29 and name CA, resi 35 and name CA
distance c6_29_34, resi 29 and name CA, resi 34 and name CA
distance c6_30_35, resi 30 and name CA, resi 35 and name CA
distance c6_29_33, resi 29 and name CA, resi 33 and name CA
color cluster6_col, c6_*
disable c6_*

distance c7_10_34, resi 10 and name CA, resi 34 and name CA
distance c7_9_32, resi 9 and name CA, resi 32 and name CA
distance c7_10_29, resi 10 and name CA, resi 29 and name CA
color cluster7_col, c7_*
disable c7_*

# Re-enable all distance objects after coloring
enable c*_*

hide labels, c*_*
#zoom contacted_residues, 10
#orient contacted_residues


# Terminal labels

select Nterm, pdb2f4k and resi 1 and name CA
select Cterm, pdb2f4k and resi 35 and name CA
label Nterm, "N-terminus"
label Cterm, "C-terminus"
set label_size, 24
set label_color, black
set label_position, [2, 2, 2]
select chain_labels, pdb2f4k and name CA and resi 1+5+10+15+20+25+30+35
label chain_labels, resi