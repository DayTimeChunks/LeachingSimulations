
#################
# 1st Pulse
#################
# Leached Metalaxyl (Observed) - Vine
# 6min, 12min, 30min, 30min @ 135, 55, 55, 30mm/h
# 0d, 10d
mx_obs_sterile_lch1 = np.array([])  # high, med-12, med-30, low
mx_obs_untreat_lch1 = np.array([])

mx_sol_leach1 = [mx_obs_sterile_lch1, mx_obs_untreat_lch1]

# Ponded Metalaxyl ( Observed) - Vine
# high, med-12, med-30, low
mx_obs_sterile_vine_roff = np.array([0., 0., 0., 0.])

mx_obs_untreat_vine_roff = np.array([0., 0., 0., 0.])

mx_sol_pond1 = [mx_obs_sterile_vine_roff, mx_obs_untreat_vine_roff]

#################
# 2nd Pulse
#################
# Observed Metalaxyl Leached - Vine
# high-6, med-12, med-30, low-30
# 0d, 10d
mx_obs_sterile_lch2 = np.array([])
mx_obs_untreat_lch2 = np.array([])

mx_sol_leach2 = [mx_obs_sterile_lch2, mx_obs_untreat_lch2]

# Observed Metalaxyl in Ponding - Vine
# high-6, med-12, med-30, low-30
# 0d, 10d
mx_obs_sterile_vine_roff2 = np.array([])

mx_obs_untreat_vine_roff2 = np.array([])

mx_sol_pond2 = [mx_obs_sterile_vine_roff2, mx_obs_untreat_vine_roff2]


########################
## DETAILED
#########################
# Organized by intensity:
# [sterile, untreat, sterile_aged, untreat_aged]
# all at 6 min, high inetnesity
leach1_mx_high_6min = np.array([])
leach2_mx_high_6min = np.array([])


# all at 12 min, med intensity
leach1_mx_med_12min = np.array([])
leach2_mx_med_12min = np.array([])

# all at 30min, med intensity
leach1_mx_med_30min = np.array([])
leach2_mx_med_30min = np.array([])

# all at 30min, low intensity
leach1_mx_low_30min = np.array([])
leach2_mx_low_30min = np.array([])

# Organized individually:
# [sterile, untreat, sterile_aged, untreat_aged]
# all at 6 min, high inetnesity

pond2_mx_high_6min = np.array([])

# all at 12 min, med intensity
pond2_mx_med_12min = np.array([])

# all at 30min, med intensity
pond2_mx_med_30min = np.array([])

# all at 30min, low intensity
pond2_mx_low_30min = np.array([])


