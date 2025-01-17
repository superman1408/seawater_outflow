import math


# ----------------------------- INPUTS DATA ----------------------------------------------
HDPE_density_rho_HDPE = 960
Outside_diameter_OD = 2.3
Concrete_Coating_thickness_t_cc = 0.55
Wall_Thickness_t_HDPE = 88.5
Volume_of_Concrete_per_meter_of_pipe_Vc = 3.5
Concrete_density_rho_c = 2400
Marine_growth_Thickness_t_mg = 0
Marine_growth_density_rho_mg = 0

Content_density_seawater_rho_cont = 0
Safety_factor_for_weight_gamma_w = 1.1
Seawater_density_rho_seawater = 1025
gravity_g = 9.807

# --------------------------------- ENVIRONMENT DATA --------------------------------------
Significant_wave_height_Hs = 3.5
Spectral_peak_period_Tp = 13
Water_depth_d = 10.97
Related_angle_btw_pipeline_and_current_direction_teta = 90
Ref_major_height_above_the_seabed_z_r = 3

#------------------------------ SOIL DATA -------------------------------------------------
Sunbmerged_unit_soil_weight_for_sand_gamma_s = 13.5




# ----------------------------------CALCULATIOn--------------------------------------------

Inside_diameter_ID = Outside_diameter_OD - 2 * Wall_Thickness_t_HDPE
print(Inside_diameter_ID)
# ---------------------Water Partical Velocity---------------------------------------------

Peak_Enhanchment_factor_phi = Spectral_peak_period_Tp/s /math.sqrt(Significant_wave_height_Hs)

if Peak_Enhanchment_factor_phi <= 3.6 :
    Peakedness_Parameter_gamma = 5
elif 3.6 < Peak_Enhanchment_factor_phi < 5.0 :
    Peakedness_Parameter_gamma = math.exp(5.75 - 1.15*Peak_Enhanchment_factor_phi)
else :
    Peakedness_Parameter_gamma = 1

print(Peakedness_Parameter_gamma)

Reference_Period_Tn = math.sqrt(Water_depth_d/gravity_g)
print(Reference_Period_Tn)

Data_extracted_for_below_figure = Reference_Period_Tn/Spectral_peak_period_Tp
print(Data_extracted_for_below_figure)





 