import math


# ----------------------------- INPUTS DATA ----------------------------------------------
HDPE_density_rho_HDPE = 960
Outside_diameter_OD = 2.3
Concrete_Coating_thickness_t_cc = 0.55
Wall_Thickness_t_HDPE = 88.5 /1000
Volume_of_Concrete_per_meter_of_pipe_Vc = 3.5
Concrete_density_rho_c = 2400
Marine_growth_Thickness_t_mg = 0
Marine_growth_density_rho_mg = 0

Content_density_seawater_rho_cont = 0
Safety_factor_for_weight_gamma_w = 1.1
Seawater_density_rho_seawater = 1025
gravity_g = 9.807

# --------------------------------------------------- ENVIRONMENT DATA --------------------------------------
Significant_wave_height_Hs = 3.5
Spectral_peak_period_Tp = 13
Water_depth_d = 10.97
Related_angle_btw_pipeline_and_current_direction_teta = 90
Ref_major_height_above_the_seabed_zr = 3

#------------------------------------------------------ SOIL DATA -------------------------------------------------
Sunbmerged_unit_soil_weight_for_sand_gamma_s = 13.5



# ----------------------------------------------------CALCULATIOn--------------------------------------------

Inside_diameter_ID = Outside_diameter_OD - 2 * Wall_Thickness_t_HDPE
print("Inside_diameter_ID",Inside_diameter_ID)

Hydodynamic_diameter_D = Outside_diameter_OD + 2 * (Marine_growth_density_rho_mg)
print(Hydodynamic_diameter_D)
# ---------------------------------------------Water Partical Velocity---------------------------------------------

Peak_Enhanchment_factor_phi = Spectral_peak_period_Tp /math.sqrt(Significant_wave_height_Hs)
print(Peak_Enhanchment_factor_phi)

if Peak_Enhanchment_factor_phi <= 3.6 :
    Peakedness_Parameter_gamma = 5
elif 3.6 < Peak_Enhanchment_factor_phi < 5.0 :
    Peakedness_Parameter_gamma = math.exp(5.75 - 1.15*Peak_Enhanchment_factor_phi)
else :
    Peakedness_Parameter_gamma = 1

print("Peakedness_Parameter_gamma",Peakedness_Parameter_gamma)

Reference_Period_Tn = math.sqrt(Water_depth_d/gravity_g)
print("Reference_Period_Tn ",Reference_Period_Tn)


Significant_wave_induced_water_particle_velocity_Us = Significant_wave_height_Hs/Reference_Period_Tn
print("Significant_wave_induced_water_particle_velocity_Us ",Significant_wave_induced_water_particle_velocity_Us)

Data_extracted_for_below_figure_2 = (Significant_wave_induced_water_particle_velocity_Us * Reference_Period_Tn)/Significant_wave_height_Hs
print("Data_extracted_for_below_figure_2 ",Data_extracted_for_below_figure_2)

# ---------------------------------------------------- Data_extracted_from_figure_3_2 ------------------------------------------------

Spectrally_Derived_mean_zero_Up_crossing_period = Reference_Period_Tn/Spectral_peak_period_Tp
print("Spectrally_Derived_mean_zero_Up_crossing_period  " , Spectrally_Derived_mean_zero_Up_crossing_period)

Mean_zero_up_crossing_period_at_seabed_level_Tu = 0.8415 * Spectral_peak_period_Tp
print("Mean_zero_up_crossing_period_at_seabed_level_Tu   ", Mean_zero_up_crossing_period_at_seabed_level_Tu)

Ratio_of_mean_zero_up_crossing_period_to_peak_wave_period = Mean_zero_up_crossing_period_at_seabed_level_Tu/Spectral_peak_period_Tp
print("Ratio_of_mean_zero_up_crossing_period_to_peak_wave_period ",Ratio_of_mean_zero_up_crossing_period_to_peak_wave_period)

# ------------------------------------------------------Input -----------------------------------------------------------------------

Reduction_factor_due_to_spectral_directionally_and_spreading_R_D = 1
Angle_btw_wave_and_pipeline_headings_teta_w = 90


# -----------------------------------------------------Calculation-------------------------------------------------------------------
Significant_wave_induced_water_particle_velocity_for_design_UsD  = Reduction_factor_due_to_spectral_directionally_and_spreading_R_D * Significant_wave_induced_water_particle_velocity_Us
print("Significant_wave_induced_water_particle_velocity_Us ",Significant_wave_induced_water_particle_velocity_Us)


# ------------------------------------------------------Current Velocity--------------------------------------------------------------
# ----------------------------------------------------input---------------------------------------------------------------------------

Current_velocity_a_distance_Zr_above_the_seabed_Vzr = 0.46
Effective_Seabed_Roughness_z0a = 4/10**5

Characteristic_current_induced_flow_velocity_at_pipeline_level_V = Current_velocity_a_distance_Zr_above_the_seabed_Vzr * (((1+Effective_Seabed_Roughness_z0a/Hydodynamic_diameter_D)* math.log(Hydodynamic_diameter_D/Effective_Seabed_Roughness_z0a +1) -1)/math.log(Ref_major_height_above_the_seabed_zr/Effective_Seabed_Roughness_z0a + 1)) * math.sin(Related_angle_btw_pipeline_and_current_direction_teta)
print("Characteristic_current_induced_flow_velocity_at_pipeline_level_V ",Characteristic_current_induced_flow_velocity_at_pipeline_level_V)

if Peakedness_Parameter_gamma == 1.0 :
    Friction_for_calculation_of_kt = 1.25

elif Peakedness_Parameter_gamma == 3.3 :
    Friction_for_calculation_of_kt = 1.21

elif Peakedness_Parameter_gamma == 5 :
    Friction_for_calculation_of_kt = 1.17

print("Friction_for_calculation_of_kt ",Friction_for_calculation_of_kt)



# Design_single_oscillation_velocity_period_at_seabed_level_T

if Reference_Period_Tn/Mean_zero_up_crossing_period_at_seabed_level_Tu <= 0.2 :
    Design_single_oscillation_velocity_period_at_seabed_level_T = Mean_zero_up_crossing_period_at_seabed_level_Tu * (Friction_for_calculation_of_kt - 5* (Friction_for_calculation_of_kt -1 )* (Reference_Period_Tn/Mean_zero_up_crossing_period_at_seabed_level_Tu))
else:
    Design_single_oscillation_velocity_period_at_seabed_level_T = 1
print("Design_single_oscillation_velocity_period_at_seabed_level_T ",Design_single_oscillation_velocity_period_at_seabed_level_T)

Current_to_Wave_Velocity_Ratio_M = Characteristic_current_induced_flow_velocity_at_pipeline_level_V/Significant_wave_induced_water_particle_velocity_Us
print("Current_to_Wave_Velocity_Ratio_M ", Current_to_Wave_Velocity_Ratio_M)

Keulegan_Carpenter_Number_K = Significant_wave_induced_water_particle_velocity_Us *Design_single_oscillation_velocity_period_at_seabed_level_T / Hydodynamic_diameter_D
print("Keulegan_Carpenter_Number_K ",Keulegan_Carpenter_Number_K)

# Weight Calculation

Area_of_pipe_AOD = math.pi * Outside_diameter_OD**2 /4
print("Area_of_pipe_AOD  " ,Area_of_pipe_AOD)

Volume_of_pipe_VOD = Area_of_pipe_AOD * 1
print("Volume_of_pipe_VOD" ,Volume_of_pipe_VOD)

Internal_Area_of_pipe_AID = math.pi * Inside_diameter_ID**2/4
print("Internal_Area_of_pipe_AID  " ,Internal_Area_of_pipe_AID)


 