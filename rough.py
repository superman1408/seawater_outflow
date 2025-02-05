import math
Current_velocity_a_distance_Zr_above_the_seabed_Vzr =0.46
Effective_Seabed_Roughness_z0a = 4/10**5
Hydodynamic_diameter_D = 3.4


Ref_major_height_above_the_seabed_zr = 3
Related_angle_btw_pipeline_and_current_direction_teta = 90


Characteristic_current_induced_flow_velocity_at_pipeline_level_V = Current_velocity_a_distance_Zr_above_the_seabed_Vzr * (((1+(Effective_Seabed_Roughness_z0a/Hydodynamic_diameter_D))* math.log((Hydodynamic_diameter_D/Effective_Seabed_Roughness_z0a) +1) -1)/math.log((Ref_major_height_above_the_seabed_zr/Effective_Seabed_Roughness_z0a) + 1)) * math.sin(Related_angle_btw_pipeline_and_current_direction_teta)
print("Characteristic_current_induced_flow_velocity_at_pipeline_level_V ",Characteristic_current_induced_flow_velocity_at_pipeline_level_V)