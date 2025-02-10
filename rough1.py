import math
Current_velocity_a_distance_Zr_above_the_seabed_Vzr =0.46
Effective_Seabed_Roughness_z0a = 4/10**5
Hydodynamic_diameter_D = 3.4


Ref_major_height_above_the_seabed_zr = 3
Related_angle_btw_pipeline_and_current_direction_teta = 90


A = 1+Effective_Seabed_Roughness_z0a/Hydodynamic_diameter_D

B = math.log(Hydodynamic_diameter_D/Effective_Seabed_Roughness_z0a+1)

C = math.log(Ref_major_height_above_the_seabed_zr/Effective_Seabed_Roughness_z0a+1)


print(A)
print(B)
print(C)


V = Current_velocity_a_distance_Zr_above_the_seabed_Vzr * (A*B-1/C)*math.sin(Related_angle_btw_pipeline_and_current_direction_teta)

print("V",V)