from siesta_package.siesta_utils import *


#-------- Main App --------#

### Select spectral range in nm
spectral_range = (350, 875)

### Select camera sensor size
# #Emergent camera
# emergent = Camera_sensor(name="Emergent", pixels_x=5120, pixels_y=5120, pixel_size=2.5)

# # Andor Zyla
# zyla = Camera_sensor(name="Andor Zyla", pixels_x=2160, pixels_y=2560, pixel_size=6.5) #mm

# ZWO ASI294MM Pro
asi294 = Camera_sensor(name="ZWO ASI294MM Pro", pixels_x=5644, pixels_y=8288, pixel_size=2.3) #mm

### Select center(s) of Echellograms, defined by the center wavelengths from spectral_range. By default, centered inside camera sensor. 

# 1 spatial element
# centers_list = [(zyla.size_x_mm/2, zyla.size_y_mm/2)]
centers_list = [(asi294.size_x_mm/2, 4/5*asi294.size_y_mm)]
# centers_list = [(5, 5)]
# centers_list = [(cmosx_max/2,0.8*cmosy_max)]

# >1 spatial elements
# shift_x = 0.1 #mm shift in x from initial position
# shift_y = 0.075 #mm shift in y from initial position
# centers_list = [(asi294.size_x_mm/2, 4/5*asi294.size_y_mm + shift_y), (asi294.size_x_mm/2, 4/5*asi294.size_y_mm - shift_y)]



# centers_list = [(cmosx_max/4,cmosy_max/2),(3/4*cmosx_max,cmosy_max/2)]
# centers_list = [(cmosx_max/2,0.505*cmosy_max),(cmosx_max/2,cmosy_max/2)]
# centers_list = [(cmosx_max/2,cmosy_max/4),(cmosx_max/2,3/4*cmosy_max)]
# centers_list = [(cmosx_max/6,0.18*cmosy_max),(cmosx_max/2,0.18*cmosy_max),(5/6*cmosx_max,0.18*cmosy_max),(cmosx_max/6,0.51*cmosy_max),(cmosx_max/2,0.51*cmosy_max),(5/6*cmosx_max,0.51*cmosy_max),(1/6*cmosx_max,0.84*cmosy_max),(1/2*cmosx_max,0.84*cmosy_max),(5/6*cmosx_max,0.84*cmosy_max)]

### Select number of spectral points to compute (for each spatial element)
spectral_res = 0.02 #nm, used for sampling spectral range.

### Select optical elements parameters ###
ech = EchelleGrating(name="thorlabs echelle", groove_density=31.6, blaze_angle=63, semi_deviation_angle_deg=3.375)
# ech = EchelleGrating(name="thorlabs echelle", groove_density=79, blaze_angle=64, semi_deviation_angle_deg=3.375)

# current disperser v0
# disp = Grating(name="disperser", groove_density=300, alpha=18)

# CaF2 prism
# disp = Prism(name="disperser prism", glass_type=["main","CaF2","Li"], beam_diameter=20, base=50, apex_angle_deg=60, input_angle_deg=25, manual=False, spectral_range_nm=spectral_range, spectral_res_nm=spectral_res) #CaF2-Li manual prism
# disp2 = Prism(name="second disperser prism", glass_type=["main","CaF2","Li"], beam_diameter=20, base=50, apex_angle_deg=60, input_angle_deg=48, prev=disp, manual=False, spectral_range_nm=spectral_range, spectral_res_nm=spectral_res) #CaF2-Li manual prism

# (N-)F2 prism
# disp = Prism(name="disperser prism", glass_type=["specs","SCHOTT-optical","N-F2"], beam_diameter=20, base=50, apex_angle_deg=45, input_angle_deg=25, manual=False, spectral_range_nm=spectral_range, spectral_res_nm=spectral_res) #sf11 manual prism
# disp2 = Prism(name="disperser prism 2", glass_type=["specs","SCHOTT-optical","F2"], beam_diameter=20, base=50, apex_angle_deg=45, input_angle_deg=18, manual=False, prev=disp, spectral_range_nm=spectral_range, spectral_res_nm=spectral_res) #sf11 manual prism

# UVFS prism (SiO2)
# disp = Prism(name="disperser prism", glass_type=["main","SiO2","Malitson"], beam_diameter=20, base=50, apex_angle_deg=60, input_angle_deg=40, manual=False, spectral_range_nm=spectral_range, spectral_res_nm=spectral_res) #sf11 manual prism
# disp2 = Prism(name="disperser prism 2", glass_type=["main","SiO2","Malitson"], beam_diameter=20, base=50, apex_angle_deg=60, input_angle_deg=26, manual=False, prev=disp, spectral_range_nm=spectral_range, spectral_res_nm=spectral_res) #sf11 manual prism

# N-BK7 prism
disp = Prism(name="disperser prism", glass_type=["specs","SCHOTT-optical","N-BK7"], beam_diameter=20, base=50,apex_angle_deg=60, input_angle_deg=38, manual=False, spectral_range_nm=spectral_range, spectral_res_nm=spectral_res)
disp2 = Prism(name="disperser prism 2", glass_type=["specs","SCHOTT-optical","N-BK7"], beam_diameter=20, base=50,apex_angle_deg=45, input_angle_deg=15, manual=False, spectral_range_nm=spectral_range, spectral_res_nm=spectral_res, prev=disp)
disp3 = Prism(name="disperser prism 3", glass_type=["specs","SCHOTT-optical","N-BK7"], beam_diameter=20, base=50,apex_angle_deg=45, input_angle_deg=46, manual=False, spectral_range_nm=spectral_range, spectral_res_nm=spectral_res, prev=disp2)
# disp2 = Prism(name="disperser prism 2", glass_type=["specs","SCHOTT-optical","N-BK7"], beam_diameter=20, base=50,apex_angle_deg=60, input_angle_deg=31, manual=False, spectral_range_nm=spectral_range, spectral_res_nm=spectral_res, prev=disp)

# ZnSe prism (infrared)
# disp = Prism(name="disperser prism", glass_type=["main","ZnSe","Marple"], beam_diameter=15, base=5,manual=False) #sf11 manual prism

# N-FK5 prism
# disp = Prism(name="disperser prism", glass_type=["specs","SCHOTT-optical","N-FK5"], beam_diameter=20, base=50, apex_angle_deg=60, input_angle_deg=40, manual=False, spectral_range_nm=spectral_range, spectral_res_nm=spectral_res) 
# disp2 = Prism(name="disperser prism 2", glass_type=["specs","SCHOTT-optical","N-FK5"], beam_diameter=20, base=50, apex_angle_deg=60, input_angle_deg=23, manual=False, prev=disp, spectral_range_nm=spectral_range, spectral_res_nm=spectral_res) 


# manual coefficients glass
# sf11 S coeffs: [2, 0.0, 1.7385, 0.0136, 0.311, 0.0616, 1.1749, 121.92] #First element 2 is to specify Sllemeir-2 formula
# N-BK7 S coeffs: [2, 0.0, 1.12735, 0.0072, 0.1244, 0.0270, 0.827, 100.38] #First element 2 is to specify Sllemeir-2 formula
# ZnSe S coeffs: [2, 3, 1.9, 0.113, 0.0, 0.0, 0.0, 0] #First element 2 is to specify Sllemeir-2 formula



camera = Lens(name="camera", focal_length=128, diameter=50)
collimator = Lens(name="collimator", focal_length=217, diameter=50)
slit = Slit(name="slit", width=0.03, height=0.2) #mm @ZimMAIN 60" = 1.629 mm and 5" = 0.136 mm, 1" = 27 um
# slit = Slit(name="slit", width=0.492, height=5.905) #mm @SST 60" = 5.905 mm and 5" = 0.492 mm, 1" = 98.4 um
if __name__ == "__main__":
    sowisp = Instrument(name="SOWISP", spectral_range=spectral_range, spatial_centers=centers_list, spectral_res_nm=spectral_res, echelle=ech, disperser=disp3, camera_lens=camera, collimator_lens=collimator, slit=slit, camera_sensor=asi294, wavelength_scan_width_nm=0.02)
    # sowisp.exportAsImage(species="Thorium", filename="Thorium_sowisp_ASI294MMPro", path="./correlation_lab_siesta", spectral_range_nm=spectral_range)
    # sowisp.exportAsFits(species=["Neon"], filename="Neon_sowisp_Zyla_slit_ESOv2", path="./correlation_lab_siesta", spectral_range_nm=spectral_range, wantSlitKernel=True)
    # sowisp.exportAsFits(species=["Thorium"], filename="Thorium_sowisp_ASI294MMPro_slit", path="./correlation_lab_siesta", spectral_range_nm=spectral_range, wantSlitKernel=True)
    # sowisp.exportDFmapping(df_mapping_list_indices=[0],filename="./exports/sowisp_mapping_asi294mm")
    sowisp.plotCD()
