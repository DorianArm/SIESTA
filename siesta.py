from siesta_package.siesta_utils import *


#-------- Main App --------#


### Select spectral range in nm
spectral_range = (507, 875)

### Select camera sensor size
#Emergent camera size
cmosx_max = 12.8
cmosy_max = 12.8

### Select center(s) of Echellograms, defined by the center wavelengths from spectral_range. By default, centered inside camera sensor. 

# 1 spatial element
centers_list = [(cmosx_max/2,cmosy_max/2)]
# centers_list = [(cmosx_max/2,0.8*cmosy_max)]

# >1 spatial elements
# centers_list = [(cmosx_max/4,cmosy_max/2),(3/4*cmosx_max,cmosy_max/2)]
# centers_list = [(cmosx_max/2,0.505*cmosy_max),(cmosx_max/2,cmosy_max/2)]
# centers_list = [(cmosx_max/2,cmosy_max/4),(cmosx_max/2,3/4*cmosy_max)]
# centers_list = [(cmosx_max/6,0.18*cmosy_max),(cmosx_max/2,0.18*cmosy_max),(5/6*cmosx_max,0.18*cmosy_max),(cmosx_max/6,0.51*cmosy_max),(cmosx_max/2,0.51*cmosy_max),(5/6*cmosx_max,0.51*cmosy_max),(1/6*cmosx_max,0.84*cmosy_max),(1/2*cmosx_max,0.84*cmosy_max),(5/6*cmosx_max,0.84*cmosy_max)]

### Select number of spectral points to compute (for each spatial element)
n_spectral = 10000

### depreciated (do not remove), was working with plotly only
save_file = False #not working with Dash, need to find a solution


### Select optical elements parameters ###
ech = EchelleGrating(name="thorlabs echelle", groove_density=31.6, blaze_angle=63, semi_deviation_angle_deg=7.5)
disp = Grating(name="disperser", groove_density=300)  
camera = Lens(name="camera", focal_length=100, diameter=50)
collimator = Lens(name="collimator", focal_length=200, diameter=50)
slit = Slit(name="slit", width=0.03, height=0.2) #mm @ZimMAIN 60" = 1.629 mm and 5" = 0.136 mm, 1" = 27 um
# slit = Slit(name="slit", width=0.492, height=5.905) #mm @SST 60" = 5.905 mm and 5" = 0.492 mm, 1" = 98.4 um
if __name__ == "__main__":
    df_mapping_list = computeCD(spatial_centers=centers_list, spectral_range=spectral_range, n_spectral=n_spectral, echelle=ech, disperser=disp, camera_lens=camera, collimator_lens=collimator, cmosx_max=cmosx_max, cmosy_max=cmosy_max, alpha_deg=15, slit=slit,write_hmtl=save_file)
# df_mapping_list = computeCD(spatial_centers=centers_list, spectral_range=spectral_range, n_spectral=n_spectral, echelle=ech, disperser=disp_prism, camera_lens=camera, collimator_lens=collimator, cmosx_max=cmosx_max, cmosy_max=cmosy_max, alpha_deg=15, slit=slit,write_hmtl=save_file)






