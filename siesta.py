import numpy as np
# import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
import dash
from dash import dcc, html
from dash.dependencies import Input, Output
from math import ceil, floor
import pandas as pd
import os.path
import refractiveindex as ri



#-------- Main App

spectral_range = (507, 875)
#Emergent camera size
# cmosx_max = 12.8
# cmosy_max = 12.8

#Arbitrary camera size
cmosx_max = 12.8
cmosy_max = 12.8

centers_list = [(cmosx_max/2,cmosy_max/2)]
# centers_list = [(cmosx_max/2,0.8*cmosy_max)]
# centers_list = [(cmosx_max/4,cmosy_max/2),(3/4*cmosx_max,cmosy_max/2)]
# centers_list = [(cmosx_max/2,0.505*cmosy_max),(cmosx_max/2,cmosy_max/2)]
# centers_list = [(cmosx_max/2,cmosy_max/4),(cmosx_max/2,3/4*cmosy_max)]
# centers_list = [(cmosx_max/6,0.18*cmosy_max),(cmosx_max/2,0.18*cmosy_max),(5/6*cmosx_max,0.18*cmosy_max),(cmosx_max/6,0.51*cmosy_max),(cmosx_max/2,0.51*cmosy_max),(5/6*cmosx_max,0.51*cmosy_max),(1/6*cmosx_max,0.84*cmosy_max),(1/2*cmosx_max,0.84*cmosy_max),(5/6*cmosx_max,0.84*cmosy_max)]
n_spectral = 10000

save_file = False #not working with Dash, need to find a solution

ech = EchelleGrating(name="thorlabs echelle", groove_density=31.6, blaze_angle=63, semi_deviation_angle_deg=7.5)
disp = Grating(name="disperser", groove_density=300)  
# disp_prism = Prism(name="disperser prism", glass_type=["main","CaF2","Li"], beam_diameter=25, base=100)
# disp_prism = Prism(name="disperser prism", glass_type=["specs","SCHOTT-optical","SF11"], beam_diameter=15, base=25,manual=True) #sf11 manual prism
# disp_prism = Prism(name="disperser prism", glass_type=["specs","SCHOTT-optical","BK7"], beam_diameter=15, base=25,manual=True) #sf11 manual prism
# disp_prism = Prism(name="disperser prism", glass_type=["main","ZnSe","Marple"], beam_diameter=15, base=5,manual=True) #sf11 manual prism
# sf11 S coeffs: [2, 0.0, 1.7385, 0.0136, 0.311, 0.0616, 1.1749, 121.92] #First element 2 is to specify Sllemeir-2 formula
# N-BK7 S coeffs: [2, 0.0, 1.12735, 0.0072, 0.1244, 0.0270, 0.827, 100.38] #First element 2 is to specify Sllemeir-2 formula
# ZnSe S coeffs: [2, 3, 1.9, 0.113, 0.0, 0.0, 0.0, 0] #First element 2 is to specify Sllemeir-2 formula

camera = Lens(name="camera", focal_length=100, diameter=50)
collimator = Lens(name="collimator", focal_length=200, diameter=50)
slit = Slit(name="slit", width=0.03, height=0.2) #mm @ZimMAIN 60" = 1.629 mm and 5" = 0.136 mm, 1" = 27 um
# slit = Slit(name="slit", width=0.492, height=5.905) #mm @SST 60" = 5.905 mm and 5" = 0.492 mm, 1" = 98.4 um

df_mapping_list = computeCD(spatial_centers=centers_list, spectral_range=spectral_range, n_spectral=n_spectral, echelle=ech, disperser=disp, camera_lens=camera, collimator_lens=collimator, cmosx_max=cmosx_max, cmosy_max=cmosy_max, alpha_deg=15, slit=slit,write_hmtl=save_file)
# df_mapping_list = computeCD(spatial_centers=centers_list, spectral_range=spectral_range, n_spectral=n_spectral, echelle=ech, disperser=disp_prism, camera_lens=camera, collimator_lens=collimator, cmosx_max=cmosx_max, cmosy_max=cmosy_max, alpha_deg=15, slit=slit,write_hmtl=save_file)
df_sample = df_mapping_list[0]
# df_sample.info()
# df_sample.loc[(df_sample["wavelengths"]>=853.2) & (df_sample["wavelengths"]<= 855.2)]





