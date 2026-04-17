###-------------------- Import Modules --------------------###
from datetime import datetime
from matplotlib import gridspec
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import numpy as np
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
import cv2

import scienceplots
import astropy.io.fits as fits
import skimage


plt.style.use('science')
plt.style.use(['science','no-latex'])
plt.rcParams["font.weight"] = "bold"


###-------------------- Classes --------------------###
# Base class for optical elements
class OpticalElement:
    def __init__(self, name: str):
        self.name = name
    
    def __str__(self):
        return f"Optical Element: {self.name}"
    

# Class for slits
class Slit(OpticalElement):
    def __init__(self, name: str, width: float, height: float):
        super().__init__(name)
        self.width = width #mm
        self.height = height #mm
    
    def __str__(self):
        return f"Slit: {self.name}, Width: {self.width} mm, Height: {self.height} mm"
    

# Class for lenses
class Lens(OpticalElement):
    def __init__(self, name: str, focal_length: float, diameter: float):
        super().__init__(name)
        self.focal_length = focal_length
        self.diameter = diameter
    
    def __str__(self):
        return f"Lens: {self.name}, Focal Length: {self.focal_length} mm, Diameter: {self.diameter} mm"


#Class for Gratings
class Grating(OpticalElement):
    def __init__(self, name: str, groove_density: float, alpha: float = 0, m_order: int = 1):
        super().__init__(name)
        self.groove_density = groove_density  # in grooves per mm
        self.alpha = alpha # in degrees
        self.diffraction_order = m_order
    
    def __str__(self):
        return f"Grating: {self.name}, Groove Density: {self.groove_density} grooves/mm"
    
    def compute_angularDisp(self,wavelength_nm: np.ndarray): #dependent on diffraction order
        if np.abs(self.groove_density*self.diffraction_order*wavelength_nm.any()/1e6 - np.sin(np.deg2rad(self.alpha))) < 1:
            A = (self.groove_density/1e3 * self.diffraction_order)/(np.cos(np.arcsin(self.groove_density/1e6*self.diffraction_order*wavelength_nm - np.sin(np.deg2rad(self.alpha))))) #mrad/nm
        else:
            print("diffraction order is not possible")
        return A

# Class for prisms
class Prism(OpticalElement):

    def __init__(self, name: str, glass_type: list, beam_diameter: float, base:float, manual: bool = False):
        super().__init__(name)
        self.glass_type = glass_type #expected as in refractiveindex database : [shelf, book, page]
        self.beam_diameter= beam_diameter  # in mm
        self.base = base # base length in mm
        self.manual = manual

    
    def __str__(self):
        return f"Prism: {self.name}, Base length [mm]: {self.base}, Beam diameter [mm]: {self.beam_diameter}, Glass Type: Shelf={self.glass_type[0]},Book={self.glass_type[1]},Page={self.glass_type[2]}"
    

    def GetScoeffs(self, verbose: bool = True, sellmeierCoeffs: list =[1,0,0,0,0,0,0]):
        if self.manual == False:
            shelf = self.glass_type[0]
            book = self.glass_type[1]
            page = self.glass_type[2]
            glass = ri.RefractiveIndexMaterial(shelf=shelf,book=book,page=page)
            coeffs = glass._coefs


            if verbose:
                print(f"Sellmeier coefficients for {shelf} {book} {page}: {coeffs}")
            
        else: #manual==True so manual input of Sellmeier coefficients
            coeffs = sellmeierCoeffs
            coeffs = coeffs.replace("[","").replace("]","").split(",")
            coeffs = [float(i) for i in coeffs]
            sellmeier_formula = coeffs[0]
            coeffs = coeffs[1:]

            if sellmeier_formula == 1.0:
                self.coeffs = coeffs
            elif sellmeier_formula == 2.0:
                coeffs[2] = np.sqrt(coeffs[2])
                coeffs[4] = np.sqrt(coeffs[4])
                coeffs[6] = np.sqrt(coeffs[6])
                self.coeffs = coeffs
            if verbose:
                print(f"Manual Sellmeier coefficients: {coeffs}\nSellmeier formula: {sellmeier_formula}")
        
        n_coeffs = len(coeffs)
        if n_coeffs != 7:
                print(f"Check Sellmeier coefficients in database or input, number of coeffs is {n_coeffs}")
        
        return coeffs
        
    def Sellmeier(self, coeffs,wavelengths_nm):
        wavelengths_um = wavelengths_nm / 1000 # nm to um
        n = np.sqrt(1 + coeffs[0] + coeffs[1]*wavelengths_um**2/(wavelengths_um**2-coeffs[2]**2) + coeffs[3]*wavelengths_um**2/(wavelengths_um**2-coeffs[4]**2)+ coeffs[5]*wavelengths_um**2/(wavelengths_um**2-coeffs[6]**2))
        return n
    

    def DerivativeNwl_per_um(self,wavelengths_nm): #returns um⁻1
        wavelengths_um = wavelengths_nm / 1000 # nm to um
        coeffs = self.GetScoeffs()
        n = self.Sellmeier(coeffs=coeffs, wavelengths_nm=wavelengths_nm)
        du_dwl = -2 * coeffs[1] * coeffs[2]**2 * np.power(wavelengths_um,-3) / (np.power(1 - np.power(coeffs[2]/wavelengths_um,2),2)) -2 * coeffs[3] * coeffs[4]**2 * np.power(wavelengths_um,-3) / (np.power(1 - np.power(coeffs[4]/wavelengths_um,2),2)) -2 * coeffs[5] * coeffs[6]**2 * np.power(wavelengths_um,-3) / (np.power(1 - np.power(coeffs[6]/wavelengths_um,2),2))
        dn_dwl = 1 / (2*np.sqrt(np.power(n,2))) * du_dwl 
        return np.abs(dn_dwl) #um⁻1, not to get inversed corss dispersion pattern
        

    def compute_angularDisp(self,wavelength_nm: np.ndarray):
        A = self.base / self.beam_diameter * self.DerivativeNwl_per_um(wavelengths_nm=wavelength_nm) #rad/um or mrad/nm
        return A 

    

# Class for Echelle gratings
class EchelleGrating(Grating):
    def __init__(self, name: str, groove_density: float, blaze_angle: float, semi_deviation_angle_deg: float):
        super().__init__(name,groove_density)
        self.blaze_angle = blaze_angle
        self.semi_deviation_angle_deg = semi_deviation_angle_deg
    
    def __str__(self):
        return f"EchelleGrating: {self.name}, Groove Density: {self.groove_density} grooves/mm, Blaze angle: {self.blaze_angle} deg, Semi-dev angle: {self.semi_deviation_angle_deg} deg"

    def compute_angularDispE(self,wavelength_nm: float):
        A = (2 * np.sin(np.deg2rad(self.blaze_angle))*np.cos(np.deg2rad(self.semi_deviation_angle_deg)))/(wavelength_nm/1e3 * np.cos(np.arcsin(self.groove_density*self.diffraction_order*wavelength_nm/1e6 - np.sin(np.deg2rad(self.blaze_angle + self.semi_deviation_angle_deg))))) #expected to be in mrad/nm (rad/um)
        # u = (self.groove_density * self.diffraction_order/1e3)/(2*np.sin(np.deg2rad(self.blaze_angle)))
        # A = (u)/(np.sqrt(1-np.power((u*wavelength_nm/1e3),2))) #expected to be in mrad/nm (rad/um)
        return A
        
    def compute_diffractionorder(self, blazewavelength: float):
        m = 2 * np.sin(np.deg2rad(self.blaze_angle)*np.cos(np.deg2rad(self.semi_deviation_angle_deg))) / (self.groove_density*blazewavelength/1e6)
        return m

    def compute_blazewavelength(self, diffraction_order: np.ndarray):
        blaze_wavelength_nm  = 2 * 1e6 * np.ones(np.shape(diffraction_order)) * np.sin(np.deg2rad(self.blaze_angle)*np.cos(np.deg2rad(self.semi_deviation_angle_deg))) / (self.groove_density*diffraction_order)
        return blaze_wavelength_nm
    
    def compute_FSR(self, blazewavelength_array: np.ndarray):
        FSR_nm = (np.ones(np.shape(blazewavelength_array)) * self.groove_density/1e6 * np.power(blazewavelength_array,2)) / (2 * np.sin(np.deg2rad(self.blaze_angle)*np.cos(np.deg2rad(self.semi_deviation_angle_deg))))
        return FSR_nm

    def compute_exitAngle(self, wavelength_nm: float | np.ndarray, m_diffraction_order_array: np.ndarray | int, index_FSR: np.ndarray):
        if np.isscalar(wavelength_nm) and np.isscalar(m_diffraction_order_array):
            exit_angle_rad = np.arcsin(self.groove_density/1e6*m_diffraction_order_array*wavelength_nm - np.sin(np.deg2rad(self.blaze_angle + self.semi_deviation_angle_deg)))
        else:
            exit_angle_rad = np.arcsin(self.groove_density/1e6*m_diffraction_order_array[np.maximum(index_FSR-1,0)]*wavelength_nm - np.sin(np.deg2rad(self.blaze_angle + self.semi_deviation_angle_deg)))
        return exit_angle_rad

class Camera_sensor(OpticalElement):
    def __init__(self, name: str, pixels_x: int, pixels_y: int, pixel_size: float):
        super().__init__(name)
        self.px_x = pixels_x
        self.px_y = pixels_y
        self.px_size = pixel_size
        self.size_x_mm, self.size_y_mm = self.getCameraSensorSize()

    def getCameraSensorSize(self):
        size_x_mm = self.px_x * self.px_size / 1000 #mm
        size_y_mm = self.px_y * self.px_size / 1000 #mm
        return size_x_mm, size_y_mm
        
    
class Instrument(OpticalElement):
    def __init__(self, name: str, slit: Slit, echelle: EchelleGrating, disperser: Grating | Prism, collimator_lens: Lens, camera_lens: Lens, camera_sensor: Camera_sensor, spatial_centers: list, wavelength_scan_width_nm: float, spectral_range: tuple, spectral_res_nm: float):
        super().__init__(name)
        self.slit = slit
        self.echelle = echelle
        self.disperser = disperser
        self.collimator_lens = collimator_lens
        self.camera_lens = camera_lens
        self.camera_sensor = camera_sensor
        self.spatial_centers = spatial_centers
        self.wavelength_scan_width_nm = wavelength_scan_width_nm
        
        # initializing methods
        self.setSpectralRange(spectral_range=spectral_range)
        self.setSpectralRes(spectral_res=spectral_res_nm)
        self.getSpectralLines()
        self.__createDFcmos()

        # computing dataset for each spatial center and creating mapping dataframes
        if hasattr(self, "spectral_range") and hasattr(self, "spectral_res"):
            df_mapping_list = [None] * len(self.spatial_centers)
            for spatial_center in self.spatial_centers:
                dataset = self.computeCD(spectral_range_nm=self.spectral_range, spectral_res_nm=self.spectral_res, spatial_center=spatial_center)
                df_mapping = self.createDFmapping(dataset=dataset, spatial_center=spatial_center)
                df_mapping_list[self.spatial_centers.index(spatial_center)] = df_mapping
            self.df_mapping_list = df_mapping_list
        else:
            raise ValueError("Spectral range and spectral resolution must be set before computing the mapping. Please set them using setSpectralRange and setSpectralRes methods.")            
    


    ###-------------------- Backend Methods --------------------###
    def getSpectralLines(self)-> dict:
        """
        getSpectralLines returns a hidden dictionary of spectral lines with their respective names/species associated with their wavelengths in nm.
        
        :param self: Instrument class object
        :return: dict(str: list) where the key is the name/species of the spectral line and the value is a list of wavelengths in nm associated with that line/species.
        """
        # If getSpectralLines was called before, return the already stored dictionary
        if hasattr(self, "_Instrument__spectral_lines"):
            return self.__spectral_lines
        
        # Neon I NIST spectral lines
        neon_lines = np.load("./NIST_Atomic-Specie/Neon.npy").squeeze().tolist()
        # Thorium I NIST spectral lines
        thorium_lines = np.load("./NIST_Atomic-Specie/Thorium0.95.npy").squeeze().tolist()
        # Argon I NIST spectral lines
        argon_lines = np.load("./NIST_Atomic-Specie/Argon0.95.npy").squeeze().tolist()


        self.__spectral_lines = {"Mg I b3 5167": [516.7] ,"Mg I b2 5172": [517.2] ,"Mg I b1 5183": [518.3] ,"Fe I 5250": [525.0],"Mn I 5399": [539.9], "He I D3 5876":[587.6], "Na I D2 5890": [589.0] ,"Na I D1 5896": [589.6],"Fe I 6173 (HMI)": [617.3], "Fe I 6301-6302": [630.15],"Ca I 6439": [643.9], "Ha": [656.3], "Ni I 6643": [664.3],"Fe I/Ca I 6718": [671.8],"K I 7699": [769.9],"Ca II 8498": [849.8], "Ca II 8542": [854.2],"Ca II 8662": [866.2], "Neon": neon_lines, "Thorium": thorium_lines, "Argon": argon_lines} #to be continued, rn only wl > 500 nm

        return self.__spectral_lines

    def setSpectralRange(self, spectral_range: tuple) -> None:
        self.spectral_range = spectral_range
        
        return None


    def setSpectralRes(self, spectral_res: float) -> None:
        self.spectral_res = spectral_res
        
        return None
    
    # NEED TO BE TESTED
    def exportDFmapping(self, df_mapping_list: list, filename: str) -> None:
        current_datetime = datetime.today().strftime('%Y-%m-%d_%H-%M-%S')
        for i, df_mapping in enumerate(df_mapping_list, start=1):
            df_mapping.to_csv(f"{filename}_spatial_{i}_{current_datetime}.csv", index=False, sep="\t")
        
        return None
    
    def exportAsImage(self, species: str, filename: str, spectral_range_nm: tuple, path: str = ".") -> None:
        if species not in self.__spectral_lines.keys():
            raise ValueError(f"Species {species} not found in spectral lines dictionary. Please check the available species and their respective wavelengths.")
            return None
        
        wavelengths = np.array([wavelength for wavelength in self.__spectral_lines[species] if spectral_range_nm[0] <= wavelength <= spectral_range_nm[1]])
        df_mapping_array = np.array([])
        # df_mapping_array = np.zeros((len(wavelengths)*len(self.spatial_centers), 2))

        for spatial_center in self.spatial_centers:
            dataset = self.computeCD(spatial_center=spatial_center, isArrayDefined=True, defined_spectral_array=np.array(wavelengths), spectral_range_nm=(None, None), spectral_res_nm=None)
            df_mapping = self.createDFmapping(dataset=dataset, spatial_center=spatial_center)
            df_mapping = df_mapping["X"].to_frame().join(df_mapping["Y"])
            df_mapping_array = np.concatenate((df_mapping_array, df_mapping.to_numpy()), axis=0) if df_mapping_array.size else df_mapping.to_numpy()
        df_mapping_array_pixels = np.round(df_mapping_array / self.camera_sensor.px_size * 1000, decimals=0).astype(int) #mm to um to pixels
        cmos_simulated_image = np.zeros((self.camera_sensor.px_y, self.camera_sensor.px_x))
        cmos_simulated_image[df_mapping_array_pixels[:,1].astype(int), df_mapping_array_pixels[:,0].astype(int)] = 1 #setting the pixels corresponding to the spectral lines to 1
        slit_size_cmos_pxl = [self.computeMagnification()[0] * self.slit.width / self.camera_sensor.px_size * 1000, self.computeMagnification()[1] * self.slit.height / self.camera_sensor.px_size * 1000] #mm to um to pixels
        kernel_slit = np.ones((int(slit_size_cmos_pxl[1]), int(slit_size_cmos_pxl[0])), dtype=np.uint8) #kernel for dilation to simulate slit width on image
        cmos_simulated_image = cv2.dilate(cmos_simulated_image, kernel_slit, iterations=1)
        cmos_simulated_image_show = Data(path="", isFits=False)
        cmos_simulated_image_show.name = filename
        cmos_simulated_image_show.data = cmos_simulated_image
        cmos_simulated_image_show.showImage(save=True, vmin=0, vmax=1, path=path)
        cv2.imwrite(os.path.join(path, filename + "_raw.png"), np.flipud(cmos_simulated_image)*255) #saving the image as png, multiplying by 255 to get values between 0 and 255 for uint8 format
        
        return None

    def exportAsFits(self, species: str, filename: str, spectral_range_nm: tuple, path: str = ".", wantSlitKernel: bool = False) -> None:
        if species not in self.__spectral_lines.keys():
            raise ValueError(f"Species {species} not found in spectral lines dictionary. Please check the available species and their respective wavelengths.")
            
        wavelengths = np.array([wavelength for wavelength in self.__spectral_lines[species] if spectral_range_nm[0] <= wavelength <= spectral_range_nm[1]])
        df_mapping_array = np.array([])
        # df_mapping_array = np.zeros((len(wavelengths)*len(self.spatial_centers), 2))

        for spatial_center in self.spatial_centers:
            dataset = self.computeCD(spatial_center=spatial_center, isArrayDefined=True, defined_spectral_array=np.array(wavelengths), spectral_range_nm=(None, None), spectral_res_nm=None)
            df_mapping = self.createDFmapping(dataset=dataset, spatial_center=spatial_center)
            df_mapping = df_mapping["X"].to_frame().join(df_mapping["Y"]).join(df_mapping["wavelengths"]).join(df_mapping["angular_dispersion_x[mrad/nm]"])
            df_mapping_array = np.concatenate((df_mapping_array, df_mapping.to_numpy()), axis=0) if df_mapping_array.size else df_mapping.to_numpy()

        # base image array
        df_mapping_array_pixels = np.concatenate((np.round(df_mapping_array[:,:2] / self.camera_sensor.px_size * 1000, decimals=0).astype(int), df_mapping_array[:, 2:4]), axis=1) #mm to um to pixels + associated wavelengths + associated angulare dispersion in mrad/nm
        cmos_simulated_image = np.zeros((self.camera_sensor.px_y, self.camera_sensor.px_x))
        cmos_simulated_image[df_mapping_array_pixels[:,1].astype(int), df_mapping_array_pixels[:,0].astype(int)] = 1 #setting the pixels corresponding to the spectral lines to 1
        
        # creation of wavelength mask from simulated image
        cmos_wavelength_mask = np.zeros((self.camera_sensor.px_y, self.camera_sensor.px_x))
        cmos_wavelength_mask[df_mapping_array_pixels[:,1].astype(int), df_mapping_array_pixels[:,0].astype(int)] = df_mapping_array_pixels[:,2] # setting the pixels corresponding to the spectral lines to their associated wavelengths

        cmos_dispersion_mask = np.zeros((self.camera_sensor.px_y, self.camera_sensor.px_x))
        cmos_dispersion_mask[df_mapping_array_pixels[:,1].astype(int), df_mapping_array_pixels[:,0].astype(int)] = df_mapping_array_pixels[:,3] # setting the pixels corresponding to the spectral lines to their associated angular dispersion in mrad/nm
        
        if wantSlitKernel:
            slit_size_cmos_pxl = [self.computeMagnification()[0] * self.slit.width / self.camera_sensor.px_size * 1000, self.computeMagnification()[1] * self.slit.height / self.camera_sensor.px_size * 1000] #mm to um to pixels
            kernel_slit = np.ones((int(slit_size_cmos_pxl[1]), int(slit_size_cmos_pxl[0])), dtype=np.uint8) #kernel for dilation to simulate slit width on image
            cmos_simulated_image = cv2.dilate(cmos_simulated_image, kernel_slit, iterations=1)
            cmos_wavelength_mask = cv2.dilate(cmos_wavelength_mask, kernel_slit, iterations=1)
            cmos_dispersion_mask = cv2.dilate(cmos_dispersion_mask, kernel_slit, iterations=1)
        # FITS file creation
        header_dict = {"target": species, "range_nm": spectral_range_nm}
        header = fits.Header()
        for key, value in header_dict.items():
            header[key] = value
        primary_hdu = fits.PrimaryHDU(data=None, header=header)
        main_image_hdu = fits.ImageHDU(data=cmos_simulated_image, name="SIMULATED_IMAGE")
        cmos_wavelength_mask_hdu = fits.ImageHDU(data=cmos_wavelength_mask, name="WAVELENGTH_MASK")
        cmos_dispersion_mask_hdu = fits.ImageHDU(data=cmos_dispersion_mask, name="DISPERSION_MASK")
        hdul = fits.HDUList([primary_hdu, main_image_hdu, cmos_wavelength_mask_hdu, cmos_dispersion_mask_hdu])
        hdul.writeto(os.path.join(path, filename + ".fits"), overwrite=True)
        
        return None
    
    def __ComputeX(self, cmosx0: float, spectral_array: np.ndarray):
        
        max_order = ceil(self.echelle.compute_diffractionorder(spectral_array[0]))
        min_order = floor(self.echelle.compute_diffractionorder(spectral_array[-1]))

        diffraction_order_array = np.arange(start=min_order, stop=max_order+1, step=1)
        blaze_wavelength_array = self.echelle.compute_blazewavelength(diffraction_order=diffraction_order_array)
        FSR_array = self.echelle.compute_FSR(blazewavelength_array=blaze_wavelength_array)
        angular_dispersion_array = self.echelle.compute_angularDispE(wavelength_nm=blaze_wavelength_array)
        FSR_bins = blaze_wavelength_array + 0.5 * FSR_array

        index_FSR = np.digitize(x=spectral_array, bins=FSR_bins)

        # cmosx = self.camera_lens.focal_length * angular_dispersion_array[np.maximum(index_FSR-1,0)]/1000 * (spectral_array - blaze_wavelength_array[np.maximum(index_FSR-1,0)]) + np.ones(np.shape(spectral_array)) * cmosx0 #mm (mrad.nm⁻1 to rad.nm⁻1 was done with /1000)
        lambda0 = np.median(spectral_array) 
        idx = np.argmin(np.abs(spectral_array - lambda0))
        m0 = index_FSR[idx]
        # beta0_rad = self.echelle.compute_exitAngle(wavelength_nm=lambda0, m_diffraction_order_array=m0, index_FSR=index_FSR) #rad])
        beta0_rad = self.echelle.compute_exitAngle(wavelength_nm=float(np.squeeze(self.echelle.compute_blazewavelength(diffraction_order=np.array([m0])))), m_diffraction_order_array=m0, index_FSR=index_FSR) #rad 
        # both betas are always of same sign
        if beta0_rad < 0:
            beta_rad_reduced = self.echelle.compute_exitAngle(wavelength_nm=spectral_array, m_diffraction_order_array=diffraction_order_array, index_FSR=index_FSR) + np.ones_like(spectral_array) * beta0_rad
        else:
            beta_rad_reduced = self.echelle.compute_exitAngle(wavelength_nm=spectral_array, m_diffraction_order_array=diffraction_order_array, index_FSR=index_FSR) - np.ones_like(spectral_array) * beta0_rad

        
        cmosx = self.camera_lens.focal_length * beta_rad_reduced + np.ones_like(spectral_array) * cmosx0 #mm 
        
        return cmosx,diffraction_order_array[np.maximum(index_FSR-1,0)],blaze_wavelength_array[np.maximum(index_FSR-1,0)], angular_dispersion_array[np.maximum(index_FSR-1,0)], FSR_array[np.maximum(index_FSR-1,0)]
    

    def __ComputeY(self, cmosy0: float, spectral_array: np.ndarray):
        
        # center wavelength of spectral_array
        lambda0 = (spectral_array[-1] + spectral_array[0]) / 2 
        # differential spectral array with respect to center wavelength
        dspectral_array = spectral_array - np.ones(np.shape(spectral_array)) * lambda0
        
        #optical system computations
        if self.disperser.__class__.__name__ == "Prism":
            Ac = self.disperser.compute_angularDisp(wavelength_nm=spectral_array) /1000 #mrad.nm⁻1 to rad.nm⁻1
            
        else:
            Ac = self.disperser.compute_angularDisp(wavelength_nm=dspectral_array) /1000 #mrad.nm⁻1 to rad.nm⁻1
        #computation y coordinate
        cmosy = np.ones(np.shape(dspectral_array))* cmosy0 + self.camera_lens.focal_length * Ac * dspectral_array
             
        return cmosy, Ac
    

    def computeMagnification(self):# r for crossdisperser ?
        r = (1 - np.tan(np.deg2rad(self.echelle.blaze_angle))*np.tan(np.deg2rad(self.echelle.semi_deviation_angle_deg)))/(1 + np.tan(np.deg2rad(self.echelle.blaze_angle))*np.tan(np.deg2rad(self.echelle.semi_deviation_angle_deg)))
        demag_y = self.camera_lens.focal_length / self.collimator_lens.focal_length
        demag_x = r * self.camera_lens.focal_length / self.collimator_lens.focal_length

        return demag_x, demag_y, r
    

    def __createDFcmos(self):
        # dataframe for camera sensor array
        df_cmos = pd.DataFrame({"X": [0, self.camera_sensor.size_x_mm, self.camera_sensor.size_x_mm, 0, 0], "Y": [0, 0, self.camera_sensor.size_y_mm, self.camera_sensor.size_y_mm, 0]})
        self.df_cmos = df_cmos

        return df_cmos
    

    def computeCD(self,spectral_range_nm: tuple, spectral_res_nm: float, spatial_center: tuple, isArrayDefined: bool = False, defined_spectral_array: np.ndarray = None):
        # definition of array of wavelengths to compute within spectral range and explicit definition of spatial center on cmos
        if isArrayDefined == False:
            spectral_array = np.arange(start=spectral_range_nm[0], stop=spectral_range_nm[1]+spectral_res_nm, step=spectral_res_nm)# in nm
        else:
            spectral_array = defined_spectral_array
        x0,y0 = spatial_center
        
        # slit size on cmos
        demag_x, demag_y, r = self.computeMagnification()
        slit_width_cmos, slit_height_cmos = self.slit.width * demag_x, self.slit.height * demag_y
        
        # compute Y coordinates (slit center)
        cmosx,diffraction_order_array,blaze_wavelength_array, angular_dispersion_x,FSR_array = self.__ComputeX(cmosx0=x0,spectral_array=spectral_array)
        
        # compute Y coordinates (slit center) 
        cmosy, angular_dispersion_y = self.__ComputeY(cmosy0=y0,spectral_array=spectral_array)

        # compute spectral resolution and theoritical R
        max_spectral_res = np.ones(np.shape(angular_dispersion_x))*((slit_width_cmos*1000) / self.camera_lens.focal_length) / angular_dispersion_x
        max_theoritical_R = spectral_array / max_spectral_res 

        return spectral_array, cmosx, cmosy, diffraction_order_array, blaze_wavelength_array, FSR_array, angular_dispersion_x, angular_dispersion_y, max_spectral_res, max_theoritical_R, slit_width_cmos, slit_height_cmos
    

    def createDFmapping(self, dataset: tuple, spatial_center: tuple):
        
        # unpacking dataset tuple created by computeCD method
        spectral_array, cmosx, cmosy, diffraction_order_array, blaze_wavelength_array, FSR_array, dispersion_x, dispersion_y, max_spectral_res, max_theoritical_R, slit_width_cmos, slit_height_cmos = dataset
        
        # assign color to each wavelength for plotting, using getColor method
        color_list = self.getColor(wavelength=spectral_array, gamma=0.8)
        
        #fitting shape of dict for centers
        spatial_centers_array = np.shape(cmosx)[0] * [spatial_center]
        slit_width_cmos_array = np.shape(cmosx)[0] * [slit_width_cmos*1000] #mm to um
        slit_height_cmos_array = np.shape(cmosx)[0] * [slit_height_cmos*1000] #mm to um
       
        #dataframes
        df_mapping = pd.DataFrame({"Center" : spatial_centers_array,"X": cmosx, "Y": cmosy, "wavelengths": spectral_array, "echelle_diffraction_order": diffraction_order_array, "echelle_blaze_wavelength[nm]": blaze_wavelength_array, "angular_dispersion_x[mrad/nm]": dispersion_x, "angular_dispersion_y[mrad/nm]": dispersion_y, "color":color_list, "Max Spectral resolution [nm]" : max_spectral_res, "Max Theoritical R" : max_theoritical_R,"Slit width [um]":slit_width_cmos_array, "Slit height [um]": slit_height_cmos_array, "Free Spectral Range [nm]":FSR_array})#, #add colors CIE in booktab

        
        return df_mapping
    

    ###-------------------- Frontend Methods --------------------###
    def getColor(self,wavelength, gamma=0.8):
        #
        #    Based on code by Dan Bruton
        #    http://www.physics.sfasu.edu/astro/color/spectra.html
        #    '''
        # Adapted from <script src="https://gist.github.com/friendly/67a7df339aa999e2bcfcfec88311abfc.js"></script>
        wl = np.asarray(wavelength, dtype=float)

        R = np.zeros_like(wl)
        G = np.zeros_like(wl)
        B = np.zeros_like(wl)

        # 380–440 nm
        m = (wl >= 380) & (wl <= 440) # m boolean mask
        att = 0.3 + 0.7 * (wl[m] - 380) / (440 - 380)
        R[m] = (-(wl[m] - 440) / (440 - 380) * att) ** gamma
        B[m] = (1.0 * att) ** gamma

        # 440–490 nm
        m = (wl > 440) & (wl <= 490)
        G[m] = ((wl[m] - 440) / (490 - 440)) ** gamma
        B[m] = 1.0

        # 490–510 nm
        m = (wl > 490) & (wl <= 510)
        G[m] = 1.0
        B[m] = (-(wl[m] - 510) / (510 - 490)) ** gamma

        # 510–580 nm
        m = (wl > 510) & (wl <= 580)
        R[m] = ((wl[m] - 510) / (580 - 510)) ** gamma
        G[m] = 1.0

        # 580–645 nm
        m = (wl > 580) & (wl <= 645)
        R[m] = 1.0
        G[m] = (-(wl[m] - 645) / (645 - 580)) ** gamma

        # 645–750 nm
        m = (wl > 645) & (wl <= 750)
        att = 0.3 + 0.7 * (750 - wl[m]) / (750 - 645)
        R[m] = (1.0 * att) ** gamma

        # Scale to [0,255]
        R = np.clip(R * 255, 0, 255).astype(np.uint8)
        G = np.clip(G * 255, 0, 255).astype(np.uint8)
        B = np.clip(B * 255, 0, 255).astype(np.uint8)

        # Convert to hex strings
        colors = [f"#{r:02x}{g:02x}{b:02x}" for r, g, b in zip(R, G, B)]
        return colors
    

    def plotCD(self):
        slit_width_cmos, slit_height_cmos = self.df_mapping_list[0]["Slit width [um]"].iloc[0]/1000, self.df_mapping_list[0]["Slit height [um]"].iloc[0]/1000 #um to mm 
        #-----Plotting Parameters------#
        lines = list(self.__spectral_lines.keys())
        selection_spectral_windows = {
                                        k: pd.concat(
                                            [
                                                self.df_mapping_list[0].loc[
                                                    (self.df_mapping_list[0]["wavelengths"] <= v + self.wavelength_scan_width_nm / 2) &
                                                    (self.df_mapping_list[0]["wavelengths"] >= v - self.wavelength_scan_width_nm / 2)
                                                ]
                                                for v in v_list
                                            ]
                                        )
                                        for k, v_list in self.__spectral_lines.items()
                                    }

        #-----Actual Plotting------#
        app = dash.Dash(__name__)
        app.title = "SIESTA - Simulating Interactive Echelle Spectrogram for Targeted Applications"

        app.layout = html.Div([
            html.H2("SIESTA - Simulating Interactive Echelle Spectrogram for Targeted Applications", style={'color': 'blue'}),
            html.Div([
            html.Button("Show All", id='show-all-button', n_clicks=0, style={'marginRight': '10px'}),
            html.Button("Clear All", id='clear-all-button', n_clicks=0)
    ], style={'marginTop': '10px',"marginBottom": '10px'}),
            dcc.Dropdown(
                id='group-selector',
                options=[{'label': k, 'value': k} for k in selection_spectral_windows.keys()],
                multi=True,
                placeholder="Select lines to highlight"
            ),

            dcc.Graph(id='grid-plot', config={'displayModeBar': False}),
            

            html.Div(id='meta-info', style={'marginTop': '20px', 'fontSize': '16px'})
        ])



         # Callback to update dropdown values
        @app.callback(
            Output('group-selector', 'value'),
            Input('show-all-button', 'n_clicks'),
            Input('clear-all-button', 'n_clicks'),
            prevent_initial_call=True
        )
        def update_dropdown(show_clicks, clear_clicks):
            # Use Dash callback context to determine which button was clicked
            ctx = dash.callback_context
            if not ctx.triggered:
                raise dash.exceptions.PreventUpdate

            triggered_id = ctx.triggered[0]['prop_id'].split('.')[0]

            if triggered_id == 'show-all-button':
                return list(self.__spectral_lines.keys())
            elif triggered_id == 'clear-all-button':
                return []
            # else:
            #     raise dash.exceptions.PreventUpdate
            
        @app.callback(
            Output('grid-plot', 'figure'),
            Input('group-selector', 'value'),
            prevent_initial_call=True
        )

        def update_plot(selected_groups):
            fig = go.Figure()
            if not selected_groups:
                # Plot all mappings
                for i, df_mapping in enumerate(self.df_mapping_list, start=1):
                    hovertexts = [
                        f"Wavelength: {wl}<br>Echelle order: {order}<br>Max resolution [A]: {10 * res}<br>Max R: {R}<br>Slit [um](x, y): ({slit_x:.0f}, {slit_y:.0f})"
                        for wl, order, res, R, slit_x, slit_y in zip(
                            df_mapping["wavelengths"],
                            df_mapping["echelle_diffraction_order"],
                            df_mapping["Max Spectral resolution [nm]"],
                            df_mapping["Max Theoritical R"],
                            df_mapping["Slit width [um]"],
                            df_mapping["Slit height [um]"]
                        )
                    ]
                    
                    fig.add_scattergl(
                        name=f"spatial {i}",
                        x=df_mapping["X"],
                        y=df_mapping["Y"],
                        mode="markers",
                        hovertext=hovertexts,
                        marker=dict(
                            size=8,
                            color=df_mapping["color"]
                        ),
                        error_x=dict(
                            type='data',
                            array=np.ones(len(df_mapping)) * slit_width_cmos / 2,
                            visible=True,
                            color="white"
                        ),
                        error_y=dict(
                            type='data',
                            array=np.ones(len(df_mapping)) * slit_height_cmos / 2,
                            visible=True,
                            color="white"
                        ),
                        opacity=0.05
                    )

                # CMOS rectangle
                fig.add_scattergl(
                    x=self.df_cmos["X"],
                    y=self.df_cmos["Y"],
                    mode="lines",
                    name=f"{self.camera_sensor.name}: {self.camera_sensor.size_x_mm}x{self.camera_sensor.size_y_mm}mm",
                    line=dict(color="red", width=2)
                )

                # Layout details
                if self.disperser.__class__.__name__ == "Prism":
                    title = f"""
                    Camera focal length : {self.camera_lens.focal_length} mm
                    Echelle groove density : {self.echelle.groove_density:.1f} mm⁻¹<br>
                    Echelle blaze angle : {self.echelle.blaze_angle:.1f}°<br>
                    Deviation angle from Littrow : {2 * self.echelle.semi_deviation_angle_deg:.1f}°<br>
                    Cross-disperser Prism: {self.disperser.base} mm
                    """
                else:
                    title = f"""
                    Camera focal length : {self.camera_lens.focal_length} mm, 
                    Echelle groove density : {self.echelle.groove_density:.1f} mm⁻¹, 
                    Echelle blaze angle : {self.echelle.blaze_angle:.1f}°, 
                    Deviation angle from Littrow : {2 * self.echelle.semi_deviation_angle_deg:.1f}°, 
                    Cross-disperser groove density : {self.disperser.groove_density:.1f} mm⁻¹
                    """#<br> was used before to \n in html

                fig.update_layout(
                    title=title,
                    width=1200,
                    height=800,
                    xaxis_title="X [mm]",
                    yaxis_title="Y [mm]",
                    xaxis=dict(scaleanchor="y", scaleratio=1),
                    showlegend=True
                )

                return fig
            else:
                
                # Plot only selected groups
                for i, df_mapping in enumerate(self.df_mapping_list, start=1):

                    hovertexts = [
                        f"Wavelength: {wl}<br>Echelle order: {order}<br>Max resolution [A]: {10 * res}<br>Max R: {R}<br>Slit [um](x, y): ({slit_x:.0f}, {slit_y:.0f})"
                        for wl, order, res, R, slit_x, slit_y in zip(
                            df_mapping["wavelengths"],
                            df_mapping["echelle_diffraction_order"],
                            df_mapping["Max Spectral resolution [nm]"],
                            df_mapping["Max Theoritical R"],
                            df_mapping["Slit width [um]"],
                            df_mapping["Slit height [um]"]
                        )
                    ]
                    
                    #background points
                    fig.add_scattergl(
                            name=f"spatial {i}",
                            x=df_mapping["X"],
                            y=df_mapping["Y"],
                            mode="markers",
                            hovertext=hovertexts,
                            marker=dict(
                                size=8,
                                color=df_mapping["color"]
                            ),
                            error_x=dict(
                                type='data',
                                array=np.ones(len(df_mapping)) * slit_width_cmos / 2,
                                visible=True,
                                color="white"
                            ),
                            error_y=dict(
                                type='data',
                                array=np.ones(len(df_mapping)) * slit_height_cmos / 2,
                                visible=True,
                                color="white"
                            ),
                            opacity=0.05
                        )
                    for group in selected_groups:
                        
                        df_sel = selection_spectral_windows[group]
                        hovertexts_hl = [
                        f"Wavelength: {wl}<br>Echelle order: {order}<br>Max resolution [A]: {10 * res}<br>Max R: {R}<br>Slit [um](x, y): ({slit_x:.0f}, {slit_y:.0f})"
                        for wl, order, res, R, slit_x, slit_y in zip(
                            df_sel["wavelengths"],
                            df_sel["echelle_diffraction_order"],
                            df_sel["Max Spectral resolution [nm]"],
                            df_sel["Max Theoritical R"],
                            df_sel["Slit width [um]"],
                            df_sel["Slit height [um]"]
                            )
                        ]
                        fig.add_scattergl(
                            x=df_sel["X"],
                            y=df_sel["Y"],
                            mode="markers",
                            marker=dict(
                            size=8,
                            color=df_sel["color"]
                            ),
                            error_x=dict(
                                type='data',
                                array=np.ones(len(df_sel)) * slit_width_cmos / 2,
                                visible=True,
                                color="white"
                            ),
                            error_y=dict(
                                type='data',
                                array=np.ones(len(df_sel)) * slit_height_cmos / 2,
                                visible=True,
                                color="white"
                            ),
                            name=f"spatial {i}: "+group,
                            hovertext=hovertexts_hl,
                            hoverinfo="text",
                            opacity=0.8
                        )

                # CMOS rectangle
                fig.add_scattergl(
                    x=self.df_cmos["X"],
                    y=self.df_cmos["Y"],
                    mode="lines",
                    name=f"{self.camera_sensor.name}: {self.camera_sensor.size_x_mm}x{self.camera_sensor.size_y_mm}mm",
                    line=dict(color="red", width=2)
                )

                # Layout details
                if self.disperser.__class__.__name__ == "Prism":
                    title = f"""
                    Camera focal length : {self.camera_lens.focal_length} mm<br>
                    Echelle groove density : {self.echelle.groove_density:.1f} mm⁻¹<br>
                    Echelle blaze angle : {self.echelle.blaze_angle:.1f}°<br>
                    Deviation angle from Littrow : {2 * self.echelle.semi_deviation_angle_deg:.1f}°<br>
                    Cross-disperser Prism: {self.disperser.base} mm
                    """
                else:
                    title = f"""
                    Camera focal length : {self.camera_lens.focal_length} mm | Echelle groove density : {self.echelle.groove_density:.1f} mm⁻¹ | Echelle blaze angle : {self.echelle.blaze_angle:.1f}° | Deviation angle from Littrow : {2 * self.echelle.semi_deviation_angle_deg:.1f}° <br>Cross-disperser groove density : {self.disperser.groove_density:.1f}mm⁻¹"""

                fig.update_layout(
                    title=title,
                    width=1200,
                    height=800,
                    xaxis_title="X [mm]",
                    yaxis_title="Y [mm]",
                    xaxis=dict(scaleanchor="y", scaleratio=1),
                    showlegend=True
                )

                return fig
        # Run the app
        app.run(debug=True,port=8050,jupyter_mode='tab') #jupyter_mode='tab' : opens automatically browser, 'external' not.
        

class Data():
    def __init__(self,path: str, isFits: bool=True) -> None:
        self.path = path
        self.name: str = path.split("/")[-1]
        self.header: dict[str, str]
        self.data: np.ndarray
        if isFits:
            self.header, self.data = self.load_data()
        return None

    def load_data(self) -> tuple[dict[str, str], np.ndarray]:
        with fits.open(self.path) as file:
            # isHduRead = False
            # for hdu in file:
            # header reading
            # if hdu.header is not None and not isHduRead:
            if file[0].header is not None:
                header: dict[str, str] = dict(file[0].header)
                # isHduRead = True
            else:
                raise ValueError("No header found in FITS file.")
            # data reading
            if file[0].data is not None:
                if file[0].data.ndim == 2:
                    data: np.ndarray = file[0].data
            else:
                raise ValueError("No data found in FITS file.")
        return header, data
    
    def showImage(self, figsize: tuple[int,int]=(8,8), fontsize: int=24, vmin: int=0, vmax: int=0, save: bool=False, path: str=None) -> tuple[plt.Figure, plt.Axes]:
        fig, ax = plt.subplots(1, 1, figsize=figsize)

        fig.suptitle(f"{self.name}", fontsize=fontsize, x=0.5, y=1.00, fontweight='bold')
        ax.set_xlabel('X [pixels]')
        ax.set_ylabel('Y [pixels]')

        im = ax.imshow(self.data, cmap='gray', vmin=vmin, vmax=vmax, origin='lower')
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = fig.colorbar(im, cax=cax)
        formatter = ticker.ScalarFormatter(useMathText=True)
        formatter.set_scientific(True)
        formatter.set_powerlimits((0, 0))
        cbar.ax.yaxis.set_major_formatter(formatter)
        cbar.set_label("Intensity [ADU]") 

        ax.tick_params(color='white', labelcolor='black', which="major", length=8,width=1)
        ax.tick_params(color='white', labelcolor='black', which="minor", length=4,width=1)
        # big ticks
        ax.xaxis.set_major_locator(MultipleLocator(1000))
        ax.xaxis.set_major_formatter('{x:.0f}')

        # For the minor ticks, use no labels; default NullFormatter.
        ax.xaxis.set_minor_locator(MultipleLocator(200))

        # Y-axis
        ax.set_ylabel("Pixels",fontsize=20)
        ax.tick_params(labelsize=20)


        # big ticks
        ax.yaxis.set_major_locator(MultipleLocator(1000))
        ax.yaxis.set_major_formatter('{x:.0f}')

        # For the minor ticks, use no labels; default NullFormatter.
        ax.yaxis.set_minor_locator(MultipleLocator(200))
        
        if save:
            fig.savefig(f"{os.path.join(path, self.name)}.png", dpi=300)
        # fig.show()

        return fig, ax
    
    def showYprofile(self, xpixel: int, figsize: tuple[int,int]=(8,8), fontsize: int=24, vmin: int=0, vmax: int=0, save: bool=False) -> None:
        grid = gridspec.GridSpec(1,2, width_ratios=[1,1])
        fig,ax = plt.figure(figsize=figsize), [plt.subplot(grid[0]), plt.subplot(grid[1])]
        fig.suptitle(f"y={xpixel} profile of {self.name}", fontsize=fontsize, x=0.5, y=1.00, fontweight='bold')

        ## Left plot: image with vertical line ##
        ax[0].imshow(self.data, cmap='gray')
        ax[0].set_xlabel('X [pixels]')
        ax[0].set_ylabel('Y [pixels]')

        im = ax[0].imshow(self.data, cmap='gray', vmin=vmin, vmax=vmax)
        ax[0].axvline(x=xpixel, color='red', linestyle='--', linewidth=1)

        ax[0].tick_params(color='white', labelcolor='black', which="major", length=8,width=1)
        ax[0].tick_params(color='white', labelcolor='black', which="minor", length=4,width=1)
        # big ticks
        ax[0].xaxis.set_major_locator(MultipleLocator(1000))
        ax[0].xaxis.set_major_formatter('{x:.0f}')

        # For the minor ticks, use no labels; default NullFormatter.
        ax[0].xaxis.set_minor_locator(MultipleLocator(200))
        # Y-axis
        ax[0].set_ylabel("Pixels",fontsize=20)
        ax[0].tick_params(labelsize=20)

        # big ticks
        ax[0].yaxis.set_major_locator(MultipleLocator(1000))
        ax[0].yaxis.set_major_formatter('{x:.0f}')
        # For the minor ticks, use no labels; default NullFormatter.
        ax[0].yaxis.set_minor_locator(MultipleLocator(200))
        


        ##-------- Right plot: Y profile at xpixel --------##
        ax[1].plot(self.data[:, xpixel], color='black')

        # ticks parameters
        ax[1].tick_params(color='black', labelcolor='black', which="major", length=8,width=1)
        ax[1].tick_params(color='black', labelcolor='black', which="minor", length=4,width=1)
        ax[1].tick_params(labelsize=20)
        
        # X-axis parameters
        ax[1].set_xlabel('Y [pixels]')
        ax[1].xaxis.set_major_locator(MultipleLocator(1000))
        ax[1].xaxis.set_major_formatter('{x:.0f}')
        ax[1].xaxis.set_minor_locator(MultipleLocator(200))

        # Y-axis parameters
        ax[1].set_ylabel('Intensity [ADU]')
        ax[1].yaxis.set_major_locator(MultipleLocator(50))
        ax[1].yaxis.set_major_formatter('{x:.0f}')
        ax[1].yaxis.set_minor_locator(MultipleLocator(10))

        if save:
            fig.savefig(f"{self.name}.png", dpi=300)
        # fig.show()
        
    
        return None    
