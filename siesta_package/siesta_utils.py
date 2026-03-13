###-------------------- Import Modules --------------------###
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
    def __init__(self, name: str, slit: Slit, echelle: EchelleGrating, disperser: Grating | Prism, collimator_lens: Lens, camera_lens: Lens, camera_sensor: Camera_sensor):
        super().__init__(name)
        self.slit = slit
        self.echelle = echelle
        self.disperser = disperser
        self.collimator_lens = collimator_lens
        self.camera_lens = camera_lens
        self.camera_sensor = camera_sensor

        # initializing methods
        self.getSpectralLines()
    
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
    
    def __ComputeX(self, cmosx0: float, spectral_array: np.ndarray):
        
        max_order = ceil(self.echelle.compute_diffractionorder(spectral_array[0]))
        min_order = floor(self.echelle.compute_diffractionorder(spectral_array[-1]))

        diffraction_order_array = np.arange(start=min_order, stop=max_order+1, step=1)
        blaze_wavelength_array = self.echelle.compute_blazewavelength(diffraction_order=diffraction_order_array)
        FSR_array = self.echelle.compute_FSR(blazewavelength_array=blaze_wavelength_array)
        angular_dispersion_array = self.echelle.compute_angularDispE(wavelength_nm=blaze_wavelength_array)
        FSR_bins = blaze_wavelength_array + 0.5 * FSR_array

        index_FSR = np.digitize(x=spectral_array, bins=FSR_bins)

        cmosx = self.camera_lens.focal_length * angular_dispersion_array[np.maximum(index_FSR-1,0)]/1000 * (spectral_array - blaze_wavelength_array[np.maximum(index_FSR-1,0)]) + np.ones(np.shape(spectral_array)) * cmosx0 #mm (mrad.nm⁻1 to rad.nm⁻1 was done with /1000)
        
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
    
    def createDFcmos(self):
        # dataframe for camera sensor array
        df_cmos = pd.DataFrame({"X": [0, self.camera_sensor.size_x_mm, self.camera_sensor.size_x_mm, 0, 0], "Y": [0, 0, self.camera_sensor.size_y_mm, self.camera_sensor.size_y_mm, 0]})

        return df_cmos
    
    def computeCD(self,spectral_range_nm: tuple, spectral_res_nm: float, spatial_centers: list):
        spectral_array = np.arange(start=spectral_range_nm[0], stop=spectral_range_nm[1]+spectral_res_nm, step=spectral_res_nm)# in nm

        color_list = self.getColor(wavelength=spectral_array, gamma=0.8)
        demag_x, demag_y, r = self.computeMagnification()
        slit_width_cmos, slit_height_cmos = self.slit.width * demag_x, self.slit.height * demag_y
        # cmosx, diffraction_order_array, blaze_wavelength_array, dispersion_x, FSR_array = self.__ComputeX(cmosx0=spatial_centers[0,], spectral_array=spectral_array)
        
        # vectorize ComputeX and ComputeY for each spatial center
        # store everything in 1 dataframe to be able to export it as pkl, with 1 column spatial element
        
        df_mapping_list = []
        for x0,y0 in spatial_centers:
            cmosx,diffraction_order_array,blaze_wavelength_array, angular_dispersion_array,FSR_array = ComputeX(cmosx0=x0,spectral_array=spectral_array)
            cmosy, Ac = self.__ComputeY(cmosy0=y0,spectral_array=spectral_array)

            df_mapping = CreateDF_mapping(spatial_centers=(x0,y0), cmosx=cmosx, cmosy=cmosy, diffraction_order_array=diffraction_order_array, blaze_wavelength_array=blaze_wavelength_array, dispersion_x=angular_dispersion_array, dispersion_y=Ac, color_array=color_list, slit_width_cmos=slit_width_cmos,slit_height_cmos=slit_height_cmos,fsr_array=FSR_array)
            df_mapping_list.append(df_mapping)
        # df_mapping = pd.concat(df_mapping_list)

        
        
        # print(demag_x)
        
        # print(slit_width_cmos, slit_height_cmos)
        # print(color_list,len(color_list))
        df_cmos = self.createDFFcmos()

        pass
    
    def createDFmapping(self, spectral_array, spatial_centers, cmosx, cmosy, diffraction_order_array, blaze_wavelength_array, dispersion_x, dispersion_y,color_array,slit_width_cmos,slit_height_cmos,fsr_array):
        #compute spectral resolution and theoritical R
        spectral_res = np.ones(np.shape(dispersion_x))*((slit_width_cmos*1000) / self.camera_lens.focal_length) / dispersion_x
        
        #fitting shape of dict for centers
        spatial_centers_array = np.shape(cmosx)[0] * [spatial_centers]
        slit_width_cmos_array = np.shape(cmosx)[0] * [slit_width_cmos*1000] #mm to um
        slit_height_cmos_array = np.shape(cmosx)[0] * [slit_height_cmos*1000] #mm to um
       
        theoritical_R = spectral_array / spectral_res 
        #dataframes
        df_mapping = pd.DataFrame({"Center" : spatial_centers_array,"X": cmosx, "Y": cmosy, "wavelengths": spectral_array, "echelle_diffraction_order": diffraction_order_array, "echelle_blaze_wavelength[nm]": blaze_wavelength_array, "angular_dispersion_x[mrad/nm]": dispersion_x, "angular_dispersion_y[mrad/nm]": dispersion_y, "color":color_array, "Spectral resolution [nm]" : spectral_res, "Theoritical R" : theoritical_R,"Slit width [um]":slit_width_cmos_array, "Slit height [um]": slit_height_cmos_array, "Free Spectral Range [nm]":fsr_array})#, #add colors CIE in booktab

        
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
    
    def plotCD(self,df_cmos, df_mapping_list, slit_width_cmos, slit_height_cmos):

        #-----Plotting Parameters------#
        lines = list(self.__spectral_lines.keys())
        selection_spectral_windows = {
                                        k: pd.concat(
                                            [
                                                df_mapping_list[0].loc[
                                                    (df_mapping_list[0]["wavelengths"] <= v + 0.05) &
                                                    (df_mapping_list[0]["wavelengths"] >= v - 0.05) #0.2 nm window
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
                for i, df_mapping in enumerate(df_mapping_list, start=1):
                    hovertexts = [
                        f"Wavelength: {wl}<br>Echelle order: {order}<br>Resolution [A]: {10 * res}<br>Theory R: {R}<br>Slit [um](x, y): ({slit_x:.0f}, {slit_y:.0f})"
                        for wl, order, res, R, slit_x, slit_y in zip(
                            df_mapping["wavelengths"],
                            df_mapping["echelle_diffraction_order"],
                            df_mapping["Spectral resolution [nm]"],
                            df_mapping["Theoritical R"],
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
                    x=df_cmos["X"],
                    y=df_cmos["Y"],
                    mode="lines",
                    name=f"Emergent cam {cmosx_max}x{cmosy_max}mm",
                    line=dict(color="red", width=2)
                )

                # Layout details
                if disperser.__class__.__name__ == "Prism":
                    title = f"""
                    Camera focal length : {camera_lens.focal_length} mm
                    Echelle groove density : {echelle.groove_density:.1f} mm⁻¹<br>
                    Echelle blaze angle : {echelle.blaze_angle:.1f}°<br>
                    Deviation angle from Littrow : {2 * echelle.semi_deviation_angle_deg:.1f}°<br>
                    Cross-disperser Prism: {disperser.base} mm
                    """
                else:
                    title = f"""
                    Camera focal length : {camera_lens.focal_length} mm, 
                    Echelle groove density : {echelle.groove_density:.1f} mm⁻¹, 
                    Echelle blaze angle : {echelle.blaze_angle:.1f}°, 
                    Deviation angle from Littrow : {2 * echelle.semi_deviation_angle_deg:.1f}°, 
                    Cross-disperser groove density : {disperser.groove_density:.1f} mm⁻¹
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
                for i, df_mapping in enumerate(df_mapping_list, start=1):

                    hovertexts = [
                        f"Wavelength: {wl}<br>Echelle order: {order}<br>Resolution [A]: {10 * res}<br>Theory R: {R}<br>Slit [um](x, y): ({slit_x:.0f}, {slit_y:.0f})"
                        for wl, order, res, R, slit_x, slit_y in zip(
                            df_mapping["wavelengths"],
                            df_mapping["echelle_diffraction_order"],
                            df_mapping["Spectral resolution [nm]"],
                            df_mapping["Theoritical R"],
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
                        f"Wavelength: {wl}<br>Echelle order: {order}<br>Resolution [A]: {10 * res}<br>Theory R: {R}<br>Slit [um](x, y): ({slit_x:.0f}, {slit_y:.0f})"
                        for wl, order, res, R, slit_x, slit_y in zip(
                            df_sel["wavelengths"],
                            df_sel["echelle_diffraction_order"],
                            df_sel["Spectral resolution [nm]"],
                            df_sel["Theoritical R"],
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
                    x=df_cmos["X"],
                    y=df_cmos["Y"],
                    mode="lines",
                    name=f"Emergent cam {cmosx_max}x{cmosy_max}mm",
                    line=dict(color="red", width=2)
                )

                # Layout details
                if disperser.__class__.__name__ == "Prism":
                    title = f"""
                    Camera focal length : {camera_lens.focal_length} mm<br>
                    Echelle groove density : {echelle.groove_density:.1f} mm⁻¹<br>
                    Echelle blaze angle : {echelle.blaze_angle:.1f}°<br>
                    Deviation angle from Littrow : {2 * echelle.semi_deviation_angle_deg:.1f}°<br>
                    Cross-disperser Prism: {disperser.base} mm
                    """
                else:
                    title = f"""
                    Camera focal length : {camera_lens.focal_length} mm | Echelle groove density : {echelle.groove_density:.1f} mm⁻¹ | Echelle blaze angle : {echelle.blaze_angle:.1f}° | Deviation angle from Littrow : {2 * echelle.semi_deviation_angle_deg:.1f}° <br>Cross-disperser groove density : {disperser.groove_density:.1f}mm⁻¹"""

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
        
    


        

###-------------------- Functions --------------------###
def computeCD(spatial_centers: list, spectral_range: tuple, spectral_res: float, echelle: EchelleGrating, disperser: Grating | Prism, collimator_lens: Lens, camera_lens: Lens, cmosx_max: float, cmosy_max: float, alpha_deg: float, slit: Slit):
    
    spectral_array = np.arange(start=spectral_range[0], stop=spectral_range[1]+spectral_res, step=spectral_res)# in nm
    # n_spectral = len(spectral_array)

    # Neon I NIST spectral lines
    neon_lines = np.load("./NIST_Atomic-Specie/Neon.npy")
    neon_lines = np.squeeze(neon_lines)
    neon_lines = neon_lines.tolist()

    # Thorium I NIST spectral lines
    thorium_lines = np.load("./NIST_Atomic-Specie/Thorium0.95.npy")
    thorium_lines = np.squeeze(thorium_lines)
    thorium_lines = thorium_lines.tolist()

    # Argon I NIST spectral lines
    argon_lines = np.load("./NIST_Atomic-Specie/Argon0.95.npy")
    argon_lines = np.squeeze(argon_lines)
    argon_lines = argon_lines.tolist()


    spectral_lines = {"Mg I b3 5167": [516.7] ,"Mg I b2 5172": [517.2] ,"Mg I b1 5183": [518.3] ,"Fe I 5250": [525.0],"Mn I 5399": [539.9], "He I D3 5876":[587.6], "Na I D2 5890": [589.0] ,"Na I D1 5896": [589.6],"Fe I 6173 (HMI)": [617.3], "Fe I 6301-6302": [630.15],"Ca I 6439": [643.9], "Ha": [656.3], "Ni I 6643": [664.3],"Fe I/Ca I 6718": [671.8],"K I 7699": [769.9],"Ca II 8498": [849.8], "Ca II 8542": [854.2],"Ca II 8662": [866.2], "Neon": neon_lines, "Thorium": thorium_lines, "Argon": argon_lines} #to be continued, rn only wl > 500 nm

    
    def ComputeX(cmosx0, spectral_array):
        #center of cmos on x axis
        # cmosx0 = 0

        max_order = ceil(echelle.compute_diffractionorder(spectral_array[0]))
        min_order = floor(echelle.compute_diffractionorder(spectral_array[-1]))

        diffraction_order_array = np.arange(start=min_order, stop=max_order+1, step=1)
        blaze_wavelength_array = echelle.compute_blazewavelength(diffraction_order=diffraction_order_array)
        FSR_array = echelle.compute_FSR(blazewavelength_array=blaze_wavelength_array)
        angular_dispersion_array = echelle.compute_angularDispE(wavelength_nm=blaze_wavelength_array)
        FSR_bins = blaze_wavelength_array + 0.5 * FSR_array

        index_FSR = np.digitize(x=spectral_array, bins=FSR_bins)

        cmosx = camera_lens.focal_length * angular_dispersion_array[np.maximum(index_FSR-1,0)]/1000 * (spectral_array - blaze_wavelength_array[np.maximum(index_FSR-1,0)]) + np.ones(np.shape(spectral_array)) * cmosx0 #mm (mrad.nm⁻1 to rad.nm⁻1 was done with /1000)
        return cmosx,diffraction_order_array[np.maximum(index_FSR-1,0)],blaze_wavelength_array[np.maximum(index_FSR-1,0)], angular_dispersion_array[np.maximum(index_FSR-1,0)], FSR_array[np.maximum(index_FSR-1,0)]

    def ComputeY(cmosy0, spectral_array):
        
        #definitions and initial position (central wavelength in y-center of cmos)
        lambda0 = (spectral_array[-1] + spectral_array[0]) / 2 
        # cmosy0 = 0
        dspectral_array = spectral_array - np.ones(np.shape(spectral_array)) * lambda0
        
        #optical system computations
        if disperser.__class__.__name__ == "Prism":
            Ac = disperser.compute_angularDisp(wavelength_nm=spectral_array) /1000 #mrad.nm⁻1 to rad.nm⁻1
            
        else:
            Ac = disperser.compute_angularDisp(wavelength_nm=dspectral_array) /1000 #mrad.nm⁻1 to rad.nm⁻1

        #computation y coordinate
        cmosy = np.ones(np.shape(dspectral_array))* cmosy0 + camera_lens.focal_length * Ac * dspectral_array
             
        return cmosy, Ac

    def ComputeMagnification(collimator_lens: Lens, camera_lens: Lens, echelle: EchelleGrating):#, disperser: Grating, alpha_deg: float): r for crossdisperser ?
        r = (1 - np.tan(np.deg2rad(echelle.blaze_angle))*np.tan(np.deg2rad(echelle.semi_deviation_angle_deg)))/(1 + np.tan(np.deg2rad(echelle.blaze_angle))*np.tan(np.deg2rad(echelle.semi_deviation_angle_deg)))
        demag_y = camera_lens.focal_length / collimator_lens.focal_length
        demag_x = r * camera_lens.focal_length / collimator_lens.focal_length

        return demag_x, demag_y, r
    
    def CreateDF_mapping(spatial_centers, cmosx, cmosy, diffraction_order_array, blaze_wavelength_array, dispersion_x, dispersion_y,color_array,slit_width_cmos,slit_height_cmos,fsr_array):
        #compute spectral resolution and theoritical R
        spectral_res = np.ones(np.shape(dispersion_x))*((slit_width_cmos*1000) / camera_lens.focal_length) / dispersion_x
        
        #fitting shape of dict for centers
        spatial_centers_array = np.shape(cmosx)[0] * [spatial_centers]
        slit_width_cmos_array = np.shape(cmosx)[0] * [slit_width_cmos*1000] #mm to um
        slit_height_cmos_array = np.shape(cmosx)[0] * [slit_height_cmos*1000] #mm to um
       
        theoritical_R = spectral_array / spectral_res 
        #dataframes
        df_mapping = pd.DataFrame({"Center" : spatial_centers_array,"X": cmosx, "Y": cmosy, "wavelengths": spectral_array, "echelle_diffraction_order": diffraction_order_array, "echelle_blaze_wavelength[nm]": blaze_wavelength_array, "angular_dispersion_x[mrad/nm]": dispersion_x, "angular_dispersion_y[mrad/nm]": dispersion_y, "color":color_array, "Spectral resolution [nm]" : spectral_res, "Theoritical R" : theoritical_R,"Slit width [um]":slit_width_cmos_array, "Slit height [um]": slit_height_cmos_array, "Free Spectral Range [nm]":fsr_array})#, #add colors CIE in booktab

        
        return df_mapping
    
    def CreateDF_cmos(cmosx_max, cmosy_max):
       
        #dataframe
        df_cmos = pd.DataFrame({"X": [0, cmosx_max, cmosx_max, 0, 0], "Y": [0, 0, cmosy_max, cmosy_max, 0]})

        
        
        return df_cmos

    def GetColor(wavelength, gamma=0.8):
        #
        #    Based on code by Dan Bruton
        #    http://www.physics.sfasu.edu/astro/color/spectra.html
        #    '''
        # Adapted from <script src="https://gist.github.com/friendly/67a7df339aa999e2bcfcfec88311abfc.js"></script>
        color_list = []
        for i in range(np.shape(wavelength)[0]):
            wl_i = wavelength[i]
            if (wl_i < 380 or wl_i > 750):
                R = 0.0
                G = 0.0
                B = 0.0
            elif (wl_i >= 380 and wl_i <= 440):
                attenuation = 0.3 + 0.7 * (wl_i - 380) / (440 - 380)
                R = ((-(wl_i - 440) / (440 - 380)) * attenuation) ** gamma
                G = 0.0
                B = (1.0 * attenuation) ** gamma
                
            elif (wl_i >= 440 and wl_i <= 490):
                R = 0.0
                G = ((wl_i - 440) / (490 - 440)) ** gamma
                B = 1.0
                
            elif (wl_i >= 490 and wl_i <= 510) :
                R = 0.0
                G = 1.0
                B = (-(wl_i - 510) / (510 - 490)) ** gamma
                
            elif (wl_i >= 510 and wl_i <= 580):
                R = ((wl_i - 510) / (580 - 510)) ** gamma
                G = 1.0
                B = 0.0
            
            elif (wl_i >= 580 and wl_i <= 645):  
                R = 1.0
                G = (-(wl_i - 645) / (645 - 580)) ** gamma
                B = 0.0
            
            elif (wl_i >= 645 and wl_i <= 750):
                attenuation = 0.3 + 0.7 * (750 - wl_i) / (750 - 645)
                R = (1.0 * attenuation) ** gamma
                G = 0.0
                B = 0.0
        
            else:
                R = 0.0
                G = 0.0
                B = 0.0
                
            R = round(R * 255)
            G = round(G * 255)
            B = round(B * 255)
            color_list.append("#%02x%02x%02x" % (R,G,B)) 
        return color_list

    def DrawGrid(df_cmos, df_mapping_list, slit_width_cmos, slit_height_cmos):

        #-----Plotting Parameters------#
        lines = list(spectral_lines.keys())
        selection_spectral_windows = {
                                        k: pd.concat(
                                            [
                                                df_mapping_list[0].loc[
                                                    (df_mapping_list[0]["wavelengths"] <= v + 0.05) &
                                                    (df_mapping_list[0]["wavelengths"] >= v - 0.05) #0.2 nm window
                                                ]
                                                for v in v_list
                                            ]
                                        )
                                        for k, v_list in spectral_lines.items()
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
                return list(spectral_lines.keys())
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
                for i, df_mapping in enumerate(df_mapping_list, start=1):
                    hovertexts = [
                        f"Wavelength: {wl}<br>Echelle order: {order}<br>Resolution [A]: {10 * res}<br>Theory R: {R}<br>Slit [um](x, y): ({slit_x:.0f}, {slit_y:.0f})"
                        for wl, order, res, R, slit_x, slit_y in zip(
                            df_mapping["wavelengths"],
                            df_mapping["echelle_diffraction_order"],
                            df_mapping["Spectral resolution [nm]"],
                            df_mapping["Theoritical R"],
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
                    x=df_cmos["X"],
                    y=df_cmos["Y"],
                    mode="lines",
                    name=f"Emergent cam {cmosx_max}x{cmosy_max}mm",
                    line=dict(color="red", width=2)
                )

                # Layout details
                if disperser.__class__.__name__ == "Prism":
                    title = f"""
                    Camera focal length : {camera_lens.focal_length} mm
                    Echelle groove density : {echelle.groove_density:.1f} mm⁻¹<br>
                    Echelle blaze angle : {echelle.blaze_angle:.1f}°<br>
                    Deviation angle from Littrow : {2 * echelle.semi_deviation_angle_deg:.1f}°<br>
                    Cross-disperser Prism: {disperser.base} mm
                    """
                else:
                    title = f"""
                    Camera focal length : {camera_lens.focal_length} mm, 
                    Echelle groove density : {echelle.groove_density:.1f} mm⁻¹, 
                    Echelle blaze angle : {echelle.blaze_angle:.1f}°, 
                    Deviation angle from Littrow : {2 * echelle.semi_deviation_angle_deg:.1f}°, 
                    Cross-disperser groove density : {disperser.groove_density:.1f} mm⁻¹
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
                for i, df_mapping in enumerate(df_mapping_list, start=1):

                    hovertexts = [
                        f"Wavelength: {wl}<br>Echelle order: {order}<br>Resolution [A]: {10 * res}<br>Theory R: {R}<br>Slit [um](x, y): ({slit_x:.0f}, {slit_y:.0f})"
                        for wl, order, res, R, slit_x, slit_y in zip(
                            df_mapping["wavelengths"],
                            df_mapping["echelle_diffraction_order"],
                            df_mapping["Spectral resolution [nm]"],
                            df_mapping["Theoritical R"],
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
                        f"Wavelength: {wl}<br>Echelle order: {order}<br>Resolution [A]: {10 * res}<br>Theory R: {R}<br>Slit [um](x, y): ({slit_x:.0f}, {slit_y:.0f})"
                        for wl, order, res, R, slit_x, slit_y in zip(
                            df_sel["wavelengths"],
                            df_sel["echelle_diffraction_order"],
                            df_sel["Spectral resolution [nm]"],
                            df_sel["Theoritical R"],
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
                    x=df_cmos["X"],
                    y=df_cmos["Y"],
                    mode="lines",
                    name=f"Emergent cam {cmosx_max}x{cmosy_max}mm",
                    line=dict(color="red", width=2)
                )

                # Layout details
                if disperser.__class__.__name__ == "Prism":
                    title = f"""
                    Camera focal length : {camera_lens.focal_length} mm<br>
                    Echelle groove density : {echelle.groove_density:.1f} mm⁻¹<br>
                    Echelle blaze angle : {echelle.blaze_angle:.1f}°<br>
                    Deviation angle from Littrow : {2 * echelle.semi_deviation_angle_deg:.1f}°<br>
                    Cross-disperser Prism: {disperser.base} mm
                    """
                else:
                    title = f"""
                    Camera focal length : {camera_lens.focal_length} mm | Echelle groove density : {echelle.groove_density:.1f} mm⁻¹ | Echelle blaze angle : {echelle.blaze_angle:.1f}° | Deviation angle from Littrow : {2 * echelle.semi_deviation_angle_deg:.1f}° <br>Cross-disperser groove density : {disperser.groove_density:.1f}mm⁻¹"""

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
        
        

    #computations
    color_list = GetColor(wavelength=spectral_array)
    demag_x, demag_y, r = ComputeMagnification(collimator_lens=collimator_lens, camera_lens=camera_lens, echelle=echelle)
    slit_width_cmos = slit.width * demag_x
    slit_height_cmos = slit.height * demag_y
    df_mapping_list = []
    for x0,y0 in spatial_centers:
        cmosx,diffraction_order_array,blaze_wavelength_array, angular_dispersion_array,FSR_array = ComputeX(cmosx0=x0,spectral_array=spectral_array)
        cmosy, Ac = ComputeY(cmosy0=y0,spectral_array=spectral_array)

        df_mapping = CreateDF_mapping(spatial_centers=(x0,y0), cmosx=cmosx, cmosy=cmosy, diffraction_order_array=diffraction_order_array, blaze_wavelength_array=blaze_wavelength_array, dispersion_x=angular_dispersion_array, dispersion_y=Ac, color_array=color_list, slit_width_cmos=slit_width_cmos,slit_height_cmos=slit_height_cmos,fsr_array=FSR_array)
        df_mapping_list.append(df_mapping)
    # df_mapping = pd.concat(df_mapping_list)

    
    
    # print(demag_x)
    
    # print(slit_width_cmos, slit_height_cmos)
    # print(color_list,len(color_list))
    df_cmos = CreateDF_cmos(cmosx_max=cmosx_max, cmosy_max=cmosy_max)
    DrawGrid(df_cmos=df_cmos, df_mapping_list=df_mapping_list,slit_height_cmos=slit_height_cmos, slit_width_cmos=slit_width_cmos)

    return df_mapping_list