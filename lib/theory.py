import numpy as np
from scipy.special import hyp2f1
from math import factorial
from lmfit import Minimizer, Parameters

class Integration():
    '''Integrates the electric field according to procedure from Sheik-Bahae, given the initial parameters'''
    def __init__(self, window, beta, n2, DPhi0, positions, d0, aperture_radius, wavelength, beamwaist, n_components, integration_steps, stype="CA"):
        super(Integration, self).__init__()
        self.window = window
        
        # Data range
        self.z = positions # evenly spaced

        # Apertures and distances
        self.d0 = d0 # distance from z=0 to the aperture plane
        self.ra = aperture_radius # in meters

        # Beam properties
        self.lda = wavelength # in meters
        self.w0 = beamwaist # in meters

        # Sample properties
        self.n2 = n2 # non-linear refractive index
        try:
            self.T = beta*self.lda/self.n2 # non-linear transmittance
        except ZeroDivisionError:
            self.T = 0
        self.DPhi0 = DPhi0 # on-axis phase shift in z=0
        
        # Integration parameters
        self.mm = n_components # number of electric field components (for Gaussian decomposition)
        self.ir = integration_steps
        
        self.stype = stype
        self.derive(self.DPhi0, self.w0, self.d0, self.ra, stype)
    
    def derive(self, DPhi0, w0, d0, ra, stype):
        '''1st method called. Called on __init__'''
        # Beam properties
        self.k = 2*np.pi/self.lda # wave vector in free space
        self.z0 = 0.5*self.k*w0**2 # diffraction length of the beam
        self.wa = w0*np.sqrt(1+d0**2/self.z0**2) # beam radius at the aperture plane

        # Aperture radius
        self.ra = ra

        # Sample properties
        self.Dphi0 = DPhi0/(1+self.z**2/self.z0**2)

        # Additional derived parameters
        self.wz = w0*np.sqrt(1+self.z**2/self.z0**2)

        self.Rz = self.z+self.z0**2/self.z
        self.d = d0-self.z
        self.g = 1+self.d/self.Rz

        match stype:
            case "CA":
                self.bigproduct()
            case "OA":
                self.calculate_Tz_for_OA(self.window.solventOA_absorptionModel_comboBox.currentText()) # CZYŻBY????????
        
    def calculate_Tz_for_OA(self, model):
        '''2nd method called. Called by derive() for stype = "OA".'''
        match self.window.fittingTabs.currentIndex():
            case 0:
                ftype = "Silica"
            case 1:
                ftype = "Solvent"
            case 2:
                ftype = "Sample"
        
        match ftype:
            case "Silica":
                return # Do not fit OA for silica
            case "Solvent":
                absorption_checkbox = self.window.solventOA_isAbsorption_checkBox.isChecked()
            case "Sample":
                absorption_checkbox = self.window.sampleOA_isAbsorption_checkBox.isChecked()
        
        if absorption_checkbox is False:
            self.Tz = 0
            # Psi_n = (n * beta_n * I0**n * L_eff)**(1/n) (n+1)-photon absorption
            # Psi1 = T*DPhi0/2/pi
            #Psi1 = self.T*self.DPhi0/(2*np.pi)
            if self.T != 0:
                #Psi1 = -(lambertw(-np.exp(-self.T)*self.T)-self.T)/self.T # Lambert W function calculates inverse of hyp2f1(1,1,2,-psi1) at self.z=0
                Psi1 = self.T
                #Psi2 = (2*window.laserI0*self.T*self.DPhi0/(2*np.pi))**(1/2)
            else:
                Psi1 = 0
            
            Psi2 = self.T*8 # 8 is a value out of the hat to obtain similar extreme point (valley or peak)
            
            match model:
                case "2PA":
                    psi1 = Psi1/(1+self.z**2/self.z0**2)
                    # self.Tz = hyp2f1(1/n,1/n,(n+1)/n,-psi**n)
                    self.Tz = hyp2f1(1,1,2,-psi1)
                    pass
                case "3PA":
                    psi2 = Psi2/(1+self.z**2/self.z0**2)
                    # self.Tz = hyp2f1(1/n,1/n,(n+1)/n,-psi**n)
                    self.Tz = hyp2f1(1/2,1/2,3/2,-(psi2)**2)
                
                case "2PA+3PA":
                    if Psi1 !=0:
                        psi1 = Psi1/(1+self.z**2/self.z0**2)
                        psi2 = Psi2/(1+self.z**2/self.z0**2)
                        f_x_psi1_psi2 = 1+psi1*(0.339*np.sin(0.498*psi2)-0.029)/(1+0.966*psi1*psi2**-0.718) # coupling function formula from OPTICS EXPRESS Vol. 13, No. 23 9231
                        self.Tz = hyp2f1(1,1,2,-psi1)*hyp2f1(1/2,1/2,3/2,-(psi2)**2)*f_x_psi1_psi2
                    
                    else:
                        self.Tz = np.ones_like(self.z)
                
                case "RSA":
                    self.Tz = np.ones_like(self.z)
                    
                case "SA":
                    self.Tz = np.ones_like(self.z)
                    
                case "2PA+SA":
                    self.Tz = np.ones_like(self.z)
                    
                case _:
                    self.Tz = np.ones_like(self.z)
            
            self.Tznorm = 2*self.Tz/(np.average(self.Tz[0:10])+np.average(self.Tz[len(self.Tz)-10:]))
        
        else:
            self.Tznorm = np.ones_like(self.z)

        return self.Tznorm

    def calculate_fm(self):
        '''3rd method called. Called by bigproduct()'''
        self.result = [(1j*self.Dphi0)**m/factorial(m)*self.product[m] for m in range(0,self.mm)]
        return self.result

    def bigproduct(self):
        '''2nd method called. Called by derive() for stype = "CA".'''
        # Big product operator (j from 1 to m)
        # if self.beta == 0:
        # self.product = np.ones(self.mm) # The result is already known. Uses numpy array for calculate_fm function to work properly

        #else: # THIS MUST BE CALCULATED TO TAKE INTO ACCOUNT TRANSMITTANCE THROUGH THE APERTURE
        self.product = []
        for m in range(0,self.mm):
            if m==0:
                self.product.append(1)
                continue # 0th value of m gives product=1 and in cumulative product it gives no contribution, so go to next iteration straight away.
            else:
                self.product.append(np.cumprod([1+1j*(j-1/2)/2/np.pi*self.T for j in range(1,m+1)])[-1]) # calculate cumulative product for given m
        #                                                                                                                          # and make 'product' contain all m products
        self.fm = self.calculate_fm()

        Tzo = self.open() # Returns final result for OA
        Tzc = self.closed() # Returns final result for CA
        
        self.Tznorm = Tzc/Tzo # This way, transmittance through aperture is taken into account
    
    def open(self):
        '''4th method called. Called by bigproduct()'''
        self.dr = 3*self.wa/self.ir # integration over radius step size
        self.open_sum = self.bigsum()
        return self.open_sum

    def closed(self):
        '''6th method called. Called by bigproduct()'''
        self.dr = self.ra/self.ir # integration over radius step size
        self.closed_sum = self.bigsum()
        return self.closed_sum
    
    def bigsum(self):
        '''5th and 7th method called. Called by open() and closed()'''
        # Big sum operator (m from 0 to "infinity")
        # integration over radius
        
        # Transmitted power
        self.Tz = 0

        for rr in range(self.ir):
            self.E = 0
            
            for m in range(0,self.mm):
                self.wm0 = self.wz/np.sqrt((2*m+1))
                self.dm = 0.5*self.k*self.wm0**2
                self.wm = self.wm0*np.sqrt(self.g**2+self.d**2/self.dm**2)
                self.tm = np.arctan(self.g/(self.d/self.dm))
                self.Rm = self.d/(1-self.g/(self.g**2+self.d**2/self.dm**2))

                self.E += self.fm[m]*np.exp(1j*self.tm)*self.wm0/self.wm/self.wz*np.exp((-1/self.wm**2+1j*np.pi/self.lda/self.Rm)*(rr*self.dr)**2)
        
            self.Tz += np.abs(self.E)**2*rr*self.dr # transmittance through aperture plane
        
        # T(z)=P_T/(S*P_i), where S=1-exp(-2*ra^2/wa^2)
        #S = 1-np.exp(-2*self.ra**2/self.wa**2)
        #I0 - instantaneous laser fluence - this is missing (we want to get it)
        #self.Tznorm = 3E8*8.854E-12*self.Tz/(S*self.w0**2/2*I0)
        
        # Normalize transmitted power before dividing CA/OA (the only way of normalization when we don't have I0)
        self.Tznorm = 2*self.Tz/(np.average(self.Tz[0:10])+np.average(self.Tz[len(self.Tz)-10:]))

        return self.Tznorm

class Fitting():
    def __init__(self, window, sample_type: Integration, amplitude, beamwaist, zero_level, centerpoint, nop, data):
        super(Fitting, self).__init__()
        self.window = window
        self.sample_type = sample_type
        self.amplitude = amplitude
        self.beamwaist = beamwaist
        self.zero_level = zero_level
        self.centerpoint = centerpoint
        self.nop = nop
        self.ydata = data

    # The parametrized function to be plotted. It is also initial guess for automatic fitting.
    def manual(self, zero_level, centerpoint, amplitude, beamwaist, z_range, stype="CA") -> list:
        '''Returns integrated field for given input parameters
        
        ONLY FOR n2 FOR NOW!!!!!!!!!!!!!'''
        self.z_range = z_range # in meters
        self.sample_type.z = np.array([self.z_range*(zz - centerpoint)/self.nop-self.z_range/2 for zz in range(self.nop)]) # in meters
        self.window.get_general_parameters()
        if stype == "CA":
            self.sample_type.derive(amplitude,beamwaist,self.window.d0,self.window.ra,stype)
            cas = self.sample_type.closed_sum # Tznorm
            oas = self.sample_type.open_sum   # Tznorm
            result = (cas/oas)+(zero_level-1)
        elif stype == "OA":
            self.sample_type.derive(amplitude,beamwaist,self.window.d0,self.window.ra,stype)
            oas = self.sample_type.Tznorm
            result = oas+(zero_level-1)
        
        return result

    # The function to be minimized in automated fitting
    def fcn2min(self,params,*weights):
        pars = params.valuesdict().values()
        ynew = self.manual(*pars)
        return [w*((yn-yi)**2) for w,yn,yi in zip(weights,ynew,self.ydata)] # SSE

    # The actual processor for automatic fitting
    def automatic(self, z_range, ftype: str, stype: str, line_xydata):
        self.z_range = z_range

        if (ftype == "Solvent" and self.window.solventCA_customBeamwaist_checkBox.isChecked() is False):#\
            #or (ftype == "Sample" and window.sampleCA_customBeamwaist_checkBox.isChecked() is False):
            vary_beamwaist = False
        else:
            vary_beamwaist = True
        
        vary_centerpoint = True
        if ftype == "Solvent" and stype == "OA":
            if self.window.solventOA_customCenterPoint_checkBox.isChecked() is False:
                vary_centerpoint = False
            else:
                vary_centerpoint = True
        
        # Apply initial values
        self.params = Parameters()
        self.params.add('Zero', value=self.zero_level, min=0.75, max=1.25)
        self.params.add('Center', value=self.centerpoint, min=-50, max=50, vary=vary_centerpoint) # in number of datapoints
        
        if stype == "CA":
            self.params.add('DPhi0', value=self.amplitude, min=-2, max=2)
            self.params.add('Beamwaist', value=self.beamwaist, min=15E-6, max=150E-6, vary=vary_beamwaist)
            self.params.add('Zrange', value=self.z_range, vary=False)
        
        elif stype == "OA":
            self.params.add('T', value=self.amplitude, min=-2, max=2)
            self.params.add('Beamwaist', value=self.beamwaist, vary=False)
            self.params.add('Zrange', value=self.z_range, vary=False)
        
        if ftype == "Silica":
            #line = window.silica_figure.axes.get_lines()[0]
            #xs, ys = line.get_data()
            xs, ys = line_xydata
            weights = np.ones(np.shape(xs))
            if hasattr(self.window, 'silica_cursorPositions'):
                if len(self.window.silica_cursorPositions) == 2:
                    # CA ranges for weighting the fit
                    x1 = self.window.silica_cursorPositions[0][0]
                    x1_index = list(xs).index(min(xs, key=lambda x:abs(x-x1)))
                    # y1 = window.silica_cursorPositions[0][1]
                    #y1_index = list(ys).index(min(ys, key=lambda x:abs(x-y1)))
                    x2 = self.window.silica_cursorPositions[1][0]
                    x2_index = list(xs).index(min(xs, key=lambda x:abs(x-x2)))
                    # y2 = window.silica_cursorPositions[1][1]
                    #y2_index = list(ys).index(min(ys, key=lambda x:abs(x-y2)))
                    
                    x_sm, x_lg = sorted([x1_index, x2_index])
                    weights = [0 if (xi < x_sm or xi > x_lg) else 1 for xi in range(len(xs))]
            
            fitter = Minimizer(self.fcn2min,self.params,fcn_args=(weights))

        elif ftype == "Solvent":
            #line = window.solventCA_figure.axes.get_lines()[0]
            #xs, ys = line.get_data()
            xs, ys = line_xydata
            weights = np.ones(np.shape(xs))
            if stype == "CA":
                attribute_name = 'solventCA_cursorPositions'
            elif stype == "OA":
                attribute_name = 'solventOA_cursorPositions'
            
            if hasattr(self.window, attribute_name):
                if stype == "CA":
                    attribute = self.window.solventCA_cursorPositions
                elif stype == "OA":
                    attribute = self.window.solventOA_cursorPositions

                if len(attribute) == 2:
                    # CA ranges for weighting the fit
                    x1 = attribute[0][0]
                    x1_index = list(xs).index(min(xs, key=lambda x:abs(x-x1)))
                    # y1 = attribute[0][1]
                    # y1_index = list(ys).index(min(ys, key=lambda x:abs(x-y1)))
                    x2 = attribute[1][0]
                    x2_index = list(xs).index(min(xs, key=lambda x:abs(x-x2)))
                    # y2 = attribute[1][1]
                    # y2_index = list(ys).index(min(ys, key=lambda x:abs(x-y2)))
                    
                    x_sm, x_lg = sorted([x1_index, x2_index])
                    weights = [0 if (xi < x_sm or xi > x_lg) else 1 for xi in range(len(xs))]
            
            fitter = Minimizer(self.fcn2min,self.params,fcn_args=(weights))
        
        result = fitter.minimize(method='least_squares', max_nfev=1000)

        res_pars = result.params.valuesdict().values()
        result_line = self.manual(*res_pars, stype) # this will have to catch two lines (n2 and OA)

        result.params.pretty_print()

        return result, result_line
