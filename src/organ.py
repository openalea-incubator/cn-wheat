# -*- coding: latin-1 -*-

'''
Created on 11 août 2014

@author: cchambon
'''

from __future__ import division # use "//" to do integer division


class Organ(object):
    
    Mstruct_axis = 2.08 # Structural mass (g) of a plant (Bertheloot, 2011)
    alpha_axis = 1
    
    def __init__(self, name):
        self.name = name


class Lamina(Organ):
    
    alpha_lamina = 1 # Proportion of leaf structral mass containing substrate
    delta_Dstorage = 0.0001
    K_storage = 20
    K_sucrose = 0.66
    Vmax_storage = 2
    Vmax_sucrose = 1
    sigma = 1.85e-07 # Conductivity
    
    def __init__(self, lamina_area, Mstruct_lamina, Assimilation, STORAGE_0, 
                 SUCROSE_lamina_0, TRIOSESP_0, name=''):
        super(Lamina, self).__init__(name)
        # parameters
        self.lamina_area = lamina_area # Leaf area (m-2) of a flag leaf (Bertheloot, 2011)
        self.Mstruct_lamina = Mstruct_lamina # Structural mass (g) of a flag leaf (Bertheloot, 2011)
        self.Assimilation = Assimilation
        
        # linear interpolation of Assimilation      
        self.An_linear_interpolation = None
        
        # initialization
        # flow to phloem
        self.Loading_sucrose = 0
        # compartments
        self.STORAGE_0 = STORAGE_0
        self.SUCROSE_lamina_0 = SUCROSE_lamina_0
        self.TRIOSESP_0 = TRIOSESP_0
        
    def calculate_Photosynthesis(self, t):
        return self.An_linear_interpolation(t) * self.lamina_area
    
    def calculate_D_storage(self, STORAGE):
        '''Flow from STORAGE to SUCROSE_lamina
        '''
        return max(0, Lamina.delta_Dstorage * (STORAGE/(self.Mstruct_lamina*Lamina.alpha_lamina)))  

    def calculate_S_storage(self, TRIOSESP):
        '''Flow from TRIOSESP to STORAGE
        '''
        return ((max(0, TRIOSESP)/(self.Mstruct_lamina*Lamina.alpha_lamina)) * Lamina.Vmax_storage) / ((max(0, TRIOSESP)/(self.Mstruct_lamina*Lamina.alpha_lamina)) +Lamina.K_storage)
        
    def calculate_S_sucrose(self, TRIOSESP):
        '''Flow from TRIOSESP to SUCROSE_lamina
        '''
        return ((max(0,TRIOSESP)/(self.Mstruct_lamina*Lamina.alpha_lamina)) * Lamina.Vmax_sucrose) / ((max(0, TRIOSESP)/(self.Mstruct_lamina*Lamina.alpha_lamina)) +Lamina.K_sucrose)
    
    def calculate_STORAGE_derivative(self, S_storage, D_storage):
        '''µmol of C'''
        return (S_storage - D_storage) * (self.Mstruct_lamina*Lamina.alpha_lamina)
    
    def calculate_SUCROSE_lamina_derivative(self, S_sucrose, D_storage, Loading_sucrose):
        '''µmol of C'''
        return (S_sucrose + D_storage - Loading_sucrose) * (self.Mstruct_lamina*Lamina.alpha_lamina)
    
    def calculate_TRIOSESP_derivative(self, Photosynthesis, S_sucrose, S_storage):
        '''µmol of C''' 
        return Photosynthesis - (S_sucrose + S_storage) * (self.Mstruct_lamina*Lamina.alpha_lamina)
    
    def calculate_Conc_TriosesP(self, TRIOSESP):
        return (TRIOSESP/self.Mstruct_lamina)/3
    
    def calculate_Conc_Storage(self, STORAGE):
        return (STORAGE/self.Mstruct_lamina)/6
    
    def calculate_Conc_Sucrose_lamina(self, SUCROSE_lamina):
        return (SUCROSE_lamina/self.Mstruct_lamina)/12
    
    def calculate_Loading_sucrose(self, SUCROSE_lamina, SUCROSE_phloem):
        return (max(0, SUCROSE_lamina)/(self.Mstruct_lamina*Lamina.alpha_lamina)) * ((max(0, SUCROSE_lamina)/(self.Mstruct_lamina*Lamina.alpha_lamina)) - (max(0, SUCROSE_phloem)/(Organ.Mstruct_axis*Organ.alpha_axis))) * (Lamina.sigma * self.Mstruct_lamina**(2/3))
    
    def get_initial_conditions(self):
        return [self.STORAGE_0, self.SUCROSE_lamina_0, self.TRIOSESP_0] # keep this order !
    
    
class Ear(Organ):
    
    Vmax_ear = 0.026
    K_ear = 2
    
    def __init__(self, Ear_value_0, name=''):
        super(Ear, self).__init__(name)
        # initialization
        # flow from phloem
        self.Unloading = 0
        # compartments
        self.Ear_value_0 = Ear_value_0
    
    def calculate_Ear_value_derivative(self, Unloading):
        return Unloading
    
    def calculate_Dry_mass_ear(self, Ear_value):
        '''g of C'''
        return Ear_value/(12000000)
    
    def calculate_Unloading(self, SUCROSE_phloem):
        return max(0, (Ear.Vmax_ear*(SUCROSE_phloem/(Organ.Mstruct_axis*Organ.alpha_axis)))/(Ear.K_ear+(SUCROSE_phloem/(Organ.Mstruct_axis*Organ.alpha_axis))))
    
    def get_initial_conditions(self):
        return [self.Ear_value_0]
    
    
class Phloem(Organ):
    
    def __init__(self, SUCROSE_phloem_0, name=''):
        super(Phloem, self).__init__(name)
        # initialization
        # variables
        self.Conc_Sucrose_phloem = 0
        # compartments
        self.SUCROSE_phloem_0 = SUCROSE_phloem_0
        
    def calculate_SUCROSE_phloem_derivative(self, organs):
        '''µmol of C'''
        SUCROSE_phloem_derivative = 0
        for organ_ in organs:
            if isinstance(organ_, Lamina):
                SUCROSE_phloem_derivative += organ_.Loading_sucrose*organ_.Mstruct_lamina*Lamina.alpha_lamina
            elif isinstance(organ_, Ear):
                SUCROSE_phloem_derivative -= organ_.Unloading*Organ.Mstruct_axis*Organ.alpha_axis
        return SUCROSE_phloem_derivative
    
    def calculate_Conc_Sucrose_phloem(self, SUCROSE_phloem):
        return (SUCROSE_phloem/Organ.Mstruct_axis)/12

    def get_initial_conditions(self):
        return [self.SUCROSE_phloem_0]

    