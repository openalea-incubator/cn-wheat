# -*- coding: latin-1 -*-

'''
Created on 11 août 2014

@author: cchambon
'''

from __future__ import division # use "//" to do integer division


class Organ(object):
    
    Mstruct_axis = 2.08 # Structural mass (g) of a plant (Bertheloot, 2011)
    alpha_axis = 1
    delta_Dstorage = 0.0001
    delta_t = 3600
    sigma = 1.85e-07 # Conductivity
    Vmax_storage = 2
    K_storage = 20
    Vmax_sucrose = 1
    K_sucrose = 0.66
    Vmax_Dfructan = 0.009
    K_Dfructan = 2000
    K_Sfructan = 20000
    Vmax_Sfructan = 0.06
    alpha = 1 # Proportion of leaf structural mass containing substrate # TODO: is it really common to all organs? 
    
    def __init__(self, name):
        self.name = name
        
    # VARIABLES
        
    def calculate_Photosynthesis(self, t):
        return self.An_linear_interpolation(t)*self.Area*Organ.delta_t
        
    def calculate_Conc_Storage(self, STORAGE):
        return (STORAGE/self.Mstruct)/6
    
    def calculate_Conc_Sucrose(self, SUCROSE):
        return (SUCROSE/self.Mstruct)/12
    
    def calculate_Conc_TriosesP(self, TRIOSESP):
        return (TRIOSESP/self.Mstruct)/3
    
    # FLOWS
    
    def calculate_D_storage(self, STORAGE):
        '''Flow from STORAGE to SUCROSE
        '''
        return max(0, Organ.delta_Dstorage * (STORAGE/(self.Mstruct*Organ.alpha)))*Organ.delta_t
    
    def calculate_S_storage(self, TRIOSESP):
        '''Flow from TRIOSESP to STORAGE
        '''
        return ((max(0, TRIOSESP)/(self.Mstruct*Organ.alpha)) * Organ.Vmax_storage) / ((max(0, TRIOSESP)/(self.Mstruct*Organ.alpha)) +Organ.K_storage)*Organ.delta_t
    
    def calculate_S_sucrose(self, TRIOSESP):
        '''Flow from TRIOSESP to SUCROSE
        '''
        return ((max(0,TRIOSESP)/(self.Mstruct*Organ.alpha)) * Organ.Vmax_sucrose) / ((max(0, TRIOSESP)/(self.Mstruct*Organ.alpha)) +Organ.K_sucrose)*Organ.delta_t

    def calculate_Loading_sucrose(self, SUCROSE, SUCROSE_phloem):
        return (max(0, SUCROSE)/(self.Mstruct*Organ.alpha)) * ((max(0, SUCROSE)/(self.Mstruct*Organ.alpha)) - (max(0, SUCROSE_phloem)/(Organ.Mstruct_axis*Organ.alpha_axis))) * (Organ.sigma * self.Mstruct**(2/3))*Organ.delta_t
    
    # COMPARTMENTS
    
    def calculate_STORAGE_derivative(self, S_storage, D_storage):
        '''µmol of C'''
        return (S_storage - D_storage) * (self.Mstruct*Organ.alpha)
    
    def calculate_SUCROSE_derivative(self, S_sucrose, D_storage, Loading_sucrose): 
        '''µmol of C'''
        return (S_sucrose + D_storage - Loading_sucrose) * (self.Mstruct*Organ.alpha)
    
    def calculate_TRIOSESP_derivative(self, Photosynthesis, S_sucrose, S_storage):
        '''µmol of C''' 
        return Photosynthesis - (S_sucrose + S_storage) * (self.Mstruct*Organ.alpha)


class Lamina(Organ):
    
    laminae_inflexion_points = {'lamina1': (600, 78.75), 
                                'lamina2': (480, 68.61), 
                                'lamina3': (360, 48.76)}
    
    def __init__(self, Area, Mstruct, Rdark, Assimilation, STORAGE_0, 
                 SUCROSE_0, TRIOSESP_0, name='lamina'):
        super(Lamina, self).__init__(name)
        # parameters
        self.Area = Area # Leaf Area (m-2) of a flag leaf (Bertheloot, 2011) TODO: if the publication is about the concept of the parameter: OK. If the publication is about the value of the parameter: move this doc to test_cn_model. 
        self.Mstruct = Mstruct # Structural mass (g) of a flag leaf (Bertheloot, 2011)
        self.Rdark = Rdark
        self.Assimilation = Assimilation
        
        # linear interpolation of Assimilation      
        self.An_linear_interpolation = None
        
        # initialization
        # flow to phloem
        self.Loading_sucrose = 0
        # compartments
        self.STORAGE_0 = STORAGE_0
        self.SUCROSE_0 = SUCROSE_0
        self.TRIOSESP_0 = TRIOSESP_0
        
    def get_initial_conditions(self):
        return [self.STORAGE_0, self.SUCROSE_0, self.TRIOSESP_0] # keep this order !
        
    # VARIABLES
    
    def calculate_Photosynthesis(self, t):
        t_inflexion, value_inflexion = Lamina.laminae_inflexion_points.get(self.name, (float("inf"), None))
        if t <= t_inflexion:
            Photosynthesis = self.An_linear_interpolation(t)*self.Area*Organ.delta_t
        else:
            Photosynthesis = max(0, self.An_linear_interpolation(t) * ((-0.0721*t + value_inflexion)/10000) * Organ.delta_t)
        return Photosynthesis
    
    def calculate_Rd(self, Photosynthesis):
        if Photosynthesis == 0:
            Rd = self.Rdark*self.Area*Organ.delta_t
        else:
            Rd = 0
        return Rd
    
    # COMPARTMENTS
    
    def calculate_SUCROSE_derivative(self, S_sucrose, D_storage, Loading_sucrose, Rd):
        '''µmol of C'''
        return (S_sucrose + D_storage - Loading_sucrose) * (self.Mstruct*Organ.alpha) - Rd
    
    
class Phloem(Organ):
    
    def __init__(self, SUCROSE_0, Respiration_0, name='phloem'):
        super(Phloem, self).__init__(name)
        # initialization
        # compartments
        self.SUCROSE_0 = SUCROSE_0
        self.Respiration_0 = Respiration_0 # TODO: self.Respiration_0 not used anywhere. Keep it anyway?
        
    def get_initial_conditions(self):
        return [self.SUCROSE_0, self.Respiration_0]
        
    # VARIABLES
    
    def calculate_Conc_Sucrose(self, SUCROSE):
        return (SUCROSE/Organ.Mstruct_axis)/12
    
    def calculate_Conc_C_Sucrose(self, SUCROSE): 
        return SUCROSE/(Organ.Mstruct_axis*Organ.alpha_axis)
    
    # FLOWS
    
    def calculate_Maintenance_respiration(self):
        return 0.004208754*Organ.delta_t
    
    # COMPARTMENTS
        
    def calculate_SUCROSE_derivative(self, organs):
        '''µmol of C'''
        SUCROSE_derivative = 0
        for organ_ in organs:
            if isinstance(organ_, Lamina):
                SUCROSE_derivative += organ_.Loading_sucrose*organ_.Mstruct*Organ.alpha
            elif isinstance(organ_, Sheath):
                SUCROSE_derivative += (organ_.Loading_sucrose + organ_.D_fructan - organ_.S_fructan) * (organ_.Mstruct*Organ.alpha)
            elif isinstance(organ_, Internode):
                SUCROSE_derivative += (organ_.Loading_sucrose + organ_.D_fructan - organ_.S_fructan) * (organ_.Mstruct*Organ.alpha)
            elif isinstance(organ_, Peduncle):
                SUCROSE_derivative += (organ_.Loading_sucrose + organ_.D_fructan - organ_.S_fructan) * (organ_.Mstruct*Organ.alpha)
            elif isinstance(organ_, Chaff):
                SUCROSE_derivative += (organ_.Loading_sucrose*organ_.Mstruct*Organ.alpha)
            elif isinstance(organ_, Grains):
                SUCROSE_derivative -= (organ_.Loading_sucrose_structure + (organ_.Loading_sucrose_storage * ((organ_.STRUCTURE/1E6)*12)))
            elif isinstance(organ_, Roots):
                SUCROSE_derivative -= (organ_.Loading_sucrose * organ_.Mstruct)
        return SUCROSE_derivative
    
    def calculate_Respiration_derivative(self, Maintenance_respiration):
        '''From Evers et al (2010)'''
        return Maintenance_respiration*Organ.Mstruct_axis
    

class Chaff(Organ):
    
    def __init__(self, Area, Mstruct, Assimilation, STORAGE_0, SUCROSE_0, 
                 TRIOSESP_0, name='chaff'):
        super(Chaff, self).__init__(name)
        # parameters
        self.Area = Area
        self.Mstruct = Mstruct
        self.Assimilation = Assimilation
        
        # linear interpolation of Assimilation      
        self.An_linear_interpolation = None
        
        # initialization
        # flow to phloem
        self.Loading_sucrose = 0
        # compartments
        self.STORAGE_0 = STORAGE_0
        self.SUCROSE_0 = SUCROSE_0
        self.TRIOSESP_0 = TRIOSESP_0
        
    def get_initial_conditions(self):
        return [self.STORAGE_0, self.SUCROSE_0, self.TRIOSESP_0] # keep this order !
        
    # COMPARTMENTS
        
    def calculate_STORAGE_derivative(self, S_storage, D_storage):
        '''µmol of C'''
        return (S_storage ) * (self.Mstruct*Organ.alpha)-D_storage
    
    def calculate_SUCROSE_derivative(self, S_sucrose, D_storage, Loading_sucrose):
        '''µmol of C'''
        return (S_sucrose  - Loading_sucrose) * (self.Mstruct*Organ.alpha)+D_storage
        

class Internode(Organ):
    
    def __init__(self, Area, Mstruct, Assimilation, FRUCTAN_0, STORAGE_0, 
                 SUCROSE_0, TRIOSESP_0, name='internode'):
        super(Internode, self).__init__(name)
        # parameters
        self.Area = Area
        self.Mstruct = Mstruct
        self.Assimilation = Assimilation
        
        # linear interpolation of Assimilation
        self.An_linear_interpolation = None
        
        # initialization
        # flow to phloem
        self.Loading_sucrose = 0
        self.D_fructan = 0
        self.S_fructan = 0
        
        # compartments
        self.FRUCTAN_0 = FRUCTAN_0
        self.STORAGE_0 = STORAGE_0
        self.SUCROSE_0 = SUCROSE_0
        self.TRIOSESP_0 = TRIOSESP_0
        
    def get_initial_conditions(self):
        return [self.FRUCTAN_0, self.STORAGE_0, self.SUCROSE_0, self.TRIOSESP_0]
        
    # VARIABLES
    
    def calculate_Conc_Fructan(self, FRUCTAN): 
        return (FRUCTAN/self.Mstruct)/6
    
    # FLOWS
    
    def calculate_D_fructan(self, SUCROSE_phloem, FRUCTAN):
        return min( (Organ.K_Dfructan * Organ.Vmax_Dfructan) / ((max(0, SUCROSE_phloem)/(Organ.Mstruct_axis*Organ.alpha_axis)) + Organ.K_Dfructan) , max(0, FRUCTAN))*Organ.delta_t
    
    def calculate_S_fructan(self, SUCROSE_phloem):
        return ((max(0, SUCROSE_phloem)/(Organ.Mstruct_axis*Organ.alpha_axis)) * Organ.Vmax_Sfructan) / ((max(0, SUCROSE_phloem)/(Organ.Mstruct_axis*Organ.alpha_axis)) + Organ.K_Sfructan)*Organ.delta_t
    
    # COMPARTMENTS

    def calculate_FRUCTAN_derivative(self, S_fructan, D_fructan): 
        return (S_fructan - D_fructan)* (self.Mstruct*Organ.alpha)
    
    
class Peduncle(Organ):
    
    def __init__(self, Area, Mstruct, Assimilation, FRUCTAN_0, STORAGE_0, 
                 SUCROSE_0, TRIOSESP_0, name='peduncle'):
        super(Peduncle, self).__init__(name)
        # parameters
        self.Area = Area
        self.Mstruct = Mstruct
        self.Assimilation = Assimilation

        # linear interpolation of Assimilation
        self.An_linear_interpolation = None
        
        # initialization
        # flow to phloem
        self.Loading_sucrose = 0
        self.D_fructan = 0
        self.S_fructan = 0
        
        # compartments
        self.FRUCTAN_0 = FRUCTAN_0
        self.STORAGE_0 = STORAGE_0
        self.SUCROSE_0 = SUCROSE_0
        self.TRIOSESP_0 = TRIOSESP_0
        
    def get_initial_conditions(self):
        return [self.FRUCTAN_0, self.STORAGE_0, self.SUCROSE_0, self.TRIOSESP_0]
        
    # VARIABLES
    
    def calculate_Conc_Fructan(self, FRUCTAN):
        return (FRUCTAN/self.Mstruct)/6
    
    # FLOWS
    
    def calculate_D_fructan(self, SUCROSE_phloem, FRUCTAN):
        return min( (Organ.K_Dfructan * Organ.Vmax_Dfructan) / ((max(0, SUCROSE_phloem)/(Organ.Mstruct_axis*Organ.alpha_axis)) + Organ.K_Dfructan) , max(0, FRUCTAN))*Organ.delta_t
    
    def calculate_S_fructan(self, SUCROSE_phloem):
        return ((max(0, SUCROSE_phloem)/(Organ.Mstruct_axis*Organ.alpha_axis)) * Organ.Vmax_Sfructan) / ((max(0, SUCROSE_phloem)/(Organ.Mstruct_axis*Organ.alpha_axis)) + Organ.K_Sfructan)*Organ.delta_t
    
    # COMPARTMENTS
    
    def calculate_FRUCTAN_derivative(self, S_fructan, D_fructan):
        return (S_fructan - D_fructan)* (self.Mstruct*Organ.alpha)
    
    
class Sheath(Organ):
    
    def __init__(self, Area, Mstruct, Assimilation, FRUCTAN_0, STORAGE_0, 
                 SUCROSE_0, TRIOSESP_0, name='sheath'):
        super(Sheath, self).__init__(name)
        
        # parameters
        self.Area = Area
        self.Mstruct = Mstruct
        self.Assimilation = Assimilation

        # linear interpolation of Assimilation
        self.An_linear_interpolation = None
        
        # initialization
        # flow to phloem
        self.Loading_sucrose = 0
        self.D_fructan = 0
        self.S_fructan = 0
        
        # compartments
        self.FRUCTAN_0 = FRUCTAN_0
        self.STORAGE_0 = STORAGE_0
        self.SUCROSE_0 = SUCROSE_0
        self.TRIOSESP_0 = TRIOSESP_0

    def get_initial_conditions(self):
        return [self.FRUCTAN_0, self.STORAGE_0, self.SUCROSE_0, self.TRIOSESP_0]
        
    # VARIABLES
    
    def calculate_Conc_Fructan(self, FRUCTAN): 
        return (FRUCTAN/self.Mstruct)/6
    
    # FLOWS
    
    def calculate_D_fructan(self, SUCROSE_phloem, FRUCTAN): 
        return min( (Organ.K_Dfructan * Organ.Vmax_Dfructan) / ((max(0, SUCROSE_phloem)/(Organ.Mstruct_axis*Organ.alpha_axis)) + Organ.K_Dfructan) , max(0, FRUCTAN))*Organ.delta_t
    
    def calculate_S_fructan(self, SUCROSE_phloem):
        return ((max(0, SUCROSE_phloem)/(Organ.Mstruct_axis*Organ.alpha_axis)) * Organ.Vmax_Sfructan) / ((max(0, SUCROSE_phloem)/(Organ.Mstruct_axis*Organ.alpha_axis)) + Organ.K_Sfructan)*Organ.delta_t
        
    # COMPARTMENTS
    
    def calculate_FRUCTAN_derivative(self, S_fructan, D_fructan): 
        return (S_fructan - D_fructan)*(self.Mstruct*Organ.alpha)
    

class Grains(Organ):
    
    Grain_storage = 0
    Grain_structure = 0
    filling_init = 360
    K_storage = 500
    K_RGR = 300
    Vmax_storage = 0.125
    Vmax_RGR = 1.9e-06
    Y_grains = 0.75
    
    def __init__(self, STORAGE_0, STRUCTURE_0, name='grains'):
        super(Grains, self).__init__(name)
        
        # initialization
        # flow to phloem
        self.Loading_sucrose_structure = 0
        self.Loading_sucrose_storage = 0
        self.STRUCTURE = 0
        
        # compartments
        self.STORAGE_0 = STORAGE_0
        self.STRUCTURE_0 = STRUCTURE_0
        
    def get_initial_conditions(self):
        return [self.STORAGE_0, self.STRUCTURE_0]
        
    # VARIABLES
    
    def calculate_Dry_mass(self, STRUCTURE, STORAGE): 
        '''g of C'''
        return ((STRUCTURE + STORAGE)/1000000)*12
    
    def calculate_RGR_structure(self, SUCROSE_phloem):
        return ((max(0, SUCROSE_phloem)/(Organ.Mstruct_axis*Organ.alpha_axis)) * Grains.Vmax_RGR) / ((max(0, SUCROSE_phloem)/(Organ.Mstruct_axis*Organ.alpha_axis)) + Grains.K_RGR)
    
    # FLOWS
    
    def calculate_Loading_sucrose_storage(self, t, SUCROSE_phloem):
        if t<=Grains.filling_init:
            Loading_sucrose_storage = 0
        else:
            Loading_sucrose_storage = ((max(0, SUCROSE_phloem)/(Organ.Mstruct_axis*Organ.alpha_axis)) * Grains.Vmax_storage) / ((max(0, SUCROSE_phloem)/(Organ.Mstruct_axis*Organ.alpha_axis)) + Grains.K_storage)*Organ.delta_t
        return Loading_sucrose_storage
    
    def calculate_Loading_sucrose_structure(self, t, STRUCTURE, RGR_structure):
        if t<=Grains.filling_init: 
            Loading_sucrose_structure = STRUCTURE * RGR_structure * Organ.delta_t
        else:
            Loading_sucrose_structure = 0
        return Loading_sucrose_structure
    
    # COMPARTMENTS
    
    def calculate_STORAGE_derivative(self, Loading_sucrose_storage, STRUCTURE):
        '''Grain filling'''
        return Loading_sucrose_storage* Grains.Y_grains * ((STRUCTURE/1E6)*12)
    
    def calculate_STRUCTURE_derivative(self, Loading_sucrose_structure): 
        '''Endosperm cell division'''
        return Loading_sucrose_structure*Grains.Y_grains
    

class Roots(Organ):
    K_roots = 100
    Vmax_roots = 0.015
    
    def __init__(self, Mstruct, Sucrose_0, name='roots'):
        super(Roots, self).__init__(name)
        
        # parameters
        self.Mstruct = Mstruct

        # initialization
        # flow to phloem
        self.Loading_sucrose = 0
        
        # compartments
        self.Sucrose_0 = Sucrose_0
        
    def get_initial_conditions(self):
        return [self.Sucrose_0] # keep this order !
    
    # VARIABLES
        
    def calculate_Dry_mass(self, Sucrose):
        '''g of C'''
        return (Sucrose*12)/1000000
    
    # FLOWS
    
    def calculate_Loading_sucrose(self, SUCROSE_phloem):
        return ((max(0, SUCROSE_phloem)/(Organ.Mstruct_axis*Organ.alpha_axis)) * Roots.Vmax_roots) / ((max(0, SUCROSE_phloem)/(Organ.Mstruct_axis*Organ.alpha_axis)) + Roots.K_roots)*Organ.delta_t
    
    # COMPARTMENTS

    def calculate_Sucrose_derivative(self, Loading_sucrose):
        return Loading_sucrose*self.Mstruct
    
    