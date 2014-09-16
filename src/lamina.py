# -*- coding: latin-1 -*-

'''
Created on 11 août 2014

@author: cchambon
'''

from __future__ import division # use "//" to do integer division

class Lamina(object):
    
    alpha_lamina = 1
    delta_Dstorage = 0.0001
    K_storage = 20
    K_sucrose = 0.66
    Vmax_storage = 2
    Vmax_sucrose = 1
    
    def __init__(self, lamina_area, Mstruct_lamina, Assimilation, STORAGE_0, 
                 SUCROSE_lamina_0, TRIOSESP_0):
        self.lamina_area = lamina_area
        self.Mstruct_lamina = Mstruct_lamina
        self.Assimilation = Assimilation
        self.STORAGE_0 = STORAGE_0
        self.SUCROSE_lamina_0 = SUCROSE_lamina_0
        self.TRIOSESP_0 = TRIOSESP_0
        
    def get_initial_conditions(self):
        return [self.STORAGE_0, self.SUCROSE_lamina_0, self.TRIOSESP_0]
        
    def calculate_Photosynthesis(self, An):
        return An * self.lamina_area
    
    def calculate_D_storage(self, STORAGE):
        '''Flow from STORAGE to SUCROSE_lamina
        '''
        return max(0, self.delta_Dstorage * (STORAGE/(self.Mstruct_lamina*self.alpha_lamina)))  

    def calculate_S_storage(self, TRIOSESP):
        '''Flow from TRIOSESP to STORAGE
        '''
        return ((max(0, TRIOSESP)/(self.Mstruct_lamina*self.alpha_lamina)) * self.Vmax_storage) / ((max(0, TRIOSESP)/(self.Mstruct_lamina*self.alpha_lamina)) +self.K_storage)
        
    def calculate_S_sucrose(self, TRIOSESP):
        '''Flow from TRIOSESP to SUCROSE_lamina
        '''
        return ((max(0,TRIOSESP)/(self.Mstruct_lamina*self.alpha_lamina)) * self.Vmax_sucrose) / ((max(0, TRIOSESP)/(self.Mstruct_lamina*self.alpha_lamina)) +self.K_sucrose)
    
    def calculate_STORAGE_derivative(self, S_storage, D_storage):
        '''µmol of C'''
        return (S_storage - D_storage) * (self.Mstruct_lamina*self.alpha_lamina)
    
    def calculate_SUCROSE_lamina_derivative(self, S_sucrose, D_storage):
        '''µmol of C'''
        return (S_sucrose + D_storage) * (self.Mstruct_lamina*self.alpha_lamina)
    
    def calculate_TRIOSESP_derivative(self, Photosynthesis, S_sucrose, S_storage):
        '''µmol of C''' 
        return Photosynthesis - (S_sucrose + S_storage) * (self.Mstruct_lamina*self.alpha_lamina)
    
    
    