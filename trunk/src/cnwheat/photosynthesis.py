#!/usr/bin/env python
# -*- coding: latin-1 -*-
"""
    cnwheat.photosynthesis
    ~~~~~~~~~~~~~~~~~~~~~~

    Model of photosynthesis based on Farquhar's approach. The model includes the dependance of photosynthesis to organ temperature and nitrogen content. Internal CO2 and organ temperature are found numerically.

    :copyright: Copyright 2014 INRA-EGC, see AUTHORS.txt.
    :license: TODO, see LICENSE.txt for details.
"""

from __future__ import division # use "//" to do integer division

from math import sqrt, log,  exp

class PhotosynthesisModel(object):

    O = 21000      #: Photosynthetic parameter: Intercellular O2 concentration, umol mol(air)-1 or Pa, from Bernacchi et al. (2001)
    Kc25 = 404     #: Photosynthetic parameter: Affinity constant of RuBisCO for C, umol mol-1 or Pa, from Bernacchi et al. (2001) (estimation in Braune et al. (2009) not enough accurate)
    Ko25 = 278.4E3 #: Photosynthetic parameter: Affinity constant of RuBisCO for O, umol mol-1 or Pa, from Bernacchi et al. (2001) (estimation in Braune et al. (2009) not enough accurate)
    Gamma25 = 39   #: Photosynthetic parameter: CO2 compensation point, umol(CO2) mol-1 (air), from Braune et al. (2009)
    theta = 0.72    #: Photosynthetic parameter: curvature parameter of J, dimensionless

    #: Nitrogen dependance of photosynthetic parameters
    #: Derived from Braune et al. (2009):
    #:     * S_Na: slope of the relation between Na and the parameters (umol g-1 s-1)
    #:     * Na_min: minimum amount of leaf nitrogen below which photosynthesis rate is zero (g (N) m-2 leaf)
    #:     * Gamma_Na1 and Gamma_Na2: parameters of alpha dependance to Na (mol mol-1 and m2 g-1 respectively)
    #:     * delta1 and delta2: parameters of m (scaling factor of gs) dependance to Na (m2 g-1 and dimensionless respectively)
    param_N = {'S_Na': {'Vc_max25': 63.2, 'Jmax25': 151, 'TPU25': 9.25, 'Rdark25': 0.493}, 'Na_min': {'Vc_max25': 0.198, 'Jmax25': 0.225, 'TPU25': 0.229, 'Rdark25': 0.118},
                'Gamma_Na1': 0.437, 'Gamma_Na2': 2.29, 'delta1': 14.7, 'delta2': -0.548}

    gsmin = 0.05            #: Stomatal conductance parameter: Minimum gs, measured in the dark (mol m-2 s-1). Braune et al. (2009).
    gb = 3.5                #: Stomatal conductance parameter: Boundary layer conductance (mol m-2 s-1). Muller et al., (2005)

    sigma = 5.6704E-8       #: Physical parameter: Stefan-Bolzmann constant (W-2 K-4)
    I0 = 1370               #: Physical parameter: Extraterrestrial solar radiation (W m-2)
    Lambda = 2260E3         #: Physical parameter: Latent heat for vaporisation of water (J kg-1)
    rhocp = 1256            #: Physical parameter: Volumetric heat capacity of air (J m-3 K-1)
    gamma = 66E-3           #: Physical parameter: Psychrometric constant (KPa K-1). Mean value
    a = 2.5                 #: Physical parameter: Attenuation coefficient of wind within a wheat canopy. From Campbell and Norman (1998), second edition. Can also be estimaed by: a = sqrt((0.2*LAI*h)/sqrt((4*width*h)/(pi*LAI))
    R = 8.3144              #: Physical parameter: Gas constant (J mol-1 K-1)
    Patm = 1.01325E5        #: Physical parameter: Atmospheric pressure (Pa)

    #TODO: trouver des bonnes valeurs de fr, ft
    fr = 0.15               #: Physical parameter: Leaf radiation reflectance
    ft = 0.15               #: Physical parameter: Leaf radiation transmittance

    #: Temperature dependance of photosynthetic parameters
    #: Parameter values derived from Braune et al. (2009) except for Kc, Ko, and Rdark (Bernacchi et al., 2001)
    #:     * deltaHa, deltaHd: enthalpie of activation and deactivation respectively (kJ mol-1)
    #:     * deltaS: entropy term (kJ mol-1 K-1)
    #:     * Tref: reference temperature (K)
    #:     * R: universal gas constant (kJ mol-1 K-1)
    param_temp = {'deltaHa': {'Vc_max': 89.7, 'Jmax': 48.9, 'TPU': 47., 'Kc': 79.43, 'Ko': 36.38, 'Gamma': 35., 'Rdark': 46.39},
                  'deltaHd': {'Vc_max': 149.3, 'Jmax': 152.3, 'TPU': 152.3},
                  'deltaS': {'Vc_max': 0.486, 'Jmax': 0.495, 'TPU': 0.495},
                  'Tref': 298.15, 'R': 8.3145E-03}

    leaf_width = 0.01 #: m
    H_canopy = 1      #: m
    H_organ = 0.5     #: m
    Na_init = 2.5     #: g m-2

    @classmethod
    def leaf_temperature(cls, leaf_width, z, H, wind0, PAR, gs, Ta, Tleaf, RH):
        """
        Energy balance for the estimation of leaf temperature
            - leaf_width (m)
            - z: height of leaf from soil (m)
            - H: canopy height (m)
            - Wind0: wind at the top of the canopy (m s-1)
            - PAR (umol m-2 s-1)
            - gs: stomatal conductance (mol m-2 s-1)
            - Ta: air temperature (degree C)
            - Tleaf: leaf temperature (degree C). Tleaf = Ta at the first iteration of the numeric resolution
            - RH: Relative humidity (decimal fraction)
        """

        # Wind speed (m s-1)
        wind = wind0 * exp(cls.a*(z/H -1))                      # From Campbell and Norman (1998), second edition. z: organ height, H: height of the canopy

        # Boundary layer restistance to heat (s m-1)
        rbh = 100*sqrt(leaf_width/wind)

        # Turbulence resistance to heat (s m-1)
        rt = (0.74 * log(((2 - 0.7*H) / (0.1*H))**2)) / (0.16*wind)

        # Net absorbed radiation Rn (PAR and NIR, J m-2 s-1)
        Iabs = (PAR*(1-cls.fr-cls.ft))/(0.55*4.55)          # Global absorbed radiation by leaf (J m-2 s-1). TODO: relation a verifier
        es_Ta = 0.611 * exp((17.4*Ta)/(239+Ta))             # Saturated vapour pressure of the air (kPa), Ta in degree Celsius
        V = RH * es_Ta                                      # Vapour pressure of the air (kPa)
        fvap = 0.56 - 0.079*sqrt(10*V)                      # Fraction of vapour pressure

        tau = Iabs/cls.I0                                   # Atmospheric transmissivity (dimensionless)
        fclear = 0.1 + 0.9*max(0, min(1, (tau-0.2)/0.5))    # Fraction sky clearness

        Rn = Iabs - cls.sigma * (Tleaf+273)**4*fvap*fclear

        # Transpiration (mm s-1), Penman-Monteith
        if Tleaf == Ta:
            Ta_K = Ta + 273.15                              # Ta in kelvin
            s = ((17.4*239)/(Ta_K + 239)**2)*es_Ta          # Slope of the curve relating saturation vapour pressure to temperature (kPa K-1)
        else:
            es_Tl = 0.611 * exp((17.4*Tleaf)/(239+Tleaf))   # Saturated vapour pressure at leaf (kPa), Tleaf in degree Celsius
            Tleaf_K, Ta_K = Tleaf + 273.15, Ta + 273.15     # Temperatures in kelvin
            s = (es_Tl - es_Ta)/(Tleaf_K - Ta_K)            # Slope of the curve relating saturation vapour pressure to temperature (kPa K-1)

        VPDa = es_Ta - V
        rbw = 0.96 * rbh                                    # Boundary layer resistance for water (s m-1)
        gsw = (1.6*gs * cls.R * (Tleaf+273.15)) / cls.Patm  # Stomatal conductance to water (m s-1). 1.6 convert gs_CO2 in gs_water. Relation given by A. Tuzet (2003)
        rswp = 1/gsw                                        # Stomatal resistance for water (s m-1)

        Ep = (s * Rn + (cls.rhocp * VPDa)/(rbh + rt)) / (cls.Lambda * (s + cls.gamma*((rbw + rt + rswp)/(rbh + rt))))

        # Leaf temperature
        Tleaf = Ta + ((rbh + rt) * (Rn - cls.Lambda*Ep)) / cls.rhocp
        return Tleaf, Ep

    @classmethod
    def stomatal_conductance(cls, Ag, An, Na, Ca, RH):
        """
        BWB model of stomatal conductance
            - Ag: global assimilation (umol m-2 s-1)
            - An: net assimilation (umol m-2 s-1)
            - Na: nitrogen content of leaf (g m-2)
            - Ca: Air CO2 (umol mol-1)
            - RH: Relative humidity (decimal fraction)
        """

        Cs = Ca - An *(1.37/(cls.gb))                            # CO2 concentration at leaf surface (umol mol-1 or Pa). From Prieto et al. (2012). gb in mol m-2 s-1
        m = cls.param_N['delta1'] * Na**cls.param_N['delta2']    # Scaling factor dependance to Na (dimensionless). This focntion is maintained although I'm not sure that it should be taken into account
        gs = (cls.gsmin + m*((Ag*RH)/(Cs)))                      # Stomatal conductance (mol m-2 s-1), from Braune et al. (2009), Muller et al. (2005): using Ag rather than An. Would be better with a function of VPD and with (Ci-gamma) instead of Cs.
        return gs

    @classmethod
    def f_temperature(cls, pname, p25, T):
        """
        Photosynthetic parameters relation to temperature
            - pname: name of parameter
            - p25: parameter value at 25 degree C
            - T: leaf temperature (degree C)
        """
        Tk = T + 273.15
        deltaHa = cls.param_temp['deltaHa'][pname]
        Tref = cls.param_temp['Tref']
        R = cls.param_temp['R']

        f_activation = exp((deltaHa * (Tk - Tref))/(R * Tref * Tk))

        if pname in ('Vc_max', 'Jmax', 'TPU'):
            deltaS = cls.param_temp['deltaS'][pname]
            deltaHd = cls.param_temp['deltaHd'][pname]
            f_deactivation = (1 + exp((Tref*deltaS - deltaHd) / (Tref*R))) / (1 + exp((Tk*deltaS - deltaHd) / (Tk*R)))
        else:
            f_deactivation = 1

        p = p25 * f_activation * f_deactivation

        return p

    @classmethod
    def photosynthesis(cls, PAR, Na, Tleaf, Ci):
        """
        Estimation of C02 assimilation. In this version, most of parameters are derived from Braune et al. (2009) on barley
            - PAR: PAR intercepted by leaf (umol m-2 s-1)
            - Na: nitrogen content of leaf (g m-2)
            - Tleaf: leaf temperature (degree C)
            - Ci: internal CO2 (umol mol-1), by default = 0.7*CO2air
        """

        ### # RuBisCO-limited carboxylation rate ###
        # RuBisCO parameters dependance to temperature
        Kc = cls.f_temperature('Kc', cls.Kc25, Tleaf)
        Ko = cls.f_temperature('Ko', cls.Ko25, Tleaf)
        Gamma = cls.f_temperature('Gamma', cls.Gamma25, Tleaf)

        # Vcmax
        Vc_max25 = 84.965 * (Na - 0.17)                                                     # Relation between Vc_max25 and Na
        Vc_max = cls.f_temperature ('Vc_max', Vc_max25, Tleaf)                              # Relation between Vc_max and temperature
        Ac = (Vc_max * (Ci-Gamma)) / (Ci + Kc * (1 + cls.O/Ko))                             # Rate of assimilation under Vc_max limitation
        ### RuBP regeneration-limited carboxylation rate via electron transport ###
        alpha = 0.0413 * Na + 0.2101                                                        # Relation between alpha and Na
        Jmax25 = 117.6 * (Na - 0.17)                                                        # Relation between Jmax25 and Na
        Jmax = cls.f_temperature('Jmax', Jmax25, Tleaf)                                     # Relation between Jmax and temperature

        # Electron transport rate
        J = ((Jmax+alpha*PAR) - sqrt((Jmax+alpha*PAR)**2 - 4*cls.theta*alpha*PAR*Jmax))/(2*cls.theta) # Muller et al. (2005), Evers et al. (2010)
        Aj = (J * (Ci-Gamma)) / (4*Ci + 8*Gamma)                                            # Rate of assimilation under RuBP regeneration limitation
        ### Triose phosphate utilisation-limited carboxylation rate ###
        TPU25 = cls.param_N['S_Na']['TPU25'] * (Na - cls.param_N['Na_min']['TPU25'])        # Relation between TPU25 and Na
        TPU = cls.f_temperature('TPU', TPU25, Tleaf)                                        # Relation between TPU and temperature
        Vomax = (Vc_max*Ko*Gamma)/(0.5*Kc*cls.O)                                            # Maximum rate of Vo (umol m-2 s-1)
        Vo = (Vomax * cls.O) / (cls.O + Ko*(1+Ci/Kc))                                       # Rate of oxygenation of RuBP (umol m-2 s-1)
        Ap = (1-Gamma/Ci)*(3*TPU) + Vo                                                      # Rate of assimilation under TPU limitation

        # Gross assimilation rate
        Ag = min(Ac, Aj, Ap)

        # Mitochondrial respiration rate of leaf in light Rd (processes other than photorespiration)
        Rdark25 = cls.param_N['S_Na']['Rdark25'] * (Na - cls.param_N['Na_min']['Rdark25'])  # Relation between Rdark25 (respiration in obscurity at 25 degree C) and Na
        Rdark = cls.f_temperature('Rdark', Rdark25, Tleaf)                                  # Relation between Rdark and temperature
        Rd = Rdark * (0.33 + (1-0.33)*(0.5)**(PAR/15))                                      # Found in Muller et al. (2005), eq. 19
        # Net C assimilation
        An = Ag - Rd
        return An, Ag, Rd


    @classmethod
    def calculate_An(cls, t, PAR, Ta, Ca, RH):
        """
        Compute CO2 assimilation following Farquhar's model and estimate internal CO2 and organ temperature numerically.

        :Parameters:

            - `PAR` (:class:`float`) - PAR

            - `Ta` (:class:`float`) - Air temperature (degree Celsius)

            - `Ca` (:class:`float`) - Air CO2 (umol mol-1)

            - `RH` (:class:`float`) - Relative humidity (decimal fraction)

        :Returns:
            An

        :Returns Type:
            :class:`float`

        """

        ### Physical parameters ###
        wind0 = 5 #: Wind speed at the top of the canopy (m s-1)

        ### Iterations to find leaf temperature and Ci ###
        Ci, Tleaf = 0.7*Ca, Ta # Initial values
        prec_Ci, prec_Tleaf = 0.1, 0.1
        count = 0
        while abs((Ci - prec_Ci)/prec_Ci) >= 0.01 or abs((Tleaf - prec_Tleaf)/prec_Tleaf) >= 0.01:
            if count >=30: # TODO: test a faire? Semble prendre du tps de calcul
                if abs((Ci - prec_Ci)/prec_Ci) >= 0.01:
                    print "Ci cannot converge at t= %s, prec_Ci= %s, Ci= %s" %(t, prec_Ci, Ci)
                else:
                    print "Tleaf cannot converge at t= %s, prec_Tleaf= %s, Tleaf= %s" %(t, prec_Tleaf, Tleaf)
                break
            else:
                prec_Ci, prec_Tleaf = Ci, Tleaf
                An, Ag, Rd = cls.photosynthesis(PAR, cls.Na_init, Tleaf, Ci)
                # Stomatal conductance
                gs = cls.stomatal_conductance(Ag, An, cls.Na_init, Ca, RH)
                # New value of Ci
                Ci = Ca - An * ((1.6/gs) + (1.37/cls.gb)) # gs and gb in mol m-2 s-1
                # New value of Tleaf
                Tleaf, E = cls.leaf_temperature(cls.leaf_width, cls.H_organ, cls.H_canopy, wind0, PAR, gs, Ta, Tleaf, RH)
                count +=1

        return An

