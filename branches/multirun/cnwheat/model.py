# -*- coding: latin-1 -*-

from __future__ import division  # use "//" to do integer division
import numpy as np
from math import exp

import parameters


"""
    cnwheat.model
    ~~~~~~~~~~~~~

    The module :mod:`cnwheat.model` defines the equations of the CN exchanges in a population of plants.

    :copyright: Copyright 2014-2017 INRA-ECOSYS, see AUTHORS.
    :license: CeCILL-C, see LICENSE for details.

    **Acknowledgments**: The research leading these results has received funding through the
    Investment for the Future programme managed by the Research National Agency
    (BreedWheat project ANR-10-BTBR-03).

    .. seealso:: Barillot et al. 2016.
"""

"""
    Information about this versioned file:
        $LastChangedBy$
        $LastChangedDate$
        $LastChangedRevision$
        $URL$
        $Id$
"""


class EcophysiologicalConstants:
    """
    Ecophysiological constants.
    """

    def __init__(self):
        pass

    C_MOLAR_MASS = 12                       #: Molar mass of carbon (g mol-1)
    NB_C_TRIOSEP = 3                        #: Number of C in 1 mol of trioseP
    NB_C_HEXOSES = 6                        #: Number of C in 1 mol of hexoses (glucose, fructose)
    NB_C_SUCROSE = 12                       #: Number of C in 1 mol of sucrose
    HEXOSE_MOLAR_MASS_C_RATIO = 0.42        #: Contribution of C in hexose mass
    TRIOSESP_MOLAR_MASS_C_RATIO = 0.21      #: Contribution of C in triosesP mass
    RATIO_C_mstruct = 0.44                 #: Mean contribution of carbon to structural dry mass (g C g-1 mstruct)

    AMINO_ACIDS_C_RATIO = 4.15              #: Mean number of mol of C in 1 mol of the major amino acids of plants (Glu, Gln, Ser, Asp, Ala, Gly)
    AMINO_ACIDS_N_RATIO = 1.25              #: Mean number of mol of N in 1 mol of the major amino acids of plants (Glu, Gln, Ser, Asp, Ala, Gly)
    PROTEINS_MOLAR_MASS_N_RATIO = 0.151     #: Mean contribution of N in protein mass (Penning De Vries 1989)
    AMINO_ACIDS_MOLAR_MASS_N_RATIO = 0.135  #: Mean contribution of N in amino acids mass of the major amino acids of plants (Glu, Gln, Ser, Asp, Ala, Gly)
    NITRATES_MOLAR_MASS_N_RATIO = 0.23      #: Contribution of N in amino acids mass
    N_MOLAR_MASS = 14                       #: Molar mass of nitrogen (g mol-1)
    AMINO_ACIDS_MOLAR_MASS_C_RATIO = 0.38   #: (Penning De Vries 1989)
    PROTEINS_MOLAR_MASS_C_RATIO = 0.38      #: As for AA


class Population(object):
    """
    The class :class:`Population` defines the CN exchanges at population scale.

    A :class:`population <Population>` must have at least one :class:`plant <Plant>`.
    """

    PARAMETERS = parameters.POPULATION_PARAMETERS  #: the internal parameters of the population

    def __init__(self, plants=None):
        if plants is None:
            plants = []
        self.plants = plants  #: the list of plants

    def calculate_integrative_variables(self):
        """Calculate the integrative variables of the population recursively.
        """
        for plant in self.plants:
            plant.calculate_integrative_variables()


class Plant(object):
    """
    The class :class:`Plant` defines the CN exchanges at plant scale.

    A :class:`plant <Plant>` must have at least one :class:`axis <Axis>`.
    """

    PARAMETERS = parameters.PLANT_PARAMETERS  #: the internal parameters of the plants

    def __init__(self, index=None, axes=None):
        self.index = index  #: the index of the plant
        if axes is None:
            axes = []
        self.axes = axes   #: the list of axes
        self.cohorts = []  #: list of cohort values - Hack to treat tillering cases : TEMPORARY

    def calculate_integrative_variables(self):
        """Calculate the integrative variables of the plant recursively.
        """
        for axis in self.axes:
            axis.calculate_integrative_variables()

    @staticmethod
    def calculate_temperature_effect_on_conductivity(Tair):
        """Effect of the temperature on phloeme translocation conductivity (Farrar 1988)
        Should multiply the rate at 20�C

        :param float Tair: Air temperature (�C)

        :return: Correction to apply to conductivity coefficients.
        :rtype: float
        """
        Q10 = 1.3
        Tref = 20.

        return Q10**((Tair-Tref)/10.)

    @staticmethod
    def calculate_temperature_effect_on_Vmax(Tair):
        """Effect of the temperature on maximal enzyme activity
        Should multiply the rate at 20�C

        :param float Tair: Air temperature (�C)

        :return: Correction to apply to enzyme activity
        :rtype: float
        """
        Tref = 20 + 273.15
        Tk = Tair + 273.15
        R = 8.3144     #: Physical parameter: Gas constant (J mol-1 K-1)
        deltaHa = 55   #: Enthalpie of activation of parameter pname (kJ mol-1)
        deltaS = 0.48  #: entropy term of parameter pname (kJ mol-1 K-1)
        deltaHd = 154  #: Enthalpie of deactivation of parameter pname (kJ mol-1)

        f_activation = np.exp((deltaHa * (Tk - Tref)) / (R * 1E-3 * Tref * Tk))  #: Energy of activation (normalized to unity)

        f_deactivation = (1 + np.exp((Tref * deltaS - deltaHd) / (Tref * R * 1E-3))) / (1 + np.exp((Tk * deltaS - deltaHd) / (Tk * R * 1E-3)))  #: Energy of deactivation (normalized to unity)

        return f_activation * f_deactivation


class Axis(object):
    """
    The class :class:`Axis` defines the CN exchanges at axis scale.

    An :class:`axis <Axis>` must have:
        * one :class:`set of roots <Roots>`,
        * one :class:`phloem <Phloem>`,
        * zero or one :class:`set of grains <Grains>`,
        * at least one :class:`phytomer<Phytomer>`.
    """

    PARAMETERS = parameters.AXIS_PARAMETERS  #: the internal parameters of the axes
    INIT_COMPARTMENTS = parameters.AXIS_INIT_COMPARTMENTS  #: the initial values of compartments and state parameters

    def __init__(self, label=None, roots=None, phloem=None, grains=None, phytomers=None, C_exudated=INIT_COMPARTMENTS.C_exudated, sum_respi_shoot=INIT_COMPARTMENTS.sum_respi_shoot,
                 sum_respi_roots=INIT_COMPARTMENTS.sum_respi_roots):

        self.label = label    #: the label of the axis
        self.roots = roots    #: the roots
        self.phloem = phloem  #: the phloem
        self.grains = grains  #: the grains
        if phytomers is None:
            phytomers = []
        self.phytomers = phytomers  #: the list of phytomers

        # state variables
        self.C_exudated = C_exudated
        self.sum_respi_shoot = sum_respi_shoot
        self.sum_respi_roots = sum_respi_roots

        # integrative variables
        self.Total_Transpiration = None  #: the total transpiration (mmol s-1)
        self.mstruct = None  #: structural mass of the axis (g)
        self.senesced_mstruct = None  #: senesced structural mass of the axis (g)

    def calculate_integrative_variables(self):
        """Calculate the integrative variables of the axis recursively.
        """
        self.mstruct = 0
        self.senesced_mstruct = 0
        if self.roots is not None:
            self.roots.calculate_integrative_variables()
            self.mstruct += self.roots.mstruct
            self.senesced_mstruct += self.roots.senesced_mstruct
        if self.phloem is not None:
            self.phloem.calculate_integrative_variables()
        if self.grains is not None:
            self.grains.calculate_integrative_variables()
            self.mstruct += self.grains.structural_dry_mass
        for phytomer in self.phytomers:
            phytomer.calculate_integrative_variables()
            self.mstruct += phytomer.mstruct * phytomer.nb_replications
            self.senesced_mstruct += phytomer.senesced_mstruct * phytomer.nb_replications

    # COMPARTMENTS

    @staticmethod
    def calculate_C_exudated(C_exudation, N_exudation, roots_mstruct):
        """delta sucrose

        :param float C_exudation: Rates of sucrose exudated (:math:`\mu mol` C g-1 mstruct h-1)
        :param float N_exudation: Rate of amino acids exudated (:math:`\mu mol` N g-1 mstruct h-1)
        :param float roots_mstruct: RStructural mass of the roots (g)

        :return: delta C loss by exudation (:math:`\mu mol` C)
        :rtype: float
        """
        return (C_exudation + N_exudation * EcophysiologicalConstants.AMINO_ACIDS_C_RATIO/EcophysiologicalConstants.AMINO_ACIDS_N_RATIO) * roots_mstruct



class Phytomer(object):
    """
    The class :class:`Phytomer` defines the CN exchanges at phytomer scale.

    A :class:`phytomer <Phytomer>` must have at least:
        * 1 photosynthetic organ: :class:`chaff <Chaff>`, :class:`peduncle <Peduncle>`,
                                  :class:`lamina <Lamina>`, :class:`internode <Internode>`,
                                  or :class:`sheath <Sheath>`.
        * or 1 :class:`hiddenzone <HiddenZone>`.
    """

    PARAMETERS = parameters.PHYTOMER_PARAMETERS  #: the internal parameters of the phytomers

    def __init__(self, index=None, chaff=None, peduncle=None, lamina=None, internode=None, sheath=None, hiddenzone=None, cohorts=[], cohorts_replications=None):
        self.index = index  #: the index of the phytomer
        self.chaff = chaff  #: the chaff
        self.peduncle = peduncle  #: the peduncle
        self.lamina = lamina  #: the lamina
        self.internode = internode  #: the internode
        self.sheath = sheath  #: the sheath
        self.hiddenzone = hiddenzone  #: the hidden zone
        self.mstruct = None  #: the structural mass of the phytomer (g)
        self.senesced_mstruct = None  #: senesced structural mass of the phytomer (g)
        self.cohorts = cohorts  #: list of cohort values - Hack to treat tillering cases : TEMPORARY. Devrait �tre port� � l'�chelle de la plante uniquement mais je ne vois pas comment faire mieux
        self.cohorts_replications = cohorts_replications  #: dictionary of number of replications per cohort rank

    def calculate_integrative_variables(self):
        """Calculate the integrative variables of the phytomer recursively.
        """
        self.mstruct = 0
        self.senesced_mstruct = 0
        for organ_ in (self.chaff, self.peduncle, self.lamina, self.internode, self.sheath, self.hiddenzone):
            if organ_ is not None:
                organ_.calculate_integrative_variables()
                self.mstruct += organ_.mstruct
                if hasattr(organ_, 'senesced_mstruct'):
                    self.senesced_mstruct += organ_.senesced_mstruct
    @property
    def nb_replications(self):
        return sum(int(v <= self.index) * self.cohorts_replications.get(v, 0) for v in self.cohorts) + 1


class Organ(object):
    """
    The class :class:`Organ` defines the CN exchanges at organ scale.

    :class:`Organ` is the base class of all organs. DO NOT INSTANTIATE IT.
    """

    def __init__(self, label):
        self.label = label  #: the label of the organ

    def initialize(self):
        """Initialize the derived attributes of the organ.
        """
        pass

    def calculate_integrative_variables(self):
        """Calculate the integrative variables of the organ recursively.
        """
        pass


class HiddenZone(Organ):
    """
    The class :class:`HiddenZone` defines the CN exchanges in an hidden zone.
    """

    PARAMETERS = parameters.HIDDEN_ZONE_PARAMETERS                #: the internal parameters of the hidden zone
    INIT_COMPARTMENTS = parameters.HIDDEN_ZONE_INIT_COMPARTMENTS  #: the initial values of compartments and state parameters

    def __init__(self, label='hiddenzone',  mstruct=INIT_COMPARTMENTS.mstruct, Nstruct=INIT_COMPARTMENTS.Nstruct,
                 sucrose=INIT_COMPARTMENTS.sucrose, fructan=INIT_COMPARTMENTS.fructan, amino_acids=INIT_COMPARTMENTS.amino_acids, proteins=INIT_COMPARTMENTS.proteins,
                 ratio_DZ=INIT_COMPARTMENTS.ratio_DZ, ratio_EOZ=INIT_COMPARTMENTS.ratio_EOZ, cohorts=[], cohorts_replications=None, index=None):

        super(HiddenZone, self).__init__(label)

        self.cohorts = cohorts  #: list of cohort values - Hack to treat tillering cases : TEMPORARY. Devrait �tre port� � l'�chelle de la plante uniquement mais je ne vois pas comment faire mieux
        self.cohorts_replications = cohorts_replications  #: dictionary of number of replications per cohort rank
        self.index = index  #: the index of the phytomer TEMPORARY

        # state parameters
        self.mstruct = mstruct                      #: g
        self.Nstruct = Nstruct                      #: g
        self.ratio_DZ = ratio_DZ
        self.ratio_EOZ = ratio_EOZ

        # state variables
        self.sucrose = sucrose                      #: :math:`\mu mol` C
        self.fructan = fructan                      #: :math:`\mu mol` C
        self.amino_acids = amino_acids              #: :math:`\mu mol` N
        self.proteins = proteins                    #: :math:`\mu mol` N

        # fluxes from phloem
        self.Unloading_Sucrose = None               #: current Unloading of sucrose from phloem to hiddenzone integrated over delta t (:math:`\mu mol` C)
        self.Unloading_Amino_Acids = None           #: current Unloading of amino acids from phloem to hiddenzone integrated over delta t (:math:`\mu mol` N)

        # other fluxes
        self.S_Proteins = None                      #: protein synthesis (:math:`\mu mol` N g-1 mstruct)
        self.S_Fructan = None                       #: fructan synthesis (:math:`\mu mol` C g-1 mstruct)
        self.D_Fructan = None                       #: fructan degradation (:math:`\mu mol` C g-1 mstruct)
        self.D_Proteins = None                      #: protein degradation (:math:`\mu mol` N g-1 mstruct)

        # intermediate variables
        self.R_residual = None  #: Residual maintenance respiration (cost from protein turn-over, cell ion gradients, futile cycles...) (:math:`\mu mol` C respired)

        # Integrated variables
        self.Total_Organic_Nitrogen = None          #: current total nitrogen amount (:math:`\mu mol` N)

    @property
    def nb_replications(self):
        return sum(int(v <= self.index) * self.cohorts_replications.get(v, 0) for v in self.cohorts) + 1

    def calculate_integrative_variables(self):
        self.Total_Organic_Nitrogen = self.calculate_Total_Organic_Nitrogen(self.amino_acids, self.proteins, self.Nstruct)

    # VARIABLES
    @staticmethod
    def calculate_Total_Organic_Nitrogen(amino_acids, proteins, Nstruct):
        """Total amount of organic N (amino acids + proteins + Nstruct).
        Used to calculate residual respiration.

        :param float amino_acids: Amount of amino acids (:math:`\mu mol` N)
        :param float proteins: Amount of proteins (:math:`\mu mol` N)
        :param float Nstruct: Structural N mass (g)

        :return: Total amount of organic N (:math:`\mu mol` N)
        :rtype: float
        """
        return amino_acids + proteins + (Nstruct / EcophysiologicalConstants.N_MOLAR_MASS)*1E6

    # FLUXES

    def calculate_Unloading_Sucrose(self, sucrose, sucrose_phloem, mstruct_axis, T_effect_conductivity):
        """Rate of sucrose Unloading from phloem to the hidden zone (:math:`\mu mol` C sucrose unloaded h-1).
        Transport-resistance equation

        :param float sucrose: Sucrose amount in the hidden zone (:math:`\mu mol` C)
        :param float sucrose_phloem: Sucrose amount in phloem (:math:`\mu mol` C)
        :param float mstruct_axis: The structural dry mass of the axis (g)
        :param float T_effect_conductivity: Effect of the temperature on the conductivity rate at 20�C (AU)

        :return: Rate of Sucrose Unloading (:math:`\mu mol` C h-1)
        :rtype: float
        """
        conc_sucrose_phloem = (sucrose_phloem / mstruct_axis)
        conc_sucrose_HZ = (sucrose / self.mstruct)
        conductance = HiddenZone.PARAMETERS.SIGMA * HiddenZone.PARAMETERS.BETA * self.mstruct**(2/3) * T_effect_conductivity  # TODO: choix valeurs paramq par rapport flux phloem-hgz

        return (conc_sucrose_phloem - conc_sucrose_HZ) * conductance * parameters.SECOND_TO_HOUR_RATE_CONVERSION

    def calculate_Unloading_Amino_Acids(self, amino_acids, amino_acids_phloem, mstruct_axis, T_effect_conductivity):
        """Rate of amino acids Unloading from phloem to the hidden zone (:math:`\mu mol` N amino acids unloaded h-1).
        Transport-resistance equation

        :param float amino_acids: Amino_acids amount in the hidden zone (:math:`\mu mol` N)
        :param float amino_acids_phloem: Amino_acids amount in phloem (:math:`\mu mol` N)
        :param float mstruct_axis: The structural dry mass of the axis (g)
        :param float T_effect_conductivity:  Effect of the temperature on the conductivity rate at 20�C (AU)

        :return: Rate of Amino_acids Unloading (:math:`\mu mol` N h-1)
        :rtype: float
        """
        conc_amino_acids_phloem = (amino_acids_phloem / mstruct_axis)
        conc_amino_acids_HZ = (amino_acids / self.mstruct)
        conductance = HiddenZone.PARAMETERS.SIGMA * HiddenZone.PARAMETERS.BETA * self.mstruct**(2/3) * T_effect_conductivity
        return (conc_amino_acids_phloem - conc_amino_acids_HZ) * conductance * parameters.SECOND_TO_HOUR_RATE_CONVERSION

    def calculate_S_proteins(self, amino_acids, T_effect_Vmax):
        """Rate of protein synthesis (:math:`\mu mol` N proteins h-1 g-1 MS).
        Michaelis-Menten function of amino acids.

        :param float amino_acids: Amino acid amount in the hidden zone (:math:`\mu mol` N)
        :param float T_effect_Vmax: Correction to apply to enzyme activity

        :return: Rate of Protein synthesis (:math:`\mu mol` N g-1 mstruct h-1)
        :rtype: float
        """
        vmax = HiddenZone.PARAMETERS.VMAX_SPROTEINS_EMZ * (1-self.ratio_DZ) + HiddenZone.PARAMETERS.VMAX_SPROTEINS_DZ * self.ratio_DZ  #: 'Mean' Vmax for the whole hidden zone
        return ((vmax * max(0, (amino_acids / self.mstruct))) / (HiddenZone.PARAMETERS.K_SPROTEINS + max(0, (amino_acids / self.mstruct)))) * parameters.SECOND_TO_HOUR_RATE_CONVERSION * T_effect_Vmax

    def calculate_D_Proteins(self, proteins, T_effect_Vmax):
        """Rate of protein degradation (:math:`\mu mol` N proteins h-1 g-1 MS).
        First order kinetic

        :param float proteins: Protein amount in the hidden zone (:math:`\mu mol` N)
        :param float T_effect_Vmax: Correction to apply to enzyme activity

        :return: Rate of Protein degradation (:math:`\mu mol` N g-1 mstruct h-1)
        :rtype: float
        """
        return max(0, (HiddenZone.PARAMETERS.delta_Dproteins * (proteins / self.mstruct))) * parameters.SECOND_TO_HOUR_RATE_CONVERSION * T_effect_Vmax

    def calculate_Regul_S_Fructan(self, Unloading_Sucrose):
        """Regulating function for fructan maximal rate of synthesis.
        Negative regulation by the loading of sucrose from the phloem ("swith-off" sigmo�dal kinetic).

        :param float Unloading_Sucrose: Sucrose unloading (:math:`\mu mol` C)

        :return: Maximal rate of fructan synthesis (:math:`\mu mol` C g-1 mstruct)
        :rtype: float
        """

        if Unloading_Sucrose >= 0:
            Vmax_Sfructans = HiddenZone.PARAMETERS.VMAX_SFRUCTAN_POT
        else:  # Regulation by sucrose unloading if hidden zone is a source for C
            rate_Loading_Sucrose_massic = -Unloading_Sucrose / self.mstruct / parameters.SECOND_TO_HOUR_RATE_CONVERSION
            Vmax_Sfructans = HiddenZone.PARAMETERS.VMAX_SFRUCTAN_POT * (HiddenZone.PARAMETERS.K_REGUL_SFRUCTAN**HiddenZone.PARAMETERS.N_REGUL_SFRUCTAN / 
                                                                        (max(0., rate_Loading_Sucrose_massic**HiddenZone.PARAMETERS.N_REGUL_SFRUCTAN) + 
                                                                         HiddenZone.PARAMETERS.K_REGUL_SFRUCTAN**HiddenZone.PARAMETERS.N_REGUL_SFRUCTAN))
        return Vmax_Sfructans

    def calculate_S_Fructan(self, sucrose, Regul_S_Fructan, T_effect_Vmax):
        """Rate of fructan synthesis (:math:`\mu mol` C fructan g-1 mstruct h-1).
        Sigmo�dal function of sucrose.

        :param float sucrose: Amount of sucrose (:math:`\mu mol` C)
        :param float Regul_S_Fructan: Maximal rate of fructan synthesis regulated by sucrose loading (:math:`\mu mol` C g-1 mstruct)
        :param float T_effect_Vmax: Correction to apply to enzyme activity

        :return: Rate of Fructan synthesis (:math:`\mu mol` C g-1 mstruct)
        :rtype: float
        """
        return ((max(0., sucrose) / self.mstruct) * HiddenZone.PARAMETERS.VMAX_SFRUCTAN_RELATIVE * Regul_S_Fructan) / \
               ((max(0., sucrose) / self.mstruct) + HiddenZone.PARAMETERS.K_SFRUCTAN) * parameters.SECOND_TO_HOUR_RATE_CONVERSION * T_effect_Vmax

    def calculate_D_Fructan(self, sucrose, fructan, T_effect_Vmax):
        """Rate of fructan degradation (:math:`\mu mol` C fructan g-1 mstruct h-1).
        Inhibition function by the end product i.e. sucrose (Bancal et al., 2012).

        :param float sucrose: Amount of sucrose (:math:`\mu mol` C)
        :param float fructan: Amount of fructan (:math:`\mu mol` C)
        :param float T_effect_Vmax: Correction to apply to enzyme activity

        :return: Rate of Fructan degradation (:math:`\mu mol` C g-1 mstruct)
        :rtype: float
        """
        d_potential = ((HiddenZone.PARAMETERS.K_DFRUCTAN * HiddenZone.PARAMETERS.VMAX_DFRUCTAN * T_effect_Vmax) /
                       ((max(0., sucrose) / self.mstruct) + HiddenZone.PARAMETERS.K_DFRUCTAN)) * parameters.SECOND_TO_HOUR_RATE_CONVERSION
        d_actual = min(d_potential, max(0., fructan))
        return d_actual

    # COMPARTMENTS

    def calculate_sucrose_derivative(self, Unloading_Sucrose, S_Fructan, D_Fructan, hiddenzone_Loading_Sucrose_contribution, R_residual):
        """delta sucrose of hidden zone.

        :param float Unloading_Sucrose: Sucrose unloaded (:math:`\mu mol` C)
        :param float S_Fructan: Fructan synthesis (:math:`\mu mol` C g-1 mstruct)
        :param float D_Fructan: Fructan degradation (:math:`\mu mol` C g-1 mstruct)
        :param float hiddenzone_Loading_Sucrose_contribution: Sucrose imported from the emerged tissues (:math:`\mu mol` C)
        :param float R_residual: Residual respiration (:math:`\mu mol` C)

        :return: delta sucrose (:math:`\mu mol` C sucrose)
        :rtype: float
        """
        return Unloading_Sucrose + (D_Fructan - S_Fructan) * self.mstruct + hiddenzone_Loading_Sucrose_contribution - R_residual

    def calculate_amino_acids_derivative(self, Unloading_Amino_Acids, S_Proteins, D_Proteins, hiddenzone_Loading_Amino_Acids_contribution):
        """delta amino acids of hidden zone.

        :param float Unloading_Amino_Acids: Amino acids unloaded (:math:`\mu mol` N)
        :param float S_Proteins: Protein synthesis (:math:`\mu mol` N g-1 mstruct)
        :param float D_Proteins: Protein degradation (:math:`\mu mol` N g-1 mstruct)
        :param float hiddenzone_Loading_Amino_Acids_contribution: Amino acids imported from the emerged tissues (:math:`\mu mol` N)

        :return: delta amino acids (:math:`\mu mol` N amino acids)
        :rtype: float
        """
        return Unloading_Amino_Acids + (D_Proteins - S_Proteins) * self.mstruct + hiddenzone_Loading_Amino_Acids_contribution

    def calculate_fructan_derivative(self, S_Fructan, D_Fructan):
        """delta fructans of hidden zone.

        :param float S_Fructan: Fructans synthesis (:math:`\mu mol` C g-1 mstruct)
        :param float D_Fructan: Fructans degradation (:math:`\mu mol` C g-1 mstruct)

        :return: delta fructans (:math:`\mu mol` C fructans)
        :rtype: float
        """
        return (S_Fructan - D_Fructan) * self.mstruct

    def calculate_proteins_derivative(self, S_Proteins, D_Proteins):
        """delta proteins of hidden zone.

        :param float S_Proteins: Protein synthesis (:math:`\mu mol` N g-1 mstruct)
        :param float D_Proteins: Protein degradation (:math:`\mu mol` N g-1 mstruct)

        :return: delta proteins (:math:`\mu mol` N proteins)
        :rtype: float
        """
        return (S_Proteins - D_Proteins) * self.mstruct


class Phloem(Organ):
    """
    The class :class:`Phloem` defines the CN exchanges in a phloem.
    """

    PARAMETERS = parameters.PHLOEM_PARAMETERS                #: the internal parameters of the phloem
    INIT_COMPARTMENTS = parameters.PHLOEM_INIT_COMPARTMENTS  #: the initial values of compartments and state parameters

    def __init__(self, label='phloem', sucrose=INIT_COMPARTMENTS.sucrose, amino_acids=INIT_COMPARTMENTS.amino_acids):

        super(Phloem, self).__init__(label)

        # state variables
        self.sucrose = sucrose          #: :math:`\mu mol` C sucrose
        self.amino_acids = amino_acids  #: :math:`\mu mol` N amino acids

    # COMPARTMENTS

    @staticmethod
    def calculate_sucrose_derivative(contributors):
        """delta sucrose

        :param list [PhotosyntheticOrganElement, Grains, Roots, HiddenZone] contributors: Organs exchanging C with the phloem

        :return: delta sucrose (:math:`\mu mol` C sucrose)
        :rtype: float
        """
        sucrose_derivative = 0
        for contributor in contributors:
            if isinstance(contributor, PhotosyntheticOrganElement):
                sucrose_derivative += contributor.Loading_Sucrose * contributor.nb_replications
            elif isinstance(contributor, Grains):
                sucrose_derivative -= contributor.S_grain_structure + (contributor.S_grain_starch * contributor.structural_dry_mass)
            elif isinstance(contributor, Roots):
                sucrose_derivative -= contributor.Unloading_Sucrose * contributor.mstruct * contributor.__class__.PARAMETERS.ALPHA
            elif isinstance(contributor, HiddenZone):
                sucrose_derivative -= contributor.Unloading_Sucrose * contributor.nb_replications

        return sucrose_derivative

    @staticmethod
    def calculate_amino_acids_derivative(contributors):
        """delta amino acids

        :param list [PhotosyntheticOrganElement, Grains, Roots, HiddenZone] contributors: Organs exchanging N with the phloem

        :return: delta amino acids (:math:`\mu mol` N amino acids)
        :rtype: float
        """
        amino_acids_derivative = 0
        for contributor in contributors:
            if isinstance(contributor, PhotosyntheticOrganElement):
                amino_acids_derivative += contributor.Loading_Amino_Acids * contributor.nb_replications
            elif isinstance(contributor, Grains):
                amino_acids_derivative -= contributor.S_Proteins
            elif isinstance(contributor, Roots):
                amino_acids_derivative -= contributor.Unloading_Amino_Acids * contributor.mstruct * contributor.__class__.PARAMETERS.ALPHA
            elif isinstance(contributor, HiddenZone):
                amino_acids_derivative -= contributor.Unloading_Amino_Acids * contributor.nb_replications

        return amino_acids_derivative


class Grains(Organ):
    """
    The class :class:`Grains` defines the CN exchanges in a set of grains.
    """

    AMINO_ACIDS_MOLAR_MASS_N_RATIO = 0.136                   #: Mean contribution of N in amino acids mass contained in gluten (Glu, Gln and Pro)

    PARAMETERS = parameters.GRAINS_PARAMETERS                #: the internal parameters of the grains
    INIT_COMPARTMENTS = parameters.GRAINS_INIT_COMPARTMENTS  #: the initial values of compartments and state parameters

    def __init__(self, label='grains', age_from_flowering=INIT_COMPARTMENTS.age_from_flowering, starch=INIT_COMPARTMENTS.starch, structure=INIT_COMPARTMENTS.structure,
                 proteins=INIT_COMPARTMENTS.proteins):

        super(Grains, self).__init__(label)

        # state variables
        self.age_from_flowering = age_from_flowering    #: seconds
        self.starch = starch                            #: :math:`\mu mol` of C starch
        self.structure = structure                      #: :math:`\mu mol` of C sucrose
        self.proteins = proteins                        #: :math:`\mu mol` of N proteins

        # derived attributes
        self.structural_dry_mass = None                 #: g of MS

        # fluxes from phloem
        self.S_grain_structure = None                   #: current synthesis of grain structure integrated over a delta t (:math:`\mu mol` C)
        self.S_grain_starch = None                      #: current synthesis of grain starch integrated over a delta t (:math:`\mu mol` C g-1 mstruct)
        self.S_Proteins = None                          #: current synthesis of grain proteins integrated over a delta t (:math:`\mu mol` N)

        # intermediate variables
        self.RGR_Structure = None                       #: RGR of grain structure (dimensionless?)
        self.R_grain_growth_struct = None               #: grain struct  respiration (:math:`\mu mol` C respired)
        self.R_grain_growth_starch = None               #: grain starch growth respiration (:math:`\mu mol` C respired)

    def initialize(self):
        """Initialize the derived attributes of the organ.
        """
        self.structural_dry_mass = self.calculate_structural_dry_mass(self.structure)

    # VARIABLES
    @staticmethod
    def calculate_structural_dry_mass(structure):
        """Grain structural dry mass.

        :param float structure: Grain structural C mass (:math:`\mu mol` C)

        :return: Grain structural dry mass (g)
        :rtype: float
        """
        return (structure*1E-6 * EcophysiologicalConstants.C_MOLAR_MASS) / EcophysiologicalConstants.RATIO_C_mstruct

    @staticmethod
    def modified_Arrhenius_equation(temperature):  # TODO: move in a seperate model
        """ Return value of equation from Johnson and Lewin (1946) for temperature. The equation is modified to return zero below zero degree.

        :param float temperature: organ temperature (degree Celsius)

        :return: Return value of Eyring equation from Johnson and Lewin (1946) for temperature (dimensionless). The equation is modified to return zero below zero degree.
        :rtype: float
        """

        # Parameters for temperature responses
        Temp_Ea_R = 8900  # Parameter Ea/R in Eyring equation from Johnson and Lewin (1946) - Parameter value fitted from Kemp and Blacklow (1982) (K)
        Temp_DS_R = 68.432  # Parameter deltaS/R in Eyring equation from Johnson and Lewin (1946) - Parameter value fitted from Kemp and Blacklow (1982) (dimensionless)
        Temp_DH_R = 20735.5  # Parameter deltaH/R in Eyring equation from Johnson and Lewin (1946) - Parameter value fitted from Kemp and Blacklow (1982) (K)
        Temp_Ttransition = 9  # Below this temperature f = linear function of temperature instead of Arrhenius-like(�C)

        def Arrhenius_equation(T):
            return T * exp(-Temp_Ea_R / T) / (1 + exp(Temp_DS_R - Temp_DH_R / T))

        temperature_K = temperature + 273.15  #: Kelvins

        if temperature < 0:
            res = 0
        elif temperature < Temp_Ttransition:
            res = temperature * Arrhenius_equation(Temp_Ttransition + 273.15) / Temp_Ttransition
        else:
            res = Arrhenius_equation(temperature_K)

        return res

    def calculate_temperature_effect_on_growth(self, Tair):
        """Effect of the temperature on elongation.
        Return value of equation from Johnson and Lewin (1946) for temperature. The equation is modified to return zero below zero degree.
        Identical to modified_Arrhenius_equation in ElongWheat.
        Should multiply the rate at 20�C

        :param float Tair: Air temperature(�C)

        :return: Correction to apply to RGR Structure of the grains (dimensionless)
        :rtype: float
        """
        return self.modified_Arrhenius_equation(Tair) / Grains.PARAMETERS.Arrhenius_ref

    @staticmethod
    def calculate_RGR_Structure(sucrose_phloem, mstruct_axis, T_effect_growth):
        """Relative Growth Rate of grain structure, regulated by sucrose concentration in phloem.

        :param float sucrose_phloem: Sucrose amount in phloem (:math:`\mu mol` C)
        :param float mstruct_axis: The structural dry mass of the axis (g)
        :param float T_effect_growth: Effect of the temperature on the growth rate at 20�C (AU)

        :return: RGR of grain structure (dimensionless?)
        :rtype: float
        """
        return ((max(0., sucrose_phloem) / (mstruct_axis * Axis.PARAMETERS.ALPHA)) * Grains.PARAMETERS.VMAX_RGR) / ((max(0., sucrose_phloem) / (mstruct_axis * Axis.PARAMETERS.ALPHA)) +
                                                                                                                    Grains.PARAMETERS.K_RGR) * T_effect_growth

    # FLUXES

    def calculate_S_grain_structure(self, prec_structure, RGR_Structure):
        """Rate of grain structure synthesis (:math:`\mu mol` C structure h-1).
        Exponential function, RGR regulated by sucrose concentration in the phloem.

        :param float prec_structure: Grain structure at t-1 (:math:`\mu mol` C)
        :param float RGR_Structure: Relative Growth Rate of grain structure (dimensionless?)

        :return: Rate of Synthesis of grain structure (:math:`\mu mol` C h-1)
        :rtype: float
        """
        if self.age_from_flowering <= Grains.PARAMETERS.FILLING_INIT:  #: Grain enlargment
            S_grain_structure = prec_structure * RGR_Structure * parameters.SECOND_TO_HOUR_RATE_CONVERSION
        else:                                                          #: Grain filling
            S_grain_structure = 0
        return S_grain_structure

    def calculate_S_grain_starch(self, sucrose_phloem, mstruct_axis, T_effect_Vmax):
        """Rate of starch synthesis in grains (i.e. grain filling) (:math:`\mu mol` C starch g-1 mstruct h-1).
        Michaelis-Menten function of sucrose concentration in the phloem.

        :param float sucrose_phloem: Sucrose concentration in phloem (:math:`\mu mol` C g-1 mstruct)
        :param float mstruct_axis: The structural dry mass of the axis (g)
        :param float T_effect_Vmax: Correction to apply to enzyme activity


        :return: Rate of Synthesis of grain starch (:math:`\mu mol` C g-1 mstruct h-1)
        :rtype: float
        """
        if self.age_from_flowering <= Grains.PARAMETERS.FILLING_INIT:   #: Grain enlargment
            S_grain_starch = 0
        elif self.age_from_flowering > Grains.PARAMETERS.FILLING_END:   #: Grain maturity
            S_grain_starch = 0
        else:                                                           #: Grain filling
            S_grain_starch = (((max(0., sucrose_phloem)/(mstruct_axis * Axis.PARAMETERS.ALPHA)) * Grains.PARAMETERS.VMAX_STARCH) /
                              ((max(0., sucrose_phloem)/(mstruct_axis * Axis.PARAMETERS.ALPHA)) + Grains.PARAMETERS.K_STARCH)) * parameters.SECOND_TO_HOUR_RATE_CONVERSION * T_effect_Vmax
        return S_grain_starch

    @staticmethod
    def calculate_S_proteins(S_grain_structure, S_grain_starch, amino_acids_phloem, sucrose_phloem, structural_dry_mass):
        """Protein synthesis in grains.
        N is assumed to be co-transported along with the unloaded sucrose from phloem (using the ratio amino acids:sucrose of phloem).

        :param float S_grain_structure: Synthesis of grain structure (:math:`\mu mol` C)
        :param float S_grain_starch: Synthesis of grain starch (:math:`\mu mol` C g-1 mstruct)
        :param float amino_acids_phloem: Amino acids concentration in phloem (:math:`\mu mol` N g-1 mstruct)
        :param float sucrose_phloem: Sucrose concentration in phloem (:math:`\mu mol` C g-1 mstruct)
        :param float structural_dry_mass: Grain structural dry mass (g)

        :return: Synthesis of grain proteins (:math:`\mu mol` N)
        :rtype: float
        """
        if sucrose_phloem > 0:
            S_Proteins = (S_grain_structure + S_grain_starch*structural_dry_mass) * (amino_acids_phloem / sucrose_phloem)
        else:
            S_Proteins = 0
        return S_Proteins

    # COMPARTMENTS
    @staticmethod
    def calculate_structure_derivative(S_grain_structure, R_growth):
        """delta grain structure.

        :param float S_grain_structure: Synthesis of grain structure (:math:`\mu mol` C)
        :param float R_growth: Grain growth respiration (:math:`\mu mol` C respired)

        :return: delta grain structure (:math:`\mu mol` C structure)
        :rtype: float
        """
        return S_grain_structure - R_growth

    @staticmethod
    def calculate_starch_derivative(S_grain_starch, structural_dry_mass, R_growth):
        """delta grain starch.

        :param float S_grain_starch: Synthesis of grain starch (:math:`\mu mol` C g-1 mstruct)
        :param float structural_dry_mass: Grain structural dry mass (g)
        :param float R_growth: Grain growth respiration (:math:`\mu mol` C respired)

        :return: delta grain starch (:math:`\mu mol` C starch)
        :rtype: float
        """
        return (S_grain_starch * structural_dry_mass) - R_growth

    @staticmethod
    def calculate_proteins_derivative(S_Proteins):
        """delta grain proteins.

        :param float S_Proteins: Synthesis of grain proteins (:math:`\mu mol` N)

        :return: delta grain proteins (:math:`\mu mol` N proteins)
        :rtype: float
        """
        return S_Proteins


class Roots(Organ):
    """
    The class :class:`Roots` defines the CN exchanges in a set of roots.
    """

    PARAMETERS = parameters.ROOTS_PARAMETERS                #: the internal parameters of the roots
    INIT_COMPARTMENTS = parameters.ROOTS_INIT_COMPARTMENTS  #: the initial values of compartments and state parameters

    def __init__(self, label='roots', mstruct=INIT_COMPARTMENTS.mstruct, senesced_mstruct=INIT_COMPARTMENTS.senesced_mstruct, Nstruct=INIT_COMPARTMENTS.Nstruct, sucrose=INIT_COMPARTMENTS.sucrose,
                 nitrates=INIT_COMPARTMENTS.nitrates,
                 amino_acids=INIT_COMPARTMENTS.amino_acids, cytokinins=INIT_COMPARTMENTS.cytokinins):

        super(Roots, self).__init__(label)

        # state parameters
        self.mstruct = mstruct                 #: Structural mass (g)
        self.senesced_mstruct = senesced_mstruct  #: Senesced structural mass (g)
        self.Nstruct = Nstruct                 #: Structural N mass (g)

        # state variables
        self.sucrose = sucrose                 #: :math:`\mu mol` C sucrose
        self.nitrates = nitrates               #: :math:`\mu mol` N nitrates
        self.amino_acids = amino_acids         #: :math:`\mu mol` N amino acids
        self.cytokinins = cytokinins           #: AU cytokinins

        # fluxes from phloem
        self.Unloading_Sucrose = None          #: current Unloading of sucrose from phloem to roots
        self.Unloading_Amino_Acids = None      #: current Unloading of amino acids from phloem to roots

        # other fluxes
        self.Export_Nitrates = None            #: Total export of nitrates from roots to shoot organs integrated over a delta t (:math:`\mu mol` N)
        self.Export_Amino_Acids = None         #: Total export of amino acids from roots to shoot organs integrated over a delta t (:math:`\mu mol` N)
        self.S_Amino_Acids = None              #: Rate of amino acid synthesis in roots integrated over a delta t (:math:`\mu mol` N g-1 mstruct)
        self.Uptake_Nitrates = None            #: Rate of nitrate uptake by roots integrated over a delta t (:math:`\mu mol` N nitrates)
        self.S_cytokinins = None               #: Rate of cytokinin synthesis integrated over a delta t (AU g-1 mstruct)
        self.Export_cytokinins = None          #: Total export of cytokinin from roots to shoot organs integrated over a delta t (AU)

        # Integrated variables
        self.Total_Organic_Nitrogen = None     #: current amount of organic N (:math:`\mu mol` N)

        # intermediate variables
        self.R_Nnit_upt = None                 #: Nitrate uptake respiration (:math:`\mu mol` C respired)
        self.R_Nnit_red = None                 #: Nitrate reduction-linked respiration (:math:`\mu mol` C respired)
        self.R_residual = None                 #: Residual maintenance respiration (cost from protein turn-over, cell ion gradients, futile cycles...) (:math:`\mu mol` C respired)
        self.C_exudation = None                #: C sucrose lost by root exudation integrated over a delta t (:math:`\mu mol` C g-1 mstruct)
        self.N_exudation = None                #: N amino acids lost by root exudation integrated over a delta t (:math:`\mu mol` N g-1 mstruct)
        self.regul_transpiration = None        #: Dimensionless regulating factor of metabolite exports from roots by shoot transpiration
        self.HATS_LATS = None                  #: Nitrate influx (:math:`\mu mol` N)
        self.sum_respi = None                  #: Sum of respirations for roots i.e. related to N uptake, amino acids synthesis and residual (:math:`\mu mol` C)

    def calculate_integrative_variables(self):
        self.Total_Organic_Nitrogen = self.calculate_Total_Organic_Nitrogen(self.amino_acids, self.Nstruct)

    # VARIABLES

    @staticmethod
    def calculate_Total_Organic_Nitrogen(amino_acids, Nstruct):
        """Total amount of organic N (amino acids + Nstruct).
        Used to calculate residual respiration.

        :param float amino_acids: Amount of amino acids (:math:`\mu mol` N)
        :param float Nstruct: Structural N mass (g)

        :return: Total amount of organic N (:math:`\mu mol` N)
        :rtype: float
        """
        return amino_acids + (Nstruct / EcophysiologicalConstants.N_MOLAR_MASS)*1E6

    @staticmethod
    def calculate_regul_transpiration(total_transpiration):
        """A function to regulate metabolite exports from roots by shoot transpiration

        :param float total_transpiration: Total transpiration (mmol s-1)

        :return:
            Dimensionless regulating factor
        :rtype: float
        """
        return total_transpiration * Roots.PARAMETERS.CST_TRANSPIRATION

    # FLUXES

    def calculate_Unloading_Sucrose(self, sucrose_roots, sucrose_phloem, mstruct_axis, T_effect_conductivity):
        """Rate of sucrose Unloading from phloem to roots (:math:`\mu mol` C sucrose unloaded g-1 mstruct h-1).
        Michaelis-Menten function of the sucrose concentration in phloem.

        :param float sucrose_roots: Amount of sucrose in roots (:math:`\mu mol` C)
        :param float sucrose_phloem: Sucrose concentration in phloem (:math:`\mu mol` C g-1 mstruct)
        :param float mstruct_axis: The structural dry mass of the axis (g)
        :param float T_effect_conductivity: Effect of the temperature on the conductivity rate at 20�C (AU)

        :return: Rate of Sucrose Unloading (:math:`\mu mol` C g-1 mstruct h-1)
        :rtype: float
        """
        conc_sucrose_roots = sucrose_roots / (self.mstruct * self.__class__.PARAMETERS.ALPHA)
        conc_sucrose_phloem = sucrose_phloem / (mstruct_axis * parameters.AXIS_PARAMETERS.ALPHA)
        #: Driving compartment (:math:`\mu mol` C g-1 mstruct)
        driving_sucrose_compartment = max(conc_sucrose_roots, conc_sucrose_phloem)
        #: Gradient of sucrose between the roots and the phloem (:math:`\mu mol` C g-1 mstruct)
        diff_sucrose = conc_sucrose_phloem - conc_sucrose_roots
        #: Conductance depending on mstruct (g2 :math:`\mu mol`-1 s-1)
        conductance = Roots.PARAMETERS.SIGMA_SUCROSE * Roots.PARAMETERS.BETA * self.mstruct ** (2/3) * T_effect_conductivity

        return driving_sucrose_compartment * diff_sucrose * conductance * parameters.SECOND_TO_HOUR_RATE_CONVERSION

    @staticmethod
    def calculate_Unloading_Amino_Acids(Unloading_Sucrose, sucrose_phloem, amino_acids_phloem):
        """Unloading of amino_acids from phloem to roots.
        Amino acids are assumed to be co-transported along with the unloaded sucrose from phloem (using the ratio amino acids:sucrose of phloem).

        :param float Unloading_Sucrose: Sucrose Unloading (:math:`\mu mol` C g-1 mstruct)
        :param float sucrose_phloem: Sucrose concentration in phloem (:math:`\mu mol` C g-1 mstruct)
        :param float amino_acids_phloem: Amino acids concentration in phloem (:math:`\mu mol` N g-1 mstruct)

        :return: Amino acids Unloading (:math:`\mu mol` N g-1 mstruct)
        :rtype: float
        """
        if amino_acids_phloem <= 0 or sucrose_phloem <= 0 or Unloading_Sucrose <= 0:
            Unloading_Amino_Acids = 0
        else:
            Unloading_Amino_Acids = Unloading_Sucrose * (amino_acids_phloem / sucrose_phloem)
        return Unloading_Amino_Acids

    def calculate_Uptake_Nitrates(self, Conc_Nitrates_Soil, nitrates_roots, sucrose_roots, T_effect_Vmax):
        """Rate of nitrate uptake by roots
            - Nitrate uptake is calculated as the sum of the 2 transport systems: HATS and LATS
            - HATS and LATS parameters are calculated as a function of root nitrate concentration (negative regulation)
            - Nitrate uptake is finally regulated by the total culm transpiration and sucrose concentration (positive regulation)

        :param float Conc_Nitrates_Soil: Soil nitrate concentration Unloading (:math:`\mu mol` N m-3 soil)
        :param float nitrates_roots: Amount of nitrates in roots (:math:`\mu mol` N)
        :param float sucrose_roots: Amount of sucrose in roots (:math:`\mu mol` C)
        :param float T_effect_Vmax: Correction to apply to enzyme activity


        :return: Nitrate uptake (:math:`\mu mol` N nitrates) and nitrate influxes HATS and LATS (:math:`\mu mol` N h-1)
        :rtype: (float, float)
        """
        conc_nitrates_roots = nitrates_roots / self.mstruct

        #: High Affinity Transport System (HATS)
        VMAX_HATS_MAX = max( 0., Roots.PARAMETERS.A_VMAX_HATS * conc_nitrates_roots + Roots.PARAMETERS.B_VMAX_HATS)  #: Maximal rate of nitrates influx at saturating soil N concentration;HATS (:math:`\mu   mol` N nitrates g-1 mstruct s-1)
        K_HATS =  max( 0., Roots.PARAMETERS.A_K_HATS * conc_nitrates_roots + Roots.PARAMETERS.B_K_HATS)   #: Affinity coefficient of nitrates influx at saturating soil N concentration;HATS (:math:`\mu mol` m-3)
        HATS = (VMAX_HATS_MAX * Conc_Nitrates_Soil) / (K_HATS + Conc_Nitrates_Soil)                         #: Rate of nitrate influx by HATS (:math:`\mu mol` N nitrates uptaked s-1 g-1 mstruct)

        #: Low Affinity Transport System (LATS)
        K_LATS = max(0., Roots.PARAMETERS.A_LATS * conc_nitrates_roots + Roots.PARAMETERS.B_LATS)  #: Rate constant for nitrates influx at low soil N concentration; LATS (m3 g-1 mstruct s-1)
        LATS = (K_LATS * Conc_Nitrates_Soil)                                                            #: Rate of nitrate influx by LATS (:math:`\mu mol` N nitrates g-1 mstruct)

        #: Nitrate influx (:math:`\mu mol` N)
        HATS_LATS = (HATS + LATS)
        nitrate_influx = HATS_LATS * parameters.SECOND_TO_HOUR_RATE_CONVERSION * T_effect_Vmax * self.mstruct

        # Regulations
        regul_C = (sucrose_roots/self.mstruct) / ((sucrose_roots/self.mstruct) + Roots.PARAMETERS.K_C)  #: Nitrate uptake regulation by root C
        if HATS_LATS < Roots.PARAMETERS.MIN_INFLUX_FOR_UPTAKE:
            net_nitrate_uptake = 0
        else:
            net_nitrate_uptake = nitrate_influx * Roots.PARAMETERS.NET_INFLUX_UPTAKE_RATIO * regul_C             #: Net nitrate uptake (:math:`\mu mol` N nitrates uptaked by roots)
        return net_nitrate_uptake, nitrate_influx

    def calculate_S_amino_acids(self, nitrates, sucrose, T_effect_Vmax):
        """Rate of amino acid synthesis in roots (:math:`\mu mol` N amino acids g-1 mstruct h-1).
        Bi-substrate Michaelis-Menten function of nitrates and sucrose.

        :param float nitrates: Amount of nitrates in roots (:math:`\mu mol` N)
        :param float sucrose: Amount of sucrose in roots (:math:`\mu mol` C)
        :param float T_effect_Vmax: Correction to apply to enzyme activity


        :return: Amino acids synthesis (:math:`\mu mol` N g-1 mstruct h-1)
        :rtype: float
        """
        return T_effect_Vmax * Roots.PARAMETERS.VMAX_AMINO_ACIDS / ((1 + Roots.PARAMETERS.K_AMINO_ACIDS_NITRATES / (nitrates / (self.mstruct * Roots.PARAMETERS.ALPHA))) *
                                                                    (1 + Roots.PARAMETERS.K_AMINO_ACIDS_SUCROSE/(sucrose / (self.mstruct * Roots.PARAMETERS.ALPHA)))
                                                                    ) * parameters.SECOND_TO_HOUR_RATE_CONVERSION

    def calculate_Export_Nitrates(self, nitrates, regul_transpiration):
        """Total export of nitrates from roots to shoot organs
        Export is calculated as a function on nitrate concentration and culm transpiration.

        :param float nitrates: Amount of nitrates in roots (:math:`\mu mol` N)
        :param float regul_transpiration: Regulating factor by transpiration (mmol H2O m-2 s-1)

        :return: Rate of Export of nitrates (:math:`\mu mol` N h-1)
        :rtype: float
        """

        f_nitrates = (nitrates / (self.mstruct * Roots.PARAMETERS.ALPHA)) * Roots.PARAMETERS.K_NITRATE_EXPORT           #: :math:`\mu mol` g-1 s-1
        Export_Nitrates = f_nitrates * self.mstruct * regul_transpiration * parameters.SECOND_TO_HOUR_RATE_CONVERSION   #: Nitrate export regulation by transpiration (:math:`\mu mol` N)
        return max(min(Export_Nitrates, nitrates), 0.)

    def calculate_Export_Amino_Acids(self, amino_acids, regul_transpiration):
        """Total export of amino acids from roots to shoot organs
        Amino acids export is calculated as a function of nitrate export using the ratio amino acids:nitrates in roots.

        :param float amino_acids: Amount of amino acids in roots (:math:`\mu mol` N)
        :param float regul_transpiration: Regulating factor by transpiration (mmol H2O m-2 s-1)

        :Returns:
            Rate of Export of amino acids (:math:`\mu mol` N h-1)
        :Returns Type:
            :class:`float`
        """
        f_amino_acids = (amino_acids/(self.mstruct * Roots.PARAMETERS.ALPHA)) * Roots.PARAMETERS.K_AMINO_ACIDS_EXPORT        #: :math:`\mu mol` g-1 s-1
        Export_Amino_Acids = f_amino_acids * self.mstruct * regul_transpiration * parameters.SECOND_TO_HOUR_RATE_CONVERSION  #: Amino acids export regulation by plant transpiration (:math:`\mu mol` N)
        return max(min(Export_Amino_Acids, amino_acids), 0.)

    @staticmethod
    def calculate_exudation(Unloading_Sucrose, sucrose_roots, amino_acids_roots, amino_acids_phloem):
        """C sucrose and N amino acids lost by root exudation (:math:`\mu mol` C or N g-1 mstruct).
            - C exudation is calculated as a fraction of C Unloading from phloem
            - N exudation is calculated from C exudation using the ratio amino acids:sucrose of the phloem

        :param float Unloading_Sucrose: Sucrose Unloading (:math:`\mu mol` C g-1 mstruct h-1)
        :param float sucrose_roots: Amount of sucrose in roots (:math:`\mu mol` C)
        :param float amino_acids_roots: Amount of amino acids in roots (:math:`\mu mol` N)
        :param float amino_acids_phloem: Amount of amino acids in phloem (:math:`\mu mol` N)

        :return: Rates of C exudated (:math:`\mu mol` C g-1 mstruct h-1) and N_exudation (:math:`\mu mol` N g-1 mstruct h-1)
        :rtype: (float, float)
        """
        if sucrose_roots <=0 or Unloading_Sucrose <=0:
            C_exudation = 0
        else:
            C_exudation = min(sucrose_roots, Unloading_Sucrose * Roots.PARAMETERS.C_EXUDATION)  #: C exudated (:math:`\mu mol` g-1 mstruct)
        if amino_acids_phloem <= 0 or amino_acids_roots <= 0 or sucrose_roots <= 0:
            N_exudation = 0
        else:
            N_exudation = min((amino_acids_roots / sucrose_roots), Roots.PARAMETERS.N_EXUDATION_MAX) * C_exudation
        return C_exudation, N_exudation #TODO: C_exudation and N_exudation should be renamed as the exudation of AA result in a loss of both C and N

    def calculate_S_cytokinins(self, sucrose_roots, nitrates_roots, T_effect_Vmax):
        """ Rate of cytokinin synthesis (AU cytokinins g-1 mstruct h-1).
        Cytokinin synthesis regulated by both root sucrose and nitrates. As a signal molecule, cytokinins are assumed have a neglected effect on sucrose.
        Thus, no cost in C is applied to the sucrose pool.

        :param float sucrose_roots: Amount of sucrose in roots (:math:`\mu mol` C)
        :param float nitrates_roots: Amount of nitrates in roots (:math:`\mu mol` N)
        :param float T_effect_Vmax: Correction to apply to enzyme activity

        :return: Rate of Cytokinin synthesis (AU g-1 mstruct h-1)
        :rtype: float
        """
        conc_sucrose = max(0, (sucrose_roots/self.mstruct))
        conc_Nitrates = max(0, (nitrates_roots/self.mstruct))

        f_sucrose = conc_sucrose**Roots.PARAMETERS.N_SUC_CYTOKININS / (conc_sucrose**Roots.PARAMETERS.N_SUC_CYTOKININS + Roots.PARAMETERS.K_SUCROSE_CYTOKININS**Roots.PARAMETERS.N_SUC_CYTOKININS)
        f_nitrates = conc_Nitrates**Roots.PARAMETERS.N_NIT_CYTOKININS / (conc_Nitrates**Roots.PARAMETERS.N_NIT_CYTOKININS + Roots.PARAMETERS.K_NITRATES_CYTOKININS**Roots.PARAMETERS.N_NIT_CYTOKININS)

        S_cytokinins = Roots.PARAMETERS.VMAX_S_CYTOKININS * f_sucrose * f_nitrates * parameters.SECOND_TO_HOUR_RATE_CONVERSION * T_effect_Vmax
        return S_cytokinins

    def calculate_Export_cytokinins(self, cytokinins, regul_transpiration):
        """Total export of cytokinin from roots to shoot organs
        Cytokinin export is calculated as a function of cytokinin concentration and culm transpiration.

        :param float cytokinins: Amount of cytokinins in roots (AU)
        :param float regul_transpiration: Regulating factor by transpiration (mmol H2O m-2 s-1)

        :return: Rate of Cytokinin export (AU h-1)
        :rtype: float
        """
        f_cytokinins = (cytokinins / (self.mstruct * Roots.PARAMETERS.ALPHA)) * Roots.PARAMETERS.K_CYTOKININS_EXPORT       #: AU g-1 s-1
        Export_cytokinins = f_cytokinins * self.mstruct * regul_transpiration * parameters.SECOND_TO_HOUR_RATE_CONVERSION  #: Cytokinin export regulation by plant transpiration (AU)

        return max(min(Export_cytokinins, cytokinins), 0.)

    # COMPARTMENTS

    def calculate_sucrose_derivative(self, Unloading_Sucrose, S_Amino_Acids, C_exudation, sum_respi):
        """delta root sucrose.

        :param float Unloading_Sucrose: Sucrose Unloading (:math:`\mu mol` C g-1 mstruct)
        :param float S_Amino_Acids: Amino acids synthesis (:math:`\mu mol` N g-1 mstruct)
        :param float C_exudation: C exudation (:math:`\mu mol` C g-1 mstruct)
        :param float sum_respi: Sum of respirations for roots i.e. related to N uptake, amino acids synthesis and residual (:math:`\mu mol` C)

        :return: delta root sucrose (:math:`\mu mol` C sucrose)
        :rtype: float
        """
        #: Contribution of sucrose to the synthesis of amino_acids
        sucrose_consumption_AA = (S_Amino_Acids / EcophysiologicalConstants.AMINO_ACIDS_N_RATIO) * EcophysiologicalConstants.AMINO_ACIDS_C_RATIO
        return (Unloading_Sucrose - sucrose_consumption_AA - C_exudation) * self.mstruct - sum_respi

    def calculate_nitrates_derivative(self, Uptake_Nitrates, Export_Nitrates, S_Amino_Acids):
        """delta root nitrates.

        :param float Uptake_Nitrates: Nitrate uptake (:math:`\mu mol` N nitrates)
        :param float Export_Nitrates: Export of nitrates (:math:`\mu mol` N)
        :param float S_Amino_Acids: Amino acids synthesis (:math:`\mu mol` N g-1 mstruct)

        :return: delta root nitrates (:math:`\mu mol` N nitrates)
        :rtype: float
        """
        import_nitrates_roots = Uptake_Nitrates
        nitrate_reduction_AA = S_Amino_Acids  #: Contribution of nitrates to the synthesis of amino_acids
        return import_nitrates_roots - Export_Nitrates - nitrate_reduction_AA * self.mstruct

    def calculate_amino_acids_derivative(self, Unloading_Amino_Acids, S_Amino_Acids, Export_Amino_Acids, N_exudation):
        """delta root amino acids.

        :param float Unloading_Amino_Acids: Amino acids Unloading (:math:`\mu mol` N g-1 mstruct)
        :param float S_Amino_Acids: Amino acids synthesis (:math:`\mu mol` N g-1 mstruct)
        :param float Export_Amino_Acids: Export of amino acids (:math:`\mu mol` N)
        :param float N_exudation: N exudated (:math:`\mu mol` g-1 mstruct)

        :return: delta root amino acids (:math:`\mu mol` N amino acids)
        :rtype: float
        """
        return (Unloading_Amino_Acids + S_Amino_Acids - N_exudation) * self.mstruct - Export_Amino_Acids

    def calculate_cytokinins_derivative(self, S_cytokinins, Export_cytokinins):
        """delta root cytokinins.

        :param float S_cytokinins: Cytokinin synthesis (AU g-1 mstruct)
        :param float Export_cytokinins: Cytokinin export (AU)

        :return: delta root cytokinins (AU cytokinins)
        :rtype: float
        """
        return S_cytokinins * self.mstruct - Export_cytokinins


class PhotosyntheticOrgan(Organ):
    """
    The class :class:`PhotosyntheticOrgan` defines the CN exchanges in a photosynthetic organ.

    A :class:`photosynthetic organ <PhotosyntheticOrgan>` must have at least 1
    :class:`photosynthetic organ element <PhotosyntheticOrganElement>`:
    :class:`chaff element <ChaffElement>`, :class:`lamina element <LaminaElement>`,
    :class:`internode element <InternodeElement>`, :class:`peduncle element <PeduncleElement>`,
    or :class:`sheath element <SheathElement>`.

    :class:`PhotosyntheticOrgan` is the base class of all photosynthetic organs. DO NOT INSTANTIATE IT.
    """

    PARAMETERS = parameters.PHOTOSYNTHETIC_ORGAN_PARAMETERS  #: the internal parameters of the photosynthetic organs

    def __init__(self, label, exposed_element, enclosed_element):

        super(PhotosyntheticOrgan, self).__init__(label)

        self.exposed_element = exposed_element    #: the exposed element
        self.enclosed_element = enclosed_element  #: the enclosed element
        self.mstruct = None                       #: the structural dry mass
        self.senesced_mstruct = None              #: senesced structural dry mass

    def calculate_integrative_variables(self):
        self.mstruct = 0
        self.senesced_mstruct = 0
        for element in (self.exposed_element, self.enclosed_element):
            if element is not None:
                element.calculate_integrative_variables()
                self.mstruct += element.mstruct
                self.senesced_mstruct += element.senesced_mstruct


class Chaff(PhotosyntheticOrgan):
    """
    The class :class:`Chaff` defines the CN exchanges in a chaff.
    """

    PARAMETERS = parameters.CHAFF_PARAMETERS  #: the internal parameters of the chaffs

    def __init__(self, label=None, exposed_element=None, enclosed_element=None):
        super(Chaff, self).__init__(label, exposed_element, enclosed_element)


class Lamina(PhotosyntheticOrgan):
    """
    The class :class:`Lamina` defines the CN exchanges in a lamina.
    """

    PARAMETERS = parameters.LAMINA_PARAMETERS  #: the internal parameters of the laminae

    def __init__(self, label=None, exposed_element=None, enclosed_element=None):
        super(Lamina, self).__init__(label, exposed_element, enclosed_element)


class Internode(PhotosyntheticOrgan):
    """
    The class :class:`Internode` defines the CN exchanges in an internode.
    """

    PARAMETERS = parameters.INTERNODE_PARAMETERS  #: the internal parameters of the internodes

    def __init__(self, label=None, exposed_element=None, enclosed_element=None):
        super(Internode, self).__init__(label, exposed_element, enclosed_element)


class Peduncle(PhotosyntheticOrgan):
    """
    The class :class:`Peduncle` defines the CN exchanges in a peduncle.
    """

    PARAMETERS = parameters.PEDUNCLE_PARAMETERS  #: the internal parameters of the peduncles

    def __init__(self, label=None, exposed_element=None, enclosed_element=None):
        super(Peduncle, self).__init__(label, exposed_element, enclosed_element)


class Sheath(PhotosyntheticOrgan):
    """
    The class :class:`Sheath` defines the CN exchanges in a sheath.
    """

    PARAMETERS = parameters.SHEATH_PARAMETERS  #: the internal parameters of the sheaths

    def __init__(self, label=None, exposed_element=None, enclosed_element=None):
        super(Sheath, self).__init__(label, exposed_element, enclosed_element)


class PhotosyntheticOrganElement(object):
    """
    The class :class:`PhotosyntheticOrganElement` defines the CN exchanges in a photosynthetic organ element.

    An element must belong to an organ of the same type (e.g. a class:`LaminaElement` must belong to a class:`Lamina`).

    :class:`PhotosyntheticOrganElement` is the base class of all photosynthetic organs elements. DO NOT INSTANTIATE IT.
    """

    PARAMETERS = parameters.PHOTOSYNTHETIC_ORGAN_ELEMENT_PARAMETERS                #: the internal parameters of the photosynthetic organs elements
    INIT_COMPARTMENTS = parameters.PHOTOSYNTHETIC_ORGAN_ELEMENT_INIT_COMPARTMENTS  #: the initial values of compartments and state parameters

    def __init__(self, label=None, green_area=INIT_COMPARTMENTS.green_area, mstruct=INIT_COMPARTMENTS.mstruct, senesced_mstruct=INIT_COMPARTMENTS.senesced_mstruct, Nstruct=INIT_COMPARTMENTS.Nstruct,
                 triosesP=INIT_COMPARTMENTS.triosesP, starch=INIT_COMPARTMENTS.starch, sucrose=INIT_COMPARTMENTS.sucrose, fructan=INIT_COMPARTMENTS.fructan,
                 nitrates=INIT_COMPARTMENTS.nitrates, amino_acids=INIT_COMPARTMENTS.amino_acids, proteins=INIT_COMPARTMENTS.proteins, cytokinins=INIT_COMPARTMENTS.cytokinins,
                 Tr=INIT_COMPARTMENTS.Tr, Ag=INIT_COMPARTMENTS.Ag, Ts=INIT_COMPARTMENTS.Ts, is_growing=INIT_COMPARTMENTS.is_growing, cohorts=[], cohorts_replications=None, index=None):

        self.label = label                   #: the label of the element
        self.cohorts = cohorts  #: list of cohort values - Hack to treat tillering cases : TEMPORARY. Devrait �tre port� � l'�chelle de la plante uniquement mais je ne vois pas comment faire mieux
        self.cohorts_replications = cohorts_replications  #: dictionary of number of replications per cohort rank
        self.index = index  #: the index of the phytomer TEMPORARY

        # state parameters
        self.mstruct = mstruct               #: Structural dry mass (g)
        self.senesced_mstruct = senesced_mstruct   #: Senesced structural dry mass (g)
        self.Nstruct = Nstruct               #: Structural N mass (g)
        self.is_growing = is_growing         #: Flag indicating if the element is growing or not (:class:`bool`)
        self.green_area = green_area         #: green area (m-2)
        self.Tr = Tr                         #: Transpiration rate (mmol m-2 s-1)
        self.Ag = Ag                         #: Gross assimilation (:math:`\mu mol` m-2 s-1)
        self.Ts = Ts                         #: Organ temperature (�C)

        # state variables
        self.triosesP = triosesP             #: :math:`\mu mol` C
        self.starch = starch                 #: :math:`\mu mol` C
        self.sucrose = sucrose               #: :math:`\mu mol` C
        self.fructan = fructan               #: :math:`\mu mol` C
        self.nitrates = nitrates             #: :math:`\mu mol` N
        self.amino_acids = amino_acids       #: :math:`\mu mol` N
        self.proteins = proteins             #: :math:`\mu mol` N
        self.cytokinins = cytokinins         #: AU

        # fluxes to phloem
        self.Loading_Sucrose = None          #: Rate of sucrose loading to phloem (:math:`\mu mol` C)
        self.Loading_Amino_Acids = None      #: Rate of amino acids loading to phloem (:math:`\mu mol` N)

        # other fluxes
        self.S_Proteins = None              #: Rate of protein synthesis (:math:`\mu mol` N g-1 mstruct)
        self.S_Amino_Acids = None           #: Rate of amino acids synthesis (:math:`\mu mol` N g-1 mstruct)
        self.Regul_S_Fructan = None         #: Maximal rate of fructan synthesis (:math:`\mu mol` C g-1 mstruct)
        self.S_Starch = None                #: Rate of starch synthesis (:math:`\mu mol` C g-1 mstruct)
        self.D_Starch = None                #: Rate of starch degradation (:math:`\mu mol` C g-1 mstruct)
        self.S_Sucrose = None               #: Rate of sucrose synthesis (:math:`\mu mol` C g-1 mstruct)
        self.S_Fructan = None               #: Rate of fructan synthesis (:math:`\mu mol` C g-1 mstruct)
        self.D_Fructan = None               #: Rate of fructan degradation ((:math:`\mu mol` C g-1 mstruct)
        self.Nitrates_import = None         #: Total nitrates imported from roots (:math:`\mu mol` N nitrates)
        self.Amino_Acids_import = None      #: Total amino acids imported from roots (:math:`\mu mol` N amino acids)
        self.k_proteins = None              #: First order kinetic regulated by cytokinins concentration
        self.D_Proteins = None              #: Rate of protein degradation (:math:`\mu mol` N g-1 mstruct)
        self.cytokinins_import = None       #: Import of cytokinins (AU)
        self.D_cytokinins = None            #: Rate of cytokinins degradation (AU g-1 mstruct)

        # Integrated variables
        self.Total_Organic_Nitrogen = None  #: current total nitrogen amount (:math:`\mu mol` N)

        # intermediate variables
        self.R_Nnit_red = None              #: Nitrate reduction-linked respiration (:math:`\mu mol` C respired)
        self.R_residual = None              #: Residual maintenance respiration (cost from protein turn-over, cell ion gradients, futile cycles...) (:math:`\mu mol` C respired)
        self.Transpiration = None           #: Surfacic transpiration rate of an element (mmol H2O s-1)
        self.R_phloem_loading = None        #: Phloem loading respiration (:math:`\mu mol` C respired)
        self.Photosynthesis = None          #: Total Photosynthesis of an element integrated over a delta t (:math:`\mu mol` C)
        self.sum_respi = None               #: Sum of respirations for the element i.e. related to C loading to phloem, amino acids synthesis and residual (:math:`\mu mol` C)

    @property
    def nb_replications(self):
        return sum(int(v <= self.index) * self.cohorts_replications.get(v, 0) for v in self.cohorts) + 1

    def calculate_integrative_variables(self):
        """Calculate the integrative variables of the element.
        """
        self.Total_Organic_Nitrogen = self.calculate_Total_Organic_Nitrogen(self.amino_acids, self.proteins, self.Nstruct)

    # VARIABLES

    @staticmethod
    def calculate_total_Photosynthesis(Ag, green_area):
        """Total Photosynthesis of an element (:math:`\mu mol` C m-2 h-1 * m2).

        :param float Ag: Gross Photosynthesis rate (:math:`\mu mol` C m-2 s-1)
        :param float green_area: Green area (m2)

        :return: Rate of Total Photosynthesis (:math:`\mu mol` C h-1)
        :rtype: float
        """
        return Ag * green_area * parameters.SECOND_TO_HOUR_RATE_CONVERSION

    @staticmethod
    def calculate_Total_Transpiration(Tr, green_area):
        """Surfacic transpiration rate of an element

        :param float Tr: Transpiration rate (mmol H2O m-2 s-1)
        :param float green_area: Green area (m2)

        :return: Total transpiration (mmol H2O s-1)
        :rtype: float
        """
        return Tr * green_area

    def calculate_Regul_S_Fructan(self, Loading_Sucrose):
        """Regulating function for fructan maximal rate of synthesis.
        Negative regulation by the loading of sucrose from the phloem ("swith-off" sigmo�dal kinetic).

        :param float Loading_Sucrose: Sucrose loading (:math:`\mu mol` C)

        :return: Maximal rate of fructan synthesis (:math:`\mu mol` C g-1 mstruct h-1)
        :rtype: float
        """
        if Loading_Sucrose <= 0:
            Vmax_Sfructans = self.__class__.PARAMETERS.VMAX_SFRUCTAN_POT
        else:  # Regulation by sucrose loading
            rate_Loading_Sucrose_massic = Loading_Sucrose/self.mstruct/parameters.SECOND_TO_HOUR_RATE_CONVERSION
            Vmax_Sfructans = ((self.__class__.PARAMETERS.VMAX_SFRUCTAN_POT * self.__class__.PARAMETERS.K_REGUL_SFRUCTAN**self.__class__.PARAMETERS.N_REGUL_SFRUCTAN) /
                              (max(0, rate_Loading_Sucrose_massic**self.__class__.PARAMETERS.N_REGUL_SFRUCTAN) +
                               self.__class__.PARAMETERS.K_REGUL_SFRUCTAN**self.__class__.PARAMETERS.N_REGUL_SFRUCTAN))
        return Vmax_Sfructans

    @staticmethod
    def calculate_Total_Organic_Nitrogen(amino_acids, proteins, Nstruct):
        """Total amount of organic N (amino acids + proteins + Nstruct).
        Used to calculate residual respiration.

        :param float amino_acids: Amount of amino acids (:math:`\mu mol` N)
        :param float proteins: Amount of proteins (:math:`\mu mol` N)
        :param float Nstruct: Structural N mass (g)

        :return: Total amount of organic N (:math:`\mu mol` N)
        :rtype: float
        """
        return amino_acids + proteins + (Nstruct / EcophysiologicalConstants.N_MOLAR_MASS)*1E6

    # FLUXES

    def calculate_S_Starch(self, triosesP, T_effect_Vmax):
        """Rate of starch synthesis (:math:`\mu mol` C starch g-1 mstruct h-1).
        Michaelis-Menten function of triose phosphates.

        :param float triosesP: Amount of triose phosphates (:math:`\mu mol` C)
        :param float T_effect_Vmax: Correction to apply to enzyme activity

        :return: Rate of Starch synthesis (:math:`\mu mol` C g-1 mstruct h-1)
        :rtype: float
        """
        if triosesP <= 0:
            S_Starch = 0
        else :
            S_Starch = (((  triosesP / (self.mstruct * self.__class__.PARAMETERS.ALPHA)) * self.__class__.PARAMETERS.VMAX_STARCH) /
                (( triosesP / (self.mstruct * self.__class__.PARAMETERS.ALPHA)) + self.__class__.PARAMETERS.K_STARCH)) * parameters.SECOND_TO_HOUR_RATE_CONVERSION * T_effect_Vmax
        return S_Starch

    def calculate_D_Starch(self, starch, T_effect_Vmax):
        """Rate of starch degradation (:math:`\mu mol` C starch g-1 mstruct h-1).
        First order kinetic.

        :param float starch: Amount of starch (:math:`\mu mol` C)
        :param float T_effect_Vmax: Correction to apply to enzyme activity

        :return: Starch degradation (:math:`\mu mol` C g-1 mstruct h-1)
        :rtype: float
        """
        return max(0, self.__class__.PARAMETERS.DELTA_DSTARCH * (starch / (self.mstruct*self.__class__.PARAMETERS.ALPHA))) * parameters.SECOND_TO_HOUR_RATE_CONVERSION * T_effect_Vmax

    def calculate_S_Sucrose(self, triosesP, T_effect_Vmax):
        """Rate of sucrose synthesis (:math:`\mu mol` C sucrose g-1 mstruct h-1).
        Michaelis-Menten function of triose phosphates.

        :param float triosesP: Amount of triose phosphates (:math:`\mu mol` C)
        :param float T_effect_Vmax: Correction to apply to enzyme activity

        :return: Rate of Sucrose synthesis (:math:`\mu mol` C g-1 mstruct h-1)
        :rtype: float
        """
        if triosesP <= 0:
            S_Sucrose = 0
        else :
            S_Sucrose =  ((( triosesP / (self.mstruct*self.__class__.PARAMETERS.ALPHA)) * self.__class__.PARAMETERS.VMAX_SUCROSE) /
                    (( triosesP / (self.mstruct*self.__class__.PARAMETERS.ALPHA)) + self.__class__.PARAMETERS.K_SUCROSE)) * parameters.SECOND_TO_HOUR_RATE_CONVERSION * T_effect_Vmax
        return S_Sucrose

    def calculate_Loading_Sucrose(self, sucrose, sucrose_phloem, mstruct_axis, T_effect_conductivity):
        """Rate of sucrose loading to phloem (:math:`\mu mol` C sucrose h-1).
        Transport-resistance model.

        :param float sucrose: Amount of sucrose in the element (:math:`\mu mol` C)
        :param float sucrose_phloem: Amount of sucrose in the phloem (:math:`\mu mol` C)
        :param float mstruct_axis: Structural dry mass of the axis (g)
        :param float T_effect_conductivity: Effect of the temperature on the conductivity rate at 20�C (AU)

        :return: Rate of Sucrose loading (:math:`\mu mol` C h-1)
        :rtype: float
        """
        conc_sucrose_element = sucrose / (self.mstruct*self.__class__.PARAMETERS.ALPHA)
        conc_sucrose_phloem = sucrose_phloem / (mstruct_axis * parameters.AXIS_PARAMETERS.ALPHA)
        #: Driving compartment (:math:`\mu mol` C g-1 mstruct)
        driving_sucrose_compartment = max(conc_sucrose_element, conc_sucrose_phloem)
        #: Gradient of sucrose between the element and the phloem (:math:`\mu mol` C g-1 mstruct)
        diff_sucrose = conc_sucrose_element - conc_sucrose_phloem
        #: Conductance depending on mstruct (g2 :math:`\mu mol`-1 s-1)
        conductance = self.__class__.PARAMETERS.SIGMA_SUCROSE * self.__class__.PARAMETERS.BETA * self.mstruct**(2/3) * T_effect_conductivity

        return driving_sucrose_compartment * diff_sucrose * conductance * parameters.SECOND_TO_HOUR_RATE_CONVERSION

    def calculate_export_sucrose(self, sucrose, sucrose_hiddenzone, mstruct_hiddenzone, T_effect_conductivity):
        """Rate of sucrose exportation to hidden zone (:math:`\mu mol` C sucrose h-1).
        Transport-resistance model.

        :param float sucrose: Amount of sucrose in the element (:math:`\mu mol` C)
        :param float sucrose_hiddenzone: Sucrose amount in the hidden zone (:math:`\mu mol` C)
        :param float mstruct_hiddenzone: mstruct of the hidden zone (g)
        :param float T_effect_conductivity: Effect of the temperature on the conductivity rate at 20�C (AU)


        :return: Rate of Sucrose export (:math:`\mu mol` C h-1)
        :rtype: float
        """
        conc_sucrose_element = sucrose / (self.mstruct*self.__class__.PARAMETERS.ALPHA)
        conc_sucrose_hiddenzone = sucrose_hiddenzone / mstruct_hiddenzone
        #: Gradient of sucrose between the element and the hidden zone (:math:`\mu mol` C g-1 mstruct)
        diff_sucrose = conc_sucrose_element - conc_sucrose_hiddenzone
        #: Conductance depending on mstruct
        conductance = HiddenZone.PARAMETERS.SIGMA * self.__class__.PARAMETERS.BETA * mstruct_hiddenzone**(2/3) * T_effect_conductivity

        return diff_sucrose * conductance * parameters.SECOND_TO_HOUR_RATE_CONVERSION

    def calculate_S_Fructan(self, sucrose, Regul_S_Fructan, T_effect_Vmax):
        """Rate of fructan synthesis (:math:`\mu mol` C fructan g-1 mstruct h-1).
        Sigmo�dal function of sucrose.

        :param float sucrose: Amount of sucrose (:math:`\mu mol` C)
        :param float Regul_S_Fructan: Maximal rate of fructan synthesis regulated by sucrose loading (:math:`\mu mol` C g-1 mstruct)
        :param float T_effect_Vmax: Correction to apply to enzyme activity

        :return: Rate of Fructan synthesis (:math:`\mu mol` C g-1 mstruct h-1)
        :rtype: float
        """
        return ((max(0., sucrose) / (self.mstruct*self.__class__.PARAMETERS.ALPHA)) * Regul_S_Fructan) /\
               ((max(0., sucrose) / (self.mstruct*self.__class__.PARAMETERS.ALPHA)) + self.__class__.PARAMETERS.K_SFRUCTAN) * parameters.SECOND_TO_HOUR_RATE_CONVERSION * T_effect_Vmax

    def calculate_D_Fructan(self, sucrose, fructan, T_effect_Vmax):
        """Rate of fructan degradation (:math:`\mu mol` C fructan g-1 mstruct h-1).
        Inhibition function by the end product i.e. sucrose (Bancal et al., 2012).

        :param float sucrose: Amount of sucrose (:math:`\mu mol` C)
        :param float fructan: Amount of fructan (:math:`\mu mol` C)
        :param float T_effect_Vmax: Correction to apply to enzyme activity

        :return: Rate of Fructan degradation (:math:`\mu mol` C g-1 mstruct h-1)
        :rtype: float
        """
        d_potential = ((self.__class__.PARAMETERS.K_DFRUCTAN * self.__class__.PARAMETERS.VMAX_DFRUCTAN) /
                       ((max(0., sucrose) / (self.mstruct * self.__class__.PARAMETERS.ALPHA)) + self.__class__.PARAMETERS.K_DFRUCTAN)) * parameters.SECOND_TO_HOUR_RATE_CONVERSION * T_effect_Vmax
        d_actual = min(d_potential, max(0., fructan))
        return d_actual

    @staticmethod
    def calculate_Nitrates_import(Export_Nitrates, element_transpiration, Total_Transpiration):
        """Total nitrates imported from roots (:math:`\mu mol` N nitrates).
        Nitrates coming from roots (fraction of uptake + direct export) are distributed according to the contribution of the element to culm transpiration.

        :param float Export_Nitrates: Exported nitrates by roots (:math:`\mu mol` N)
        :param float element_transpiration: Element transpiration (mmol H2O s-1)
        :param float Total_Transpiration: Culm transpiration (mmol H2O s-1)

        :return: Total nitrates import (:math:`\mu mol` N nitrates)
        :rtype: float
        """
        if Total_Transpiration > 0:
            Nitrates_import = Export_Nitrates * (element_transpiration / Total_Transpiration)      #: Proportion of exported nitrates from roots to element
        else:  # Avoids further float division by zero error
            Nitrates_import = 0
        return Nitrates_import

    @staticmethod
    def calculate_Amino_Acids_import(roots_exported_amino_acids, element_transpiration, Total_Transpiration):
        """Total amino acids imported from roots  (:math:`\mu mol` N amino acids).
        Amino acids exported by roots are distributed according to the contribution of the element to culm transpiration.

        :param float roots_exported_amino_acids: Exported amino acids by roots (:math:`\mu mol` N)
        :param float element_transpiration: Element transpiration (mmol H2O s-1)
        :param float Total_Transpiration: Culm transpiration (mmol H2O s-1)

        :return: Total amino acids import (:math:`\mu mol` N amino acids)
        :rtype: float
        """
        if Total_Transpiration > 0:
            Amino_Acids_import = roots_exported_amino_acids * (element_transpiration / Total_Transpiration)  #: Proportion of exported amino acids from roots to organ
        else:
            Amino_Acids_import = 0
        return Amino_Acids_import

    def calculate_S_amino_acids(self, nitrates, triosesP, T_effect_Vmax):
        """Rate of amino acids synthesis (:math:`\mu mol` N amino acids h-1 g-1 MS).
        Bi-substrate Michaelis-Menten function of nitrates and triose phosphates.

        :param float nitrates: Amount of nitrates (:math:`\mu mol` N)
        :param float triosesP: Amount of triosesP (:math:`\mu mol` C)
        :param float T_effect_Vmax: Correction to apply to enzyme activity

        :return: Rate of Amino acids synthesis (:math:`\mu mol` N h-1 g-1 mstruct)
        :rtype: float
        """
        if nitrates <= 0 or triosesP <= 0:
            calculate_S_amino_acids = 0
        else:
            calculate_S_amino_acids = self.__class__.PARAMETERS.VMAX_AMINO_ACIDS /\
                                      ((1 + self.__class__.PARAMETERS.K_AMINO_ACIDS_NITRATES / (nitrates / (self.mstruct*self.__class__.PARAMETERS.ALPHA))) *
                                       (1 + self.__class__.PARAMETERS.K_AMINO_ACIDS_TRIOSESP / (triosesP / (self.mstruct*self.__class__.PARAMETERS.ALPHA)))) * \
                                      parameters.SECOND_TO_HOUR_RATE_CONVERSION * T_effect_Vmax
        return calculate_S_amino_acids

    def calculate_S_proteins(self, amino_acids, T_effect_Vmax):
        """Rate of protein synthesis (:math:`\mu mol` N proteins h-1 g-1 MS).
        Michaelis-Menten function of amino acids.

        :param float amino_acids: Amount of amino acids (:math:`\mu mol` N)
        :param float T_effect_Vmax: Correction to apply to enzyme activity

        :return: Protein synthesis (:math:`\mu mol` N h-1 g-1 mstruct)
        :rtype: float
        """
        calculate_S_proteins = (((max(0., amino_acids) / (self.mstruct * self.__class__.PARAMETERS.ALPHA)) * self.__class__.PARAMETERS.VMAX_SPROTEINS) /
                                ((max(0., amino_acids) / (self.mstruct * self.__class__.PARAMETERS.ALPHA)) + self.__class__.PARAMETERS.K_SPROTEINS)
                                ) * parameters.SECOND_TO_HOUR_RATE_CONVERSION * T_effect_Vmax
        return calculate_S_proteins

    def calculate_D_Proteins(self, proteins, cytokinins, T_effect_Vmax):
        """Rate of protein degradation (:math:`\mu mol` N proteins s-1 g-1 MS h-1).
        First order kinetic regulated by cytokinins concentration.

        :param float proteins: Amount of proteins (:math:`\mu mol` N)
        :param float cytokinins: Amount of cytokinins (AU)
        :param float T_effect_Vmax: Correction to apply to enzyme activity

        :return: Rate of protein degradation (:math:`\mu mol` N g-1 mstruct)
        :rtype: float
        """
        conc_proteins = proteins / (self.mstruct*self.__class__.PARAMETERS.ALPHA)
        conc_cytokinins = max(0, cytokinins / self.mstruct)

        regul_cytokinins = (self.__class__.PARAMETERS.VMAX_DPROTEINS_CYTOK * self.__class__.PARAMETERS.K_DPROTEINS_CYTOK**self.__class__.PARAMETERS.N_DPROTEINS) /\
                           (conc_cytokinins**self.__class__.PARAMETERS.N_DPROTEINS + self.__class__.PARAMETERS.K_DPROTEINS_CYTOK**self.__class__.PARAMETERS.N_DPROTEINS)

        return max(0, (conc_proteins * self.__class__.PARAMETERS.VMAX_DPROTEINS / (conc_proteins + self.__class__.PARAMETERS.K_DPROTEINS)) *
                   parameters.SECOND_TO_HOUR_RATE_CONVERSION * regul_cytokinins * T_effect_Vmax)

    def calculate_Loading_Amino_Acids(self, amino_acids, amino_acids_phloem, mstruct_axis, T_effect_conductivity):
        """Rate of amino acids loading to phloem (:math:`\mu mol` N amino acids h-1).
        Transport-resistance model.

        :param float amino_acids: Amount of amino acids in the element (:math:`\mu mol` N)
        :param float amino_acids_phloem: Amount of amino acids in the phloem (:math:`\mu mol` N)
        :param float mstruct_axis: Structural dry mass of the axis (g)
        :param float T_effect_conductivity: Effect of the temperature on the conductivity rate at 20�C (AU)

        :return: Amino acids loading (:math:`\mu mol` N h-1)
        :rtype: float
        """
        Conc_Amino_Acids_element = amino_acids / (self.mstruct * self.__class__.PARAMETERS.ALPHA)
        Conc_Amino_Acids_phloem = amino_acids_phloem / (mstruct_axis * parameters.AXIS_PARAMETERS.ALPHA)
        #: Driving compartment (:math:`\mu mol` N g-1 mstruct)
        driving_amino_acids_compartment = max(Conc_Amino_Acids_element, Conc_Amino_Acids_phloem)
        #: Gradient of amino acids between the element and the phloem (:math:`\mu mol` N g-1 mstruct)
        diff_amino_acids = Conc_Amino_Acids_element - Conc_Amino_Acids_phloem
        #: Conductance depending on mstruct (g2 :math:`\mu mol`-1 s-1)
        conductance = self.__class__.PARAMETERS.SIGMA_AMINO_ACIDS * self.__class__.PARAMETERS.BETA * self.mstruct**(2/3) * T_effect_conductivity

        return driving_amino_acids_compartment * diff_amino_acids * conductance * parameters.SECOND_TO_HOUR_RATE_CONVERSION

    def calculate_Export_Amino_Acids(self, amino_acids, amino_acids_hiddenzone, mstruct_hiddenzone, T_effect_conductivity):
        """Rate of amino acids exportation to hidden zone (:math:`\mu mol` N amino acids h-1).
        Transport-resistance model.

        :param float amino_acids: Amount of amino acids in the element (:math:`\mu mol` N)
        :param float amino_acids_hiddenzone: Amino acids amount in the hidden zone (:math:`\mu mol` N)
        :param float mstruct_hiddenzone: mstruct of the hidden zone (g)
        :param float T_effect_conductivity: Effect of the temperature on the conductivity rate at 20�C (AU)

        :return: Rate of Amino acids export (:math:`\mu mol` N h-1)
        :rtype: float
        """
        Conc_Amino_Acids_element = amino_acids / (self.mstruct * self.__class__.PARAMETERS.ALPHA)
        Conc_Amino_Acids_hiddenzone = amino_acids_hiddenzone / mstruct_hiddenzone
        #: Gradient of amino acids between the element and the hidden zone (:math:`\mu mol` N g-1 mstruct)
        diff_amino_acids = Conc_Amino_Acids_element - Conc_Amino_Acids_hiddenzone
        #: Conductance depending on mstruct
        conductance = HiddenZone.PARAMETERS.SIGMA * self.__class__.PARAMETERS.BETA * mstruct_hiddenzone**(2/3) * T_effect_conductivity

        return diff_amino_acids * conductance * parameters.SECOND_TO_HOUR_RATE_CONVERSION

    @staticmethod
    def calculate_cytokinins_import(roots_exporteD_cytokinins, element_transpiration, Total_Transpiration):
        """Import of cytokinins (AU).
        Cytokinin exported by roots are distributed according to the contribution of the element to culm transpiration.

        :param float roots_exporteD_cytokinins: Exported cytokinins from roots (AU)
        :param float element_transpiration: Element transpiration (mmol H2O s-1)
        :param float Total_Transpiration: Culm transpiration (mmol H2O s-1)

        :return: Cytokinin import (AU)
        :rtype: float
        """
        if Total_Transpiration > 0:
            cytokinins_import = roots_exporteD_cytokinins * (element_transpiration / Total_Transpiration)
        else:
            cytokinins_import = 0
        return cytokinins_import

    def calculate_D_cytokinins(self, cytokinins, T_effect_Vmax):
        """Rate of cytokinins degradation (AU g-1 mstruct h-1).
        First order kinetic. Vary with organ temperature.

        :param float cytokinins: Amount of cytokinins (AU)
        :param float T_effect_Vmax: Correction to apply to enzyme activity

        :return: Rate of Cytokinin degradation (AU g-1 mstruct h-1)
        :rtype: float
        """
        return max(0, self.__class__.PARAMETERS.DELTA_D_CYTOKININS * (cytokinins / (self.mstruct * self.__class__.PARAMETERS.ALPHA))) * parameters.SECOND_TO_HOUR_RATE_CONVERSION * T_effect_Vmax

    # COMPARTMENTS

    def calculate_triosesP_derivative(self, Photosynthesis, S_Sucrose, S_Starch, S_Amino_Acids):
        """ delta triose phosphates of element.

        :param float Photosynthesis: Total gross Photosynthesis (:math:`\mu mol` C)
        :param float S_Sucrose: Sucrose synthesis (:math:`\mu mol` C g-1 mstruct)
        :param float S_Starch: Starch synthesis (:math:`\mu mol` C g-1 mstruct)
        :param float S_Amino_Acids: Amino acids synthesis (:math:`\mu mol` N g-1 mstruct)

        :return: delta triose phosphates (:math:`\mu mol` C triose phosphates)
        :rtype: float
        """
        #: Contribution of triosesP to the synthesis of amino_acids
        triosesP_consumption_AA = (S_Amino_Acids / EcophysiologicalConstants.AMINO_ACIDS_N_RATIO) * EcophysiologicalConstants.AMINO_ACIDS_C_RATIO
        return Photosynthesis - (S_Sucrose + S_Starch + triosesP_consumption_AA) * (self.mstruct * self.__class__.PARAMETERS.ALPHA)

    def calculate_starch_derivative(self, S_Starch, D_Starch):
        """delta starch of element.

        :param float S_Starch: Starch synthesis (:math:`\mu mol` C g-1 mstruct)
        :param float D_Starch: Starch degradation (:math:`\mu mol` C g-1 mstruct)

        :return: delta starch (:math:`\mu mol` C starch)
        :rtype: float
        """
        return (S_Starch - D_Starch) * (self.mstruct*self.__class__.PARAMETERS.ALPHA)

    def calculate_sucrose_derivative(self, S_Sucrose, D_Starch, Loading_Sucrose, S_Fructan, D_Fructan, sum_respi):
        """delta sucrose of element.

        :param float S_Sucrose: Sucrose synthesis (:math:`\mu mol` C g-1 mstruct)
        :param float D_Starch: Starch degradation (:math:`\mu mol` C g-1 mstruct)
        :param float Loading_Sucrose: Sucrose loading (:math:`\mu mol` C)
        :param float S_Fructan: Fructan synthesis (:math:`\mu mol` C g-1 mstruct)
        :param float D_Fructan: Fructan degradation (:math:`\mu mol` C g-1 mstruct)
        :param float sum_respi: Sum of respirations for the element i.e. related to C loading to phloem, amino acids synthesis and residual (:math:`\mu mol` C)

        :return: delta sucrose (:math:`\mu mol` C sucrose)
        :rtype: float
        """
        return (S_Sucrose + D_Starch + D_Fructan - S_Fructan) * self.mstruct - sum_respi - Loading_Sucrose

    def calculate_fructan_derivative(self, S_Fructan, D_Fructan):
        """delta fructan of element.

        :param float S_Fructan: Fructan synthesis (:math:`\mu mol` C g-1 mstruct)
        :param float D_Fructan: Fructan degradation (:math:`\mu mol` C g-1 mstruct)

        :return: delta fructan (:math:`\mu mol` C fructan)
        :rtype: float
        """
        return (S_Fructan - D_Fructan) * (self.mstruct*self.__class__.PARAMETERS.ALPHA)

    def calculate_nitrates_derivative(self, Nitrates_import, S_Amino_Acids):
        """delta nitrates of element.

        :param float Nitrates_import: Nitrate import from roots (:math:`\mu mol` N)
        :param float S_Amino_Acids: Amino acids synthesis (:math:`\mu mol` N g-1 mstruct)

        :return: delta nitrates (:math:`\mu mol` N nitrates)
        :rtype: float
        """
        nitrate_reduction_AA = S_Amino_Acids  #: Contribution of nitrates to the synthesis of amino_acids
        return Nitrates_import - (nitrate_reduction_AA * self.mstruct*self.__class__.PARAMETERS.ALPHA)

    def calculate_amino_acids_derivative(self, Amino_Acids_import, S_Amino_Acids, S_Proteins, D_Proteins, Loading_Amino_Acids):
        """delta amino acids of element.

        :param float Amino_Acids_import: Amino acids import from roots (:math:`\mu mol` N)
        :param float S_Amino_Acids: Amino acids synthesis (:math:`\mu mol` N g-1 mstruct)
        :param float S_Proteins: Protein synthesis (:math:`\mu mol` N g-1 mstruct)
        :param float D_Proteins: Protein degradation (:math:`\mu mol` N g-1 mstruct)
        :param float Loading_Amino_Acids: Amino acids loading (:math:`\mu mol` N)

        :return: delta amino acids (:math:`\mu mol` N amino acids)
        :rtype: float
        """
        return Amino_Acids_import - Loading_Amino_Acids + (S_Amino_Acids + D_Proteins - S_Proteins) * (self.mstruct * self.__class__.PARAMETERS.ALPHA)

    def calculate_proteins_derivative(self, S_Proteins, D_Proteins):
        """delta proteins of element.

        :param float S_Proteins: Protein synthesis (:math:`\mu mol` N g-1 mstruct)
        :param float D_Proteins: Protein degradation (:math:`\mu mol` N g-1 mstruct)

        :return: delta proteins (:math:`\mu mol` N proteins)
        :rtype: float
        """
        return (S_Proteins - D_Proteins) * (self.mstruct * self.__class__.PARAMETERS.ALPHA)

    def calculate_cytokinins_derivative(self, import_cytokinins, D_cytokinins):
        """delta cytokinins of element.

        :param float import_cytokinins: Cytokinin import (AU)
        :param float D_cytokinins: Cytokinin degradation (AU g-1 mstruct)

        :return: delta cytokinins (AU cytokinins)
        :rtype: float
        """
        return import_cytokinins - D_cytokinins * (self.mstruct * self.__class__.PARAMETERS.ALPHA)


class ChaffElement(PhotosyntheticOrganElement):
    """
    The class :class:`ChaffElement` defines the CN exchanges in a chaff element.
    """

    PARAMETERS = parameters.CHAFF_ELEMENT_PARAMETERS  #: the internal parameters of the chaffs elements


class LaminaElement(PhotosyntheticOrganElement):
    """
    The class :class:`LaminaElement` defines the CN exchanges in a lamina element.
    """

    PARAMETERS = parameters.LAMINA_ELEMENT_PARAMETERS  #: the internal parameters of the laminae elements


class InternodeElement(PhotosyntheticOrganElement):
    """
    The class :class:`InternodeElement` defines the CN exchanges in an internode element.
    """

    PARAMETERS = parameters.INTERNODE_ELEMENT_PARAMETERS  #: the internal parameters of the internodes elements


class PeduncleElement(PhotosyntheticOrganElement):
    """
    The class :class:`PeduncleElement` defines the CN exchanges in a peduncle element.
    """

    PARAMETERS = parameters.PEDUNCLE_ELEMENT_PARAMETERS  #: the internal parameters of the peduncles elements


class SheathElement(PhotosyntheticOrganElement):
    """
    The class :class:`SheathElement` defines the CN exchanges in a sheath element.
    """

    PARAMETERS = parameters.SHEATH_ELEMENT_PARAMETERS  #: the internal parameters of the sheaths elements


class Soil(object):
    """
    The class :class:`Soil` defines the amount of nitrogen in the volume of soil explored by roots.
    """

    PARAMETERS = parameters.SOIL_PARAMETERS  #: the internal parameters of the soil

    def __init__(self, volume=None, nitrates=None, Tsoil=None):

        # state parameters
        self.volume = volume                   #: volume of soil explored by roots (m3)
        self.Tsoil = Tsoil                     #: soil temperature (�C)
        self.constant_Conc_Nitrates = False    #: If True, the model run with a constant soil nitrate concentration (bool)

        # state variables
        self.nitrates = nitrates               #: :math:`\mu mol` N nitrates

        # intermediate variables
        self.Conc_Nitrates_Soil = None         #: soil nitrate concentration Unloading (:math:`\mu mol` N m-3 soil)
        self.mineralisation = None             #: mineralisation on organic N into nitrates in soil (:math:`\mu mol`)

    @staticmethod
    def calculate_temperature_effect_on_Vmax(Tsoil):
        """Effect of the temperature on maximal enzyme activity
        Should multiply the rate at 20�C

        :param float Tsoil: Soil temperature (�C)

        :return: Correction to apply to enzyme activity
        :rtype: float
        """
        Tref = 20 + 273.15
        Tk = Tsoil + 273.15
        R = 8.3144  #: Physical parameter: Gas constant (J mol-1 K-1)
        deltaHa = 55  # 89.7  #: Enthalpie of activation of parameter pname (kJ mol-1)
        deltaS = 0.48  # 0.486  #: entropy term of parameter pname (kJ mol-1 K-1)
        deltaHd = 154   # 149.3 #: Enthalpie of deactivation of parameter pname (kJ mol-1)

        f_activation = np.exp((deltaHa * (Tk - Tref)) / (R * 1E-3 * Tref * Tk))  #: Energy of activation (normalized to unity)

        f_deactivation = (1 + np.exp((Tref * deltaS - deltaHd) / (Tref * R * 1E-3))) / (1 + np.exp((Tk * deltaS - deltaHd) / (Tk * R * 1E-3)))  #: Energy of deactivation (normalized to unity)

        return f_activation * f_deactivation

    # VARIABLES

    def calculate_Conc_Nitrates(self, nitrates):
        """Nitrate concentration in soil.

        :param float nitrates: Amount of nitrates (:math:`\mu mol` N)

        :return: Nitrate concentration (:math:`\mu mol` nitrates m-3)
        :rtype: float
        """
        return max(0, (nitrates / self.volume))

    # FLUX
    @staticmethod
    def calculate_mineralisation(T_effect_Vmax):
        """Mineralisation on organic N into nitrates in soil.

        :param float T_effect_Vmax: Correction to apply to enzyme activity

        :return: Rate of Nitrate mineralisation (:math:`\mu mol` h-1)
        :rtype: float
        """
        return parameters.SOIL_PARAMETERS.MINERALISATION_RATE * parameters.SECOND_TO_HOUR_RATE_CONVERSION * T_effect_Vmax

    # COMPARTMENTS

    @staticmethod
    def calculate_nitrates_derivative(mineralisation, soil_contributors, culm_density, constant_Conc_Nitrates):
        """delta soil nitrates.

        :param float mineralisation: N mineralisation in soil (:math:`\mu mol` m-2 N nitrates)
        :param (float, int) soil_contributors: A tuple with (Nitrate uptake per axis (:math:`\mu mol` N nitrates), the plant id)
        :param dict [plant_id, culm_density] culm_density: A dictionary of culm density (culm_density = {plant_id: culm_density, ...})
        :param bool constant_Conc_Nitrates: If True, the model run with a constant soil nitrate concentration.

        :return: delta nitrates (:math:`\mu mol` N nitrates)
        :rtype: float
        """
        delta_Nitrates = 0
        if not constant_Conc_Nitrates:
            Uptake_Nitrates = 0
            for root_uptake, plant_id in soil_contributors:
                Uptake_Nitrates += root_uptake * culm_density[plant_id]  # TODO: temporary, will be removed in next version
            delta_Nitrates = mineralisation - Uptake_Nitrates
        return delta_Nitrates
