# Following https://github.com/PyORBIT-Collaboration/py-orbit/blob/2d4f5a2062849b1e2ec237107f6371adae6bf2d7/src/orbit/MaterialInteractions/cross_sections.cc
import numpy as np

# Energy in GeV
energy = np.array([0.0005, 0.001, 0.0015, 0.002, 0.0025, 0.003, 0.0035, 0.004, 0.0045, 0.005, 0.0055, 0.006, 0.0065, 0.007, 0.0075, 0.008, 0.009, 0.01, 0.011, 0.012, 0.013, 0.014, 0.015, 0.0175, 0.02, 0.0225, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.051, 0.055, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.14, 0.16, 0.18, 0.2, 0.225, 0.25, 0.275, 0.3, 0.325, 0.35, 0.375, 0.4, 0.5, 0.7, 1.0, 1.5, 2.0, 2.5])

# Elastic cross sections for atomic number Z=6
z6_elastic = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.321, 0.335, 0.349, 0.363, 0.377, 0.386, 0.395, 0.404, 0.413, 0.423, 0.434, 0.456, 0.479, 0.509, 0.539, 0.569, 0.598, 0.628, 0.685, 0.743, 0.783, 0.822, 0.878, 0.915, 0.938, 0.813, 0.688, 0.678, 0.641, 0.594, 0.515, 0.448, 0.392, 0.345, 0.305, 0.272, 0.214, 0.167, 0.138, 0.117, 0.098, 0.085, 0.077, 0.072, 0.068, 0.067, 0.066, 0.067, 0.077, 0.09, 0.102, 0.108, 0.114, 0.113])

# Inelastic cross sections for atomic number Z=6
z6_inelastic = np.array([0., 0., 0., 0., 0., 0., 0., 0., 0., 0.013, 0.078, 0.131, 0.175, 0.212, 0.244, 0.271, 0.314, 0.346, 0.370, 0.389, 0.403, 0.413, 0.421, 0.432, 0.434, 0.432, 0.427, 0.412, 0.394, 0.376, 0.359, 0.344, 0.288, 0.281, 0.272, 0.257, 0.245, 0.235, 0.228, 0.223, 0.220, 0.223, 0.213, 0.212, 0.212, 0.212, 0.213, 0.214, 0.216, 0.217, 0.219, 0.221, 0.223, 0.229, 0.246, 0.261, 0.261, 0.257, 0.255])

# Elastic cross sections for atomic number Z=13
z13_elastic = np.array([0., 0., 0., 0., 0., .091, .164, .237, .311, .384, .428, .473, .518, .562, .607, .623, .654, .686, .692, .699, .701, .700, .700, .707, .714, .737, .761, .832, .905, .976, .993, 1.010, 1.007, .994, .979, .921, .854, .786, .720, .658, .601, .490, .382, .320, .274, .231, .202, .181, .167, .159, .153,.151, .151, .170, .198, .220, .232, .243, .242])

# Inelastic cross sections for atomic number Z=13
z13_inelastic = np.array([0., 0., 0., 0., .123, .280, .391, .474, .537, .587, .626, .658, .684, .705, .722, .736, .757, .771, .780, .785, .787, .787, .786, .777, .764, .749, .732, .699, .667, .638, .611, .587, .530, .518, .503, .478, .456, .437, .422, .410, .402, .397, .394, .393, .392, .393, .395, .397, .399, .402, .405, .408, .411, .415, .438, .457, .460, .454, .450])

# Elastic cross sections for atomic number Z=26
z26_elastic = np.array([0., 0., 0., 0., 0., .001, .015, .030, .045, .059, .109, .159, .209, .259, .309, .364, .474, .584, .668, .753, .819, .867, .915, .964, 1.013, 1.030, 1.047, 1.063, 1.075, 1.098, 1.121, 1.143, 1.153, 1.192, 1.240, 1.293, 1.306, 1.289, 1.250, 1.198, 1.137, .984, .777, .669, .582, .499, .438, .393, .361, .340, .325, .317, .314, .346, .389, .425, .444, .461, .460])

# Inelastic cross sections for atomic number Z=26
z26_inelastic = np.array([0., 0., 0., 0., 0., 0., 0., .017, .174, .301, .404, .490, .562, .623, .676, .721, .795, .852, .896, .931, .959, .981, .999, 1.028, 1.044, 1.05, 1.05, 1.04, 1.021, .999, .976, .952, .920, .899, .873, .834, .790, .753, .719, .698, .685, .672, .677, .674, .673, .674, .675, .678, .681, .685, .689, .693,.697, .679, .710, .735, .744, .740, .739])

# Elastic cross sections for atomic number Z=29
z29_elastic = np.array([0., 0., 0., 0., 0., 0., .017, .025, .033, .077, .121, .165, .208, .252, .310, .426,.541, .634, .726, .800, .855, .910, .971, 1.033, 1.058, 1.083, 1.113, 1.132, 1.156, 1.161, 1.165, 1.175, 1.216, 1.268, 1.337, 1.368, 1.368, 1.343, 1.300, 1.245, 1.092, .869, .752, .657, .565, .496, .446, .410, .384, .368, .358, .353, .385, .429, .468, .488, .505, .504])

# Inelastic cross sections for atomic number Z=29
z29_inelastic = np.array([0., 0., 0., 0., 0., 0., 0., .063, .225, .354, .460, .548, .622, .685, .739, .786, .862, .920, .966, 1.003, 1.032, 1.055, 1.073, 1.104, 1.121, 1.128, 1.128, 1.118, 1.100, 1.077, 1.053, 1.029, .986, .966, .942, .902, .859, .822, .788, .765, .749, .732, .736, .733, .732, .732, .734, .737, .740, .744,.748, .752, .756, .736, .768, .794, .805, .804, .806])

# Elastic cross sections for atomic number Z=73
z73_elastic = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, .004, .011, .019, .072, .126, .207, .316, .425, .704, .982, 1.124, 1.266, 1.401, 1.506, 1.586, 1.593, 1.600, 1.608, 1.639, 1.679, 1.773, 1.884, 2.001, 2.109, 2.198, 2.261, 2.255, 1.984, 1.831, 1.679, 1.506, 1.359, 1.240, 1.146, 1.07, 1.022, .985, .963, .997, 1.053, 1.119, 1.186, 1.193, 1.192])

# Inelastic cross sections for atomic number Z=73
z73_inelastic = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, .012, .286, .505, .685, .833, .958, 1.065, 1.156, 1.334, 1.462, 1.557, 1.628, 1.723, 1.779, 1.810, 1.825, 1.831, 1.824, 1.823, 1.821, 1.807, 1.784, 1.761, 1.736, 1.705, 1.676, 1.629, 1.604, 1.591, 1.578, 1.574, 1.568, 1.564, 1.560, 1.555, 1.551, 1.546, 1.541, 1.549, 1.599, 1.641, 1.666, 1.669, 1.674])

# Elastic cross sections for atomic number Z=74
z74_elastic = np.array([0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., .003, .010, .017, .067, .118, .197, .305, .412, .694, .975, 1.121, 1.267, 1.404, 1.511, 1.591, 1.600, 1.610, 1.618, 1.649, 1.689, 1.782, 1.893, 2.010, 2.120, 2.211, 2.276, 2.275, 2.006, 1.852, 1.700, 1.526, 1.378, 1.258, 1.163, 1.090, 1.037, 1.000, .977, 1.011, 1.067, 1.134, 1.171, 1.209, 1.207])

# Inelastic cross sections for atomic number Z=74
z74_inelastic = np.array([0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., .265, .489, .672, .823, .951, 1.059, 1.152, 1.334, 1.465, 1.562, 1.635, 1.733, 1.790, 1.823, 1.840, 1.846, 1.840, 1.839, 1.838, 1.825, 1.804, 1.782, 1.759, 1.728, 1.700, 1.653, 1.626, 1.612, 1.599, 1.594, 1.588, 1.583, 1.578, 1.573, 1.568, 1.562, 1.557, 1.568, 1.618, 1.659, 1.686, 1.689, 1.693])

# Elastic cross sections for atomic number Z=78 (newer table)
z78_elastic = np.array([0. ,0. ,0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.002, 0.007, 0.011, 0.053, 0.095, 0.166, 0.267, 0.368, 0.660, 0.951, 1.110, 1.268, 1.414, 1.527, 1.611, 1.628, 1.644, 1.652, 1.684, 1.724, 1.815, 1.926, 2.044, 2.159, 2.257, 2.330, 2.344, 2.082, 1.930, 1.776, 1.599, 1.447, 1.323, 1.224, 1.148, 1.093, 1.053, 1.029, 1.063, 1.120, 1.188, 1.226, 1.265, 1.263])

# Inelastic cross sections for atomic number Z=78 (newer table)
z78_inelastic = np.array([0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.176, 0.415, 0.611, 0.774, 0.911, 1.028, 1.128, 1.325, 1.467, 1.573, 1.653, 1.762, 1.827, 1.866, 1.888, 1.898, 1.947, 1.935, 1.920, 1.896, 1.867, 1.833, 1.799, 1.762, 1.729, 1.678, 1.658, 1.649, 1.641, 1.639, 1.638, 1.637, 1.637, 1.637, 1.637, 1.637, 1.638, 1.633, 1.684, 1.727, 1.755, 1.758, 1.763])

# Elastic cross sections for atomic number Z=82
z82_elastic = np.array([0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., .001, .004, .007, .040, .074, .138, .231, .325, .625, .925, 1.097, 1.270, 1.426, 1.545, 1.633, 1.659, 1.684, 1.692, 1.725, 1.765, 1.853, 1.962, 2.083, 2.203, 2.309, 2.392, 2.423, 2.170, 2.020, 1.865, 1.685, 1.528, 1.399, 1.296, 1.216, 1.158, 1.116, 1.090, 1.124, 1.181, 1.251, 1.291, 1.330, 1.328])

# Inelastic cross sections for atomic number Z=82
z82_inelastic = np.array([0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., .090, .346, .556, .731, .878, 1.004, 1.112, 1.325, 1.479, 1.594, 1.682, 1.803, 1.877, 1.922, 1.948, 1.962, 2.074, 2.048, 2.015, 1.979, 1.939, 1.892, 1.845, 1.801, 1.762, 1.706, 1.695, 1.691, 1.690, 1.691, 1.695, 1.699, 1.705, 1.712, 1.718, 1.725, 1.733, 1.711, 1.763, 1.807, 1.836, 1.840,1.845])


class MP:
    def __init__(self, z, a, rho, radlength, R, ecross, icross):
        self.z = z # atomic number
        self.a = a # atomic weight in amu
        self.rho = rho # density in kg/m^3
        self.radlength = radlength # in meters
        self.R = R # ?
        self.ecross = ecross # in barn (1e-24 cm^2)
        self.icross = icross # in barn (1e-24 cm^2)


carbon = MP(6.0, 12.01, 2.27e+03, 18.8/100.0, 0.94*pow(12.01, 1./3.), z6_elastic, z6_inelastic)
aluminum = MP(13.0, 26.92, 2.7e+03, 8.9/100.0, 0.94*pow(26.92, 1./3.), z13_elastic, z13_inelastic)
iron = MP(26, 55.85, 7.87e+03, 1.76/100.0, 0.94*pow(55.85, 1./3.), z26_elastic, z26_inelastic)
copper = MP(29, 63.546, 8.96e+03, 1.43/100.0, 0.94*pow(63.55, 1./3.), z29_elastic, z29_inelastic)
tantalum = MP(73, 180.95, 16.6e+03, 0.411/100.0, 0.94*pow(180.95, 1./3.), z73_elastic, z73_inelastic)
tungstun = MP(74, 183.84, 19.3e+03, 0.35/100.0, 0.94*pow(183.84, 1./3.), z74_elastic, z74_inelastic)
platinum = MP(78, 195.08, 21.45e+03, 0.305/100.0, 0.94*pow(195.08, 1./3.), z78_elastic, z78_inelastic)
lead = MP(82, 207.2, 11.35e+03, 0.56/100.0, 0.94*pow(207.2, 1./3.), z82_elastic, z82_inelastic)


def get_z(ma_index):
    if ma_index == 0:
        z = carbon.z
    elif ma_index == 1:
        z = aluminum.z
    elif ma_index == 2:
        z = iron.z
    elif ma_index == 3:
        z = copper.z
    elif ma_index == 4:
        z = tantalum.z
    elif ma_index == 5:
        z = tungstun.z
    elif ma_index == 6:
        z = platinum.z
    else:
        z = lead.z

    return z


def get_a(ma_index):
    if ma_index == 0:
        a = 12.01
    elif ma_index == 1:
        a = 26.92
    elif ma_index == 2:
        a = 55.85
    elif ma_index == 3:
        a = 63.546
    elif ma_index == 4:
        a = 180.95
    elif ma_index == 5:
        a = 183.84
    elif ma_index == 6:
        a = 195.08
    else:
        a = 207.2
        
    return a


def get_rho(ma_index):
    if ma_index == 0:
        rho = carbon.rho
    elif ma_index == 1:
        rho = aluminum.rho
    elif ma_index == 2:
        rho = iron.rho
    elif ma_index == 3:
        rho = copper.rho
    elif ma_index == 4:
        rho = tantalum.rho
    elif ma_index == 5:
        rho = tungstun.rho
    elif ma_index == 6:
        rho = platinum.rho
    else:
        rho = lead.rho
        
    return rho


def get_radlength(ma_index):
    if ma_index == 0:
        radlength = carbon.radlength
    elif ma_index == 1:
        radlength = aluminum.radlength
    elif ma_index == 2:
        radlength = iron.radlength
    elif ma_index == 3:
        radlength = copper.radlength
    elif ma_index == 4:
        radlength = tantalum.radlength
    elif ma_index == 5:
        radlength = tungstun.radlength
    elif ma_index == 6:
        radlength = platinum.radlength
    else:
        radlength = lead.radlength

    return radlength


def get_elastic_crosssection(energyrequest, ma_index, energy=energy):
    ec = None
    
    if ma_index == 0:
        ec = carbon.ecross
    elif ma_index == 1:
        ec = aluminum.ecross
    elif ma_index == 2:
        ec = iron.ecross
    elif ma_index == 3:
        ec = copper.ecross
    elif ma_index == 4:
        ec = tantalum.ecross
    elif ma_index == 5:
        ec = tungstun.ecross
    elif ma_index == 6:
        ec = platinum.ecross
    else:
        ec = lead.ecross
    
    energytracker = 0.0

    ilow = 57
    for i in range(58):
        if energyrequest <= energy[i+1]:
            ilow = i
            break
            
    efrac = (energyrequest - energy[ilow]) / (energy[ilow + 1] - energy[ilow])
    
    cross_section = ec[ilow] + efrac*(ec[ilow + 1] - ec[ilow])
    #print("efrac, ilow, ec ", energyrequest, ilow, ec[ilow])
    return cross_section


def get_inelastic_crosssection(energyrequest, ma_index, energy=energy):
    
    if ma_index == 0:
        ic = carbon.icross
    elif ma_index == 1:
        ic = aluminum.icross
    elif ma_index == 2:
        ic = iron.icross
    elif ma_index == 3:
        ic = copper.icross
    elif ma_index == 4:
        ic = tantalum.icross
    elif ma_index == 5:
        ic = tungstun.icross
    elif ma_index == 6:
        ic = platinum.icross
    else:
        ic = lead.icross
        
    energytracker = 0.0
    
    ilow = 57
    for i in range(58):
        if energyrequest <= energy[i+1]:
            ilow = i
            break
            
    efrac = (energyrequest - energy[ilow]) / (energy[ilow + 1] - energy[ilow])
    
    cross_section = ic[ilow] + efrac*(ic[ilow + 1] - ic[ilow])
    
    return cross_section
