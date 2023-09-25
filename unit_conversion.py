import os, sys
import numpy as np

# three constant values
#time   = 1 # au
#length = 1 # au ie bohr
#energy = 1 # au ie Hartree
FS   = 41.341374575751 # au
BOHR = 0.529177249 # AA
H2EV = 27.21140795 # hartree
EV2J = 1.602176634*1e-19 # ev to J=C*V
Mole = 6.022*1e23
Cal2J = 4.184
EV2kJM = EV2J*Mole*1e-3

Boltzmann = 1.380649*1e-23 # J/K
Planck = 6.62607015*1e-34 # Js
C = 299792458 # speed of light m/s

WN = EV2J/C/Planck*1e-2 # ev to wavenumber cm^-1
EV2K = EV2J/Boltzmann # ev to temperature K


# units of these qualities
# milli(m), micro(u), nano(n), pico(p), femto(f), atto(a)
# deci(d), centi(c), angstrom(aa)
units_long = {
        'time':   ['day', 'hour', 'minute', 'second', 'millisecond', \
                   'microsecond', 'nanosecond', 'picosecond', \
                   'femtosecond', 'atomicunit', 'attosecond'],
        'length': ['meter', 'decimeter', 'centimeter', 'millimeter', \
                   'micrometer', 'nanometer', 'angstrom', 'bohr', \
                   'picometer', 'femtometer', 'attometer'],
        'energy': ['hartree', 'electronvolt', 'milliev', 'kcal/mol', 'kj/mol'],
        'frequency': ['terahertz', 'gigahertz', 'megahertz', 'kilohertz', 'hertz', 'cm^-1'],
        'temperature': ['kelvin'],
        }

units_short = {
        'time':   ['d', 'h', 'm', 's', 'ms', 'us', 'ns', 'ps', 'fs', 'au', 'as'],
        'length': ['m', 'dm', 'cm', 'mm', 'um', 'nm', 'aa', 'bohr', 'pm', 'fm', 'am'],
        'energy': ['h', 'ev', 'mev', 'kcal', 'kj'],
        'frequency': ['thz', 'ghz', 'mhz', 'khz', 'hz', 'cm-1'],
        'temperature': ['k'],
        }
properties = list(units_long.keys())
# indices of atomic units for converting different properties
Idx = ['ns', 'nm', 'ev', 'cm-1', 'k']
for i, s in enumerate(properties):
    Idx[i] = units_short[s].index(Idx[i])
#print(Idx)

units_conversion = {
        'time':        np.array([8.64*1e12, 3.6*1e11, 6.*1e10, 1e9, 1e6, 1e3, 1., 1e-3, 1e-6, 1e-6/FS, 1e-9]),
        'length':      np.array([1e9, 1e8, 1e7, 1e6, 1e3, 1., .1, BOHR/10, 1e-3, 1e-6, 1e-9]),
        'energy':      np.array([H2EV, 1., 1e-3, Cal2J/EV2kJM, 1./EV2kJM]),
        'frequency':   np.array([1e12, 1e9, 1e6, 1e3, 1., C*1e2]),
        'temperature': np.array([1.]),
        'time_to_length':        C,   # ns to nm
        'energy_to_frequency':   WN,  # ev to cm-1
        'energy_to_temperature': EV2K, # ev to K
        'frequency_to_length':   1e7, # cm-1 to nm
        'frequency_to_time':     1e7/C,
        'energy_to_length':      1e7*WN,
        'energy_to_time':        1e7*WN/C,
        'frequency_to_temperature': EV2K/WN,
        }


def find_properties(unit0='au', unit1='fs'):
    prop, index = [], []
    for u in [unit0, unit1]:
        for s in properties:
            for name in [units_short, units_long]:
                if u in name[s]:
                    prop.append( s)
                    index.append( name[s].index(u))
    if len(prop) > 2:
        print(prop, index)
        raise ValueError('please specify units with full names:')
    return prop, index


def convert_same_units(value, prop, index):
    conv = units_conversion[prop]
    i, j = index
    return value * conv[i] / conv[j]


def convert_different_units(value, prop, index):
    c = [0]*2
    for k, p in enumerate(prop):
        conv = units_conversion[p]
        c[k] = conv[index[k]] / conv[Idx[properties.index(p)]]

    c2 = 0
    if prop[0]+'_to_'+prop[1] in units_conversion:
        c2 = units_conversion[prop[0]+'_to_'+prop[1]]
    elif prop[1]+'_to_'+prop[0] in units_conversion:
        c2 = 1./units_conversion[prop[1]+'_to_'+prop[0]]
    else:
        raise ValueError('uneccepted property conversion.')

    if 'length' in prop and ('energy' in prop or 'frequency' in prop):
        value = 1./value

    return value * c[0] * c2 / c[1]


def convert_units(value, unit0='au', unit1='fs'):
    prop, index = find_properties(unit0, unit1)
    if prop[0] == prop[1]:
        return convert_same_units(value, prop[0], index)
    else:
        return convert_different_units(value, prop, index)


if __name__ == '__main__':
    value0 = float(sys.argv[1])
    unit0, unit1 = sys.argv[2].lower(), sys.argv[3].lower()
    value1 = convert_units(value0, unit0, unit1)
    print(str(value0)+' '+unit0+' = '+str(value1)+' '+unit1)

