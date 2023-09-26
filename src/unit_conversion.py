import os, sys
import numpy as np

PI2    = 2.*np.pi
FS     = 41.341374575751   # fs to atomic unit time
BOHR   = 0.529177249       # bohr to angstrom AA
H2EV   = 27.21140795       # hartree to ev
EV2J   = 1.602176634*1e-19 # ev to J=C*V
Mole   = 6.022*1e23
Cal2J  = 4.184             # cal to J
EV2kJM = EV2J*Mole*1e-3    # ev to kcal/mol

C         = 299792458        # speed of light m/s
Boltzmann = 1.380649*1e-23   # J/K
Planck    = 6.62607015*1e-34 # J/Hz
PlanckBar = Planck/PI2       # Js

WN    = EV2J/C/Planck*1e-2    # ev to wavenumber cm^-1
EV2ns = Planck/EV2J*1e9       # ev to ns ### wikipedia use planckbar!
EV2nm = EV2ns*C               # ev to ns
EV2K  = EV2J/Boltzmann        # ev to temperature K
D2kg  = 1.4924180856045*1e-10 # Dalton to kg

# units of these qualities
# milli(m), micro(u), nano(n), pico(p), femto(f), atto(a)
# deci(d), centi(c), angstrom(aa)
units_long = {
        'time':        ['day', 'hour', 'minute', 'second', 'millisecond', \
                        'microsecond', 'nanosecond', 'picosecond', \
                        'femtosecond', 'atomicunit', 'attosecond'],
        'length':      ['meter', 'decimeter', 'centimeter', 'millimeter', \
                        'micrometer', 'nanometer', 'angstrom', 'bohr', \
                        'picometer', 'femtometer', 'attometer'],
        'energy':      ['hartree', 'electronvolt', 'milliev', 'kcal/mol', 'kj/mol'],
        'frequency':   ['terahertz', 'gigahertz', 'megahertz', 'kilohertz', 'hertz', 'cm^-1'],
        'temperature': ['kelvin'],
        'mass':        ['kilogram', 'gram', 'dalton'],
        }

units_short = {
        'time':        ['d', 'h', 'm', 's', 'ms', 'us', 'ns', 'ps', 'fs', 'au', 'as'],
        'length':      ['m', 'dm', 'cm', 'mm', 'um', 'nm', 'aa', 'b', 'pm', 'fm', 'am'],
        'energy':      ['h', 'ev', 'mev', 'kcal', 'kj'],
        'frequency':   ['thz', 'ghz', 'mhz', 'khz', 'hz', 'cm-1'],
        'temperature': ['k'],
        'mass':        ['kg', 'g', 'u'],
        }
properties = list(units_long.keys())
# indices of atomic units for converting different properties
Idx_name = ['ns', 'nm', 'ev', 'cm-1', 'k', 'kg']
Idx = [0]*len(Idx_name)
for i, s in enumerate(properties):
    Idx[i] = units_short[s].index(Idx_name[i])
#print(Idx)

units_conversion = {
        'time':        np.array([8.64*1e12, 3.6*1e11, 6.*1e10, 1e9, 1e6, 1e3, 1., 1e-3, 1e-6, 1e-6/FS, 1e-9]),
        'length':      np.array([1e9, 1e8, 1e7, 1e6, 1e3, 1., .1, BOHR/10, 1e-3, 1e-6, 1e-9]),
        'energy':      np.array([H2EV, 1., 1e-3, Cal2J/EV2kJM, 1./EV2kJM]),
        'frequency':   np.array([1e12, 1e9, 1e6, 1e3, 1., C*1e2]),
        'temperature': np.array([1.]),
        'mass':        np.array([1e3, 1., D2kg]),
        'time_to_length':           C,          # ns to nm
        'energy_to_time':           EV2ns,      # ev to ns
        'energy_to_length':         EV2nm,      # ev to nm
        'energy_to_frequency':      WN,         # ev to cm-1
        'energy_to_temperature':    EV2K,       # ev to K
        'energy_to_mass':           EV2J/C**2,  # ev to kg
        'frequency_to_time':        EV2ns/WN,   # hz to ns
        'frequency_to_length':      EV2nm/WN,   # hz to nm
        'frequency_to_temperature': EV2K/WN,    # hz to K
        'time_to_temperature':      EV2K/EV2ns, # ns to K
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
        raise ValueError('unaccepted property conversion.')

    if ('time' in prop or 'length' in prop) and ('energy' in prop or 'frequency' in prop):
        value = 1./value
        if prop[1]+'_to_'+prop[0] in units_conversion: c2 = 1./c2

    #print(c[0], c2, c[1])
    return value * c[0] * c2 / c[1]


def convert_units(value, unit0='au', unit1='fs'):
    prop, index = find_properties(unit0, unit1)
    if prop[0] == prop[1]:
        return convert_same_units(value, prop[0], index)
    else:
        return convert_different_units(value, prop, index)


def convert_other_property(value, unit0='au'):
    prop, index = find_properties(unit0, None)
    prop, index = prop[0], index[0]
    i = properties.index(prop)

    value1 = np.zeros(len(properties))
    for j, p in enumerate(properties):
        k = Idx[j]
        if j == i:
            value1[j] = convert_same_units(value, prop, [index, k])
        elif j != i:
            value1[j] = convert_different_units(value, [prop, p], [index, k])

        print(value1[j], Idx_name[j], end='  ')
    print('')
    return value1


if __name__ == '__main__':
    value0 = float(sys.argv[1])
    unit0, unit1 = sys.argv[2].lower(), sys.argv[3].lower()
    value1 = convert_units(value0, unit0, unit1)
    print(str(value0)+' '+unit0+' = '+str(value1)+' '+unit1)

    #convert_other_property(value0, unit0)
