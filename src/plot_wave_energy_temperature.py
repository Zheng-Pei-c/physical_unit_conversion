import os, sys
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import ticker
import matplotlib.colors as mcolors

from unit_conversion import convert_units

def plot_wave_energy_frequency():
    # 1nm to 1mm (1e6nm)
    wave_names = {'EUV':    [1,     100], # nm extreme UV       1
                  'VUV':    [100,   190], # nm vacuum UV
                  'DUV':    [190,   280], # nm deep UV          3
                  'UV_B':   [280,   315], # nm mid UV
                  'UV_A':   [315,   380], # nm near UV          5 #
                  'violet': [380,   435], # nm violet           6
                  'blue':   [435,   500], # nm blue
                  'cyan':   [500,   520], # nm cyan
                  'green':  [520,   565], # nm green
                  'yellow': [565,   590], # nm yellow
                  'orange': [590,   625], # nm orange
                  'red':    [625,   780], # nm red             12 #
                  'NIR_A':  [780,  1400], # nm near infrared I 13
                  'NIR_B':  [1400, 3000], # nm near infrared II   #
                  'NIR_C':  [3,    50  ], # um mid infrared    15
                  'FIR':    [50,   1000], # um far infrared    16
                  }
    wavelengths = []
    for key, [v0, v1] in wave_names.items():
        wavelengths.append(np.linspace(v0, v1, (v1-v0)//2))
    #wavelengths = np.hstack(wavelengths)
    wavelengths = [np.hstack(wavelengths[0]), np.hstack(wavelengths[1]),
                   np.hstack(wavelengths[2:12]), np.hstack(wavelengths[12:13]),
                   np.hstack(wavelengths[14]), np.hstack(wavelengths[15]),
                   ]
    #print(wavelengths.shape)

    n = len(wavelengths)
    energy, frequency, temperature, time = [0]*n, [0]*n, [0]*n, [0]*n

    units = {'l': ['nm']*n,
             'e': ['eV']*n,
             'f': ['THz']*n,
             'T': ['K']*n,
             't': ['fs']*n,
             }
    units['l'][-2:] = ['um']*2
    units['e'][-3]  = 'kcal'
    units['e'][-2:] = ['meV']*2
    units['f'][-3:] = ['cm^-1']*3
    units['t'][0]   = 'as'
    units['t'][1]   = 'au'


    for i, length in enumerate(wavelengths):
        e = convert_units(length, units['l'][i], units['e'][i])
        energy[i] = np.copy(e)
        frequency[i] = convert_units(e, units['e'][i], units['f'][i])
        temperature[i] = convert_units(e, units['e'][i], units['T'][i])
        time[i] = convert_units(e, units['e'][i], units['t'][i])



    fig_name = 'spectrum_reference'
    fig = plt.figure(figsize=(12, 6), dpi=300)
    #fig = plt.figure()
    nrow, ncol = 2, 3
    colors = list(wave_names.keys())

    for i, length in enumerate(wavelengths):
        ax1 = plt.subplot(nrow, ncol, i+1)
        ax2 = ax1.twinx()
        ax3 = ax1.twiny()

        ax1.plot(length, energy[i], color='b')
        ax2.plot(length, temperature[i], color='b') # has save curve as energy, just for yticks

        unit = '$\mu$m' if 'um' in units['l'][i] else units['l'][i]
        ax1.set_xlabel('Wavelength ('+unit+')')
        ax1.set_ylabel('Energy ('+units['e'][i]+')')
        ax2.set_ylabel('Temperature ('+units['T'][i]+')')

        if i < 3:
            if i < 2:
                ax3.plot(time[i], energy[i], color='b')
            else:
                ax3.plot(time[i], energy[i], color='g')
            ax3.set_xlabel('Time ('+units['t'][i]+')')
        else:
            ax3.plot(frequency[i], energy[i], color='g')
            unit = 'cm$^{-1}$' if '-1' in units['f'][i] else units['f'][i]
            ax3.set_xlabel('Frequency ('+unit+')')

        alpha = .15 if i==2 else .4
        ax1.grid(ls='--', alpha=alpha)

        if i == 2:
            hight = ax1.get_ylim()
            for key, (v0, v1) in list(wave_names.items())[5:12]:
                width = v1 - v0
                ax1.bar(v0, hight, width, color=key, align='edge', alpha=.8)
                ax1.bar(v0, hight, width, color='none', edgecolor=key, align='edge', alpha=.85, hatch='|', zorder=0)


    plt.tight_layout()

#    plt.show()
    plt.savefig(fig_name+'.png')


if __name__ == '__main__':
    plot_wave_energy_frequency()
