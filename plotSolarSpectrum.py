"""
Plot Solar Spectrum
Written By: Dean Keithly
Written On: 6/12/2019
"""
#Note: the military likely has access to low latency solar spectrum data

import pyspectral.blackbody as BB
import pyspectral.solar as SOLAR
from pyspectral.solar import TOTAL_IRRADIANCE_SPECTRUM_2000ASTM
import numpy as np
import matplotlib.pyplot as plt

wavelengths = np.linspace(start=250.*10.**-9.,stop=2500.*10.**-9., num=200) #in nm
#Inputs wavelength in meters
sunTemp = 5800. # this was basiclly from wikipedia. TODO:UPDATE
BBirradiance = BB.blackbody(wavel=wavelengths,temp=sunTemp)/10**9. # in W/m^2/nm
SOLARirradiance = SOLAR.SolarIrradianceSpectrum(SOLAR.TOTAL_IRRADIANCE_SPECTRUM_2000ASTM, dlambda=0.005)
#SOLARirradiance.wavelength is in um

fig = plt.figure(3654862)
plt.rc('axes',linewidth=2)
plt.rc('lines',linewidth=2)
plt.rcParams['axes.linewidth']=2
plt.rc('font',weight='bold')
ax1 = fig.add_subplot(1,1,1)
ax1.plot(wavelengths*10**9.,BBirradiance, color='k', linewidth=2, label='Black Body Irradiance @ 5800K')
ax1.set_xlabel('Wavelength (nm)', weight='bold')
ax1.set_ylabel(r'Black Body Radiance $(W/m^2/Sr/nm)$', weight='bold')
ax1.set_xlim(left=200.,right=1200.)
ax1.set_ylim(bottom=0.,top=28000.)
ax2 = ax1.twinx()
ax2.plot(SOLARirradiance.wavelength*10**3.,SOLARirradiance.irradiance*(4.*np.pi), color='r', linewidth=2, label='Solar Irradiance')
ax2.set_ylim(bottom=0.,top=28000.)
ax2.spines['right'].set_color('red') # setting the right side axis to red
ax2.xaxis.label.set_color('red')
ax2.tick_params(axis='y', colors='red')
ax2.set_ylabel(r'Solar Irradiance $(W/m^2/Sr/nm)$', weight='bold', color='red')

plt.subplots_adjust(left=0.15,right=0.85)


plt.show(block=False)
