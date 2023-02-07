import json
import numpy as np
import matplotlib.pyplot as plt
from astropy import constants as con

"""
Author: Zhenbo SU @ USTC
Mail: zbsu@mail.ustc.edu.cn
"""

con_hp = 6.62606957e-27
con_kb = 1.3806488e-16
con_c = 29979245800.0
G = con.G.cgs.value
mdot_ucvt = con.M_sun.cgs.value / (365*24*60*60)

w_start = 1
w_end = 10000
wavelength = np.linspace(w_start, w_end, w_end - w_start + 1)
freq = con_c / (wavelength * 1.E-8)

plt.rcParams.update({
    'font.size': 12, 
    'font.family': 'Times',
    'axes.linewidth': 0.8,
    'axes.spines.bottom': True,
    'axes.spines.left': True,
    'axes.spines.right': True,
    'axes.spines.top': True,
    'xtick.direction': 'in',
    'xtick.top': True,
    'ytick.right': True,
    'ytick.direction': 'in',
    'xtick.major.size': 6,
    'xtick.major.width': 0.8,
    'ytick.major.size': 6,
    'ytick.major.width': 0.8,
    'xtick.minor.size': 3,
    'xtick.minor.width': 0.6,
    'ytick.minor.size': 3,
    'ytick.minor.width': 0.6,
    "text.usetex": True,
    "axes.titlesize":"medium",
    "figure.dpi" : 150
})

def blackbody_cgs(wavelength_in_bb, temp):
    # wavelength_in_bb: AA
    # temp : K
    wavelength_in_bb = wavelength_in_bb * 1.E-8
    tmp = con_hp * con_c / (con_kb * temp)
    with np.errstate(over='ignore'):
        sed = 2 * con_hp * con_c ** 2 / (wavelength_in_bb ** 5) / (
            np.expm1(tmp / wavelength_in_bb))  # [erg/s/cm^2/cm]
    return sed * 1.E-8  # [erg/s/cm^2/A]

def rms(aspin):
    z1 = 1 + np.power(1 - aspin ** 2, 1/3) * (np.power(1 + aspin, 1/3) + np.power(1 - aspin, 1/3))
    z2 = np.sqrt(3 * np.power(aspin, 2) + np.power(z1, 2))
    rms = 0.5 * (3 + z2 - np.sqrt((3-z1)*(3+z1+2*z2))) # in rg.
    return rms

class WindDisk:
    
    def __init__(self, mbh, mdot, rin, fracring, rratio, eps=0, grav=False) -> None:
        
        self.mbh = mbh # in solar mass
        self.disk_mdot = mdot # in Mun/yr 
        self.rin = rin # in Rg
        self.fracring = fracring
        self.rratio = rratio
        self.eps = eps # TODO
        self.grav = grav # if grav, \sum(F, g).
        
        self.bhmass = self.mbh * con.M_sun.to('g')
        self.bhmass = self.bhmass.value
        self.rgravcm = con.G.cgs.value * self.bhmass / (con.c.cgs.value ** 2)
        self.wavelength = np.linspace(w_start, w_end, w_end - w_start + 1)
        self.freq = con_c / (self.wavelength * 1.E-8)
        
        # plt.plot(self.radius, teff)
        # plt.xscale('log')
        # plt.yscale('log')
        # plt.show()

    def winddisk(self, disk_mdot_in, calsed=False, plot=False):
        
        # radius = np.empty((self.fracring, 4))
        # radius.fill(0.0)
        
        rin_zones = np.zeros(self.fracring)
        rout_zones = np.zeros(self.fracring)
        mdot = np.zeros(self.fracring)
        wrf = np.zeros(self.fracring)
        
        flux = np.zeros(self.fracring)
        wrf[0] = 0
        mdot[0] = disk_mdot_in


        rin_zones[0] = self.rin
        
        rout_zones[0] = self.rin * self.rratio

        for i in range(self.fracring - 1):
            rin_zones[i+1] = rin_zones[i] * self.rratio
            rout_zones[i+1] = rout_zones[i] * self.rratio
        
        rmid_zones = np.sqrt(rin_zones * rout_zones)
        # print(rmid_zones)

        area_zones = np.pi * (rout_zones ** 2 - rin_zones ** 2) * (self.rgravcm ** 2)
        
        self.radius = rmid_zones
        radius_cm = rmid_zones * self.rgravcm
        # print(radius_cm)
        if self.grav:
            for i in range(self.fracring - 1):
                
                flux[i] = 0.75 * wrf[i] * (G * self.bhmass / radius_cm[i] ** 3) ** 0.5

                mdot[i+1] = mdot[i] + (radius_cm[i+1] - radius_cm[i]) * (4 * np.pi * radius_cm[i]) * 3.05e-9 * (flux[i] / 6.3e10) ** 1.21 / mdot_ucvt

                wrf[i+1] = wrf[i] + (radius_cm[i+1] - radius_cm[i]) * (mdot[i] * mdot_ucvt / (4 * np.pi * radius_cm[i] ** 2) * (G * self.bhmass / radius_cm[i]) ** 0.5 - 2 * wrf[i] / radius_cm[i])

        else:
            for i in range(self.fracring - 1):
                flux[i] = 0.75 * wrf[i] * (G * self.bhmass / radius_cm[i] ** 3) ** 0.5
                
                mdot[i+1] = mdot[i] + (radius_cm[i+1] - radius_cm[i]) * (4 * np.pi * radius_cm[i]) * 2.6e-12 * (flux[i] / 6.3e10) ** 1.9 / mdot_ucvt
                
                wrf[i+1] = wrf[i] + (radius_cm[i+1] - radius_cm[i]) * (mdot[i] * mdot_ucvt / (4 * np.pi * radius_cm[i] ** 2) * (G * self.bhmass / radius_cm[i]) ** 0.5 - 2 * wrf[i] / radius_cm[i])     
        
        teff = (flux / 5.67e-5) ** 0.25
        mdot[np.isnan(mdot)] = 9999.0
        
        if plot:
            fig = plt.figure(figsize=(12, 4))
            plt.plot(self.radius, mdot)
            # plt.plot()
            plt.xscale('log')
            plt.yscale('log')
            plt.show()

        if calsed:
            self.sed = 0 * self.wavelength
            for i_area in range(np.size(area_zones)):
                tmp_sed = area_zones[i_area] * blackbody_cgs(self.wavelength, teff[i_area])
                self.sed += tmp_sed

            return [self.radius, mdot, wrf, flux, teff, area_zones, wavelength, self.sed]
        else:
            return [self.radius, mdot, wrf, flux, teff, area_zones]


        # return mdot[-1]


    def bisection(self, delta = 1e-8):
        
        init_mdot = [0, self.disk_mdot]
        left = init_mdot[0]
        right = init_mdot[1]
        mid = (left + right) / 2
        count = 0
        
        while (np.abs(self.winddisk(mid)[1][-1] - self.disk_mdot) > delta):
            if self.winddisk(mid)[1][-1] < self.disk_mdot:
                left = mid
            elif self.winddisk(mid)[1][-1] > self.disk_mdot:
                right = mid
            elif self.winddisk(mid)[1][-1] == self.disk_mdot:
                break

            count+=1
        
            mid = (right + left) / 2
        
            print("iteration: ", count, "accretion rate: ", mid)

        self.disk_mdot_in = mid

        return self.disk_mdot_in

class ThinDisk:
    
    def __init__(self, mbh, mdot, rin, fracring, rratio, eps=0, grav=False) -> None:
        
        self.mbh = mbh # in solar mass
        self.disk_mdot = mdot # in Mun/yr 
        self.rin = rin # in Rg
        self.fracring = fracring
        self.rratio = rratio
        self.eps = eps # TODO
        self.grav = grav # if grav, \sum(F, g).
        
        self.bhmass = self.mbh * con.M_sun.to('g')
        self.bhmass = self.bhmass.value
        self.rgravcm = con.G.cgs.value * self.bhmass / (con.c.cgs.value ** 2)
        self.wavelength = np.linspace(w_start, w_end, w_end - w_start + 1)
        self.freq = con_c / (self.wavelength * 1.E-8)

    def thindisk(self, calsed=True):
        rin_zones = np.zeros(self.fracring)
        rout_zones = np.zeros(self.fracring)
        wrf = np.zeros(self.fracring)
        
        flux = np.zeros(self.fracring)
        wrf[0] = 0

        rin_zones[0] = self.rin
        
        rout_zones[0] = self.rin * self.rratio
        for i in range(self.fracring - 1):
            rin_zones[i+1] = rin_zones[i] * self.rratio
            rout_zones[i+1] = rout_zones[i] * self.rratio
        
        rmid_zones = np.sqrt(rin_zones * rout_zones) 
        
        area_zones = np.pi * (rout_zones ** 2 - rin_zones ** 2) * (self.rgravcm **2 )
        
        self.radius = rmid_zones
        radius_cm = rmid_zones * self.rgravcm
        # print(radius_cm)
        self.area = area_zones

        for i in range(self.fracring - 1):
            
            flux[i] = 0.75 * wrf[i] * (G * self.bhmass / radius_cm[i] ** 3) ** 0.5
            wrf[i+1] = wrf[i] + (radius_cm[i+1] - radius_cm[i]) * (self.disk_mdot * mdot_ucvt / (4 * np.pi * radius_cm[i] ** 2) * (G * self.bhmass / radius_cm[i]) ** 0.5 - 2 * wrf[i] / radius_cm[i]) 
        
        self.teff = (flux / 5.67e-5) ** 0.25
        
        if calsed:

            self.sed = 0 * self.wavelength
            for i_area in range(np.size(area_zones)):
                tmp_sed = area_zones[i_area] * blackbody_cgs(self.wavelength, self.teff[i_area])
                self.sed += tmp_sed

        return self.teff


    def grthindisk(self, aspin=0, calsed=True):
        rin_zones = np.zeros(self.fracring)
        rout_zones = np.zeros(self.fracring)
        wrf = np.zeros(self.fracring)
        
        flux = np.zeros(self.fracring)
        wrf[0] = 0

        rin_zones[0] = self.rin
        
        rout_zones[0] = self.rin * self.rratio
        for i in range(self.fracring - 1):
            rin_zones[i+1] = rin_zones[i] * self.rratio
            rout_zones[i+1] = rout_zones[i] * self.rratio
        
        rmid_zones = np.sqrt(rin_zones * rout_zones) 
        
        area_zones = np.pi * (rout_zones ** 2 - rin_zones ** 2) * (self.rgravcm **2)
        
        self.radius = rmid_zones
        radius_cm = rmid_zones * self.rgravcm
        # print(radius_cm)
        self.area = area_zones
        
        
        # aspin = 0.2
        rms = self.rin
        
        C = 1 - 3/2/self.radius + 2 * aspin * np.sqrt(1 / (8 * np.power(self.radius, 3)))
        
        r1 = np.sqrt(2) * np.cos(1 / 3 * np.power(np.cos(aspin), -1) * np.pi/180 - np.pi / 3)
        r2 = np.sqrt(2) * np.cos(1 / 3 * np.power(np.cos(aspin), -1) * np.pi/180 + np.pi / 3)
        r3 = - np.sqrt(2) * np.cos(1 / 3 * np.power(np.cos(aspin), -1) * np.pi/180)
    
        # fr = 1 - np.sqrt(self.rin / self.radius)
        # fr = 1 / (C * np.sqrt(self.radius)) * (
        #     np.sqrt(self.radius) - np.sqrt(rms) - 3/2/np.sqrt(2) * aspin * np.log(np.sqrt(self.radius) / np.sqrt(rms)) - 3/2 * ((np.sqrt(r1) - aspin) ** 2 / (np.sqrt(r1) * (np.sqrt(r1) - np.sqrt(r2)) * (np.sqrt(r1) - np.sqrt(r3))) * np.log((np.sqrt(self.radius) - np.sqrt(r1)) / (np.sqrt(rms) - np.sqrt(r1)))) -     3/2 * ((np.sqrt(r2) - aspin) ** 2 / (np.sqrt(r2)*(np.sqrt(r2) - np.sqrt(r1))*(np.sqrt(r2)-np.sqrt(r3)))) * np.log((np.sqrt(self.radius) - np.sqrt(r2)) / (np.sqrt(rms) - np.sqrt(r2))) - 3/2 * ((np.sqrt(r3) - aspin) ** 2 / (np.sqrt(r3) * (np.sqrt(r3) - np.sqrt(r1)) * (np.sqrt(r3) - np.sqrt(r2)))) * np.log((np.sqrt(self.radius) - np.sqrt(r3)) / (np.sqrt(rms) - np.sqrt(r3)))
        # )

        self.fr =1 / (C * np.sqrt(self.radius)) * (
            np.sqrt(self.radius) - np.sqrt(rms) - 3/2/np.sqrt(2) * aspin * np.log(np.sqrt(self.radius/rms)) - 3/2 * (r1 - aspin) ** 2 / (r1 * (r1-r2) * (r1-r3)) * np.log((np.sqrt(self.radius) - r1)/(np.sqrt(rms) - r1)) - 3/2 * (r2-aspin) ** 2 / (r2 * (r2 - r1) * (r2 - r3)) * np.log((np.sqrt(self.radius)-r2) / (np.sqrt(rms) - r2)) - 3/2 * (r3-aspin) ** 2 / (r3 * (r3-r1) * (r3-r2)) * np.log((np.sqrt(self.radius) - r3) / (np.sqrt(rms) - r3))
        )

        flux = 3 / 8 / np.pi * (G * self.bhmass * self.disk_mdot * mdot_ucvt / radius_cm ** 3) * (self.fr)

        self.teff = (flux / 5.67e-5) ** 0.25
        
        if calsed:

            self.sed = 0 * self.wavelength
            for i_area in range(np.size(area_zones)):
                tmp_sed = area_zones[i_area] * blackbody_cgs(self.wavelength, self.teff[i_area])
                self.sed += tmp_sed

        return self.teff
