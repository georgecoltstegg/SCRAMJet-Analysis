# Project 2 for AERO 351 Fall 2021
# Texas A&M University
# George J. Colts-Tegg
# UIN: 627005116
# Professor Paul G. Cizmas

from math import *
import matplotlib.pyplot as plt
import numpy as np

# Problem Statement:
# A scramjet engine flies at altitude H = 10, 000 + A m. The inlet diffuser losses are
# σ0d = 70 + B%, the combustor losses are σ0c = 90 + C%, the combustor efficiency is
# ξcomb = 96%, and the nozzle losses are σ0n = 92 + D%. One can assume the engine is
# using standard fuel.

# Part A
# Assuming the ratio of specific heat capacities is constant, γ = 1.4, calculate the
# specific thrust, the thrust specific fuel consumption, the propulsion, thermal and
# overall efficiencies, and the fuel to air ratio of the engine operating between Mach
# 2 and 5, at T03 = 2400 K. Plot the variation of these quantities vs. Mach number
# using a Mach number increment of 0.2 or less.

# UIN parameters
a = 5
b = 1
c = 1
d = 6

A = 100*a
B = 0.5*b
C = 0.4*c
D = 0.5*d

M = np.arange(2, 5.05, 0.05)
thrust_array = []
TSFC_array = []
propEff_array = []
thermEff_array = []
overEff_array = []
f_array = []
# assuming p_a / p_5 = 1


def Scramjet(diff_loss, comb_loss, nozz_loss, T_03, T_a, gam, r, lhv, comb_eff, c_p, p_a, a):
    for i in range(len(M)):
        # specific thrust
        p0a = p_a*(1+((gam - 1)/2)*(M[i]**2))**(gam/(gam-1))
        T_0a = T_a*(1+((gam-1)/2)*(M[i]**2))
        u_inlet = M[i]*a
        beta = (1 + ((gam - 1)/2)*(M[i]**2))*(((diff_loss*comb_loss*nozz_loss)*1)**((gam-1)/gam))
        # M_e = sqrt((2/(gam-1))*(beta - 1))
        u_e = sqrt(((2*gam)/(gam-1))*((beta-1)/beta)*r*T_03)
        f = ((T_03/T_0a) - 1)/((comb_eff*lhv)/(c_p*T_0a - T_03/T_0a))
        thrust = (1+f)*u_e*(1+(((gam-1)/(2*gam))*(1/(beta-1))*(1-1))) - M[i]*sqrt(gam*r*T_a)
        thrust_array.append(thrust)
        f_array.append(f)
        # TSFC
        TSFC = f/thrust
        TSFC_array.append(TSFC)
        # propulsion efficiency
        prop_eff = (thrust*u_inlet)/((1+f)*((u_e**2)/2)-((u_inlet**2)/2))
        propEff_array.append(prop_eff)
        # thermal efficiency
        therm_eff = (((1+f)*((u_e**2)/2))-((u_inlet**2)/2))/(f*lhv)
        thermEff_array.append(therm_eff)
        # overall efficiency
        over_eff = therm_eff * prop_eff
        overEff_array.append(over_eff)
    return thrust_array, TSFC_array, propEff_array, thermEff_array, overEff_array, f_array


diffLoss = (70+B)/100
combLoss = (90+C)/100
combEff = 96/100
nozzLoss = (92+D)/100
T03 = 2400
Ta = 220.67  # at altitude of 10500m
pa = 24540  # std atmosphere
a = 297.4  # speed of sound standard atmosphere
Scramjet(diffLoss, combLoss, nozzLoss, T03, Ta, 1.4, 287.16, 43500000, combEff, 1004, pa, a)

# Plotting
figure, axis = plt.subplots(2, 2)
# specific thrust vs mach
axis[0, 0].plot(M, thrust_array, color='blue')
axis[0, 0].set_xlabel('Mach Number, M')
axis[0, 0].set_ylabel('Specific Thrust [N-s/kg]')
axis[0, 0].set_title('Specific Thrust vs. Mach Number')
# TSFC vs mach
axis[0, 1].plot(M, TSFC_array, color='blue')
axis[0, 1].set_xlabel('Mach Number, M')
axis[0, 1].set_ylabel('TSFC [kg/(N-s)]')
axis[0, 1].set_title('TSFC vs. Mach Number')
# Efficiencies vs mach
axis[1, 0].plot(M, propEff_array, color='blue', label='Propulsion Efficiency')
axis[1, 0].plot(M, thermEff_array, color='red', label='Thermal Efficiency')
axis[1, 0].plot(M, overEff_array, color='orange', label='Overall Efficiency')
axis[1, 0].set_xlabel('Mach Number, M')
axis[1, 0].set_ylabel('Efficiencies')
axis[1, 0].set_title('Efficiencies vs. Mach Number')
axis[1, 0].legend()
# fuel to air ratio vs mach
axis[1, 1].plot(M, f_array, color='blue')
axis[1, 1].set_xlabel('Mach Number, M')
axis[1, 1].set_ylabel('Fuel to Air Ratio, f')
axis[1, 1].set_title('Fuel to Air Ratio vs. Mach Number')
plt.show()

# Part B
# Calculate the specific thrust, the thrust specific fuel consumption, the propulsion,
# thermal and overall efficiencies, and the fuel to air ratio of the engine operating
# at Mach 5 with T03 = 2400 K using variable γ.

# givens
diffLoss = (70+B)/100
combLoss = (90+C)/100
combEff = 96/100
nozzLoss = (92+D)/100
T03 = 2400
R = 287.16
minL = 14.66
lhv = 43500000
M = 5

# state a (inlet)
Ta = 220.67  # std atmosphere at 10500 m
pa = .24540  # bar std atmosphere at 10500 m
a = 297.4  # std atmosphere at 10500 m
ha = 220  # air tables
sa = 6.3927 - .28716*log(pa)
u = M * a

# state 0a
h0a = ha + (u**2)/(2*1000)
s0a = sa
s0a_prime = 8.219
p0a = exp(-(s0a - s0a_prime)/.28716)

# state 02
p02 = p0a * diffLoss
h02 = h0a
s02 = s0a - .28716*log(diffLoss)

# state 03
p03 = p02 * combLoss
h03_air = 2784.792  # interpolated from tables
h03_lam = 3031.476  # interpolated from tables
lam = (((h03_lam*(1+minL)) - (combEff*lhv) - (h03_air*minL))/(minL*(h02-h03_air)))/1000
r = (1+minL)/(1+lam*minL)
q = ((lam - 1)*minL)/(1+lam*minL)
h03 = (r*h03_lam + q*h03_air)

# state 05
h05 = h03
p05 = diffLoss*combLoss*nozzLoss*p0a
s05_lam = 9.423  # interpolated from tables
s05_air = 9.159  # interpolated from tables
s05_comb = (r*s05_lam + q*s05_air)
s05 = s05_comb - .28716*log(p05)

# state 5
s5 = s05
p5 = pa
s5_prime = s5 + .28716*log(p5)
h5 = 744.6
u_e = sqrt((h05-h5)*2000)

# fuel to air ratio
f = 1/(lam*minL)
# specific thrust
T_specific = u_e*(1+f) - u
# thrust specific fuel consumption (TSFC)
TSFC_b = f/T_specific
# thermal efficiency
n_th = (((1+f)*((u_e**2)/2))-((u**2)/2))/(f*lhv)
# propulsion efficiency
n_prop = (T_specific*u)/((1+f)*((u_e**2)/2)-((u**2)/2))
# overall efficiency
n_o = n_th * n_prop

print('Part B\nFuel to air ratio:', f, '\nSpecific Thrust:', T_specific, '\nTSFC:', TSFC_b, '\nThermal Efficiency:', n_th
      , '\nPropulsion Efficiency:', n_prop, '\nOverall Efficiency:', n_o)
