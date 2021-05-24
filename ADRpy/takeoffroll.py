#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
.. _takeoffroll_module:

Takeoff roll module
--------------------

This module contains tools for calculating the takeoff roll distance of an aircraft
both on land and water. It also contains tools to help find float displacements given
a simple float design.

Original code seed by Tom Smith

"""

__author__ = "Zach Tait"

import matplotlib.pyplot as plt
import numpy as np


class Floatplane:

    def __init__(self, MTOW, wing_A, concept, floats):
        self.m = MTOW  # Mass
        self.W = self.m * 9.8  # Weight
        self.B = floats['beam']  # Float Beam
        self.Len = floats['length']  # Float Length
        self.rho = 1.2  # Air Density
        self.rho_w = 1025  # Water Density
        self.C_l = concept.performance['CLmaxTO']
        self.C_d = concept.performance['CDTO']
        self.S_T = wing_A
        self.a_t1 = floats['trim_start']  # Starting Trim Angle
        self.a_t2 = floats['trim_end']  # End Trim Angle
        self.AoI = 4.5  # Wing angle of incidence with Keel

    def TakeOffRun(self, t_int, Thrust):
        v = 0
        C_buoy0 = self.W / (self.rho_w * self.B ** 3)  # Stationary Coefficient of Buoyancy
        S_wet = self.Len * self.B  # Wet Surface Area
        d = 0  # Starting Distance
        onwater = True
        t = 0  # starting time
        t_s, d_s, v_s, L_s, D_s, R_s, R_frs = [], [], [], [], [], [], []
        while onwater is True:
            # Calculate all the forces
            t_s.append(t)
            d_s.append(d)
            v_s.append(v)
            L = 0.5 * self.rho * v * v * self.S_T * self.C_l
            L_s.append(L)
            D = 0.5 * self.rho * v * v * self.S_T * self.C_d
            D_s.append(D)
            C_v = v / (9.8 * self.B) ** 0.5
            C_r = 0.0011 * C_v ** 3 - 0.0221 * C_v ** 2 + 0.1062 * C_v - 0.0149
            C_v2 = C_v * (self.B / self.Len) ** (1 / 3)
            if C_v == 0:
                a_trim = self.a_t1
            else:
                A_t = 5.294 / (C_v2 - C_v)
                B_t = -(2.647 + A_t * C_v)
                a_trim = self.a_t1 + ((self.a_t1 + self.a_t2) / 2) * (1 + np.tanh(A_t * C_v + B_t))
            # AoA = a_trim + AoI
            # Can do lift curve slope stuff here
            Buoy = self.W - L
            C_buoy = Buoy / (self.rho_w * self.B ** 3)
            # C_rbuoy = C_r * (C_buoy / C_buoy0)
            # R = 2*C_rbuoy*(rho_w*B**3)
            # if C_v>7.5:
            #    R = 0
            # R_fr = 2*0.012*S_wet*v**2
            # R_frs.append(R_fr)
            if v < 30:
                R = 3000
            else:
                R = 0.0103 * v ** 3 - 0.0686 * v ** 2 + 3.101 * v - 2.392
            R = 2 * R
            R_s.append(R)
            F_net = Thrust - D - R  # - R_fr
            a = F_net / self.m
            d = d + v * t_int + a * t_int ** 2
            v = v + a * t_int
            t = t + t_int
            if t > 100:
                break
            if Buoy < 0:
                onwater = False

        W_s = []
        for i in range(len(t_s)):
            W_s.append(self.W)
        water_to = {'time': t_s,
                    'velocity': v_s,
                    'lift': L_s,
                    'drag': D_s,
                    'resistance': R_s,
                    'weight': W_s}
        return water_to

    # Flt_Geo = {'L': 0.7, 'B': 0.2, 'Dead': 40, 'Req_mass': 5}
    # Flt_Geo['L'] / Flt_Geo['B']
    #
    # Flt_disp = (Flt_Geo['Req_mass'] * 9.81) / (9.81 * 1025)
    # A_req = Flt_disp / Flt_Geo['L']
    # Wtr_line = ((A_req) / np.tan(np.pi / 2 - (Flt_Geo['Dead'] * np.pi) / 180)) ** 0.5
    # Wtr_b = Wtr_line * np.tan(np.pi / 2 - (Flt_Geo['Dead'] * np.pi) / 180)
    # H = Flt_Geo['B'] / (2 * np.tan(np.pi / 2 - (Flt_Geo['Dead'] * np.pi) / 180))
    # Flt_x = [0, -1 * Flt_Geo['B'] / 2, Flt_Geo['B'] / 2, 0]
    # Flt_y = [0, H, H, 0]
    # Wtr_x = [0, -1 * Wtr_b, Wtr_b, 0]
    # Wtr_y = [0, Wtr_line, Wtr_line, 0]
    # plt.plot(Flt_x, Flt_y, color='red', label='Float Hull')
    # plt.plot(Wtr_x, Wtr_y, color='blue', label='Water Line')
    # plt.legend(loc='lower right')
    # plt.grid()
    # print(H, Flt_x[2])


class Wheels:
    def TakeOffRun(t_int, T):
        m = 7  # Mass
        W = m * 9.8  # Weight
        t = 0  # Starting Time
        v = 0  # Starting Velocity
        rho = 1.225  # Air Density
        rho_w = 1025  # Water Density
        wind_speed = 0  # Wind Speed
        AR = 6.686
        e = 0.961
        S = 1
        a_t1 = 6 * (np.pi / 180)  # Starting Trim Angle

        def C_l(a):
            return 4.745 * a + 0.33

        C_l_0 = C_l(a_t1)
        C_D0 = 0.02
        C_Di = C_l_0 ** 2 / (np.pi * AR * e)
        d = 0  # Starting Distance
        t_s, d_s, v_s, L_s, D_s, F_R_s = [], [], [], [], [], [],
        while v < 15:
            # Calculate all the forces
            t_s.append(t)
            d_s.append(d)
            v_s.append(v)
            C_l = 4.745 * a_t1 + 0.33
            L = 0.5 * rho * v * v * S * C_l
            L_s.append(L)
            C_Di = C_l ** 2 / (np.pi * AR * e)
            C_d = C_Di + C_D0
            D = 0.5 * rho * v * v * S * C_d
            D_s.append(D)
            F_R = 0.1 * (W - L)  # Rolling Resistance
            if F_R < 0:
                F_R = 0
            F_R_s.append(F_R)
            F_net = T - D - F_R
            a = F_net / m
            d = d + v * t_int + a * t_int ** 2
            v = v + a * t_int
            t = t + t_int
            if t > 100:
                break

        W_s = []
        for i in range(len(t_s)):
            W_s.append(W)

        plt.plot(t_s, v_s, label='Velocity [m/s]')
        plt.plot(t_s, L_s, label='Lift [N]')
        plt.plot(t_s, D_s, label='Drag[N]')
        plt.plot(t_s, F_R_s, label='Rolling Resistance Resistance[N]')
        plt.plot(t_s, W_s, label='Weight[N]')
        plt.xlabel('Time[s]')

        # plt.plot(v_s, D_s, label='Drag[N]')
        # plt.plot(v_s, R_s, label='Water Resistance[N]')
        # plt.plot(v_s, R_frs, label='Hydroplaning Resistance [N]')
        # plt.xlabel('Velocity[m/s]')

        plt.legend()
        plt.grid()
        plt.title('Take Off Run with a Thrust of ' + str(T) + 'N')
        plt.show()
