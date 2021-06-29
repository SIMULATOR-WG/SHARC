#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 14 15:04:04 2021

Parâmetros orbitais de simulação para o SHARC-LEO
(variam conforme o âgulo de elevação)
Saídas: imt_long_diff_deg, intersite_distance e fator de segmento S

Ref.:
    M. Richharia
    Satellite Communication Systems - Design Principles
    Macmillan Education UK (1999)

    Luciano Camilo Alexandre
    Maximum line of sigh view (horizon view)
    between satellite/pseudo-satellite and Earth considering Earth curvature
    (Python script)

@author: leandropires
"""

import numpy as np

from obspy.geodetics import degrees2kilometers, kilometers2degrees

## Parametros gerais

h = 590.0 # Altitude do satelite em km
eta = 20 # Angulo de elevacao em graus (modificar para diferentes parâmetros)
gamma = 127.13/2 # Angulo de nadir = Angulo de visibilidade sobre 2
theta_e = 0.0 # Latitude da estacao em graus
theta_s = 0.0 # Latitude do satelite em graus
# theta_e = -15.809422 # Latitude da estacao em graus
# theta_s = -15.809422 # Latitude do satelite em graus
R = 6378.0 # Raio da Terra em km

# Area de servico maxima em km^2
    # Inflexão em 53°, para a mesma visibilidade na qual
    # A = área do Brasil quando eta = 0
A_max = 22249828

# Diferenca de longitude maxima ate o horizonte
    # script longDiff.py (Luciano)
longDiffMax = 25

## Parametros de lancamento da malha IMT
BS_density_Urban = 30 # [BS/km**2]
BS_density_Subrban = 10 # [BS/km**2]
Ra_Urban = 0.07
Ra_Suburban = 0.03
Rb = 0.05

n_clusters = 7
n_sectors = 57
n_hotspots_sector = 3

# Area de servico maxima em km^2
    # Inflexão em 53°, para a mesma visibilidade na qual
    # A = área do Brasil quando eta = 0
if eta > 53:
    ## Calculos angulares
    H = h / np.sin(np.deg2rad(eta)) # altitude corrigida pela elevacao
    SEC = 180 - np.rad2deg(np.arcsin((H+R)*np.sin(np.deg2rad(gamma))/R)) #Satelite-Earth-Centre
    ECS = 180 - SEC - gamma #Earth-Centre-Satelite
    A_s = 2 * np.pi * R**2 * (1-np.cos(np.deg2rad(ECS)))
else:
    A_s = A_max

## Distância entre sítios (centros de dois hexagonos adjacentes)
A_sector = (A_s / n_clusters) / n_sectors # Area do setor em km**2
ISD = 3 * np.sqrt(2 * A_sector / (3 * np.sqrt(3))) # Distancia em km


## Fator de segmento
# deploymentRatio em [BS/km**2]
deploymentRatio = Rb * (BS_density_Urban * Ra_Urban + BS_density_Subrban * Ra_Suburban)
N_BS = deploymentRatio * A_s # numero de estacoes base
# numero de estacoes de base por lancamento da simulacao de Monte-Carlo
n_hotspots_snapshot = n_clusters * n_sectors * n_hotspots_sector
# Fator de segmento
S_factor = int(np.ceil(N_BS / n_hotspots_snapshot))


# Calculo da constante sigma (Eq. 217.a)
sigma = R/(R + h)

# Angulo de cobertura psi (Eq. 2.17b)
psi = 90 - eta - np.rad2deg(np.arcsin(np.cos(np.deg2rad(eta)) / sigma**(-1)))

# Diferenca de longitude \phi_{es} (Eq. 2.17c adaptada)
# Se theta_e = 0, longDiff = psi
longDiff = np.rad2deg((np.arccos(np.cos(np.deg2rad(psi)))/(np.cos(np.deg2rad(theta_e)))))

# Cálculo das distâncias: raio IMT, distancia estacao-satelite e linha do horizonte
# Mesma latitude satélite e estação terrena
sys_lat_rad = np.deg2rad(theta_e)
sat_lat_rad = sys_lat_rad
distancia = R * np.arccos(np.sin(sat_lat_rad) * np.sin(sys_lat_rad ) + np.cos(sat_lat_rad) * np.cos(sys_lat_rad ) * np.cos(np.deg2rad(longDiff)))
# ou
distancia = degrees2kilometers(longDiff, radius=R)
#print(distancia)
raioIMT = np.sqrt(A_s/np.pi)
horizonte = R * np.arccos(np.sin(sat_lat_rad) * np.sin(sys_lat_rad ) + np.cos(sat_lat_rad) * np.cos(sys_lat_rad ) * np.cos(np.deg2rad(longDiffMax)))
# ou
horizonte = degrees2kilometers(longDiffMax, radius = R)

if distancia + raioIMT > horizonte:
    longDiffdesc = longDiff - kilometers2degrees((distancia + raioIMT - horizonte), radius=R)
else:
    longDiffdesc = longDiff

# Print dos parâmetros orbitais do SHARC
print('Parametros orbitais do SHARC SS_MLEO para h = ' + str(h) + ' km e elevação = ' + str(eta) + ' graus')
print(f'Parametros orbitais do SHARC SS_MLEO para h = {h} km e elevação = {eta} graus')
print('imt_long_diff_deg = ' + str(longDiffdesc) + ' graus')
print('intersite_distance = ' + str(ISD*1000) + ' metros')
print('fator de segmento S = ' + str(S_factor))








