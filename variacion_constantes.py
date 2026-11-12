import argparse
import math
import os

import numpy as np

import tp_el_nino_parte2a as tp

import matplotlib.pyplot as plt

# --------------------------------------------------------------------------
# Datos del problema (valores nominales)
# --------------------------------------------------------------------------
PARAMS_NOMINALES = dict(
    Tmax=4.0,      # anomalía máxima sostenible           [°C]
    Roce=0.30,     # tasa intrínseca de calentamiento     [1/mes]
    alpha=0.05,    # enfriamiento por evaporación
    beta=0.20,     # eficiencia de evaporación
    gamma=0.10,    # consumo de humedad por precipitación
    delta=0.40,    # conversión de humedad a lluvia
    sigma=0.70,    # descarga / disipación                [1/mes]
    A=0.45,        # amplitud de la modulación estacional [1/mes]
)
OMEGA = 2.0 * math.pi / 12.0          # rad/mes (ciclo anual)
Y0 = (1.0, 0.0, 0.0)                  # T(0), H(0), P(0)
T_FIN = 60.0                          # meses
T_ULTIMO_CICLO = 48.0                 # el último ciclo es 48 <= t <= 60
PERIODO = 12.0                        # meses

H_DEFECTO = 0.001                     # paso de Euler usado en la Parte 2 [mes]

# --------------------------------------------------------------------------
# Modelo e integradores
# --------------------------------------------------------------------------


def rhs(t, T, H, P, p):
    """Lado derecho f(t, y) del sistema; devuelve (dT/dt, dH/dt, dP/dt)."""
    R = p["Roce"] + p["A"] * math.cos(OMEGA * t - math.pi)
    dT = R * T * (1.0 - T / p["Tmax"]) - p["alpha"] * H
    dH = p["beta"] * T * T - p["gamma"] * H * P
    dP = p["delta"] * H * H - p["sigma"] * P
    return dT, dH, dP


def euler(p, h, t_fin=T_FIN, y0=Y0):
    """Euler explícito: y_{n+1} = y_n + h f(t_n, y_n). Devuelve (t, y[n,3])."""
    n = int(round(t_fin / h))
    y = np.empty((n + 1, 3))
    T, H, P = y0
    y[0] = y0
    for k in range(n):
        dT, dH, dP = rhs(k * h, T, H, P, p)
        T, H, P = T + h * dT, H + h * dH, P + h * dP
        y[k + 1] = (T, H, P)
    return np.arange(n + 1) * h, y

# --------------------------------------------------------------------------
# ESTIMACION MAXIMO DE SEVERIDAD
# --------------------------------------------------------------------------

"""SP = []
colores_plt = []
colores = [["red", "green", "blue"], ["orange", "purple", "brown"], ["cyan", "magenta", "olive"]]
valores = []

for i in range(7):
    PARAMS_NOMINALES["A"] = 0.25 + i * 0.05833334
    for j in range(3):
        PARAMS_NOMINALES["sigma"] = 0.4 + j * 0.3
        for k in range(3):
            PARAMS_NOMINALES["gamma"] = 0.05 + k * 0.075
            for l in range(3):
                PARAMS_NOMINALES["beta"] = 0.1 + l * 0.1

                t, y = euler(PARAMS_NOMINALES, H_DEFECTO)

                indicadores = tp.indicadores(t, y)
                SP.append(indicadores["SP"])
                
                colores_plt.append(colores[j][k])

                if colores_plt[-1] == "cyan":
                    valores.append([f"{PARAMS_NOMINALES['beta']:.2f}", f"{SP[-1]:.2f}"])

print(valores)


from matplotlib.collections import LineCollection

fig1 = plt.figure()

x = np.arange(len(SP))
puntos = np.array([x, SP]).T.reshape(-1, 1, 2)
segmentos = np.concatenate([puntos[:-1], puntos[1:]], axis=1)

lc = LineCollection(segmentos, colors=colores_plt[:-1], linewidths=1.5)

ax = plt.gca()
ax.add_collection(lc)

ax.set_xlim(x.min(), x.max())
ax.set_ylim(min(SP), max(SP))

plt.grid()
plt.savefig("grafico_Tot.png", dpi=300, bbox_inches="tight")

"""


# --------------------------------------------------------------------------
# CALCULO SP Y dtTP VARIANDO LAS CONSTANTES A sigma, gamma y beta
# --------------------------------------------------------------------------


SP_A = []
dtTP_A = []
SP_sigma = []
dtTP_sigma = []
SP_gamma = []
dtTP_gamma = []
SP_beta = []
dtTP_beta = []

for i in range(70):
    PARAMS_NOMINALES["A"] = 0.25 + i * 0.005
    t, y = euler(PARAMS_NOMINALES, H_DEFECTO)

    indicadores = tp.indicadores(t, y)
    SP_A.append(indicadores["SP"])
    dtTP_A.append(indicadores["dtTP"])

for i in range(60):
    PARAMS_NOMINALES["sigma"] = 0.4 + i * 0.01
    t, y = euler(PARAMS_NOMINALES, H_DEFECTO)

    indicadores = tp.indicadores(t, y)
    SP_sigma.append(indicadores["SP"])
    dtTP_sigma.append(indicadores["dtTP"])

for i in range(150):
    PARAMS_NOMINALES["gamma"] = 0.05 + i * 0.001
    t, y = euler(PARAMS_NOMINALES, H_DEFECTO)

    indicadores = tp.indicadores(t, y)
    SP_gamma.append(indicadores["SP"])
    dtTP_gamma.append(indicadores["dtTP"])


for i in range(40):
    PARAMS_NOMINALES["beta"] = 0.1 + i * 0.005
    t, y = euler(PARAMS_NOMINALES, H_DEFECTO)

    indicadores = tp.indicadores(t, y)
    SP_beta.append(indicadores["SP"])
    dtTP_beta.append(indicadores["dtTP"])


fig1 = plt.figure()
plt.plot(SP_A, label="SP")
plt.plot(dtTP_A, label="dtTP")
plt.legend()
plt.grid()
plt.savefig("grafico_A.png", dpi=300, bbox_inches="tight")

fig2 = plt.figure()
plt.plot(SP_sigma, label="SP")
plt.plot(dtTP_sigma, label="dtTP")
plt.legend()
plt.grid()
plt.savefig("grafico_sigma.png", dpi=300, bbox_inches="tight")

fig3 = plt.figure()
plt.plot(SP_gamma, label="SP")
plt.plot(dtTP_gamma, label="dtTP")
plt.legend()
plt.grid()
plt.savefig("grafico_gamma.png", dpi=300, bbox_inches="tight")

fig4 = plt.figure()
plt.plot(SP_beta, label="SP")
plt.plot(dtTP_beta, label="dtTP")
plt.legend()
plt.grid()
plt.savefig("grafico_beta.png", dpi=300, bbox_inches="tight")

plt.show()