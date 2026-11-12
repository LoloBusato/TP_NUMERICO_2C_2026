#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TP "El superniño" - CB051 Modelación Numérica (2do cuatrimestre 2026)

Núcleo común a todos los scripts: parámetros, modelo, integradores (Euler
explícito y RK4), indicadores del último ciclo anual y jacobiano.

Modelo:
    dT/dt = R(t) T (1 - T/Tmax) - alpha H
    dH/dt = beta T^2 - gamma H P
    dP/dt = delta H^2 - sigma P
    R(t)  = Roce + A cos(w t - pi),   w = 2 pi / 12
"""

import math
import os

import numpy as np

# np.trapezoid existe desde numpy 2.0; np.trapz en versiones anteriores
_trapecio = getattr(np, "trapezoid", None) or np.trapz

# --------------------------------------------------------------------------
# Datos del problema
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
OMEGA = 2.0 * math.pi / 12.0     # rad/mes (ciclo anual)
Y0 = (1.0, 0.0, 0.0)             # T(0), H(0), P(0)
T_FIN = 60.0                     # meses
T_ULTIMO_CICLO = 48.0            # último ciclo anual: 48 <= t <= 60
PERIODO = 12.0                   # meses

# Carpeta de salida común (figuras, tablas, resultados.json)
DIR_RESULTADOS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resultados")


def carpeta_resultados():
    os.makedirs(DIR_RESULTADOS, exist_ok=True)
    return DIR_RESULTADOS


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


def euler(p, h, t_fin=T_FIN, y0=Y0, limite=1e6):
    """
    Euler explícito: y_{n+1} = y_n + h f(t_n, y_n).

    Devuelve (t, y) con y de forma (n+1, 3) = columnas T, H, P.
    Si h no divide a t_fin, el último paso es más corto para terminar
    exactamente en t_fin. Si la solución supera `limite` en módulo
    (divergencia; los valores físicos son de orden 1-10) se corta la
    integración y el resto de la trayectoria se completa con NaN.
    """
    n = int(math.ceil(t_fin / h - 1e-9))
    h_ult = t_fin - (n - 1) * h
    y = np.full((n + 1, 3), np.nan)
    T, H, P = y0
    y[0] = y0
    Tmax, Roce, al, be, ga, de, si, A = (p["Tmax"], p["Roce"], p["alpha"],
                                         p["beta"], p["gamma"], p["delta"],
                                         p["sigma"], p["A"])
    cos = math.cos
    for k in range(n):
        hk = h_ult if k == n - 1 else h
        R = Roce + A * cos(OMEGA * (k * h) - math.pi)
        dT = R * T * (1.0 - T / Tmax) - al * H
        dH = be * T * T - ga * H * P
        dP = de * H * H - si * P
        T, H, P = T + hk * dT, H + hk * dH, P + hk * dP
        if not (abs(T) < limite and abs(H) < limite and abs(P) < limite):
            break
        y[k + 1] = (T, H, P)
    t = np.arange(n + 1) * h
    t[-1] = t_fin
    return t, y


def diverge(y):
    """True si la trayectoria tiene NaN (el integrador cortó por divergencia)."""
    return bool(np.isnan(y).any())


def rk4(p, h, t_fin=T_FIN, y0=Y0):
    """Runge-Kutta de orden 4. Solo se usa como verificación independiente."""
    n = int(round(t_fin / h))
    y = np.empty((n + 1, 3))
    y[0] = y0
    T, H, P = y0
    h2, h6 = h / 2.0, h / 6.0
    for k in range(n):
        t = k * h
        a1, b1, c1 = rhs(t, T, H, P, p)
        a2, b2, c2 = rhs(t + h2, T + h2 * a1, H + h2 * b1, P + h2 * c1, p)
        a3, b3, c3 = rhs(t + h2, T + h2 * a2, H + h2 * b2, P + h2 * c2, p)
        a4, b4, c4 = rhs(t + h, T + h * a3, H + h * b3, P + h * c3, p)
        T += h6 * (a1 + 2 * a2 + 2 * a3 + a4)
        H += h6 * (b1 + 2 * b2 + 2 * b3 + b4)
        P += h6 * (c1 + 2 * c2 + 2 * c3 + c4)
        y[k + 1] = (T, H, P)
    return np.arange(n + 1) * h, y


def jacobiano(t, T, H, P, p):
    """Matriz jacobiana df/dy del sistema, evaluada en (t, T, H, P)."""
    R = p["Roce"] + p["A"] * math.cos(OMEGA * t - math.pi)
    return np.array([
        [R * (1.0 - 2.0 * T / p["Tmax"]), -p["alpha"], 0.0],
        [2.0 * p["beta"] * T, -p["gamma"] * P, -p["gamma"] * H],
        [0.0, 2.0 * p["delta"] * H, -p["sigma"]],
    ])


# --------------------------------------------------------------------------
# Indicadores del último ciclo anual
# --------------------------------------------------------------------------
def tiempo_del_maximo(tc, v):
    """
    Instante del máximo de v sobre la grilla tc. Si el máximo es interior se
    refina con una parábola por 3 puntos (error mucho menor que h).
    Devuelve (t_max, es_interior).
    """
    i = int(np.argmax(v))
    interior = 0 < i < len(v) - 1
    if interior:
        a, b, c = v[i - 1], v[i], v[i + 1]
        den = a - 2.0 * b + c
        if den < 0.0:
            return tc[i] + 0.5 * (a - c) / den * (tc[1] - tc[0]), True
    return tc[i], interior


def indicadores(t, y, t_ini=T_ULTIMO_CICLO, periodo=PERIODO):
    """SP, dtTP y estadísticos de P usando únicamente t_ini <= t <= t_fin."""
    i0 = int(np.searchsorted(t, t_ini - 1e-9))
    tc, T, P = t[i0:], y[i0:, 0], y[i0:, 2]

    Pmax, Pmin = float(P.max()), float(P.min())
    # Promedio temporal (integral trapezoidal / duración): no cuenta doble
    # los extremos de la grilla.
    Pmed = float(_trapecio(P, tc) / (tc[-1] - tc[0]))
    SP = (Pmax - Pmin) / Pmed

    tT, ok_T = tiempo_del_maximo(tc, T)
    tP, ok_P = tiempo_del_maximo(tc, P)

    # "Máximo posterior de P": si dentro de la ventana el máximo de P ocurre
    # ANTES que el de T, el máximo posterior es el del ciclo siguiente
    # (régimen periódico), es decir t_P + periodo.
    dt = tP - tT
    siguiente_ciclo = dt < 0.0
    if siguiente_ciclo:
        dt += periodo

    # Chequeo de régimen anual repetitivo: último ciclo vs. el anterior.
    i_prev = int(np.searchsorted(t, t_ini - periodo - 1e-9))
    prev = y[i_prev:i0 + 1]
    last = y[i0:]
    desvio = np.abs(last - prev).max(axis=0) if len(prev) == len(last) else None

    return dict(
        Pmax=Pmax, Pmin=Pmin, Pmed=Pmed, SP=SP,
        tT=tT, Tpico=float(T.max()), tP=tP, dtTP=dt,
        siguiente_ciclo=siguiente_ciclo,
        borde=(not ok_T) or (not ok_P),
        desvio_periodicidad=desvio,
    )


# --------------------------------------------------------------------------
# Formato de números (el informe usa coma decimal)
# --------------------------------------------------------------------------
def fmt(x, nd=4):
    """Número con coma decimal y nd decimales (para los textos del informe)."""
    return f"{x:.{nd}f}".replace(".", ",")


# --------------------------------------------------------------------------
# Estilo de los gráficos
# --------------------------------------------------------------------------
# Paleta: azul = T, naranja = H, verde = P
C_T, C_H, C_P = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
GRID, SURFACE, BAND = "#e1e0d9", "#fcfcfb", "#ecebe6"


def estilo_graficos(mostrar=False):
    """Configura matplotlib (backend sin ventana salvo que se pida mostrar)."""
    import matplotlib
    if not mostrar:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "axes.edgecolor": MUTED, "axes.labelcolor": INK2,
        "xtick.color": INK2, "ytick.color": INK2, "text.color": INK,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
        "axes.axisbelow": True, "font.size": 10,
        "axes.titlesize": 11, "axes.titleweight": "regular",
    })
    return plt
