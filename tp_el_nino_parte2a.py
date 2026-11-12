#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CB051 Modelación Numérica / 9504 Análisis Numérico I - 2do cuatrimestre 2026
Modelo simplificado de El Niño - PARTE 2, punto (a)

Qué hace
--------
1. Integra el sistema no lineal
       dT/dt = R(t) T (1 - T/Tmax) - alpha H
       dH/dt = beta T^2 - gamma H P
       dP/dt = delta H^2 - sigma P
   con R(t) = Roce + A cos(w t - pi), w = 2 pi / 12, usando Euler explícito
   con un paso h fino (elegido según la Parte 1), durante 60 meses.
2. Calcula, SOLO sobre el último ciclo anual (48 <= t <= 60 meses):
       SP    = (Pmax - Pmin) / Pmed
       dtTP  = t(Pmax) - t(Tmax)
3. Muestra una tabla de verificación del paso h elegido (Euler con distintos h
   contra una referencia RK4 muy fina) para justificarlo en el informe.
4. Genera dos figuras:
       fig_parte2a_60meses.png   T(t), H(t), P(t) en 0..60 meses
       fig_parte2a_ultimo_ciclo.png   último ciclo anual con los indicadores

Uso
---
    python tp_el_nino_parte2a.py                 # h = 0.001 mes
    python tp_el_nino_parte2a.py --h 0.005       # otro paso
    python tp_el_nino_parte2a.py --outdir figs   # carpeta de salida
    python tp_el_nino_parte2a.py --show          # además muestra las figuras

Requiere: numpy y matplotlib.
"""

import argparse
import math
import os

import numpy as np
import matplotlib

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
    Pmed = float(np.trapezoid(P, tc) / (tc[-1] - tc[0]))
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
# Verificación del paso h (Parte 1 aplicada a esta parte)
# --------------------------------------------------------------------------
def tabla_verificacion(p, h_elegido, h_ref=0.0005):
    """Euler con varios h vs. RK4 fino: error máx. en T, SP y dtTP."""
    t_ref, y_ref = rk4(p, h_ref)
    ind_ref = indicadores(t_ref, y_ref)
    hs = [0.5, 0.25, 0.1, 0.05, 0.025, 0.01, 0.005, 0.0025, 0.001]
    if h_elegido not in hs:
        hs.append(h_elegido)
    hs = sorted(set(hs), reverse=True)

    filas = []
    for h in hs:
        t, y = euler(p, h)
        paso = int(round(h / h_ref))
        if abs(paso * h_ref - h) > 1e-12:
            err = float("nan")
        else:
            err = float(np.abs(y[:, 0] - y_ref[::paso, 0]).max())
        ind = indicadores(t, y)
        filas.append((h, err, ind["SP"], ind["dtTP"]))
    return ind_ref, filas


# --------------------------------------------------------------------------
# Gráficos
# --------------------------------------------------------------------------
# Paleta categórica (slots 1-3) y tintas: azul=T, naranja=H, aqua=P
C_T, C_H, C_P = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
GRID, SURFACE, BAND = "#e1e0d9", "#fcfcfb", "#ecebe6"


def estilo():
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


SERIES = [("T", r"$T$ [°C]", C_T), ("H", r"$H$ [u.h.]", C_H),
          ("P", r"$P$ [u.p.]", C_P)]


def grafico_60_meses(t, y, ind, h, ruta):
    fig, axs = plt.subplots(3, 1, figsize=(9, 7.5), sharex=True)
    for ax, (nom, etiq, col), j in zip(axs, SERIES, range(3)):
        ax.axvspan(T_ULTIMO_CICLO, T_FIN, color=BAND, zorder=0)
        ax.plot(t, y[:, j], color=col, lw=1.8)
        ax.set_ylabel(etiq)
        ax.set_xlim(0, T_FIN)
    axs[0].text(T_ULTIMO_CICLO + 0.3, axs[0].get_ylim()[1], "último ciclo",
                color=INK2, va="top", fontsize=9)
    axs[-1].set_xlabel("t [meses]")
    axs[-1].set_xticks(range(0, 61, 12))
    fig.suptitle("Modelo de El Niño - parámetros nominales (Euler, "
                 f"h = {h:g} mes)", x=0.07, ha="left", fontsize=12)
    fig.tight_layout()
    fig.savefig(ruta, dpi=160)
    return fig


def grafico_ultimo_ciclo(t, y, ind, h, ruta):
    i0 = int(np.searchsorted(t, T_ULTIMO_CICLO - 1e-9))
    tc = t[i0:]
    fig, axs = plt.subplots(3, 1, figsize=(9, 7.5), sharex=True)
    for ax, (nom, etiq, col), j in zip(axs, SERIES, range(3)):
        ax.plot(tc, y[i0:, j], color=col, lw=1.8)
        ax.set_ylabel(etiq)
        ax.set_xlim(T_ULTIMO_CICLO, T_FIN)
        # marcas de los instantes de máximo de T y de P
        ax.axvline(ind["tT"], color=C_T, lw=1.0, ls="--", alpha=0.9)
        ax.axvline(ind["tP"], color=C_P, lw=1.0, ls="--", alpha=0.9)

    # Eje de P: líneas de Pmax, Pmin y Pmed con etiqueta a la derecha
    axP = axs[2]
    for valor, nom in [(ind["Pmax"], r"$P_{max}$"), (ind["Pmin"], r"$P_{min}$"),
                       (ind["Pmed"], r"$P_{med}$")]:
        axP.axhline(valor, color=INK2, lw=0.9, ls=":")
        axP.text(T_FIN + 0.12, valor, f"{nom} = {valor:.3f}", va="center",
                 fontsize=9, color=INK2, clip_on=False)

    # Etiquetas de los dos máximos
    axs[0].plot([ind["tT"]], [ind["Tpico"]], "o", color=C_T, ms=7,
                mec=SURFACE, mew=1.5)
    axP.plot([ind["tP"]], [ind["Pmax"]], "o", color=C_P, ms=7,
             mec=SURFACE, mew=1.5)
    # Rótulos en el margen libre sobre el primer panel, hacia afuera de cada línea
    lo, hi = axs[0].get_ylim()
    axs[0].set_ylim(lo, hi + 0.22 * (hi - lo))
    ytxt = axs[0].get_ylim()[1]
    axs[0].text(ind["tT"] - 0.12, ytxt, rf"$t_{{T,max}}$ = {ind['tT']:.2f}",
                ha="right", va="top", fontsize=9, color=INK2)
    axs[0].text(ind["tP"] + 0.12, ytxt, rf"$t_{{P,max}}$ = {ind['tP']:.2f}",
                ha="left", va="top", fontsize=9, color=INK2)

    axs[-1].set_xlabel("t [meses]")
    axs[-1].set_xticks(range(48, 61, 2))
    nota = ("  (máximo de P del ciclo siguiente: $t_P$ + 12)"
            if ind["siguiente_ciclo"] else "")
    fig.suptitle(
        rf"Último ciclo anual (48 ≤ t ≤ 60):  $S_P$ = {ind['SP']:.3f},  "
        rf"$\Delta t_{{TP}}$ = {ind['dtTP']:.2f} meses{nota}",
        x=0.07, ha="left", fontsize=12)
    fig.tight_layout()
    fig.subplots_adjust(right=0.82)
    fig.savefig(ruta, dpi=160)
    return fig


# --------------------------------------------------------------------------
# Programa principal
# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--h", type=float, default=H_DEFECTO,
                    help="paso de Euler en meses (default 0.001)")
    ap.add_argument("--outdir", default=".", help="carpeta de salida de las figuras")
    ap.add_argument("--show", action="store_true", help="mostrar las figuras en pantalla")
    args = ap.parse_args()

    p = PARAMS_NOMINALES
    h = args.h
    n_ciclo = (T_ULTIMO_CICLO / h)
    if abs(n_ciclo - round(n_ciclo)) > 1e-9:
        raise SystemExit(f"h = {h} no divide a 48: t = 48 no cae sobre la grilla.")

    if not args.show:
        matplotlib.use("Agg")
    os.makedirs(args.outdir, exist_ok=True)

    # --- Solución con Euler -------------------------------------------------
    t, y = euler(p, h)
    ind = indicadores(t, y)

    print("=" * 66)
    print(f" Parte 2a - parámetros nominales - Euler explícito, h = {h:g} mes")
    print("=" * 66)
    print(f" Último ciclo anual: 48 <= t <= 60 meses")
    print(f"   Pmax = {ind['Pmax']:.5f} u.p.")
    print(f"   Pmin = {ind['Pmin']:.5f} u.p.")
    print(f"   Pmed = {ind['Pmed']:.5f} u.p.   (promedio temporal, trapecios)")
    print(f"   SP   = (Pmax - Pmin)/Pmed = {ind['SP']:.5f}")
    print(f"   t(Tmax) = {ind['tT']:.4f} meses   (Tmax del ciclo = {ind['Tpico']:.4f} °C)")
    print(f"   t(Pmax) = {ind['tP']:.4f} meses")
    print(f"   dtTP    = {ind['dtTP']:.4f} meses", end="")
    print("   [máx. de P del ciclo siguiente: t(Pmax)+12]" if ind["siguiente_ciclo"] else "")
    if ind["borde"]:
        print("   AVISO: algún máximo cae en el borde de la ventana; "
              "revisar el valor de dtTP.")
    if ind["desvio_periodicidad"] is not None:
        d = ind["desvio_periodicidad"]
        print(" Régimen repetitivo (último ciclo vs. ciclo 36-48), desvío máx.:")
        print(f"   |dT| = {d[0]:.2e}   |dH| = {d[1]:.2e}   |dP| = {d[2]:.2e}")

    # --- Verificación del paso ---------------------------------------------
    print()
    print(" Verificación del paso (Euler vs. RK4 con h = 0.0005):")
    ind_ref, filas = tabla_verificacion(p, h)
    print(f"   {'h [mes]':>9} {'max|T - Tref|':>15} {'SP':>10} {'dtTP':>9}")
    for hh, err, sp, dtp in filas:
        marca = "  <- elegido" if abs(hh - h) < 1e-12 else ""
        print(f"   {hh:>9g} {err:>15.3e} {sp:>10.5f} {dtp:>9.4f}{marca}")
    print(f"   {'RK4 ref.':>9} {'-':>15} {ind_ref['SP']:>10.5f} {ind_ref['dtTP']:>9.4f}")
    print("=" * 66)

    # --- Figuras -----------------------------------------------------------
    estilo()
    f1 = os.path.join(args.outdir, "fig_parte2a_60meses.png")
    f2 = os.path.join(args.outdir, "fig_parte2a_ultimo_ciclo.png")
    grafico_60_meses(t, y, ind, h, f1)
    grafico_ultimo_ciclo(t, y, ind, h, f2)
    print(f" Figuras guardadas en:\n   {f1}\n   {f2}")
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
