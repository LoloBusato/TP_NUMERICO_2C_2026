#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PARTE 2 - Severidad estacional, desfase y sensibilidad paramétrica.

Este script es el informe (secciones 3.1 a 3.4) ejecutado en orden:

    3.1  Parámetros nominales: SP y dtTP, gráficos de 60 meses y último ciclo
    3.2  Influencia de A, beta, gamma y sigma (un parámetro por vez)
    3.3  Búsqueda en grilla de la combinación de mayor severidad SP
    3.4  Comparación con el caso nominal y análisis de signos

Todos los indicadores se calculan SOLO sobre el último ciclo (48 <= t <= 60).

Uso:   python parte2_severidad.py [--h 0.001] [--rapido] [--show]
       --rapido usa h = 0.01 en la grilla de 3.3 (≈10 veces más rápido)
"""

import argparse
import csv
import itertools
import json
import os

import numpy as np

import tp_comun as tp

H_PARTE2 = 0.001          # paso de Euler elegido a partir de la Parte 1 [mes]
H_GRILLA_RAPIDA = 0.01
N_BARRIDO = 15            # puntos por parámetro en 3.2
N_GRILLA = 7              # valores por parámetro en la grilla de 3.3
TOL_PERIODICO = 0.01      # desvío relativo máximo entre los dos últimos ciclos

RANGOS = {"A": (0.25, 0.60), "beta": (0.10, 0.30),
          "gamma": (0.05, 0.20), "sigma": (0.40, 1.00)}
ETIQ = {"A": "A [1/mes]", "beta": "β", "gamma": "γ", "sigma": "σ [1/mes]"}

RES = {}


def encabezado(txt):
    print("\n" + "=" * 72 + f"\n {txt}\n" + "=" * 72)


def con(**cambios):
    """Parámetros nominales con los cambios indicados (no modifica los nominales)."""
    p = dict(tp.PARAMS_NOMINALES)
    p.update(cambios)
    return p


def evaluar(p, h):
    """Simula y devuelve los indicadores + chequeo de régimen repetitivo."""
    t, y = tp.euler(p, h)
    if tp.diverge(y):
        return None
    ind = tp.indicadores(t, y)
    i0 = int(np.searchsorted(t, tp.T_ULTIMO_CICLO - 1e-9))
    esc = np.abs(y[i0:]).max(axis=0)
    d = ind["desvio_periodicidad"]
    ind["desvio_rel"] = float(np.max(d / esc)) if d is not None else float("nan")
    ind["amplitud"] = ind["Pmax"] - ind["Pmin"]
    ind["Tmin_total"] = float(y[:, 0].min())
    ind["_t"], ind["_y"] = t, y
    return ind


def sin_series(ind):
    return {k: v for k, v in ind.items() if not k.startswith("_") and k != "desvio_periodicidad"}


# --------------------------------------------------------------------------
# Gráficos (estilo del script original de tp_el_nino_parte2a)
# --------------------------------------------------------------------------
SERIES = [("T", r"$T$ [°C]", tp.C_T), ("H", r"$H$ [u.h.]", tp.C_H),
          ("P", r"$P$ [u.p.]", tp.C_P)]


def grafico_60_meses(plt, t, y, h, ruta):
    fig, axs = plt.subplots(3, 1, figsize=(9, 7.5), sharex=True)
    for ax, (_, etiq, col), j in zip(axs, SERIES, range(3)):
        ax.axvspan(tp.T_ULTIMO_CICLO, tp.T_FIN, color=tp.BAND, zorder=0)
        ax.plot(t, y[:, j], color=col, lw=1.8)
        ax.set_ylabel(etiq)
        ax.set_xlim(0, tp.T_FIN)
    axs[0].text(tp.T_ULTIMO_CICLO + 0.3, axs[0].get_ylim()[1], "último ciclo",
                color=tp.INK2, va="top", fontsize=9)
    axs[-1].set_xlabel("t [meses]")
    axs[-1].set_xticks(range(0, 61, 12))
    fig.suptitle(f"Modelo de El Niño - parámetros nominales (Euler, h = {h:g} mes)",
                 x=0.07, ha="left", fontsize=12)
    fig.tight_layout()
    fig.savefig(ruta, dpi=160)
    plt.close(fig)


def grafico_ultimo_ciclo(plt, t, y, ind, ruta):
    i0 = int(np.searchsorted(t, tp.T_ULTIMO_CICLO - 1e-9))
    tc = t[i0:]
    fig, axs = plt.subplots(3, 1, figsize=(9, 7.5), sharex=True)
    for ax, (_, etiq, col), j in zip(axs, SERIES, range(3)):
        ax.plot(tc, y[i0:, j], color=col, lw=1.8)
        ax.set_ylabel(etiq)
        ax.set_xlim(tp.T_ULTIMO_CICLO, tp.T_FIN)
        ax.axvline(ind["tT"], color=tp.C_T, lw=1.0, ls="--", alpha=0.9)
        ax.axvline(ind["tP"], color=tp.C_P, lw=1.0, ls="--", alpha=0.9)
    axP = axs[2]
    for valor, nom in [(ind["Pmax"], r"$P_{max}$"), (ind["Pmin"], r"$P_{min}$"),
                       (ind["Pmed"], r"$P_{med}$")]:
        axP.axhline(valor, color=tp.INK2, lw=0.9, ls=":")
        axP.text(tp.T_FIN + 0.12, valor, f"{nom} = {valor:.3f}", va="center",
                 fontsize=9, color=tp.INK2, clip_on=False)
    axs[0].plot([ind["tT"]], [ind["Tpico"]], "o", color=tp.C_T, ms=7, mec=tp.SURFACE, mew=1.5)
    axP.plot([ind["tP"]], [ind["Pmax"]], "o", color=tp.C_P, ms=7, mec=tp.SURFACE, mew=1.5)
    lo, hi = axs[0].get_ylim()
    axs[0].set_ylim(lo, hi + 0.22 * (hi - lo))
    ytxt = axs[0].get_ylim()[1]
    axs[0].text(ind["tT"] - 0.12, ytxt, rf"$t_{{T,max}}$ = {ind['tT']:.2f}",
                ha="right", va="top", fontsize=9, color=tp.INK2)
    axs[0].text(ind["tP"] + 0.12, ytxt, rf"$t_{{P,max}}$ = {ind['tP']:.2f}",
                ha="left", va="top", fontsize=9, color=tp.INK2)
    axs[-1].set_xlabel("t [meses]")
    axs[-1].set_xticks(range(48, 61, 2))
    fig.suptitle(rf"Último ciclo anual (48 ≤ t ≤ 60):  $S_P$ = {ind['SP']:.3f},  "
                 rf"$\Delta t_{{TP}}$ = {ind['dtTP']:.2f} meses",
                 x=0.07, ha="left", fontsize=12)
    fig.tight_layout()
    fig.subplots_adjust(right=0.82)
    fig.savefig(ruta, dpi=160)
    plt.close(fig)


# --------------------------------------------------------------------------
# 3.1 Parámetros nominales
# --------------------------------------------------------------------------
def seccion_3_1(plt, h):
    encabezado("3.1  Parámetros nominales")
    nom = evaluar(con(), h)
    print(f" Euler explícito, h = {h:g} mes, último ciclo 48 <= t <= 60")
    print(f"   Pmax = {nom['Pmax']:.5f}   Pmin = {nom['Pmin']:.5f}   Pmed = {nom['Pmed']:.5f}")
    print(f"   SP   = (Pmax - Pmin)/Pmed = {nom['SP']:.5f}")
    print(f"   t(Tmax) = {nom['tT']:.4f}   t(Pmax) = {nom['tP']:.4f}   dtTP = {nom['dtTP']:.4f} meses")
    print(f"   Régimen repetitivo: desvío relativo último ciclo vs 36-48 = {nom['desvio_rel']:.2e}")

    # Justificación del paso: Euler con varios h vs. RK4 muy fino
    t_r, y_r = tp.rk4(con(), 0.0005)
    ref = tp.indicadores(t_r, y_r)
    print("\n Verificación del paso (Euler vs RK4 con h = 0.0005):")
    print(f"   {'h [mes]':>9} {'max|T - Tref|':>15} {'SP':>10} {'dtTP':>9}")
    filas = []
    for hh in (0.5, 0.25, 0.1, 0.05, 0.025, 0.01, 0.005, 0.0025, 0.001):
        t, y = tp.euler(con(), hh)
        r = int(round(hh / 0.0005))
        err = float(np.abs(y[:, 0] - y_r[::r, 0]).max())
        ind = tp.indicadores(t, y)
        filas.append(dict(h=hh, err=err, SP=ind["SP"], dtTP=ind["dtTP"]))
        marca = "  <- elegido" if abs(hh - h) < 1e-12 else ""
        print(f"   {hh:>9g} {err:>15.3e} {ind['SP']:>10.5f} {ind['dtTP']:>9.4f}{marca}")
    print(f"   {'RK4 ref.':>9} {'-':>15} {ref['SP']:>10.5f} {ref['dtTP']:>9.4f}")

    d = os.path.join(tp.carpeta_resultados())
    grafico_60_meses(plt, nom["_t"], nom["_y"], h, os.path.join(d, "fig_3_1_60meses.png"))
    grafico_ultimo_ciclo(plt, nom["_t"], nom["_y"], nom, os.path.join(d, "fig_3_1_ultimo_ciclo.png"))
    RES["h"] = h
    RES["nominal"] = sin_series(nom)
    RES["verificacion_h"] = filas
    RES["rk4_ref"] = dict(SP=ref["SP"], dtTP=ref["dtTP"])
    return nom


# --------------------------------------------------------------------------
# 3.2 Influencia de cada parámetro
# --------------------------------------------------------------------------
def seccion_3_2(plt, h, nom):
    encabezado("3.2  Influencia de A, beta, gamma y sigma (uno por vez)")
    barrido = {}
    filas_csv = []
    for k, (a, b) in RANGOS.items():
        q = np.linspace(a, b, N_BARRIDO)
        inds = [evaluar(con(**{k: v}), h) for v in q]
        barrido[k] = (q, inds)
        for v, i in zip(q, inds):
            filas_csv.append([k, f"{v:.5f}", f"{i['SP']:.5f}", f"{i['dtTP']:.4f}",
                              f"{i['Pmax']:.5f}", f"{i['Pmin']:.5f}", f"{i['Pmed']:.5f}"])
        i0, i1 = inds[0], inds[-1]
        print(f" {k:>5}: [{a:g}, {b:g}]  SP {i0['SP']:.4f} -> {i1['SP']:.4f}   "
              f"dtTP {i0['dtTP']:.3f} -> {i1['dtTP']:.3f}   "
              f"Pmed {i0['Pmed']:.3f} -> {i1['Pmed']:.3f}   "
              f"(Pmax-Pmin) {i0['amplitud']:.3f} -> {i1['amplitud']:.3f}")
    with open(os.path.join(tp.carpeta_resultados(), "tabla_3_2_sensibilidad.csv"),
              "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["parametro", "valor", "SP", "dtTP", "Pmax", "Pmin", "Pmed"])
        w.writerows(filas_csv)

    fig, ax = plt.subplots(2, 4, figsize=(15, 6))
    for j, (k, (q, inds)) in enumerate(barrido.items()):
        ax[0, j].plot(q, [i["SP"] for i in inds], "o-", color=tp.C_T, ms=4)
        ax[1, j].plot(q, [i["dtTP"] for i in inds], "o-", color="#d6336c", ms=4)
        for r in (0, 1):
            ax[r, j].axvline(tp.PARAMS_NOMINALES[k], color="k", ls=":", lw=1)
            ax[r, j].set_xlabel(ETIQ[k])
        ax[0, j].set_ylabel("SP"); ax[1, j].set_ylabel("Δt_TP [meses]")
    fig.suptitle("Sensibilidad (línea punteada: valor nominal)", x=0.07, ha="left")
    fig.tight_layout()
    fig.savefig(os.path.join(tp.carpeta_resultados(), "fig_3_2_sensibilidad.png"), dpi=160)
    plt.close(fig)

    RES["barrido"] = {
        k: dict(rango=list(RANGOS[k]),
                q=[float(v) for v in q],
                SP=[i["SP"] for i in inds], dtTP=[i["dtTP"] for i in inds],
                Pmed=[i["Pmed"] for i in inds], amplitud=[i["amplitud"] for i in inds])
        for k, (q, inds) in barrido.items()}


# --------------------------------------------------------------------------
# 3.3 Búsqueda de la máxima severidad
# --------------------------------------------------------------------------
def seccion_3_3(plt, h):
    encabezado("3.3  Grilla 7x7x7x7 de A, beta, gamma, sigma")
    vals = {k: np.linspace(a, b, N_GRILLA) for k, (a, b) in RANGOS.items()}
    print(f" {N_GRILLA ** 4} combinaciones, h = {h:g} mes (esto puede tardar un poco)...")
    filas, invalidas = [], 0
    for A, be, ga, si in itertools.product(*vals.values()):
        ind = evaluar(con(A=A, beta=be, gamma=ga, sigma=si), h)
        if ind is None or ind["desvio_rel"] > TOL_PERIODICO or not np.isfinite(ind["SP"]):
            invalidas += 1
            continue
        filas.append((ind["SP"], A, be, ga, si, ind["dtTP"], ind["desvio_rel"]))
    filas.sort(reverse=True)
    print(f" Válidas (acotadas y con régimen anual repetitivo, desvío < {TOL_PERIODICO:.0%}): "
          f"{len(filas)} de {N_GRILLA ** 4}  (descartadas: {invalidas})")
    print(f" {'SP':>8} {'A':>6} {'beta':>6} {'gamma':>6} {'sigma':>6} {'dtTP':>7} {'desvío':>9}")
    for f in filas[:10]:
        print(f" {f[0]:8.4f} {f[1]:6.3f} {f[2]:6.3f} {f[3]:6.3f} {f[4]:6.3f} {f[5]:7.3f} {f[6]:9.2e}")
    print(f" SP mínimo de la grilla: {filas[-1][0]:.4f}   (rango SP: {filas[-1][0]:.3f} - {filas[0][0]:.3f})")

    with open(os.path.join(tp.carpeta_resultados(), "tabla_3_3_grilla.csv"),
              "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["SP", "A", "beta", "gamma", "sigma", "dtTP", "desvio_rel"])
        w.writerows([[f"{x:.6g}" for x in fila] for fila in filas])

    # Mapa de SP en el plano (gamma, A) con beta y sigma en su valor óptimo
    mejor = filas[0]
    mapa = np.full((N_GRILLA, N_GRILLA), np.nan)
    for sp, A, be, ga, si, _, _ in filas:
        if abs(be - mejor[2]) < 1e-9 and abs(si - mejor[4]) < 1e-9:
            mapa[int(np.argmin(abs(vals["A"] - A))), int(np.argmin(abs(vals["gamma"] - ga)))] = sp
    fig, ax = plt.subplots(figsize=(5.6, 4.4))
    mg = (RANGOS["gamma"][1] - RANGOS["gamma"][0]) / (N_GRILLA - 1) / 2   # media celda
    mA = (RANGOS["A"][1] - RANGOS["A"][0]) / (N_GRILLA - 1) / 2
    im = ax.imshow(mapa, origin="lower", aspect="auto", cmap="viridis",
                   extent=[RANGOS["gamma"][0] - mg, RANGOS["gamma"][1] + mg,
                           RANGOS["A"][0] - mA, RANGOS["A"][1] + mA])
    ax.grid(False)
    ax.set_xlabel("γ"); ax.set_ylabel("A [1/mes]")
    ax.set_title(f"SP en el plano (γ, A)  –  β = {mejor[2]:g}, σ = {mejor[4]:g}")
    fig.colorbar(im, ax=ax, label="SP")
    fig.tight_layout()
    fig.savefig(os.path.join(tp.carpeta_resultados(), "fig_3_3_mapa_SP.png"), dpi=160)
    plt.close(fig)

    RES["grilla"] = dict(n_validas=len(filas), n_total=N_GRILLA ** 4,
                         h_grilla=h, top=[dict(SP=f[0], A=f[1], beta=f[2], gamma=f[3],
                                               sigma=f[4], dtTP=f[5]) for f in filas[:10]],
                         SP_min=filas[-1][0])
    return dict(A=mejor[1], beta=mejor[2], gamma=mejor[3], sigma=mejor[4])


# --------------------------------------------------------------------------
# 3.4 Comparación con el nominal
# --------------------------------------------------------------------------
def seccion_3_4(plt, h, nom, mejor_p):
    encabezado("3.4  Comparación nominal vs. mayor severidad")
    opt = evaluar(con(**mejor_p), h)
    print(f" {'':>10} {'nominal':>10} {'mayor SP':>10}")
    for k in ("A", "beta", "gamma", "sigma"):
        print(f" {k:>10} {tp.PARAMS_NOMINALES[k]:>10.3f} {mejor_p[k]:>10.3f}")
    for k, n in (("SP", "SP"), ("dtTP", "dtTP"), ("Pmax", "Pmax"), ("Pmin", "Pmin"), ("Pmed", "Pmed")):
        print(f" {n:>10} {nom[k]:>10.4f} {opt[k]:>10.4f}")
    print(f" Régimen repetitivo del escenario severo: desvío relativo = {opt['desvio_rel']:.2e}")

    # ¿Es "todo al máximo"?
    maximos = {k: RANGOS[k][1] for k in RANGOS}
    todo_max = evaluar(con(**maximos), h)
    print(f"\n Todo en su valor máximo (A=0,60; β=0,30; γ=0,20; σ=1,00): "
          f"SP = {todo_max['SP']:.4f}  (escenario severo: {opt['SP']:.4f})")

    # Desde el óptimo, llevar cada parámetro al extremo opuesto
    print("\n Desde el escenario severo, un parámetro al extremo opuesto:")
    flips = []
    for k in RANGOS:
        otro = RANGOS[k][0] if abs(mejor_p[k] - RANGOS[k][1]) < 1e-9 else RANGOS[k][1]
        i = evaluar(con(**{**mejor_p, k: otro}), h)
        flips.append(dict(parametro=k, valor=otro, SP=i["SP"]))
        print(f"   {k:>5} = {otro:<5g} -> SP = {i['SP']:.4f}  (cambia {i['SP'] - opt['SP']:+.4f})")

    d = tp.carpeta_resultados()
    fig, ax = plt.subplots(2, 3, figsize=(14, 6.2))
    for fila, (a, b, tit) in enumerate([(0, 60, "60 meses"), (48, 60, "último ciclo (48–60)")]):
        for col, (j, (_, etiq, _)) in enumerate(zip(range(3), SERIES)):
            m = (nom["_t"] >= a - 1e-9) & (nom["_t"] <= b + 1e-9)
            ax[fila, col].plot(nom["_t"][m], nom["_y"][m, j], color=tp.C_T, lw=1.6, label="nominal")
            ax[fila, col].plot(opt["_t"][m], opt["_y"][m, j], color="#d6336c", lw=1.6, label="mayor SP")
            ax[fila, col].set_ylabel(etiq); ax[fila, col].set_xlabel("t [meses]")
        ax[fila, 1].set_title(tit)
    ax[0, 0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(d, "fig_3_4_comparacion.png"), dpi=160)
    plt.close(fig)

    grafico_ultimo_ciclo(plt, opt["_t"], opt["_y"], opt, os.path.join(d, "fig_3_4_ultimo_ciclo_severo.png"))
    RES["severo"] = dict(parametros=mejor_p, **sin_series(opt))
    RES["todo_maximo"] = sin_series(todo_max)
    RES["flips"] = flips


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--h", type=float, default=H_PARTE2, help="paso de Euler (default 0.001)")
    ap.add_argument("--rapido", action="store_true", help="grilla de 3.3 con h = 0.01")
    ap.add_argument("--show", action="store_true", help="mostrar las figuras")
    args = ap.parse_args()
    if abs(tp.T_ULTIMO_CICLO / args.h - round(tp.T_ULTIMO_CICLO / args.h)) > 1e-9:
        raise SystemExit(f"h = {args.h} no divide a 48: t = 48 no cae sobre la grilla.")
    plt = tp.estilo_graficos(args.show)

    nom = seccion_3_1(plt, args.h)
    seccion_3_2(plt, args.h, nom)
    mejor_p = seccion_3_3(plt, H_GRILLA_RAPIDA if args.rapido else args.h)
    seccion_3_4(plt, args.h, nom, mejor_p)

    ruta = os.path.join(tp.carpeta_resultados(), "resultados_parte2.json")
    with open(ruta, "w", encoding="utf-8") as f:
        json.dump(RES, f, indent=2, ensure_ascii=False, default=float)
    print(f"\n Resultados guardados en {ruta}")
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
