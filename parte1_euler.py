#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PARTE 1 - Convergencia, orden y estabilidad del método de Euler explícito.

Este script es el informe (secciones 2.1 a 2.4) ejecutado en orden:

    2.1  Solución de referencia (Euler con h_ref muy chico + control)
    2.2  Resolución con sucesivos pasos h y error máximo en T
    2.3  Verificación experimental del orden: p = log(Err(2h)/Err(h)) / log 2
    2.4  Paso temporal aumentado: inestabilidad, umbral y análisis teórico
         con el jacobiano (autovalores y multiplicadores de Floquet)

Genera figuras y tablas en ./resultados y guarda los números en
resultados/resultados_parte1.json (los usa generar_informe.py).

Uso:   python parte1_euler.py [--show]
"""

import argparse
import json
import math
import os

import numpy as np

import tp_comun as tp

P = tp.PARAMS_NOMINALES
H_REF = 0.0001                                    # paso de la referencia [mes]
HS_ERROR = [0.025, 0.05, 0.1, 0.2, 0.4, 0.8]      # cada h es el doble del anterior
HS_ESTAB = [0.8, 1.0, 1.2, 1.25, 1.26, 1.4, 1.6, 3.2, 6.4]

RES = {}                                          # resultados que van al informe


def encabezado(txt):
    print("\n" + "=" * 72 + f"\n {txt}\n" + "=" * 72)


def error_maximo(h, t_ref, y_ref):
    """Err(h) = max |T_h - T_ref| en los nodos de h (todos son nodos de la ref)."""
    t, y = tp.euler(P, h)
    if tp.diverge(y):
        return None, None, t, y
    r = int(round(h / H_REF))
    e = np.abs(y[:, 0] - y_ref[::r, 0][:len(t)])
    return float(e.max()), float(t[int(np.argmax(e))]), t, y


# --------------------------------------------------------------------------
# 2.1 Solución de referencia
# --------------------------------------------------------------------------
def seccion_2_1():
    encabezado("2.1  Solución de referencia")
    t_ref, y_ref = tp.euler(P, H_REF)
    t_ref2, y_ref2 = tp.euler(P, H_REF / 2)
    ctrl = float(np.abs(y_ref[:, 0] - y_ref2[::2, 0]).max())
    print(f" h_ref = {H_REF:g} mes  ->  {len(t_ref) - 1} pasos")
    print(f" Control: max|T(h_ref) - T(h_ref/2)| = {ctrl:.3e} °C")
    print(" (error estimado de la referencia ~ ese valor; debe ser mucho menor")
    print("  que el error de las discretizaciones que se comparan)")
    RES["h_ref"] = H_REF
    RES["pasos_ref"] = len(t_ref) - 1
    RES["control_ref"] = ctrl
    return t_ref, y_ref


# --------------------------------------------------------------------------
# 2.2 Sucesivos pasos h
# --------------------------------------------------------------------------
def seccion_2_2(t_ref, y_ref, plt):
    encabezado("2.2  Error máximo en T para h, 2h, 4h, ...")
    filas, trayect = [], {}
    print(f" {'h':>7} {'pasos':>7} {'Err(h) [°C]':>13} {'t del error máx':>16}")
    for h in HS_ERROR:
        err, t_err, t, y = error_maximo(h, t_ref, y_ref)
        filas.append(dict(h=h, pasos=len(t) - 1, err=err, t_err=t_err))
        trayect[h] = (t, y)
        print(f" {h:>7g} {len(t) - 1:>7d} {err:>13.6f} {t_err:>16.1f}")
    RES["errores"] = filas

    # Por qué NO usar h = 0,1 como referencia (lo que figuraba en el borrador)
    _, yr01 = trayect[0.1]
    mal = []
    for h in (0.2, 0.4, 0.8):
        t, y = trayect[h]
        r = int(round(h / 0.1))
        mal.append(float(np.abs(y[:, 0] - yr01[::r, 0]).max()))
    print("\n Contraste: Err(h) usando como referencia Euler con h = 0,1")
    print("   h = 0,2 / 0,4 / 0,8 ->", ", ".join(f"{e:.6f}" for e in mal))
    RES["errores_ref_h01"] = dict(zip(["0.2", "0.4", "0.8"], mal))
    err_propio = [f["err"] for f in filas if f["h"] == 0.1][0]
    print(f"   El propio Err(0,1) respecto de la referencia fina es {err_propio:.6f}:")
    print("   del mismo orden que Err(0,2), por eso h = 0,1 no sirve de referencia.")

    # Figura: T(t) y error en función de t
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    ax[0].plot(t_ref, y_ref[:, 0], "k", lw=2.2, label=f"referencia (h = {H_REF:g})")
    for h, col in zip((0.1, 0.4, 0.8), ("#1baf7a", "#2a78d6", "#eb6834")):
        t, y = trayect[h]
        ax[0].plot(t, y[:, 0], color=col, lw=1.3, label=f"h = {h:g}")
    ax[0].set_xlim(0, 24)
    ax[0].set_xlabel("t [meses]"); ax[0].set_ylabel("T [°C]")
    ax[0].set_title("Temperatura (primeros 24 meses)"); ax[0].legend(fontsize=8)
    for h in HS_ERROR[2:]:
        t, y = trayect[h]
        r = int(round(h / H_REF))
        ax[1].plot(t, np.abs(y[:, 0] - y_ref[::r, 0][:len(t)]), lw=1.3, label=f"h = {h:g}")
    ax[1].set_xlabel("t [meses]"); ax[1].set_ylabel("|T_h − T_ref| [°C]")
    ax[1].set_title("Error punto a punto"); ax[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(tp.carpeta_resultados(), "fig_1_2_errores.png"), dpi=160)
    return filas


# --------------------------------------------------------------------------
# 2.3 Orden observado
# --------------------------------------------------------------------------
def seccion_2_3(filas, plt):
    encabezado("2.3  Orden observado  p = log(Err(2h)/Err(h)) / log(2)")
    ords = []
    print(f" {'h -> 2h':>14} {'Err(h)':>11} {'Err(2h)':>11} {'p':>8}")
    for a, b in zip(filas[:-1], filas[1:]):
        p = math.log(b["err"] / a["err"]) / math.log(2.0)
        ords.append(dict(h=a["h"], h2=b["h"], err_h=a["err"], err_2h=b["err"], p=p))
        print(f" {a['h']:>6g} -> {b['h']:<5g} {a['err']:>11.6f} {b['err']:>11.6f} {p:>8.4f}")
    pm = float(np.mean([o["p"] for o in ords]))
    print(f" p medio = {pm:.4f}   (esperado para Euler: p = 1)")
    RES["ordenes"] = ords
    RES["p_medio"] = pm

    hs = np.array([f["h"] for f in filas])
    es = np.array([f["err"] for f in filas])
    plt.figure(figsize=(5.8, 4.2))
    plt.loglog(hs, es, "o-", color="#2a78d6", label="Err(h) medido")
    plt.loglog(hs, es[0] * hs / hs[0], "k--", lw=1, label="pendiente 1 (teórica)")
    plt.xlabel("h [meses]"); plt.ylabel("max |T − T_ref| [°C]")
    plt.title("Convergencia de Euler (escala log-log)"); plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(tp.carpeta_resultados(), "fig_1_3_orden.png"), dpi=160)


# --------------------------------------------------------------------------
# 2.4 Paso aumentado: inestabilidad
# --------------------------------------------------------------------------
def divergencia_en(h):
    t, y = tp.euler(P, h)
    return tp.diverge(y)


def seccion_2_4(t_ref, y_ref, plt):
    encabezado("2.4  Paso temporal aumentado: estabilidad")
    filas = []
    print(f" {'h':>6} {'estado':>10} {'Err [°C]':>12} {'detalle':>34}")
    for h in HS_ESTAB:
        t, y = tp.euler(P, h)
        if tp.diverge(y):
            i = int(np.argmax(np.isnan(y[:, 0])))
            det = f"supera 1e6 en t ≈ {t[i]:.1f} meses"
            filas.append(dict(h=h, estable=False, err=None, detalle=det, t_div=float(t[i])))
            print(f" {h:>6g} {'DIVERGE':>10} {'-':>12} {det:>34}")
        else:
            r = int(round(h / H_REF))
            e = float(np.abs(y[:, 0] - y_ref[::r, 0][:len(t)]).max())
            det = f"T entre {y[:, 0].min():.2f} y {y[:, 0].max():.2f}"
            filas.append(dict(h=h, estable=True, err=e, detalle=det))
            print(f" {h:>6g} {'estable':>10} {e:>12.4f} {det:>34}")
    RES["estabilidad"] = filas

    # Umbral por bisección (h arbitrario; el último paso se acorta si hace falta)
    lo, hi = 1.25, 1.26
    assert not divergencia_en(lo) and divergencia_en(hi)
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        if divergencia_en(mid):
            hi = mid
        else:
            lo = mid
    scan = [(round(h, 3), divergencia_en(h)) for h in np.arange(1.25, 1.401, 0.01)]
    monotono = all(d for h, d in scan if h >= 1.26) and not scan[0][1]
    print(f"\n Umbral de divergencia (bisección en [1,25; 1,26]): h* ≈ {lo:.5f}")
    print(f" Barrido 1,25 -> 1,40 de a 0,01: estable solo en 1,25; "
          f"diverge desde 1,26 en adelante: {monotono}")
    RES["umbral"] = lo
    RES["umbral_monotono"] = monotono

    # Pérdida de la periodicidad anual de la solución numérica (T en t = 12k)
    per = []
    print("\n T(t = 12k) de la solución de Euler (¿se repite cada año?):")
    for h in (0.2, 0.8, 1.0, 12 / 11, 1.2, 1.25):
        t, y = tp.euler(P, h)
        Tk = [float(np.interp(12 * k, t, y[:, 0])) for k in (1, 2, 3, 4, 5)]
        per.append(dict(h=h, T_12k=Tk))
        print(f"   h = {h:7.4f}: " + "  ".join(f"{v:6.3f}" for v in Tk))
    RES["periodicidad"] = per

    # Figura de estabilidad
    fig, ax = plt.subplots(1, 3, figsize=(13.5, 4))
    ax[0].plot(t_ref, y_ref[:, 0], "k", lw=2, label="referencia")
    for h, col in zip((0.8, 1.0, 1.25), ("#1baf7a", "#2a78d6", "#eb6834")):
        t, y = tp.euler(P, h)
        ax[0].plot(t, y[:, 0], color=col, lw=1.3, label=f"h = {h:g}")
    ax[0].set_xlabel("t [meses]"); ax[0].set_ylabel("T [°C]")
    ax[0].set_title("Pasos estables"); ax[0].legend(fontsize=8)
    for h, col in zip((1.26, 1.6), ("#d6336c", "#7048e8")):
        t, y = tp.euler(P, h, limite=1e12)
        ax[1].plot(t, y[:, 0], color=col, lw=1.3, marker=".", label=f"h = {h:g}")
    ax[1].set_ylim(-200, 200); ax[1].set_xlabel("t [meses]")
    ax[1].set_title("Pasos inestables (eje recortado)"); ax[1].legend(fontsize=8)
    ok = [f for f in filas if f["estable"]]
    hh = [f["h"] for f in ok] + [1.26]
    ee = [f["err"] for f in ok]
    ax[2].plot(hh[:-1], ee, "o-", color="#2a78d6")
    ax[2].axvline(RES["umbral"], color="#d6336c", ls="--", lw=1.2,
                  label=f"umbral ≈ {RES['umbral']:.3f}")
    ax[2].set_xlabel("h [meses]"); ax[2].set_ylabel("Err(h) [°C]")
    ax[2].set_title("Error en el borde de estabilidad"); ax[2].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(tp.carpeta_resultados(), "fig_1_4_estabilidad.png"), dpi=160)


# --------------------------------------------------------------------------
# 2.4 (cont.) Referencia teórica con el jacobiano
# --------------------------------------------------------------------------
def h_critico_puntual(autovalores):
    """Euler: |1 + h*lam| < 1  ->  h < -2 Re(lam) / |lam|^2, solo si Re(lam) < 0."""
    v = [-2.0 * l.real / abs(l) ** 2 for l in autovalores if l.real < 0.0]
    return min(v) if v else math.inf


def seccion_2_4_jacobiano(plt):
    encabezado("2.4  Referencia teórica: jacobiano y multiplicadores de Floquet")
    # Órbita "exacta" (RK4 fino) para evaluar el jacobiano
    h_o = 0.005
    t_o, y_o = tp.rk4(P, h_o)
    orbita = lambda t: [float(np.interp(t, t_o, y_o[:, j])) for j in range(3)]

    print(" Autovalores del jacobiano sobre la solución:")
    aut = []
    for tt in (0.0, 6.0, 12.0, 60.0):
        T, H, Pp = orbita(tt)
        ev = np.linalg.eigvals(tp.jacobiano(tt, T, H, Pp, P))
        ev = sorted(ev, key=lambda z: (abs(z.imag), z.real))
        txt = ";  ".join(f"{z.real:.4f}" + (f"{z.imag:+.4f}i" if abs(z.imag) > 1e-9 else "") for z in ev)
        aut.append(dict(t=tt, autovalores=txt,
                        eig=[(float(z.real), float(z.imag)) for z in ev]))
        print(f"   t = {tt:>4g}:  {txt}")
    RES["autovalores"] = aut

    T, H, Pp = orbita(60.0)
    ev60 = np.linalg.eigvals(tp.jacobiano(60.0, T, H, Pp, P))
    comp = [z for z in ev60 if abs(z.imag) > 1e-9][0]
    hc60 = -2 * comp.real / abs(comp) ** 2
    print(f"\n y(60) = ({T:.6f}; {H:.6f}; {Pp:.6f})")
    print(f" h_crítico puntual en t = 60 (par complejo {comp.real:.5f} ± {abs(comp.imag):.5f}i): "
          f"{hc60:.4f} meses")
    print(" (el autovalor real positivo es crecimiento real del sistema, no")
    print("  inestabilidad numérica: el criterio solo aplica a los Re(lambda) < 0)")

    tt = np.arange(0.0, 60.0 + 1e-9, 0.05)
    hc = np.array([h_critico_puntual(np.linalg.eigvals(tp.jacobiano(t, *orbita(t), P))) for t in tt])
    ult = tt >= 48.0
    i_min = int(np.argmin(hc))
    print(f" h_crítico puntual sobre toda la trayectoria: mínimo {hc.min():.3f} en t = {tt[i_min]:.2f}"
          " (transitorio; el límite local tiende a 0 donde un par complejo cruza el eje imaginario)")
    print(f" h_crítico puntual en el último ciclo (48-60): entre {hc[ult].min():.3f} y {hc[ult].max():.3f}")
    RES["hcrit_60"] = float(hc60)
    RES["hcrit_min_total"] = float(hc.min())
    RES["hcrit_min_total_t"] = float(tt[i_min])
    RES["hcrit_ult_min"] = float(hc[ult].min())
    RES["hcrit_ult_max"] = float(hc[ult].max())

    # Multiplicadores de Floquet del mapa de Euler linealizado sobre la órbita,
    # un ciclo anual [48, 60]:  M = prod (I + h J_n)
    def rho(N):
        hh = 12.0 / N
        M = np.eye(3)
        for n in range(N):
            t = 48.0 + n * hh
            M = (np.eye(3) + hh * tp.jacobiano(t, *orbita(t), P)) @ M
        return float(max(abs(np.linalg.eigvals(M)))), hh

    rho_exacto = rho(2400)[0]
    print(f"\n Multiplicador de Floquet exacto (h -> 0): {rho_exacto:.4f}  (<1: órbita anual estable)")
    print(f" {'N':>4} {'h':>8} {'rho(Euler lineal)':>19}")
    fl = []
    for N in range(20, 6, -1):
        r, hh = rho(N)
        fl.append(dict(N=N, h=hh, rho=r))
        print(f" {N:>4d} {hh:>8.4f} {r:>19.4f}")
    RES["floquet"] = fl
    RES["floquet_exacto"] = rho_exacto
    cruce = [f for f in fl if f["rho"] >= 1.0]
    h_inf = max(f["h"] for f in fl if f["rho"] < 1.0 and f["h"] < (cruce[-1]["h"] if cruce else 99))
    RES["floquet_h_estable_max"] = h_inf
    RES["floquet_h_inestable_min"] = min(f["h"] for f in cruce) if cruce else None
    print(f" rho < 1 hasta h = {h_inf:.4f}; rho >= 1 desde h = {RES['floquet_h_inestable_min']:.4f}")

    # Figura: h_crítico puntual en el último ciclo y rho(h)
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    ax[0].plot(tt[ult], hc[ult], color="#2a78d6", lw=1.6)
    ax[0].axhline(RES["umbral"], color="#d6336c", ls="--", lw=1.2, label=f"umbral experimental ({RES['umbral']:.3f})")
    ax[0].axhline(hc60, color="#52514e", ls=":", lw=1.2, label=f"h_crítico(60) = {hc60:.3f}")
    ax[0].set_xlabel("t [meses]"); ax[0].set_ylabel("h crítico puntual [meses]")
    ax[0].set_title("Criterio |1 + hλ| < 1 en el último ciclo"); ax[0].legend(fontsize=8)
    ax[1].semilogy([f["h"] for f in fl], [f["rho"] for f in fl], "o-", color="#eb6834")
    ax[1].axhline(1.0, color="k", lw=1)
    ax[1].axvline(RES["umbral"], color="#d6336c", ls="--", lw=1.2)
    ax[1].set_xlabel("h [meses]"); ax[1].set_ylabel("radio espectral de M (Floquet)")
    ax[1].set_title("Estabilidad de la órbita anual del mapa de Euler"); ax[1].set_ylim(0.01, 100)
    fig.tight_layout()
    fig.savefig(os.path.join(tp.carpeta_resultados(), "fig_1_4_jacobiano.png"), dpi=160)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--show", action="store_true", help="mostrar las figuras")
    args = ap.parse_args()
    plt = tp.estilo_graficos(args.show)

    t_ref, y_ref = seccion_2_1()
    filas = seccion_2_2(t_ref, y_ref, plt)
    seccion_2_3(filas, plt)
    seccion_2_4(t_ref, y_ref, plt)
    seccion_2_4_jacobiano(plt)

    ruta = os.path.join(tp.carpeta_resultados(), "resultados_parte1.json")
    with open(ruta, "w", encoding="utf-8") as f:
        json.dump(RES, f, indent=2, ensure_ascii=False, default=float)
    print(f"\n Resultados guardados en {ruta}")
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
