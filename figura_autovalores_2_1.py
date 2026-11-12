#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figura de la sección 2.1: autovalores del jacobiano sobre la solución (Euler,
h = 0,01) y cotas de paso asociadas. Escribe resultados/fig_2_1_autovalores.png.
"""
import os

import numpy as np

import tp_comun as tp

p = tp.PARAMS_NOMINALES
t, y = tp.euler(p, 0.01)
paso = 10
tt = t[::paso]
ev = np.array([np.linalg.eigvals(tp.jacobiano(t[k], *y[k], p))
               for k in range(0, len(t), paso)])

mod_max = np.abs(ev).max(axis=1)
h_2_sobre_lmax = 2.0 / mod_max

h_crit = np.full(len(tt), np.inf)
for i, e in enumerate(ev):
    for lam in e:
        if lam.real < 0:
            h_crit[i] = min(h_crit[i], -2.0 * lam.real / abs(lam) ** 2)

plt = tp.estilo_graficos()
fig, ax = plt.subplots(1, 2, figsize=(12, 4.2))

real_ = np.array([e[np.argmin(np.abs(e.imag))].real for e in ev])
cplx = np.array([e[np.argmax(e.imag)] for e in ev])
ax[0].plot(tt, real_, color=tp.C_T, lw=1.6, label="autovalor real")
ax[0].plot(tt, cplx.real, color=tp.C_H, lw=1.6, label="Re λ (par complejo)")
ax[0].plot(tt, cplx.imag, color=tp.C_P, lw=1.6, label="Im λ (par complejo)")
ax[0].axhline(0, color=tp.MUTED, lw=0.8)
ax[0].set_xlabel("t [meses]")
ax[0].set_ylabel("autovalores de J [1/mes]")
ax[0].set_title("Autovalores del jacobiano sobre la solución")
ax[0].set_ylim(-0.85, 1.65)
ax[0].legend(frameon=False, loc="upper center", ncol=3, fontsize=8)

ax[1].plot(tt, h_2_sobre_lmax, color=tp.C_T, lw=1.6, label=r"$2/|\lambda_{max}|$")
ax[1].plot(tt, np.minimum(h_crit, 6.5), color=tp.C_H, lw=1.6,
           label=r"$-2\,\mathrm{Re}(\lambda)/|\lambda|^2$ (modos con Re<0)")
ax[1].axhline(2.35, color=tp.INK2, ls=":", lw=1.2, label="2,35 (valor del texto)")
ax[1].axhline(1.25, color=tp.C_P, ls="--", lw=1.2, label="último h estable (1,25)")
ax[1].set_ylim(0, 6.6)
ax[1].set_xlabel("t [meses]")
ax[1].set_ylabel("h [meses]")
ax[1].set_title("Cotas locales de paso (autovalores congelados)")
ax[1].legend(frameon=False, fontsize=8, loc="upper right")

fig.tight_layout()
ruta = os.path.join(tp.carpeta_resultados(), "fig_2_1_autovalores.png")
fig.savefig(ruta, dpi=150)
print("OK", ruta)
print("2/|lmax|: min %.3f, max %.3f" % (h_2_sobre_lmax.min(), h_2_sobre_lmax.max()))
print("h_crit complejo: min %.3f, max %.3f" % (h_crit.min(), h_crit[np.isfinite(h_crit)].max()))
