#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Ejecuta el informe completo en orden:

    1. parte1_euler.py            -> secciones 2.1 a 2.4
    2. parte2_severidad.py        -> secciones 3.1 a 3.4
    3. generar_excels_parte1.py   -> Excel de la Parte 1 (carpeta excel/)
    4. generar_informe.py         -> INFORME_autogenerado.md con los números calculados

Uso:   python ejecutar_todo.py            # completo (~5 min)
       python ejecutar_todo.py --rapido   # grilla de 3.3 con h = 0.01 (~1 min)
"""

import subprocess
import sys
import time
from pathlib import Path

RAIZ = Path(__file__).resolve().parent


def correr(script, *args):
    print("\n" + "#" * 78 + f"\n# {script} {' '.join(args)}\n" + "#" * 78, flush=True)
    t0 = time.time()
    r = subprocess.run([sys.executable, str(RAIZ / script), *args], cwd=RAIZ)
    if r.returncode != 0:
        sys.exit(f"\n{script} terminó con error (código {r.returncode}).")
    print(f"\n[{script}: {time.time() - t0:.0f} s]", flush=True)


def main():
    rapido = ["--rapido"] if "--rapido" in sys.argv[1:] else []
    t0 = time.time()
    correr("parte1_euler.py")
    correr("parte2_severidad.py", *rapido)
    correr("generar_excels_parte1.py")
    correr("generar_informe.py")
    print(f"\nListo en {time.time() - t0:.0f} s. Resultados en resultados/, excel/ e INFORME_autogenerado.md")


if __name__ == "__main__":
    main()
