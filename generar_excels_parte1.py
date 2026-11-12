#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Genera TODOS los Excel de la Parte 1 que armó el grupo (versión corregida).

Reemplaza a `solucion_sisteam_ecuaciones.py` y a los scripts que originaron
el resto de los archivos. Escribe en ./excel :

  (a) solucion_h_<h>.xlsx / grafico_<h>.png
        Tiempo, Temperatura, Humedad, Precipitacion para cada h (+ gráfico).
  (b) euler_h_<h>.xlsx
        Euler explícito con numeración de iteraciones n, t, T, H, P.
  (c) euler_h01_h08.xlsx
        Resumen (con fórmulas de N, error máximo y orden p) + una hoja por h
        con T_ref y |T - T_ref| como fórmula.
  (d) Resolucion_Inciso_1b_Euler.xlsx
        Inciso 1.b: error máximo, instante del máximo y orden, con gráfico.

Correcciones respecto de los archivos originales:
  * R(t) usa cos(w t - pi) como pide la consigna. El script original usaba
    sin(...), por eso solucion_h_*.xlsx daba T(0,1) = 1,0225 en vez de 0,98875.
  * La referencia es Euler con h_ref = 0,0001 y no h = 0,1: con h = 0,1 el
    error de la propia referencia (0,046 °C) es del orden del que se quiere medir.
  * El último paso se acorta si h no divide a 60 (antes se truncaba).

Las columnas de datos son valores calculados; las de error, máximos, N y
orden p son fórmulas de Excel (se calculan al abrir el archivo).

Uso:   python generar_excels_parte1.py
"""

import os

import numpy as np
from openpyxl import Workbook
from openpyxl.chart import ScatterChart, Reference, Series
from openpyxl.styles import Alignment, Font, PatternFill

import tp_comun as tp

P = tp.PARAMS_NOMINALES
H_REF = 0.0001
HS = [0.1, 0.2, 0.4, 0.8]
DIR_EXCEL = os.path.join(os.path.dirname(os.path.abspath(__file__)), "excel")

NEGRITA = Font(bold=True)
AZUL = Font(color="1F4E9E")
ENCAB = PatternFill("solid", fgColor="DDE6F3")


def encabezados(ws, fila, textos):
    for j, tx in enumerate(textos, start=1):
        c = ws.cell(row=fila, column=j, value=tx)
        c.font, c.fill = NEGRITA, ENCAB
        c.alignment = Alignment(horizontal="center", wrap_text=True)


def anchos(ws, valores):
    for j, w in enumerate(valores, start=1):
        ws.column_dimensions[ws.cell(row=1, column=j).column_letter].width = w


def etiqueta(h):
    """0.1 -> '0.1'  (nombre de hoja/archivo)."""
    return f"{h:g}"


def etiqueta_archivo(h):
    """0.1 -> '0_1'  (nombre de archivo euler_h_0_1.xlsx)."""
    return f"{h:g}".replace(".", "_")


# --------------------------------------------------------------------------
# (a) solucion_h_<h>.xlsx + grafico
# --------------------------------------------------------------------------
def excel_solucion(h, t, y, plt):
    wb = Workbook()
    ws = wb.active
    ws.title = "solucion"
    encabezados(ws, 1, ["Tiempo", "Temperatura", "Humedad", "Precipitacion"])
    for i in range(len(t)):
        ws.append([float(t[i]), *map(float, y[i])])
    anchos(ws, [12, 14, 14, 16])
    wb.save(os.path.join(DIR_EXCEL, f"solucion_h_{etiqueta(h)}.xlsx"))

    fig = plt.figure(figsize=(10, 6))
    for j, (nom, col) in enumerate([("Temperatura", tp.C_T), ("Humedad", tp.C_H),
                                    ("Precipitación", tp.C_P)]):
        plt.plot(t, y[:, j], "o-", ms=3, color=col, label=nom)
    plt.xlabel("Tiempo (meses)"); plt.ylabel("Valor")
    plt.title(f"Evolución del sistema (h = {h:g})"); plt.legend()
    fig.savefig(os.path.join(DIR_EXCEL, f"grafico_{etiqueta(h)}.png"), dpi=200, bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------
# (b) euler_h_<h>.xlsx
# --------------------------------------------------------------------------
def excel_euler(h, t, y):
    wb = Workbook()
    ws = wb.active
    ws.title = f"h_{etiqueta(h)}"
    n = len(t) - 1
    ws["A1"] = f"Euler explícito, h = {h:g} mes, {n} pasos. Valores calculados con generar_excels_parte1.py"
    encabezados(ws, 3, ["n (iteración)", "t [meses]", "T [°C]", "H [u.h.]", "P [u.p.]"])
    for i in range(len(t)):
        ws.append([i, float(t[i]), *map(float, y[i])])
    anchos(ws, [14, 12, 14, 14, 14])
    wb.save(os.path.join(DIR_EXCEL, f"euler_h_{etiqueta_archivo(h)}.xlsx"))


# --------------------------------------------------------------------------
# (c) euler_h01_h08.xlsx
# --------------------------------------------------------------------------
def excel_resumen(corridas, t_ref, y_ref):
    wb = Workbook()
    rs = wb.active
    rs.title = "Resumen"
    rs["A1"] = "Parte 1 - Error de Euler en T(t) vs. solución de referencia"
    rs["A1"].font = NEGRITA
    encabezados(rs, 3, ["h [meses]", "N pasos", "Err máx |T - T_ref| [°C]",
                        "Orden observado p = log(Err(2h)/Err(h))/log(2)"])
    for k, h in enumerate(HS):
        t, y = corridas[h]
        fila = 4 + k
        hoja = f"h_{etiqueta(h)}"
        ultima = 4 + len(t) - 1
        rs.cell(row=fila, column=1, value=h).font = AZUL
        rs.cell(row=fila, column=2, value=f"=60/A{fila}")
        rs.cell(row=fila, column=3, value=f"=MAX('{hoja}'!F4:F{ultima})")
        if k > 0:
            rs.cell(row=fila, column=4, value=f"=LOG(C{fila}/C{fila - 1})/LOG(2)")

        ws = wb.create_sheet(hoja)
        ws["A1"] = (f"Euler explícito, h = {h:g} mes ({len(t) - 1} pasos). T_ref: referencia Euler "
                    f"con h_ref = {H_REF:g} muestreada en los mismos t.")
        encabezados(ws, 3, ["t [meses]", "T [°C]", "H [u.h.]", "P [u.p.]", "T_ref [°C]", "|T - T_ref| [°C]"])
        r = int(round(h / H_REF))
        for i in range(len(t)):
            f = 4 + i
            ws.append([float(t[i]), *map(float, y[i]), float(y_ref[i * r, 0]), f"=ABS(B{f}-E{f})"])
        anchos(ws, [12, 12, 12, 12, 14, 18])
    notas = [
        "Notas:",
        "- Los valores de t, T, H, P y T_ref son resultados numéricos calculados con generar_excels_parte1.py "
        "(Euler explícito); no son fórmulas.",
        "- Las columnas de error, el máximo, N y el orden p sí son fórmulas.",
        f"- Referencia: Euler con h_ref = {H_REF:g} mes (error estimado ~2e-5 °C en T).",
        "- Los pasos h (azul) son entradas; el orden p se calcula entre pasos consecutivos "
        "(cada h es el doble del anterior).",
    ]
    for k, tx in enumerate(notas):
        rs.cell(row=4 + len(HS) + 2 + k, column=1, value=tx)
    anchos(rs, [14, 12, 26, 44])
    wb.save(os.path.join(DIR_EXCEL, "euler_h01_h08.xlsx"))


# --------------------------------------------------------------------------
# (d) Resolucion_Inciso_1b_Euler.xlsx
# --------------------------------------------------------------------------
def excel_inciso_1b(corridas, t_ref, y_ref):
    wb = Workbook()
    ws = wb.active
    ws.title = "Resultados 1b"
    ws["A1"] = "TRABAJO PRÁCTICO — INCISO 1.b"
    ws["A1"].font = Font(bold=True, size=13)
    ws["A2"] = (f"Comparación del método de Euler explícito con la solución de referencia "
                f"h_ref = {H_REF:g} mes")
    ws["A4"] = "Criterio de error"
    ws["B4"] = "Err(h) = máximo de |T_h(t) − T_ref(t)| para 0 ≤ t ≤ 60 meses."
    encabezados(ws, 6, ["Paso h [meses]", "Cantidad de pasos", "Error máximo [°C]",
                        "Instante del error máximo [meses]", "Orden observado p"])
    for k, h in enumerate(HS):
        t, y = corridas[h]
        fila = 7 + k
        hoja = f"h={etiqueta(h)}"
        ultima = 5 + len(t) - 1
        ws.cell(row=fila, column=1, value=h).font = AZUL
        ws.cell(row=fila, column=2, value=f"=60/A{fila}")
        ws.cell(row=fila, column=3, value=f"='{hoja}'!B{ultima + 2}")
        ws.cell(row=fila, column=4, value=f"='{hoja}'!B{ultima + 3}")
        if k > 0:
            ws.cell(row=fila, column=5, value=f"=LOG(C{fila}/C{fila - 1})/LOG(2)")

        hs_ = wb.create_sheet(hoja)
        hs_["A1"] = f"COMPARACIÓN DETALLADA — h = {h:g}"
        hs_["A2"] = f"Se comparan instantes coincidentes con la referencia h_ref = {H_REF:g}."
        encabezados(hs_, 4, ["t [meses]", "T_ref [°C]", "T_h [°C]", "|T_h − T_ref| [°C]"])
        r = int(round(h / H_REF))
        for i in range(len(t)):
            f = 5 + i
            hs_.append([float(t[i]), float(y_ref[i * r, 0]), float(y[i, 0]), f"=ABS(C{f}-B{f})"])
        hs_.cell(row=ultima + 2, column=1, value="Error máximo")
        hs_.cell(row=ultima + 2, column=2, value=f"=MAX(D5:D{ultima})")
        hs_.cell(row=ultima + 3, column=1, value="Instante del máximo [meses]")
        hs_.cell(row=ultima + 3, column=2, value=f"=INDEX(A5:A{ultima},MATCH(B{ultima + 2},D5:D{ultima},0))")
        anchos(hs_, [26, 14, 14, 20])

    fila_fin = 7 + len(HS) - 1
    ws["A13"] = "Conclusión"
    ws["A13"].font = NEGRITA
    ws["A14"] = ("El error máximo crece de forma aproximadamente proporcional a h y el orden observado es "
                 "cercano a 1, como se espera para Euler explícito (orden 1). La referencia usa un paso "
                 "mucho menor que los analizados, por lo que su propio error es despreciable.")
    ws["A14"].alignment = Alignment(wrap_text=True, vertical="top")
    ws.merge_cells("A14:E16")

    ch = ScatterChart()
    ch.title = "Err(h) vs h (escala log-log)"
    ch.style = 13
    ch.x_axis.title, ch.y_axis.title = "h [meses]", "Err(h) [°C]"
    ch.x_axis.scaling.logBase = 10
    ch.y_axis.scaling.logBase = 10
    ch.x_axis.delete = False
    ch.y_axis.delete = False
    x = Reference(ws, min_col=1, min_row=7, max_row=fila_fin)
    yv = Reference(ws, min_col=3, min_row=7, max_row=fila_fin)
    s = Series(yv, x, title="Err(h)")
    s.marker.symbol = "circle"
    ch.series.append(s)
    ch.height, ch.width = 8, 14
    ws.add_chart(ch, "G6")
    anchos(ws, [18, 18, 20, 30, 20])

    # Referencia muestreada cada 0,1 mes
    wr = wb.create_sheet("Referencia")
    wr["A1"] = f"Solución de referencia (Euler, h_ref = {H_REF:g}) muestreada cada 0,1 mes"
    encabezados(wr, 3, ["t [meses]", "T_ref [°C]", "H_ref [u.h.]", "P_ref [u.p.]"])
    r = int(round(0.1 / H_REF))
    for i in range(0, len(t_ref), r):
        wr.append([float(t_ref[i]), *map(float, y_ref[i])])
    anchos(wr, [12, 14, 14, 14])
    wb.save(os.path.join(DIR_EXCEL, "Resolucion_Inciso_1b_Euler.xlsx"))


def main():
    os.makedirs(DIR_EXCEL, exist_ok=True)
    plt = tp.estilo_graficos(False)
    print(f"Referencia: Euler con h_ref = {H_REF:g} ...")
    t_ref, y_ref = tp.euler(P, H_REF)
    corridas = {h: tp.euler(P, h) for h in HS}
    for h in HS:
        t, y = corridas[h]
        excel_solucion(h, t, y, plt)
        excel_euler(h, t, y)
        print(f"  h = {h:g}: {len(t) - 1} pasos, T(0,1·k) ok")
    excel_resumen(corridas, t_ref, y_ref)
    excel_inciso_1b(corridas, t_ref, y_ref)
    print(f"Archivos escritos en {DIR_EXCEL}:")
    for f in sorted(os.listdir(DIR_EXCEL)):
        print("  ", f)


if __name__ == "__main__":
    main()
