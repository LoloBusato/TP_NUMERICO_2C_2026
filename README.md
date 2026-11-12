# TP_NUMERICO_2C_2026

TP "El superniño" – CB051 Modelación Numérica / 9504 Análisis Numérico I (2do cuatrimestre 2026, FIUBA).
Modelo simplificado de El Niño (T, H, P) resuelto con Euler explícito.

**El informe está en [INFORME.md](INFORME.md)** (con figuras y fórmulas). Los scripts son ese informe ejecutado en orden: cada número y cada figura del informe sale de ellos.

## Cómo correrlo

```
pip install -r requirements.txt
python ejecutar_todo.py            # completo, ~5 min
python ejecutar_todo.py --rapido   # grilla de 3.3 con h = 0.01, ~1 min
```

Cada script también se puede correr por separado (en este orden).

| Script | Qué hace | Sección del informe |
|---|---|---|
| `tp_comun.py` | Parámetros, modelo, Euler explícito, RK4 (solo control), indicadores SP y ΔtTP, jacobiano | — |
| `parte1_euler.py` | Referencia, errores, orden observado, estabilidad, umbral y análisis con el jacobiano / Floquet | 2.1 – 2.4 |
| `parte2_severidad.py` | Caso nominal, sensibilidad a A, β, γ, σ, grilla de 2401 combinaciones, comparación | 3.1 – 3.4 |
| `generar_excels_parte1.py` | Todos los Excel de la Parte 1 (en `excel/`) | 2.2, 2.3 |
| `generar_informe.py` | Escribe `INFORME.md` a partir de `resultados/*.json` | — |
| `ejecutar_todo.py` | Corre todo en orden | — |

Salidas: `resultados/` (figuras PNG, tablas CSV, JSON con los números) y `excel/`.

## Resultados principales

- **Parte 1:** con referencia `h_ref = 0,0001` el orden observado de Euler es ≈ 1,02 (esperado 1). Es estable hasta `h = 1,25`; con `h = 1,26` diverge dentro de los 60 meses. El análisis lineal (jacobiano / Floquet) predice la pérdida de estabilidad de la solución anual cerca de `h ≈ 1,0`.
- **Parte 2** (`h = 0,001`): nominal `SP = 0,4285`, `Δt = 0,83 meses`. Mayor severidad con `A = 0,60`, `β = 0,30`, `γ = 0,05`, `σ = 1,00` (`SP = 0,9908`). No es "todo al máximo": `γ` va en su mínimo.

## Qué se corrigió respecto de los scripts anteriores

Los tres scripts originales (`solucion_sisteam_ecuaciones.py`, `tp_el_nino_parte2a.py`, `variacion_constantes.py`) se reemplazaron por los de arriba (siguen en el historial de git). Cambios:

- `solucion_sisteam_ecuaciones.py` usaba `sin(ωt − π)` en `R(t)`; la consigna pide `cos(ωt − π)`. Por eso `solucion_h_*.xlsx` daba `T(0,1) = 1,0225` en vez de `0,98875`. Además truncaba el último paso cuando `h` no divide a 60.
- `variacion_constantes.py` iba modificando el diccionario de parámetros sin restaurarlo: el barrido de `σ` se hacía con `A = 0,595`, el de `γ` con `A` y `σ` ya cambiados, etc., en vez de "un parámetro por vez con los demás nominales". Tampoco llegaba a los extremos de los intervalos (`A` hasta 0,595; `σ` hasta 0,99; `γ` hasta 0,199; `β` hasta 0,295) y los gráficos no tenían el valor del parámetro en el eje x.
- `tp_el_nino_parte2a.py` usa `np.trapezoid`, que no existe en numpy < 2.0; ahora hay compatibilidad con ambas versiones. Sus indicadores (promedio temporal con trapecios y máximos refinados con parábola) se mantuvieron en `tp_comun.py`.
- La solución de referencia de la Parte 1 era Euler con `h = 0,1`, cuyo error propio (0,046 °C) es del orden del que se quiere medir; daba órdenes de 1,61 y 1,28. Ahora la referencia es Euler con `h = 0,0001`.
- El "paso crítico teórico 2,35" salía de `h ≤ 2/|λmax|`, válido solo para autovalores reales negativos; acá los autovalores relevantes son complejos y el criterio correcto da ≈ 1,0.
