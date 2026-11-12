#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Escribe INFORME_autogenerado.md (no pisa INFORME.md, que se edita a mano) a partir de los resultados calculados por parte1_euler.py
y parte2_severidad.py (resultados/*.json y resultados/*.png).

Todos los números del texto salen de los JSON: nada se copia a mano.
El texto va en plantillas con marcadores @clave@ que se reemplazan por los
valores calculados (así el LaTeX queda escrito tal cual).

Uso:   python generar_informe.py        (después de correr las Partes 1 y 2)
"""

import inspect
import json
import math
import os
import re

import tp_comun as tp

DIR = tp.DIR_RESULTADOS
RAIZ = os.path.dirname(os.path.abspath(__file__))
REPO_URL = "https://github.com/agustinlucchelli/TP_NUMERICO_2C_2026"

V = {}          # clave -> texto que reemplaza a @clave@
CONT = {"fig": 0, "tab": 0}


def c(x, nd=4):
    """Número con coma decimal."""
    return tp.fmt(x, nd)


def cs(x, nd=2):
    """Notación científica con coma decimal."""
    return ("%.*e" % (nd, x)).replace(".", ",")


def cargar(nombre):
    with open(os.path.join(DIR, nombre), encoding="utf-8") as f:
        return json.load(f)


def tabla_md(enc, filas):
    s = "| " + " | ".join(enc) + " |\n"
    s += "|" + "|".join(["---"] * len(enc)) + "|\n"
    for f in filas:
        s += "| " + " | ".join(str(x) for x in f) + " |\n"
    return s


def tabla(clave, titulo, enc, filas):
    """Tabla numerada; guarda V[clave] (título + tabla) y V[clave+'_n'] (número)."""
    CONT["tab"] += 1
    V[clave + "_n"] = str(CONT["tab"])
    V[clave] = "**Tabla %d.** %s\n\n%s" % (CONT["tab"], titulo, tabla_md(enc, filas))


def tabla_sin_titulo(clave, enc, filas):
    V[clave] = tabla_md(enc, filas)


def figura(clave, archivo, texto):
    CONT["fig"] += 1
    V[clave] = "![Figura %d](resultados/%s)\n\n*Figura %d: %s*\n" % (CONT["fig"], archivo, CONT["fig"], texto)


def eig_texto(eig):
    partes, hecho = [], False
    for re, im in sorted(eig, key=lambda z: (abs(z[1]), z[0])):
        if abs(im) < 1e-9:
            partes.append(c(re))
        elif im > 0 and not hecho:
            partes.append("%s ± %si" % (c(re), c(im)))
            hecho = True
    return "; ".join(partes)


def render(plantilla):
    for k, v in V.items():
        plantilla = plantilla.replace("@" + k + "@", str(v))
    sin = set(re.findall(r"@([A-Za-z0-9_]+)@", plantilla))
    if sin:
        raise KeyError("marcadores sin reemplazar: %s" % sorted(sin))
    return plantilla


# ==========================================================================
CARATULA = r"""# TRABAJO PRÁCTICO "EL SUPER NIÑO"

**UNIVERSIDAD DE BUENOS AIRES**
**FACULTAD DE INGENIERÍA – DPTO. DE COMPUTACIÓN**

| | |
|---|---|
| **T.P. N°** | 1 |
| **Asignatura** | Modelación numérica (CB-051) / Análisis Numérico I (9504) |
| **Cuatrimestre** | 2do. Cuatrimestre 2026 |
| **Grupo de trabajo práctico** | N° 1 |
| **Fecha de realización** | 1/11/2026 |
| **Fecha de entrega** | 13/11/2026 |

**Integrantes del trabajo práctico**

| Nombre y apellido | Padrón |
|---|---|
| Ezequiel Ariel Miño | 107646 |
| Agustín Lucchelli | 109337 |
| Lorenzo Busato | 110490 |

---

## Índice

1. [Introducción al trabajo práctico](#1-introducción-al-trabajo-práctico)
2. [Parte 1: Convergencia, orden y estabilidad del método de Euler explícito](#2-parte-1-convergencia-orden-y-estabilidad-del-método-de-euler-explícito)
   - 2.1. Solución de referencia
   - 2.2. Resolución con sucesivos pasos h
   - 2.3. Verificación experimental
   - 2.4. Resolución con el paso temporal aumentado
3. [Parte 2: Severidad estacional, desfase y sensibilidad paramétrica](#3-parte-2-severidad-estacional-desfase-y-sensibilidad-paramétrica)
   - 3.1. Determinación de parámetros nominales
   - 3.2. Análisis e influencia de los factores
   - 3.3. Combinaciones posibles de los factores
   - 3.4. Comparación de casos
4. [Código de Python](#4-código-de-python)
5. [Bibliografía y fuentes](#5-bibliografía-y-fuentes)

---
"""

INTRO = r"""## 1. Introducción al trabajo práctico

El fenómeno El Niño – Oscilación del Sur es un ejemplo de sistema dinámico acoplado: una perturbación térmica del océano puede producir, con cierto retardo, respuestas en la humedad atmosférica y en la precipitación. En este trabajo se utiliza un modelo temporal simplificado de tres variables para estudiar, por un lado, propiedades numéricas del método de Euler explícito aplicado a un problema de valores iniciales (convergencia, orden y estabilidad) y, por otro, la respuesta del modelo frente a cambios en sus parámetros.

Las variables del modelo son $T(t)$, la anomalía de temperatura superficial del mar [°C]; $H(t)$, la anomalía de humedad en la atmósfera baja [u.h.]; y $P(t)$, una variable simplificada asociada a la precipitación o nubosidad [u.p.]. Se resuelve el sistema de ecuaciones diferenciales ordinarias no lineales

$$\frac{dT}{dt} = R(t)\,T\left(1-\frac{T}{T_{max}}\right) - \alpha H$$

$$\frac{dH}{dt} = \beta T^2 - \gamma H P$$

$$\frac{dP}{dt} = \delta H^2 - \sigma P$$

con una excitación estacional de período anual en la tasa de crecimiento térmico:

$$R(t) = R_{oce} + A\cos(\omega t - \pi), \qquad \omega = \frac{2\pi}{12}$$

Cada término tiene un papel físico: $R(t)T(1-T/T_{max})$ es un crecimiento logístico de la temperatura modulado por las estaciones; $-\alpha H$ es el enfriamiento por evaporación; $\beta T^2$ es la fuente de humedad (un mar más caliente evapora más); $-\gamma HP$ es el consumo de humedad por la lluvia; $\delta H^2$ es la conversión de humedad en lluvia y $-\sigma P$ es la descarga o disipación de la precipitación sobre el continente.

**Parámetros nominales**

| Parámetro | Significado físico | Valor nominal |
|---|---|---|
| $T_{max}$ | Anomalía máxima térmicamente sostenible | 4,0 °C |
| $R_{oce}$ | Tasa intrínseca de calentamiento oceánico | 0,30 mes⁻¹ |
| $\alpha$ | Coeficiente de enfriamiento por evaporación | 0,05 |
| $\beta$ | Eficiencia de evaporación | 0,20 |
| $\gamma$ | Tasa de consumo de humedad por precipitación | 0,10 |
| $\delta$ | Eficiencia de conversión de humedad a lluvia | 0,40 |
| $\sigma$ | Tasa de descarga/disipación sobre el continente | 0,70 mes⁻¹ |
| $A$ | Amplitud de la modulación estacional | 0,45 mes⁻¹ |
| $\omega$ | Frecuencia angular de un ciclo anual | $2\pi/12$ rad/mes |

Condiciones iniciales: $T(0)=1{,}0$ °C, $H(0)=0$ u.h., $P(0)=0$ u.p. Intervalo de integración: $0 \le t \le 60$ meses.

El trabajo se divide en dos partes. En la Parte 1 se analiza el comportamiento numérico del método de Euler explícito. En la Parte 2 se usa una discretización suficientemente precisa y estable para estudiar la severidad estacional de las precipitaciones, el desfase entre temperatura y precipitación y la sensibilidad a los parámetros $A$, $\beta$, $\gamma$ y $\sigma$. Todo el código está en Python (sección 4): cada script reproduce, en el mismo orden, lo que se presenta en este informe.

---
"""

P1_21 = r"""## 2. Parte 1: Convergencia, orden y estabilidad del método de Euler explícito

### 2.1. Solución de referencia

Se resuelve el sistema con el método de Euler explícito,

$$y_{n+1} = y_n + h\,f(t_n, y_n), \qquad y = [T, H, P]^T$$

que, escrito por componentes, es

$$T_{n+1} = T_n + h\left[R(t_n)\,T_n\left(1-\frac{T_n}{T_{max}}\right) - \alpha H_n\right]$$

$$H_{n+1} = H_n + h\left[\beta T_n^2 - \gamma H_n P_n\right]$$

$$P_{n+1} = P_n + h\left[\delta H_n^2 - \sigma P_n\right]$$

con $R(t_n) = R_{oce} + A\cos(\omega t_n - \pi)$.

La consigna pide una solución de referencia, también obtenida con Euler, pero con un paso $h_{ref}$ tan pequeño que su error sea despreciable frente al de las discretizaciones que se analizan. Se adoptó

$$h_{ref} = @h_ref@\ \text{meses}$$

lo que equivale a @pasos_ref@ pasos. Para verificar que este paso alcanza, se repitió el cálculo con $h_{ref}/2$ y se comparó: la diferencia máxima en $T$ es de @control_ref@ °C. Como el menor error que se mide más adelante es @err_min@ °C, el error de la referencia es del orden del @pct_ref@ % del menor error analizado, por lo que es despreciable.

"""

P1_22 = r"""### 2.2. Resolución con sucesivos pasos h

Se resolvió nuevamente el problema con los pasos $h = 0{,}025;\ 0{,}05;\ 0{,}1;\ 0{,}2;\ 0{,}4;\ 0{,}8$ meses, donde cada valor es el doble del anterior. Todos son múltiplos de $h_{ref}$, de modo que las soluciones se comparan en los mismos instantes. El error se cuantifica con el máximo error absoluto de $T$ sobre todo el intervalo:

$$Err(h) = \max_{0 \le t \le 60} \left| T_h(t) - T_{ref}(t) \right|$$

@tab_errores@
@fig_errores@
El error crece al aumentar el paso y lo hace de forma aproximadamente proporcional a $h$: cuando $h$ se duplica, el error se duplica aproximadamente. El máximo error se produce siempre alrededor de $t \approx$ @t_err_min@–@t_err_max@ meses, que es donde $T$ crece más rápido (interpretación: Euler utiliza la pendiente al inicio de cada paso y se retrasa donde la pendiente cambia más).

"""

P1_23 = r"""### 2.3. Verificación experimental

En el rango convergente se calcula el orden observado entre cada par de pasos consecutivos:

$$p = \frac{\log\left(Err(2h)/Err(h)\right)}{\log 2}$$

@tab_orden@
El orden observado es $p \approx$ @p_medio@ en promedio (entre @p_min@ y @p_max@), prácticamente igual al valor teórico $p = 1$ del método de Euler explícito (error global $O(h)$). El orden tiende a 1 a medida que $h$ disminuye, que es el comportamiento asintótico esperado.

@fig_orden@
**Observación sobre la elección de la referencia.** Si se tomara como referencia la solución de Euler con $h = 0{,}1$, los errores de $h = 0{,}2;\ 0{,}4;\ 0{,}8$ serían @e01_02@; @e01_04@ y @e01_08@ °C, y los órdenes resultarían superiores a 1 (aproximadamente @p01_a@ y @p01_b@). Eso no indica un problema del método sino de la referencia: Euler con $h = 0{,}1$ tiene un error propio de @err_01@ °C respecto de la referencia fina, comparable al error que se quiere medir en $h = 0{,}2$ (@err_02@ °C). Por eso la consigna exige una referencia cuyo error sea despreciable.

"""

P1_24 = r"""### 2.4. Resolución con el paso temporal aumentado

#### Resultados experimentales

Se continuó aumentando el paso. Se consideró que una simulación **diverge** cuando alguna variable supera $10^6$ en módulo (los valores físicos son de orden 1 a 10).

@tab_estab@
@fig_estab@
Hasta $h = 1{,}25$ la simulación completa los 60 meses y el error sigue creciendo de manera suave; con $h = 1{,}26$ la solución explota (cruza $10^6$ en $t \approx$ @t_div_126@ meses). Un refinamiento por bisección ubica el umbral en $h^* \approx$ @umbral@ meses, y barriendo de a 0,01 entre 1,25 y 1,40 se observa que solo $h = 1{,}25$ es estable. Entre los pasos de la sucesión doble (0,8 → 1,6) se encuentra entonces el inicio de la inestabilidad.

Además, antes de divergir, la solución numérica pierde el carácter periódico anual. La siguiente tabla muestra $T(t = 12k)$, que en la solución de referencia se repite año a año:

@tab_period@
Con $h \le 0{,}8$ los valores se estabilizan (para $h = 0{,}8$: @T_08_estab@), mientras que desde $h \approx 1{,}0$ ya no se repiten de un año al otro aunque la solución siga acotada.

#### Error de discretización vs. inestabilidad numérica

- **Error de discretización**: aparece porque Euler aproxima la solución continua con pasos finitos. Para $h \le 1{,}25$ la solución sigue la forma de la referencia y el error crece de forma suave y aproximadamente proporcional a $h$ (orden 1). Se reduce achicando el paso.
- **Inestabilidad numérica**: aparece cuando el esquema amplifica las perturbaciones en cada paso. La solución numérica deja de parecerse a la del problema continuo: cambia de signo, oscila con amplitud creciente y toma valores de orden $10^6$ o más. No es "más error": es otro comportamiento, y el error deja de ser proporcional a $h$.

#### Referencia teórica con el jacobiano

Para tener un valor teórico con el cual comparar el umbral experimental se analiza el jacobiano del sistema,

$$J(t,T,H,P) = \begin{pmatrix} R(t)\left(1-\dfrac{2T}{T_{max}}\right) & -\alpha & 0 \\ 2\beta T & -\gamma P & -\gamma H \\ 0 & 2\delta H & -\sigma \end{pmatrix}$$

evaluado sobre la solución. Para Euler explícito, cada autovalor $\lambda$ con $\mathrm{Re}(\lambda) < 0$ exige $|1 + h\lambda| < 1$, es decir

$$h < h_{crit}(\lambda) = \frac{-2\,\mathrm{Re}(\lambda)}{|\lambda|^2}$$

(para autovalores reales negativos se reduce a $h < 2/|\lambda|$). Un autovalor con parte real positiva representa crecimiento real del sistema, no inestabilidad numérica, por lo que el criterio solo se aplica a los modos con $\mathrm{Re}(\lambda) < 0$.

@tab_autoval@
En $t = 60$ el par complejo da $h_{crit}(60) =$ @hcrit_60@ meses. Sin embargo, este valor es solo local (con coeficientes congelados) y el sistema es no lineal y no autónomo, por lo que no garantiza estabilidad. Sobre todo el intervalo el mínimo puntual es @hcrit_min@ meses en $t \approx$ @hcrit_min_t@, en el transitorio, donde un par de autovalores se acerca al eje imaginario y el límite local tiende a cero; ese valor no es representativo del régimen periódico. En el último ciclo ($48 \le t \le 60$) el límite puntual varía entre @hcrit_ult_min@ y @hcrit_ult_max@ meses.

Una referencia más apropiada para un régimen periódico es linealizar el mapa de Euler a lo largo de un ciclo anual $[48, 60]$, $M = \prod_n (I + h\,J_n)$, y estudiar su radio espectral $\rho(M)$ (multiplicadores de Floquet). La órbita anual es numéricamente estable mientras $\rho(M) < 1$. Con la solución exacta ($h \to 0$) resulta $\rho =$ @floquet_exacto@.

@tab_floquet@
@fig_jac@
El radio espectral cruza 1 entre $h =$ @fl_estable@ y $h =$ @fl_inestable@ meses, de modo que la teoría lineal predice la pérdida de estabilidad de la solución anual alrededor de **$h \approx 1{,}0$**, coherente con $h_{crit}(60) =$ @hcrit_60_3@ y con la pérdida de periodicidad observada.

#### Conclusiones de la Parte 1

- El error de Euler explícito es proporcional a $h$ (orden observado ≈ @p_medio_2@, esperado 1).
- La estabilidad numérica se pierde en dos etapas: alrededor de $h \approx 1{,}0$ la solución deja de repetirse cada año (concordante con el análisis del jacobiano) y para $h > 1{,}25$ explota dentro de los 60 meses. Es una interpretación, no verificada en detalle, que los términos no lineales retardan la explosión entre ambos valores.
- Se debe elegir un paso muy por debajo de estos límites. Para la Parte 2 se adopta $h =$ @h2@ meses, unas @razon_h@ veces menor que el umbral experimental.

---
"""

P2_31 = r"""## 3. Parte 2: Severidad estacional, desfase y sensibilidad paramétrica

### 3.1. Determinación de parámetros nominales

Debido al régimen transitorio asociado a las condiciones iniciales, los indicadores se calculan únicamente con el último ciclo anual, $48 \le t \le 60$ meses. La severidad estacional de las precipitaciones y el desfase temperatura–precipitación se definen como

$$S_P = \frac{P_{max} - P_{min}}{P_{med}}, \qquad \Delta t_{TP} = t_{P,max} - t_{T,max}$$

donde $P_{max}$, $P_{min}$ y $P_{med}$ se calculan en el último ciclo ($P_{med}$ es el promedio temporal, por integral trapezoidal dividida por la duración), $t_{T,max}$ es el instante del máximo de $T$ y $t_{P,max}$ el del máximo de $P$. Los máximos se localizan sobre la grilla y se refinan con una parábola por tres puntos.

**Elección del paso.** Según la Parte 1, Euler es estable hasta $h \approx 1{,}25$ y su error es proporcional a $h$. Para esta parte se verificó cómo cambian los indicadores al achicar el paso, comparando contra una solución de Runge-Kutta de orden 4 con $h = 0{,}0005$ (solo como control independiente):

@tab_verif@
Con $h =$ @h2@ el error máximo en $T$ es de @err_h2@ °C y los indicadores coinciden con la referencia en @dif_SP@ para $S_P$ y en @dif_dt@ meses para $\Delta t_{TP}$. Se adopta entonces $h =$ @h2@ meses.

**Resultados con parámetros nominales** (Euler explícito, $h =$ @h2@):

@tab_nominal@
@fig_60m@
@fig_ult@
Luego de un transitorio de unos 2 a 3 años, las tres variables se repiten cada 12 meses. Se verificó que el último ciclo coincide con el anterior (diferencia relativa máxima de @desvio_nom@), por lo que usar el último ciclo es válido. El máximo de precipitación llega @dt_nom_2@ meses después del máximo de temperatura.

"""

P2_32 = r"""### 3.2. Análisis e influencia de los factores

Se modificó un parámetro por vez, manteniendo los restantes en su valor nominal, en los intervalos de la consigna (15 valores equiespaciados por parámetro).

@tab_sens@
@fig_sens@
**Tendencias observadas (datos calculados).** $S_P$ crece con $A$, $\beta$ y $\sigma$, y decrece con $\gamma$. $\Delta t_{TP}$ disminuye al aumentar cualquiera de los cuatro parámetros: la precipitación responde antes a la temperatura. Como $S_P$ es un cociente, conviene mirar también la amplitud $P_{max} - P_{min}$ y el promedio $P_{med}$ (Tabla @tab_sens_n@): con $A$ y $\beta$ crece la amplitud; con $\gamma$ cae la amplitud; con $\sigma$ la amplitud casi no cambia (@amp_sig_0@ → @amp_sig_1@) y lo que hace subir $S_P$ es que $P_{med}$ baja de @pmed_sig_0@ a @pmed_sig_1@.

**Interpretación a partir de los términos de las ecuaciones.** Estas lecturas son consistentes con los datos, pero los datos muestran la tendencia, no demuestran el mecanismo:

- **$A$** está en $R(t)$: una modulación estacional más fuerte hace oscilar más a $T$ y esa oscilación se transmite a $H$ y $P$. Por eso $S_P$ crece de forma casi lineal con $A$.
- **$\beta$** multiplica la fuente de humedad $\beta T^2$ (signo $+$): más humedad produce más lluvia vía $\delta H^2$. Crecen $P_{med}$ y la amplitud.
- **$\gamma$** multiplica el sumidero $-\gamma HP$ (signo $-$): cuanto mayor es, más humedad consume la lluvia, y $H$ y $P$ quedan más chicos y con menos contraste. Por eso $S_P$ disminuye con $\gamma$.
- **$\sigma$** multiplica el sumidero $-\sigma P$ (signo $-$): con $\sigma$ grande, $P$ decae rápido y sigue más de cerca a $H$ (menos inercia), lo que adelanta su máximo ($\Delta t_{TP}$ disminuye) y reduce su promedio.

"""

P2_33 = r"""### 3.3. Combinaciones posibles de los factores

Se exploró una grilla discreta de @n_total@ combinaciones ($7 \times 7 \times 7 \times 7$ valores equiespaciados de $A$, $\beta$, $\gamma$ y $\sigma$ dentro de sus intervalos), simulando cada una con $h =$ @h_grilla@ meses. Una combinación se considera válida si la solución es acotada y si al final de la simulación el régimen anual se repite (diferencia relativa máxima entre los dos últimos ciclos menor al 1 %). @validas_txt@

@tab_top@
La combinación de mayor severidad es **$A =$ @ps_A@, $\beta =$ @ps_beta@, $\gamma =$ @ps_gamma@, $\sigma =$ @ps_sigma@**, con $S_P =$ @sp_sev@. Las demás combinaciones del podio repiten $A = 0{,}60$ y $\gamma = 0{,}05$, en línea con las tendencias de 3.2. El óptimo queda en un vértice del dominio de búsqueda; esto es coherente con que $S_P$ sea monótono en cada parámetro en los barridos individuales, aunque la grilla no demuestra que no exista un máximo mayor entre nodos.

@fig_mapa@
"""

P2_34 = r"""### 3.4. Comparación de casos

@tab_comp@
@fig_comp@
@fig_sev@
La severidad se multiplica por @razon_sp@ (de @sp_nom_3@ a @sp_sev_3@): $P_{max}$ pasa de @pmax_nom@ a @pmax_sev@ y $P_{min}$ casi no cambia (@pmin_nom@ → @pmin_sev@). El desfase pasa de @dt_nom_2@ a @dt_sev_2@ meses, es decir que casi no se modifica.

**¿El escenario de mayor severidad corresponde a tomar todos los parámetros en su valor máximo? No.** $A$, $\beta$ y $\sigma$ están en su máximo, pero $\gamma$ está en su **mínimo**. Si se toman los cuatro en su máximo ($A = 0{,}60$, $\beta = 0{,}30$, $\gamma = 0{,}20$, $\sigma = 1{,}00$) resulta $S_P =$ @sp_todo_max@, mucho menor que @sp_sev@. Partiendo del escenario severo y llevando de a un parámetro al extremo opuesto, $S_P$ siempre empeora:

@tab_flips@
**Justificación por el signo y la función de cada término:**

- $A$ (en $R(t)$): refuerza la modulación estacional de la temperatura y por lo tanto la oscilación anual de toda la cadena $T \to H \to P$. Máximo $A$ favorece la severidad.
- $\beta$ (término $+\beta T^2$ en $dH/dt$): es la fuente de humedad que alimenta a la lluvia a través de $+\delta H^2$. Máximo $\beta$ favorece la severidad.
- $\sigma$ (término $-\sigma P$ en $dP/dt$): es un sumidero de $P$, pero con $\sigma$ grande $P$ responde rápido a $H$ y su promedio baja, lo que aumenta el cociente $(P_{max} - P_{min})/P_{med}$ (en los barridos casi no cambia la amplitud). Máximo $\sigma$ favorece la severidad.
- $\gamma$ (término $-\gamma HP$ en $dH/dt$): es un sumidero de $H$ que actúa cuando hay lluvia y atenúa los picos. Su efecto sobre $S_P$ tiene signo opuesto al de los otros tres, por lo que la severidad máxima se logra con $\gamma$ en su **mínimo**.

Entonces, no basta con "subir todo": hay que distinguir entre parámetros que **alimentan o amplifican** la cadena ($A$, $\beta$), parámetros que **disminuyen el promedio** de lluvia ($\sigma$) y parámetros que **amortiguan** ($\gamma$), cuyo efecto sobre la severidad es el opuesto.

---
"""

CODIGO = r"""## 4. Código de Python

Todo el código está en el repositorio del grupo: @repo@. Cada script es el informe ejecutado en orden y escribe sus figuras y tablas en `resultados/`.

| Archivo | Contenido | Sección del informe |
|---|---|---|
| `tp_comun.py` | Parámetros, modelo, Euler explícito, RK4 (control), indicadores, jacobiano | 2.1, 3.1 |
| `parte1_euler.py` | Referencia, errores, orden, estabilidad y análisis con el jacobiano | 2.1 a 2.4 |
| `parte2_severidad.py` | Nominal, barridos, grilla, comparación | 3.1 a 3.4 |
| `generar_excels_parte1.py` | Todos los Excel de la Parte 1 (en `excel/`) | 2.2, 2.3 |
| `generar_informe.py` | Escribe este informe (`INFORME.md`) con los números calculados | — |
| `ejecutar_todo.py` | Corre todo lo anterior en orden | — |

Ejecución completa (el paso de la grilla de 3.3 tarda unos minutos; `--rapido` lo acelera usando $h = 0{,}01$ solo en la grilla):

```
pip install -r requirements.txt
python ejecutar_todo.py
```

Función central, Euler explícito (`tp_comun.py`):

```python
@fuente_euler@
```

---

## 5. Bibliografía y fuentes

- Consigna del trabajo práctico, *CB051 Modelación Numérica / 9504 Análisis Numérico I*, 2do cuatrimestre 2026, FIUBA.
- Burden, R. L. y Faires, J. D., *Análisis numérico*, Cengage Learning.
- Chapra, S. C. y Canale, R. P., *Métodos numéricos para ingenieros*, McGraw-Hill.
- Documentación de NumPy y Matplotlib (https://numpy.org, https://matplotlib.org).
- Repositorio del grupo con todo el código: @repo@.
- Herramientas: Python 3, NumPy, Matplotlib, openpyxl. Se utilizó una herramienta de inteligencia artificial como apoyo para el desarrollo del código y la redacción; los resultados se generan con los scripts del repositorio y fueron contrastados con una solución de referencia independiente (RK4).
"""


def main():
    r1, r2 = cargar("resultados_parte1.json"), cargar("resultados_parte2.json")
    nom, sev, tmax = r2["nominal"], r2["severo"], r2["todo_maximo"]
    pn, ps = tp.PARAMS_NOMINALES, r2["severo"]["parametros"]
    h2 = r2["h"]
    err = {e["h"]: e for e in r1["errores"]}
    est = {e["h"]: e for e in r1["estabilidad"]}
    o = r1["ordenes"]
    bar = r2["barrido"]
    g = r2["grilla"]
    vh = {v["h"]: v for v in r2["verificacion_h"]}
    nombres = {"A": "A", "beta": "β", "gamma": "γ", "sigma": "σ"}
    V["repo"] = REPO_URL
    secciones = [CARATULA, INTRO]

    # ---------------- 2.1
    V["h_ref"] = c(r1["h_ref"], 4).replace(",", "{,}")
    V["pasos_ref"] = ("{:,}".format(r1["pasos_ref"])).replace(",", ".")
    V["control_ref"] = cs(r1["control_ref"])
    V["err_min"] = c(err[0.025]["err"])
    V["pct_ref"] = c(100 * r1["control_ref"] / err[0.025]["err"], 1)
    secciones.append(render(P1_21))

    # ---------------- 2.2
    tabla("tab_errores", "Error máximo en T respecto de la referencia.",
          ["h [meses]", "Pasos", "Err(h) [°C]", "Instante del error máximo [meses]"],
          [[c(e["h"], 3), e["pasos"], c(e["err"], 6), c(e["t_err"], 1)] for e in r1["errores"]])
    figura("fig_errores", "fig_1_2_errores.png",
           "Temperatura para distintos pasos frente a la referencia (izquierda) y error punto a punto (derecha).")
    V["t_err_min"] = c(min(e["t_err"] for e in r1["errores"]), 1)
    V["t_err_max"] = c(max(e["t_err"] for e in r1["errores"]), 1)
    secciones.append(render(P1_22))

    # ---------------- 2.3
    tabla("tab_orden", "Orden observado de convergencia.",
          ["h → 2h", "Err(h) [°C]", "Err(2h) [°C]", "p"],
          [["%s → %s" % (c(x["h"], 3), c(x["h2"], 3)), c(x["err_h"], 6), c(x["err_2h"], 6), c(x["p"], 4)] for x in o])
    figura("fig_orden", "fig_1_3_orden.png", "Err(h) en escala logarítmica: la pendiente medida coincide con la teórica (1).")
    V["p_medio"] = c(r1["p_medio"], 3)
    V["p_min"] = c(min(x["p"] for x in o), 3)
    V["p_max"] = c(max(x["p"] for x in o), 3)
    e01 = r1["errores_ref_h01"]
    V["e01_02"], V["e01_04"], V["e01_08"] = c(e01["0.2"], 6), c(e01["0.4"], 6), c(e01["0.8"], 6)
    V["p01_a"] = c(math.log(e01["0.4"] / e01["0.2"]) / math.log(2), 2)
    V["p01_b"] = c(math.log(e01["0.8"] / e01["0.4"]) / math.log(2), 2)
    V["err_01"], V["err_02"] = c(err[0.1]["err"], 4), c(err[0.2]["err"], 4)
    secciones.append(render(P1_23))

    # ---------------- 2.4
    filas_est = []
    for e in r1["estabilidad"]:
        det = e["detalle"].replace(".", ",")
        if e["estable"]:
            filas_est.append([c(e["h"], 2), "estable", c(e["err"], 4), det])
        else:
            filas_est.append([c(e["h"], 2), "**diverge**", "—", det])
    tabla("tab_estab", "Comportamiento de Euler al aumentar h.",
          ["h [meses]", "Estado", "Err(h) [°C]", "Observación"], filas_est)
    figura("fig_estab", "fig_1_4_estabilidad.png",
           "Izquierda: soluciones con pasos estables. Centro: soluciones inestables (eje recortado). Derecha: error en el borde de estabilidad.")
    V["t_div_126"] = c(est[1.26]["t_div"], 1)
    V["umbral"] = c(r1["umbral"], 4)
    tabla("tab_period", "Temperatura de Euler en t = 12, 24, 36, 48 y 60 meses.",
          ["h [meses]", "T(12)", "T(24)", "T(36)", "T(48)", "T(60)"],
          [[c(p["h"], 4)] + [c(v, 3) for v in p["T_12k"]] for p in r1["periodicidad"]])
    p08 = [p for p in r1["periodicidad"] if abs(p["h"] - 0.8) < 1e-9][0]["T_12k"]
    V["T_08_estab"] = "; ".join(c(v, 3) for v in p08[1:])
    tabla("tab_autoval", "Autovalores del jacobiano sobre la solución.",
          ["t [meses]", "Autovalores"], [[c(a["t"], 0), eig_texto(a["eig"])] for a in r1["autovalores"]])
    V["hcrit_60"] = c(r1["hcrit_60"], 4)
    V["hcrit_60_3"] = c(r1["hcrit_60"], 3)
    V["hcrit_min"] = c(r1["hcrit_min_total"], 3)
    V["hcrit_min_t"] = c(r1["hcrit_min_total_t"], 1)
    V["hcrit_ult_min"] = c(r1["hcrit_ult_min"], 2)
    V["hcrit_ult_max"] = c(r1["hcrit_ult_max"], 2)
    V["floquet_exacto"] = c(r1["floquet_exacto"], 4)
    fl = {f["N"]: f for f in r1["floquet"]}
    tabla("tab_floquet", "Radio espectral de M (Euler linealizado sobre un ciclo anual).",
          ["h [meses]", "ρ(M)"], [[c(fl[N]["h"], 4), c(fl[N]["rho"], 4)] for N in (16, 15, 14, 13, 12, 11, 10)])
    figura("fig_jac", "fig_1_4_jacobiano.png",
           "Izquierda: límite local h_crit(t) en el último ciclo. Derecha: radio espectral de Floquet frente a h.")
    V["fl_estable"] = c(r1["floquet_h_estable_max"], 3)
    V["fl_inestable"] = c(r1["floquet_h_inestable_min"], 3)
    V["p_medio_2"] = c(r1["p_medio"], 2)
    V["h2"] = c(h2, 3)
    V["razon_h"] = c(r1["umbral"] / h2, 0)
    secciones.append(render(P1_24))

    # ---------------- 3.1
    tabla("tab_verif", "Indicadores nominales según el paso utilizado.",
          ["h [meses]", "Err máx. en T [°C]", "S_P", "Δt_TP [meses]"],
          [[c(v["h"], 4), cs(v["err"], 2), c(v["SP"], 5), c(v["dtTP"], 4)] for v in r2["verificacion_h"]] +
          [["RK4 (referencia)", "—", c(r2["rk4_ref"]["SP"], 5), c(r2["rk4_ref"]["dtTP"], 4)]])
    V["err_h2"] = cs(vh[h2]["err"], 1)
    V["dif_SP"] = cs(abs(vh[h2]["SP"] - r2["rk4_ref"]["SP"]), 1)
    V["dif_dt"] = c(abs(vh[h2]["dtTP"] - r2["rk4_ref"]["dtTP"]), 4)
    tabla_sin_titulo("tab_nominal", ["Magnitud", "Valor"], [
        ["$P_{max}$", c(nom["Pmax"], 4)], ["$P_{min}$", c(nom["Pmin"], 4)], ["$P_{med}$", c(nom["Pmed"], 4)],
        ["$S_P$", "**%s**" % c(nom["SP"], 4)],
        ["$t_{T,max}$", "%s meses" % c(nom["tT"], 2)], ["$t_{P,max}$", "%s meses" % c(nom["tP"], 2)],
        ["$\\Delta t_{TP}$", "**%s meses**" % c(nom["dtTP"], 2)]])
    figura("fig_60m", "fig_3_1_60meses.png", "T(t), H(t) y P(t) para los 60 meses con parámetros nominales (sombreado: último ciclo).")
    figura("fig_ult", "fig_3_1_ultimo_ciclo.png", "Último ciclo anual usado para calcular S_P y Δt_TP, con los instantes de máximo de T y de P.")
    V["desvio_nom"] = cs(nom["desvio_rel"], 1)
    V["dt_nom_2"] = c(nom["dtTP"], 2)
    secciones.append(render(P2_31))

    # ---------------- 3.2
    filas32 = []
    for k in ("A", "beta", "gamma", "sigma"):
        b = bar[k]
        for lab, idx in (("mín.", 0), ("máx.", -1)):
            filas32.append([nombres[k], "%s (%s)" % (c(b["q"][idx], 3), lab), c(b["SP"][idx], 4),
                            c(b["dtTP"][idx], 3), c(b["Pmed"][idx], 3), c(b["amplitud"][idx], 3)])
    tabla("tab_sens", "Efecto de cada parámetro en los extremos de su intervalo.",
          ["Parámetro", "Valor", "S_P", "Δt_TP [meses]", "P_med", "P_max − P_min"], filas32)
    figura("fig_sens", "fig_3_2_sensibilidad.png",
           "S_P(q) (arriba) y Δt_TP(q) (abajo) para cada parámetro; la línea punteada marca el valor nominal.")
    V["amp_sig_0"], V["amp_sig_1"] = c(bar["sigma"]["amplitud"][0], 2), c(bar["sigma"]["amplitud"][-1], 2)
    V["pmed_sig_0"], V["pmed_sig_1"] = c(bar["sigma"]["Pmed"][0], 2), c(bar["sigma"]["Pmed"][-1], 2)
    secciones.append(render(P2_32))

    # ---------------- 3.3
    V["n_total"] = ("{:,}".format(g["n_total"])).replace(",", ".")
    V["h_grilla"] = c(g["h_grilla"], 3)
    if g["n_validas"] == g["n_total"]:
        V["validas_txt"] = "Las %s combinaciones cumplieron ambos requisitos." % V["n_total"]
    else:
        V["validas_txt"] = "%d de las %s combinaciones cumplieron ambos requisitos." % (g["n_validas"], V["n_total"])
    tabla("tab_top", "Las cinco combinaciones de mayor severidad.",
          ["S_P", "A", "β", "γ", "σ", "Δt_TP [meses]"],
          [[c(t["SP"], 4), c(t["A"], 3), c(t["beta"], 3), c(t["gamma"], 3), c(t["sigma"], 3), c(t["dtTP"], 3)]
           for t in g["top"][:5]])
    V["ps_A"], V["ps_beta"] = c(ps["A"], 2), c(ps["beta"], 2)
    V["ps_gamma"], V["ps_sigma"] = c(ps["gamma"], 2), c(ps["sigma"], 2)
    V["sp_sev"] = c(sev["SP"], 4)
    figura("fig_mapa", "fig_3_3_mapa_SP.png", "S_P en el plano (γ, A) con β y σ en su valor óptimo.")
    secciones.append(render(P2_33))

    # ---------------- 3.4
    tabla("tab_comp", "Caso nominal frente al escenario de mayor severidad (último ciclo anual).",
          ["", "Nominal", "Mayor severidad"], [
              ["$A$", c(pn["A"], 2), c(ps["A"], 2)], ["$\\beta$", c(pn["beta"], 2), c(ps["beta"], 2)],
              ["$\\gamma$", c(pn["gamma"], 2), c(ps["gamma"], 2)], ["$\\sigma$", c(pn["sigma"], 2), c(ps["sigma"], 2)],
              ["**$S_P$**", "**%s**" % c(nom["SP"], 4), "**%s**" % c(sev["SP"], 4)],
              ["**$\\Delta t_{TP}$ [meses]**", "**%s**" % c(nom["dtTP"], 3), "**%s**" % c(sev["dtTP"], 3)],
              ["$P_{max}$", c(nom["Pmax"], 3), c(sev["Pmax"], 3)], ["$P_{min}$", c(nom["Pmin"], 3), c(sev["Pmin"], 3)],
              ["$P_{med}$", c(nom["Pmed"], 3), c(sev["Pmed"], 3)]])
    figura("fig_comp", "fig_3_4_comparacion.png", "Caso nominal y escenario de mayor severidad: 60 meses (arriba) y último ciclo (abajo).")
    figura("fig_sev", "fig_3_4_ultimo_ciclo_severo.png", "Último ciclo del escenario de mayor severidad.")
    V["razon_sp"] = c(sev["SP"] / nom["SP"], 1)
    V["sp_nom_3"], V["sp_sev_3"] = c(nom["SP"], 3), c(sev["SP"], 3)
    V["pmax_nom"], V["pmax_sev"] = c(nom["Pmax"], 2), c(sev["Pmax"], 2)
    V["pmin_nom"], V["pmin_sev"] = c(nom["Pmin"], 2), c(sev["Pmin"], 2)
    V["dt_sev_2"] = c(sev["dtTP"], 2)
    V["sp_todo_max"] = c(tmax["SP"], 4)
    fl_ = {f["parametro"]: f for f in r2["flips"]}
    tabla("tab_flips", "S_P al llevar un parámetro al extremo opuesto del escenario severo.",
          ["Parámetro modificado", "Nuevo valor", "S_P", "Variación"],
          [[nombres[k], c(fl_[k]["valor"], 2), c(fl_[k]["SP"], 4), ("%+.4f" % (fl_[k]["SP"] - sev["SP"])).replace(".", ",")]
           for k in ("A", "beta", "gamma", "sigma")])
    secciones.append(render(P2_34))

    # ---------------- 4 y 5
    V["fuente_euler"] = inspect.getsource(tp.euler).rstrip()
    secciones.append(render(CODIGO))

    ruta = os.path.join(RAIZ, "INFORME_autogenerado.md")
    with open(ruta, "w", encoding="utf-8") as f:
        f.write("\n".join(secciones))
    print("Informe escrito en %s" % ruta)


if __name__ == "__main__":
    main()
