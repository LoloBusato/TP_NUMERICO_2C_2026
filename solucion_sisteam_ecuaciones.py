
import matplotlib.pyplot as plt
import pandas as pd
import math

#                  PARAMETROS DEL SISTEMA

T_max = 4.0
Roce = 0.30
alfa = 0.05
beta = 0.20
gamma = 0.10
delta = 0.40
sigma = 0.70
A = 0.45
omega = 2 * math.pi / 12

h_inicial = 0.9
particiones = int(60 / h_inicial)

#                  CONDICIONES INICIALES

temperatura_n = [1.0]
humedad_n = [0]
precipitacion_n = [0]
tiempos = [0]

#                  FUNCIONES DEL SISTEMA

def temperatura(R_t, T_max, T, alfa, humedad):
    return R_t*T*(1-T/T_max) - alfa*humedad

def humedad(beta, T, gamma, humedad, precipitacion):
    return beta*T**2 - gamma*humedad*precipitacion

def precipitacion(delta, humedad, sigma, precipitacion):
    return delta*humedad**2 - sigma*precipitacion

#                  METODO DE EULER
   

for i in range(1, particiones + 1):

    roce_t = Roce + A * math.sin(omega * tiempos[-1] - math.pi)

    h_i = h_inicial

    T_nueva = temperatura_n[-1] + h_i * temperatura(
        roce_t, T_max, temperatura_n[-1], alfa, humedad_n[-1])

    H_nueva = humedad_n[-1] + h_i * humedad(
        beta, temperatura_n[-1], gamma, humedad_n[-1], precipitacion_n[-1])

    P_nueva = precipitacion_n[-1] + h_i * precipitacion(
        delta, humedad_n[-1], sigma, precipitacion_n[-1])

    temperatura_n.append(T_nueva)
    humedad_n.append(H_nueva)
    precipitacion_n.append(P_nueva)

    tiempos.append(tiempos[-1] + h_i)

#                  GUARDAR DATOS EN EXCEL

datos = {
    "Tiempo": tiempos,
    "Temperatura": temperatura_n,
    "Humedad": humedad_n,
    "Precipitacion": precipitacion_n
}

df = pd.DataFrame(datos)

df.to_excel(f"solucion_h_{h_inicial}.xlsx", index=False)

#                  GRAFICAR RESULTADOS

plt.figure(figsize=(10, 6))

plt.plot(tiempos, temperatura_n, 'o-', markersize=3, label="Temperatura")
plt.plot(tiempos, humedad_n, 'o-', markersize=3, label="Humedad")
plt.plot(tiempos, precipitacion_n, 'o-', markersize=3, label="Precipitación")

plt.xlabel("Tiempo (meses)")
plt.ylabel("Valor")
plt.title(f"Evolución del sistema (h = {h_inicial})")
plt.legend()
plt.grid()

plt.savefig(f"grafico_{h_inicial}.png", dpi=300, bbox_inches="tight")

plt.show()
