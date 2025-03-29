import matplotlib.pylab as plt
import numpy as np
from scipy import signal #Para funcion triangular, sierra y cuadrada
import matplotlib.pyplot as plt #Para la visualizacion de tt y xx.

'''mi_funcion_sen'''
# Produce una funcion seno de los parametros solicitados
# Retorna 2 matrices Nx1 correspondientes a los valores de la variable y su imagen
#------------------------------------------------------
'''
Vmax la amplitud máxima de la senoidal (volts)
dc su valor medio (volts)
ff la frecuencia (Hz)
ph la fase (radianes)
nn la cantidad de muestras digitalizada por el ADC (# muestras)
fs la frecuencia de muestreo del ADC. '''

#Si delta_f es normalizado, fs = N.
def mi_funcion_sen(vmax=1,dc=0,ff=1,ph=0,nn=1000,fs=1000):
    t = np.linspace(0, nn/fs, nn, endpoint=False) #False evita que se incluya valor limite como ultimo valor (Recomendado por chatGPT)
    #Utilizo linspace en vez de arange por recomendacion de numpy (Referirse a la warning en https://numpy.org/doc/2.1/reference/generated/numpy.arange.html)
    sin = vmax * np.sin(t*2*np.pi*ff + ph) + dc
    return t, sin

def mi_funcion_triangular(vmax=1,dc=0,ff=1,ph=0,nn=1000,fs=1000):
    t = np.linspace(0, nn/fs, nn, endpoint=False)
    sierra = vmax * signal.sawtooth((t*2*np.pi*ff + ph), 0.5) + dc
    #Segundo parametro = 0,5 para obtener una funcion triangular
    return t, sierra

def mi_funcion_sierra(vmax=1,dc=0,ff=1,ph=0,nn=1000,fs=1000):
    t = np.linspace(0, nn/fs, nn, endpoint=False)
    sierra = vmax * signal.sawtooth(t*2*np.pi*ff + ph) + dc
    return t, sierra

def mi_funcion_cuadrada(vmax=1,dc=0,ff=1,ph=0,nn=1000,fs=1000):
    t = np.linspace(0, nn/fs, nn, endpoint=False)
    sierra = vmax * signal.square(t*2*np.pi*ff + ph) + dc
    return t, sierra

#------------
#Testeo
#------------
def graficar_funcion(abscisas,ordenadas,label_y="",label_x="",titulo=""):
    plt.plot(abscisas, ordenadas)
    plt.xlabel(label_x)
    plt.ylabel(label_y)
    plt.title(titulo)
    plt.grid()
    plt.axis('tight')
    #plt.show() #Comentar si se usan subplots
    return

tt, xx = mi_funcion_sen(vmax=2.5,ff=2)
plt.subplot(2, 2, 1)
graficar_funcion(tt,xx,"Seno[V]","Tiempo[s]")

tt, xx = mi_funcion_cuadrada(ff=1,nn=5000,dc=1,vmax=1)
plt.subplot(2, 2, 2)
graficar_funcion(tt,xx,"Cuadrada[V]","Tiempo[s]")

tt, xx = mi_funcion_triangular(vmax=2.5,ff=5)
plt.subplot(2, 2, 3)
graficar_funcion(tt,xx,"Triangular[V]","Tiempo[s]")

tt, xx = mi_funcion_sierra(vmax=1,ff=3,nn=2000,ph=np.pi)
plt.subplot(2, 2, 4)
graficar_funcion(tt,xx,"Sierra[V]","Tiempo[s]")

plt.show()

# Testeo valores default
'''
tt, xx = mi_funcion_sen()
plt.subplot(2, 2, 1)
graficar_funcion(tt,xx,"Seno[V]","Tiempo[s]")

tt, xx = mi_funcion_cuadrada()
plt.subplot(2, 2, 2)
graficar_funcion(tt,xx,"Cuadrada[V]","Tiempo[s]")

tt, xx = mi_funcion_triangular()
plt.subplot(2, 2, 3)
graficar_funcion(tt,xx,"Triangular[V]","Tiempo[s]")

tt, xx = mi_funcion_sierra()
plt.subplot(2, 2, 4)
graficar_funcion(tt,xx,"Sierra[V]","Tiempo[s]")

plt.show()
'''