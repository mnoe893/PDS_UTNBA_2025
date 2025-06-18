import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema


def grafico_plantilla(H, w, wb, wp, a_max, a_min, etiquetas):
    # Calculo magnitud en escala logaritmica
    mag_db = 20 * np.log10(np.abs(H)) 

    max_indices = argrelextrema(mag_db, np.greater)[0]
    min_indices = argrelextrema(mag_db, np.less)[0]

    # Grafico
    plt.figure()
    plt.semilogx(w, mag_db, linewidth=2)

    # Marcar los extremos
    plt.plot(w[max_indices], mag_db[max_indices], 'ro', label=f'Máximos locales: (~{int(w[max_indices])} rad/s;{mag_db[max_indices]} dB)')
    plt.plot(w[min_indices], mag_db[min_indices], 'bo', label=f'Mínimos locales: (~{int(w[min_indices])} rad/s;{mag_db[min_indices]} dB)')

    # Lineas limite
    plt.axvline(wp, color='green', linestyle='--')
    if(etiquetas):
        plt.text(wp*(0.9), -50, f'ωp = {wp}', color='green',rotation=270, va='center', ha='right', fontsize=10, backgroundcolor='white')

    plt.axhline(-a_max, color='green', linestyle='--')
    if(etiquetas):
        plt.text(100, -5, f'α_max = {a_max}', color='green', va='center', ha='right', fontsize=10, backgroundcolor='white')

    plt.axvline(wr, color='red', linestyle='--')
    if(etiquetas):
        plt.text(wr*(0.9), -50, f'ωr = {wr}', color='red',rotation=270, va='center', ha='right', fontsize=10, backgroundcolor='white')

    plt.axhline(-a_min, color='red', linestyle='--')
    if(etiquetas):    
        plt.text(9000, -27, f'α_min = {a_min}', color='red', va='center', ha='right', fontsize=10, backgroundcolor='white')

    # Areas pintadas
    ymin, ymax = plt.ylim()
    xmin = w[0]
    plt.fill_betweenx([-a_max, ymax], xmin, wp, color='green', alpha=0.2)
    plt.ylim(ymin, ymax)

    xmax = w[-1]
    plt.fill_betweenx([ymin, -a_min], wr, xmax, color='red', alpha=0.2)
    plt.ylim(ymin, ymax)

    plt.title('Plantilla de atenuación')
    plt.xlabel('Frecuencia ω [rad/s]')
    plt.ylabel('Atenuación α [dB]')
    plt.grid(True, which='both', linestyle='--')
    if(etiquetas):  
        plt.legend()
    plt.tight_layout()
    plt.show()

# Parámetros del filtro
ftran = 0.1
fstop = min(1 / fs_os + ftran / 2, 1 / fs_os * 5 / 4)
fpass = max(fstop - ftran / 2, fstop * 3 / 4)
ripple = 0.5         # dB
attenuation = 40     # dB

# Ganancias en dB y lineal
gains_db = [0, -ripple/2, -attenuation/2, -attenuation/2]
gains = 10 ** (np.array(gains_db) / 20)


# Filtro Butterworth
orden_filtro = 4
fc = (fs / 2) * 0.95                   # Hz
f_nyquist = fs_os / 2                 # Hz
Wn = fc / f_nyquist                   # Normalizado para butter()
b, a = butter(orden_filtro, Wn, btype='low')

# Frecuencias para análisis
w = np.logspace(0, 5, 1000)           # rad/s
w_Hz = w / (2 * np.pi)                # Hz para freqs()
_, H = freqs(b, a, worN=w)

# Graficar
grafico_plantilla(H=H, w=w, wp=wp, wr=wr, a_max=a_max, a_min=a_min, etiquetas=0)



