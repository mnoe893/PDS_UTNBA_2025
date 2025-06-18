import numpy as np
import matplotlib.pyplot as plt
from numpy.fft import fft
from scipy.signal import windows

fs = 1000       # Frecuencia de muestreo [Hz]
N = fs          # cantidad de muestras
fo = fs/4       # [Hz]

# Declaro arr de variables aleatorias
fr = np.random.uniform(-2, 2, N) #Incertidumbre de frecuencia
sigma = 1
na = np.random.normal(0, sigma, N) #Incertidumbre de amplitud

a0 = 2 #Amplitud señal
omega0 = np.pi/2 #Frecuencia señal
omega1 = omega0 + fr * ((2*np.pi) / N) #Frecuencia con incertidumbre

k = np.arange(N) #Array de enteros de 0 a N-1
x = a0 * np.sin(omega1*k) + na

Wrectangular = windows.boxcar(N)                # Arr de 1 constantes, sin ventana
Wflattop = windows.flattop(N)                   # Ventana flattop (sombrero de bruja)
Wblackmanharris = windows.blackmanharris(N)     # Ventana blackmanharris (campana)
Wcosine = windows.cosine(N)                   # Venana cosine (medio coseno, loma)
arr_W = [Wrectangular, Wflattop, Wblackmanharris, Wcosine]
nombres_ventanas = ['Rectangular', 'Flattop', 'Blackman-Harris', 'Cosine']

def estimar_amplitud(x, w):
    xw = x * w # Ventaneo la señal
    Xw = np.abs(fft(xw))
    omega0_ind = int(N // 4) # Índice de omega0 = pi/2 → frecuencia fo = fs/4 → índice = N/4
    return Xw[omega0_ind] # Retorno el valor modulo en pi/2 (pico de señal seno)

def estimar_frecuencia(x, w):
    xw = x * w
    Xw = np.abs(fft(xw))
    f_ind = np.argmax(Xw) #Retorna argumento con el maximo pico, que deberia corresponder a la frecuencia fundamental
    omega = 2 * np.pi * f_ind / N # Paso frecuencia a radianes/s
    return omega


cant_realizaciones = 200
SNR_dB = [3, 10]
SNR_lin = [10**(SNR_dB[0]/10), 10**(SNR_dB[1]/10)]

for idx_snr, snr in enumerate(SNR_lin):
    print(f"\nSNR = {SNR_dB[idx_snr]} dB")

    # Ajuste de ruido para que la potencia de la señal senoidal sea 1 W
    P_signal = 1
    P_noise = P_signal / snr
    sigma = np.sqrt(P_noise)  # sigma del ruido blanco gaussiano

    resultados_amp = []
    resultados_freq = []

    for w in arr_W: # Recorro los 4 tipos de ventana (rectangular, flattop, blackmanharris, cosine)
        estimaciones_amp = []  # Arr de 200 valores de amplitud estimada
        estimaciones_freq = [] # Arr de 200 valores de frecuencia estimada

        for _ in range(cant_realizaciones):
            # Generación de ruido y señal, igual que ya hacés
            fr = np.random.uniform(-2, 2, N)
            omega1 = omega0 + fr * ((2*np.pi) / N)
            na = np.random.normal(0, sigma, N)
            x = a0 * np.sin(omega1*k) + na

            for idx_w, w in enumerate(arr_W):
                a_hat = estimar_amplitud(x, w)
                f_hat = estimar_frecuencia(x, w)

                estimaciones_amp.append(a_hat)
                estimaciones_freq.append(f_hat)

        # Convertimos a arrays para facilitar cálculo
        estimaciones_amp = np.array(estimaciones_amp)
        estimaciones_freq = np.array(estimaciones_freq)

        # Sesgo y varianza de amplitud
        mu_amp = np.mean(estimaciones_amp)
        sesgo_amp = mu_amp - a0
        varianza_amp = np.mean((estimaciones_amp - mu_amp)**2)

        # Sesgo y varianza de frecuencia
        mu_freq = np.mean(estimaciones_freq)
        sesgo_freq = mu_freq - omega0
        varianza_freq = np.mean((estimaciones_freq - mu_freq)**2)

        print(f"[SNR {snr} - ventana {idx_w}]")
        print(f"Sesgo amplitud: {sesgo_amp:.4f}, Varianza amplitud: {varianza_amp:.4f}")
        print(f"Sesgo frecuencia: {sesgo_freq:.4f}, Varianza frecuencia: {varianza_freq:.4f}")
'''
import numpy as np
import matplotlib.pyplot as plt
from numpy.fft import fft

fs = 1000       # Frecuencia de muestreo [Hz]
N = fs          # cantidad de muestras
fo = fs/4       # [Hz]

Vmax = np.sqrt(2) #[Volts]
df = fs/N # resolución espectral

bin = fs/N
def sen(t):
    rand_var = np.random.normal(1, 0.2)
    s = Vmax * np.sin(t*2*np.pi*fo*rand_var)

    return s

tt = np.arange(N) / fs
s_puro = Vmax * np.sin(tt*2*np.pi*fo)
s_var = Vmax * np.sin(tt*2*np.pi*(fo-(bin/2)))
sen_arr = 0
for i in range(20):
    sen_arr += Vmax * np.sin(tt*2*np.pi*(fo-(i*bin/2)))

s_var = np.concatenate((np.zeros(N*9), s_var))
ff = np.arange(0,fs,df/10) # 10 porque aumente por 10 el N
#ff = np.arange(0,fs,df)
fft_s = fft(s_puro) / N
fft_svar = fft(s_var) / N
fft_arr = fft(sen_arr) / N

plt.figure(figsize=(18, 4))
plt.plot(tt, s_puro, label='Seno fo')
plt.plot(tt, s_var, label='Seno fo-bin/2')
plt.legend(); plt.xlim(0, 0.2); plt.grid();plt.title("Señal en tiempo")
plt.show()

plt.figure(figsize=(18, 4))
#plt.plot(ff, 20*np.log10(np.abs(fft_s)*np.sqrt(2)), ':x', label='Seno en fo',)
plt.plot(ff, 20*np.log10(np.abs(fft_svar)*np.sqrt(2)), ':o', label='Seno en fo-bin/2')
#plt.plot(ff, 20*np.log10(np.abs(fft_arr)*np.sqrt(2)), label='Sen arr')
plt.legend(); plt.xlim(fo*0.97,fo*1.03);plt.grid();plt.title("Potencia señal pura y corrida bin/2 = {:3.0f}".format(bin))
plt.show()
'''
