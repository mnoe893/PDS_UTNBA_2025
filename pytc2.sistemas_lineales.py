from scipy.signal import butter, lfilter, freqz
import scipy.signal as sp
import numpy as np
import scipy.signal as sig
import matplotlib.pyplot as plt
from numbers import Integral, Real

# Gráficos interactivos
#%matplotlib ipympl
# Gráficos estáticos
#%matplotlib inline

#from pytc2.sistemas_lineales import plot_plantilla
def plot_plantilla(filter_type='', fpass=0.25, ripple=0.5, fstop=0.6, attenuation=40, fs=2):
    """
    Plotea una plantilla de diseño de filtro digital.

    Parameters
    -----------
    filter_type : str, opcional
        Tipo de filtro ('lowpass', 'highpass', 'bandpass', 'bandstop'). Por defecto es 'lowpass'.
    fpass : float o tupla, opcional
        Frecuencia de paso o tupla de frecuencias de paso para los filtros 'bandpass' o 'bandstop'.
    ripple : float, opcional
        Máxima ondulación en la banda de paso (en dB). Por defecto es 0.5 dB.
    fstop : float o tupla, opcional
        Frecuencia de detención o tupla de frecuencias de detención para los filtros 'bandpass' o 'bandstop'.
    attenuation : float, opcional
        Atenuación mínima en la banda de detención (en dB). Por defecto es 40 dB.
    fs : float, opcional
        Frecuencia de muestreo. Por defecto es 2.
        
    Returns
    --------
    None


    Raises
    ------
    ValueError
        Si los argumentos no son del tipo o valor correcto.
    

    See Also
    -----------
    :func:`analyze_sys`

    
    Example
    --------
    >>> # Analiza un sistema con w0 = 1 rad/s y Q = sqrt(2)/2
    >>> import numpy as np
    >>> from scipy import signal as sig
    >>> import matplotlib.pyplot as plt
    >>> from pytc2.sistemas_lineales import bodePlot, plot_plantilla
    >>> Q = np.sqrt(2)/2
    >>> w0 = 1
    >>> num = np.array([w0**2])
    >>> den = np.array([1., w0 / Q, w0**2])
    >>> H1 = sig.TransferFunction(num, den)
    >>> fig_id, axes_hdl = bodePlot(H1, fig_id=1, axes_hdl='none', filter_description='Filtro pasa bajos', worN=1000, digital=False, xaxis='omega', fs=2*np.pi)
    >>> plt.sca(axes_hdl[0])
    >>> plot_plantilla(filter_type='lowpass', fpass=1.0, ripple=3, fstop=3.0, attenuation=20, fs=2)

    """

    if not isinstance(fpass, (tuple, np.ndarray, Integral, Real)):
        raise ValueError("fpass debe ser un float o una tupla de frecuencias de paso.")
    
    if not isinstance(fstop, (tuple, np.ndarray, Integral, Real)):
        raise ValueError("fstop debe ser un float o una tupla de frecuencias de detención.")
    
    if not isinstance(attenuation, (tuple, np.ndarray, Integral, Real)):
        try:
            attenuation = np.float64(attenuation)
        except ValueError:
            raise ValueError("attenuation debe ser un valor numérico o convertible a float.")    
    
    if not isinstance(ripple, (tuple, np.ndarray, Integral, Real)):
        try:
            ripple = np.float64(ripple)
        except ValueError:
            raise ValueError("attenuation debe ser un valor numérico o convertible a float.")    
    
    if not isinstance(fs, (Integral, Real)):
        try:
            fs = np.float64(fs)
        except ValueError:
            raise ValueError("fs debe ser un valor numérico o convertible a float.")    

    # Obtener los límites actuales de los ejes
    xmin, xmax, ymin, ymax = plt.axis()

    # Banda de paso digital
    plt.fill([xmin, xmin, fs / 2, fs / 2], [ymin, ymax, ymax, ymin], 'lightgreen', alpha=0.2, lw=1, label='Banda de paso digital')

    # analizar los valores de la plantilla para ver qué tipo de plantilla es
    tipos_permitidos = ['lowpass', 'highpass', 'bandpass', 'bandstop']
    if filter_type not in tipos_permitidos:
            
        if isinstance(fpass, (tuple, np.ndarray)) and isinstance(fstop, (tuple, np.ndarray)):
            if fstop[0] < fpass[0]:
                filter_type = 'bandpass'
            else:                
                filter_type = 'bandstop'
        else:
            if fstop < fpass:
                filter_type = 'highpass'
            else:                
                filter_type = 'lowpass'

    if filter_type == 'lowpass':
        # Definir regiones de banda de detención y banda de paso para el filtro pasa bajos
        fstop_start = fstop
        fstop_end = xmax
        fpass_start = xmin
        fpass_end = fpass
    elif filter_type == 'highpass':
        # Definir regiones de banda de detención y banda de paso para el filtro pasa altos
        fstop_start = xmin
        fstop_end = fstop
        fpass_start = fpass
        fpass_end = xmax
    elif filter_type == 'bandpass':
        if len(fpass) == 2 and len(fstop) == 2:
            fstop_start = xmin
            fstop_end = fstop[0]
            fpass_start = fpass[0]
            fpass_end = fpass[1]
            fstop2_start = fstop[1]
            fstop2_end = xmax
        else:
            raise ValueError("En modo bandpass, fpass y fstop deben ser tuplas con 2 valores.")
    elif filter_type == 'bandstop':
        if len(fpass) == 2 and len(fstop) == 2:
            fpass_start = xmin
            fpass_end = fpass[0]
            fstop_start = fstop[0]
            fstop_end = fstop[1]
            fpass2_start = fpass[1]
            fpass2_end = xmax
        else:
            raise ValueError("En modo bandstop, fpass y fstop deben ser tuplas con 2 valores.")
    else:
        raise ValueError("filtro_type debe ser 'lowpass', 'highpass', 'bandpass', o 'bandstop'.")

    # Plotea regiones de banda de detención y banda de paso
    plt.fill([fstop_start, fstop_end, fstop_end, fstop_start], [-attenuation, -attenuation, ymax, ymax], 'lightgrey', alpha=0.4, hatch='x', lw=1, ls='--', ec='k')
    plt.fill([fpass_start, fpass_start, fpass_end, fpass_end], [ymin, -ripple, -ripple, ymin], 'lightgrey', alpha=0.4, hatch='x', lw=1, ls='--', ec='k', label='Plantilla')
    
    # Plotea región adicional de banda de detención para filtro pasa banda
    if filter_type == 'bandpass':
        plt.fill([fstop2_start, fstop2_end, fstop2_end, fstop2_start], [-attenuation, -attenuation, ymax, ymax], 'lightgrey', alpha=0.4, hatch='x', lw=1, ls='--', ec='k')
    
    # Plotea región adicional de banda de paso para filtro rechaza banda
    if filter_type == 'bandstop':
        plt.fill([fpass2_start, fpass2_start, fpass2_end, fpass2_end], [ymin, -ripple, -ripple, ymin], 'lightgrey', alpha=0.4, hatch='x', lw=1, ls='--', ec='k')
    
    # Establece los límites de los ejes
    plt.axis([xmin, xmax, np.max([ymin, -100]), np.max([ymax, 1])])
# Tipo de aproximación.
        
aprox_name = 'butter'
# aprox_name = 'cheby1'
# aprox_name = 'cheby2'
# aprox_name = 'ellip'

# Por qué no hay bessel ?
#aprox_name = 'bessel'

# Requerimientos de plantilla

filter_type = 'lowpass'
# filter_type = 'highpass'
# filter_type = 'bandpass'
# filter_type = 'bandstop'


# plantillas normalizadas a Nyquist y en dB

if filter_type == 'lowpass':

    # fpass = 1/2/np.pi # 
    fpass = 0.25 # 
    ripple = 0.5 # dB
    fstop = 0.6 # Hz
    attenuation = 40 # dB

elif filter_type == 'highpass':

    fpass = 0.6 
    ripple = 0.5 # dB
    fstop = 0.25
    attenuation = 40 # dB

elif filter_type == 'bandpass':

    fpass = np.array( [0.4, 0.6] ) 
    ripple = 0.5 # dB
    fstop = np.array( [0.25, 0.75] ) 
    attenuation = 40 # dB
    
else:

    # bandstop
    fpass = np.array( [0.25, 0.75] ) 
    ripple = 0.5 # dB
    fstop = np.array( [0.4, 0.6] ) 
    attenuation = 40 # dB


    # Cálculo del filtro

# frecuencia de muestreo normalizada (Nyquist = 1)
fs = 2

if aprox_name == 'butter':

    order, wcutof = sig.buttord( 2*np.pi*fpass*fs/2, 2*np.pi*fstop*fs/2, ripple, attenuation, analog=True)
    orderz, wcutofz = sig.buttord( fpass, fstop, ripple, attenuation, analog=False)

elif aprox_name == 'cheby1':

    order, wcutof = sig.cheb1ord( 2*np.pi*fpass*fs/2, 2*np.pi*fstop*fs/2, ripple, attenuation, analog=True)
    orderz, wcutofz = sig.cheb1ord( fpass, fstop, ripple, attenuation, analog=False)
    
elif aprox_name == 'cheby2':

    order, wcutof = sig.cheb2ord( 2*np.pi*fpass*fs/2, 2*np.pi*fstop*fs/2, ripple, attenuation, analog=True)
    orderz, wcutofz = sig.cheb2ord( fpass, fstop, ripple, attenuation, analog=False)
    
elif aprox_name == 'ellip':
   
    order, wcutof = sig.ellipord( 2*np.pi*fpass*fs/2, 2*np.pi*fstop*fs/2, ripple, attenuation, analog=True)
    orderz, wcutofz = sig.ellipord( fpass, fstop, ripple, attenuation, analog=False)


# Diseño del filtro analógico

num, den = sig.iirfilter(order, wcutof, rp=ripple, rs=attenuation, btype=filter_type, analog=True, ftype=aprox_name)

my_analog_filter = sig.TransferFunction(num,den)
my_analog_filter_desc = aprox_name + '_ord_' + str(order) + '_analog'

# Diseño del filtro digital

numz, denz = sig.iirfilter(orderz, wcutofz, rp=ripple, rs=attenuation, btype=filter_type, analog=False, ftype=aprox_name)

my_digital_filter = sig.TransferFunction(numz, denz, dt=1/fs)
my_digital_filter_desc = aprox_name + '_ord_' + str(orderz) + '_digital'



# Plantilla de diseño

plt.figure(1)
plt.cla()

npoints = 1000
w_nyq = 2*np.pi*fs/2

w, mag, _ = my_analog_filter.bode(npoints)
plt.plot(w/w_nyq, mag, label=my_analog_filter_desc)

w, mag, _ = my_digital_filter.bode(npoints)
plt.plot(w/w_nyq, mag, label=my_digital_filter_desc)

plt.title('Plantilla de diseño')
plt.xlabel('Frecuencia normalizada a Nyq [#]')
plt.ylabel('Amplitud [dB]')
plt.grid(which='both', axis='both')

plt.gca().set_xlim([0, 1])

plot_plantilla(filter_type = filter_type , fpass = fpass, ripple = ripple , fstop = fstop, attenuation = attenuation, fs = fs)

plt.legend()




    