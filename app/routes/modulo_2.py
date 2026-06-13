import numpy as np
import time
from app.routes.fr_daily_rad import solar_angles_day
from scipy.integrate import trapezoid, dblquad, simpson
import pandas as pd

# --- CÁLCULOS SOLARES Y DE IRRADIANCIA ---

def calculate_radiation_components(latitude, longitude, day_of_year, freq_min, S_gd_measured=None, clear_sky=True):
    """
    Calcula la evolución diurna de la radiación global y sus componentes
    PAR directa y difusa según Spitters et al. (1986).

    Parámetros
    ----------
    latitude, longitude : float
        Coordenadas geográficas del sitio (grados).
    day_of_year : int
        Día del año (1–365).
    freq_min : str
        Resolución temporal (p. ej. "10min").
    S_gd_measured : float, opcional
        Radiación global diaria medida (MJ/m²/día). Si no se proporciona,
        se asume cielo claro estándar.
    clear_sky : bool, opcional
        Si True, aplica la corrección de cielo claro para la fracción difusa
        en PAR.

    Retorna
    -------
    df : pandas.DataFrame
        Serie temporal con ángulos solares y flujos radiativos instantáneos:
        Sg, PAR, PAR directa y difusa.
    PAR_dia : float
        PAR diaria total (mol PAR/m²).
    PAR_directa_dia : float
        PAR directa diaria (mol PAR/m²).
    PAR_difusa_dia : float
        PAR difusa diaria (mol PAR/m²).
    """
    
    df = solar_angles_day(latitude, longitude, day_of_year, freq_min)
    df.rename(columns={'elevation': 'theta', 'azimuth': 'phi'}, inplace=True)
    theta_rad = np.deg2rad(df['theta']) # Convertir theta (elevación) a radianes para cálculos
    
    #Calcular S_o_d (radiación extraterrestre diaria)
    S_sc = 1370  # W/m²
    S_sc_MW = S_sc / 1e6  
    
    # Día del año en radianes para corrección orbital
    t_d = day_of_year
    orbital_correction = 1 + 0.033 * np.cos(2 * np.pi * t_d / 365)
    
    day_mask = df['theta'] > 0
    sin_theta_day = np.sin(theta_rad[day_mask])
    
    time_interval = pd.Timedelta(freq_min).total_seconds()  # segundos
    I_sin_theta = trapezoid(sin_theta_day, dx=time_interval)
    
    S_o_d = S_sc_MW * orbital_correction * I_sin_theta # Radiación extraterrestre diaria MJ/m^2
    
    if S_gd_measured is not None:
        S_gd = S_gd_measured
        K_t = S_gd / S_o_d #indice de claridad
    else:
        # Para cielo claro estándar: S_gd = 0.76 * S_o_d
        S_gd = 0.76 * S_o_d # Calcular S_gd (radiación global diaria) 
        K_t = 0.76

    #print(f"Indice de claridad: {K_t}")
    
    # Calcular fracción difusa diaria según Spitters
    if K_t <= 0.07:
        f_diff_daily = 1.0
    elif K_t <= 0.35:
        f_diff_daily = 1 - 2.3 * (K_t - 0.07)**2
    elif K_t <= 0.75:
        f_diff_daily = 1.33 - 1.46 * K_t
    else:
        f_diff_daily = 0.23
    
    if clear_sky:
        f_diff_PAR = 1.4 * f_diff_daily # Para cielo claro: fracción difusa en PAR = 1.4 * fracción difusa global
    else:
        f_diff_PAR = f_diff_daily  # O usar eq. 10 de Spitters para precisión
    
    f_diff_PAR = min(f_diff_PAR, 0.95) # Asegurar que f_diff_PAR no exceda 1
    f_direct_PAR = 1 - f_diff_PAR
    
    b_over_a = 0.4 # Coeficiente b/a para corrección diurna (Spitters usa ~0.4)
    sin_theta = np.sin(theta_rad)
    sin_theta_corrected = sin_theta * (1 + b_over_a * sin_theta)
    I_sin_theta_corrected = trapezoid(sin_theta_corrected[day_mask], dx=time_interval) 
    
    Sg_t = np.zeros(len(df))
    Sg_t[day_mask] = (sin_theta_corrected[day_mask] * 
                      S_gd / I_sin_theta_corrected)
    
    PAR_t = (0.5/0.217) * Sg_t  # Convertir a PAR (el 50% es PAR y 1 mol de PAR es aprox 0.217 MJ)
    PARdf_t = f_diff_PAR * PAR_t
    PARdr_t = f_direct_PAR * PAR_t

    df['Sg (MW/m^2)'] = Sg_t
    df['PAR (mol PAR/m^2/s)'] = PAR_t
    df['PARdf (mol PAR/m^2/s)'] = PARdf_t
    df['PARdr (mol PAR/m^2/s)'] = PARdr_t

    PAR_dia = (0.5/0.217)*S_gd
    return df, PAR_dia, PAR_dia*f_direct_PAR, PAR_dia*f_diff_PAR

# --- CÁLCULO DE CAMINO ÓPTICO EN SETO TRAPEZOIDAL ---

def longitud_rayo_en_trapecio(
    x0, y0, dx, dy,
    w, t, h,
    eps=1e-9
):
    """
    Longitud geométrica del rayo dentro de un trapecio centrado en x=0.
    dx, dy pueden ser arrays.
    """

    # Inicializar intervalo [s_min, s_max]
    s_min = np.zeros_like(dx)
    s_max = np.full_like(dx, np.inf)

    # ---------
    # y >= 0
    # ---------
    mask = np.abs(dy) > eps
    s1 = (0 - y0) / dy
    s2 = (h - y0) / dy
    s_low = np.minimum(s1, s2)
    s_high = np.maximum(s1, s2)

    s_min = np.where(mask, np.maximum(s_min, s_low), s_min)
    s_max = np.where(mask, np.minimum(s_max, s_high), s_max)

    # Pendiente de los lados
    m = (w - t) / (2 * h)

    # -----------------
    # Lado izquierdo
    # x >= -w/2 + m y
    # -----------------
    A = dx - m * dy
    B = (-w/2 + m * y0) - x0

    mask = np.abs(A) > eps
    s = B / A
    s_min = np.where(mask & (A > 0), np.maximum(s_min, s), s_min)
    s_max = np.where(mask & (A < 0), np.minimum(s_max, s), s_max)

    # -----------------
    # Lado derecho
    # x <= w/2 - m y
    # -----------------
    A = dx + m * dy
    B = (w/2 - m * y0) - x0

    mask = np.abs(A) > eps
    s = B / A
    s_min = np.where(mask & (A < 0), np.maximum(s_min, s), s_min)
    s_max = np.where(mask & (A > 0), np.minimum(s_max, s), s_max)

    # Longitud geométrica
    length = np.maximum(0.0, s_max - s_min)

    return length

def calcular_Pr_traps_vec(
    w, t, h, d,
    alpha_c, beta_c,
    x0, y0,
    vecinos=2,
    eps=1e-9
):
    """
    Calcula la distancia efectiva recorrida por un rayo al atravesar un conjunto
    periódico de setos de sección trapezoidal.
    
    Para cada dirección definida por (alpha_c, beta_c), se determina de forma
    analítica la longitud total del rayo contenida en el follaje, sumando las
    contribuciones de varios setos vecinos. Para generalizar fuera dl plano 
    transversal dividimos por cos(beta_c).
    
    La función es vectorizable y no utiliza geometría discreta.
    """


    alpha_c = np.asarray(alpha_c)
    beta_c  = np.asarray(beta_c)
    x0      = np.asarray(x0)
    y0      = np.asarray(y0)

    theta = np.pi - alpha_c
    dx = np.sin(theta)
    dy = -np.cos(theta)
    dy = np.where(np.abs(dy) < eps, -eps, dy)

    cos_beta = np.cos(beta_c)
    cos_beta = np.where(np.abs(cos_beta) < eps, eps, cos_beta)

    P = w + d

    out_shape = np.broadcast(x0, dx).shape
    total_length = np.zeros(out_shape)

    for k in range(-vecinos, vecinos + 1):
        x_shift = x0 - k * P

        Lk = longitud_rayo_en_trapecio(
            x_shift, y0,
            dx, dy,
            w, t, h,
            eps
        )

        total_length += Lk

    Pr = total_length / cos_beta
 

    return Pr

# --- CÁLCULOS DE TRNSMISIÓN DE LA LUZ (SETO DE SECCIÓN TRAPEZOIDAL) ---

def transformar_angulos(beta, phi_sun, eta):
    alpha = phi_sun-eta
    alpha_prime = alpha % (2*np.pi) # entre 0 y 2pi
    sin_beta_c = np.cos(alpha) * np.cos(beta)
    beta_c = np.arcsin(sin_beta_c) 
    
    cos_beta_c = np.cos(beta_c)
    if abs(cos_beta_c) < 1e-10:
        cos_beta_c = 1e-10  # Evitar división por cero
    cos_alpha_c = np.sin(beta)/cos_beta_c
    alpha_c = np.arccos(cos_alpha_c)
    if alpha_prime>np.pi:
        alpha_c = -alpha_c # criterio de signos: si eta=0,  alpha_c>0 en los cuadrantes NE y SE y alpha<0 en los cuadrnates NO y SO
    return alpha_c, beta_c

def calcular_PAR_directa_traps(param_arquitectura, param_modelo_luz, num_celdas_x, num_celdas_y, df_datos_doy, freq_min):
    """
    Calcula la PAR directa diaria acumulada por celda en una hilera
    de sección trapezoidal, integrando la radiación directa a lo
    largo del día según la geometría solar y la arquitectura del seto.

    Parámetros
    ----------
    param_arquitectura : dict
        Geometría del seto: w, t, h, d (m) y orientación eta (°).
    param_modelo_luz : dict
        Parámetros ópticos: O_av, LAI y sigma.
    num_celdas_x, num_celdas_y : int
        Número de celdas en las direcciones horizontal y vertical.
    df_datos_doy : pandas.DataFrame
        Datos diarios con columnas: theta (°), phi (°) y
        'PARdr (mol PAR/m^2/s)'.
    freq_min : str
        Resolución temporal de los datos (p. ej. "10min").

    Retorna
    -------
    PAR_directa_acum : ndarray
        Matriz (num_celdas_y × num_celdas_x) con la PAR directa diaria
        acumulada por celda (mol PAR/m²/día).
    """
    w = param_arquitectura["w"]
    t = param_arquitectura["t"]
    h = param_arquitectura["h"]
    d = param_arquitectura["d"]
    eta = np.radians(param_arquitectura["eta"])
    W = max(t, w)
    
    sigma = param_modelo_luz["O_av"]
    O_av = param_modelo_luz["O_av"]
    LAI = param_modelo_luz["LAI"]
    Ld = LAI / h

    x_edges = np.linspace(-W/2, W/2, num_celdas_x + 1)
    y_edges = np.linspace(h, 0, num_celdas_y + 1)
    x_centers = (x_edges[:-1] + x_edges[1:]) / 2
    y_centers = (y_edges[:-1] + y_edges[1:]) / 2
    x_grid, y_grid = np.meshgrid(x_centers, y_centers, indexing="xy")
     
    PAR_directa_acum = np.zeros((num_celdas_y, num_celdas_x))
    delta_t = int(pd.Timedelta(freq_min).total_seconds()) 
    for idx, row in df_datos_doy.iterrows():
        # Saltar si no hay radiación directa
        if row['PARdr (mol PAR/m^2/s)'] <= 0:
            continue
        
        beta = np.radians(row['theta'])  # Elevación
        phi_sun = np.radians(row['phi'])  # Azimut solar
        PAR_dr_inst = row['PARdr (mol PAR/m^2/s)'] # PAR directa incidente en este instante
        alpha_c, beta_c = transformar_angulos(beta, phi_sun, eta)
        
        K_dir_bl = O_av / np.sin(beta) if np.sin(beta) > 0 else 100  # Evitar división por cero
        K_dir = K_dir_bl * np.sqrt(1 - sigma)
        
        Pr = calcular_Pr_traps_vec(w, t, h, d, alpha_c, beta_c, x_grid, y_grid)
        Lt = Pr*Ld
        I_dir = PAR_dr_inst * np.exp(-K_dir * Lt) # / np.sin(beta) (Beer-Lambert). 
                                                  # No dividir por sin(beta) -> incidente sobre el plano horizontal
                                          
        PAR_directa_acum += I_dir*delta_t 
    return PAR_directa_acum

def calcular_PAR_difusa_traps(param_arquitectura, param_modelo_luz, num_celdas_x, num_celdas_y, PARdf):
    """
    Calcula la PAR difusa diaria distribuida por celda en una hilera
    de sección trapezoidal, integrando la radiación difusa del cielo
    y la componente reflejada por el suelo.

    Parámetros
    ----------
    param_arquitectura : dict
        Geometría del seto: w, t, h y d (m).
    param_modelo_luz : dict
        Parámetros ópticos: O_av, LAI, sigma, rho y albedo.
    num_celdas_x, num_celdas_y : int
        Número de celdas en las direcciones horizontal y vertical.
    PARdf : float
        PAR difusa diaria incidente (mol PAR/m²).

    Retorna
    -------
    PAR_difusa_celda : ndarray
        Matriz (num_celdas_y × num_celdas_x) con la PAR difusa diaria
        total por celda (mol PAR/m²/día).
    """
    w = param_arquitectura["w"]
    t = param_arquitectura["t"]
    h = param_arquitectura["h"]
    d = param_arquitectura["d"]
    W = max(t, w)
    sigma = param_modelo_luz["O_av"]
    O_av = param_modelo_luz["O_av"]
    LAI = param_modelo_luz["LAI"]
    rho = param_modelo_luz["rho"]
    albedo = param_modelo_luz["albedo"]
    Ld = LAI / h
     
    x_edges = np.linspace(-W/2, W/2, num_celdas_x + 1)
    y_edges = np.linspace(h, 0, num_celdas_y + 1)
    x_centers = (x_edges[:-1] + x_edges[1:]) / 2
    y_centers = (y_edges[:-1] + y_edges[1:]) / 2
    x_grid, y_grid = np.meshgrid(x_centers, y_centers, indexing="xy")
    
    num_divisiones = 20
    angulos = np.linspace(-np.pi/2, np.pi/2, num_divisiones)
    alpha_c_grid, beta_c_grid = np.meshgrid(angulos, angulos)

    PAR_difusa_celda_arriba = np.zeros((num_celdas_y, num_celdas_x))
    PAR_difusa_celda_abajo = np.zeros((num_celdas_y, num_celdas_x))

    def integrando_difusa(beta_c, alpha_c, x, y, arriba = True):
        if arriba:
            Pr = calcular_Pr_traps_vec(w, t, h, d, alpha_c, beta_c, x, y)
        else:
            Pr = calcular_Pr_traps_vec(t, w, h, w-t+d, alpha_c, beta_c, x, y)
        Lt = Pr*Ld
        return 1/np.pi*(1-rho)*np.exp(-O_av*Lt*np.sqrt(1-sigma))*np.cos(beta_c)**2*np.cos(alpha_c)
    
    def integrar_angulos(x, y):
        integrando_eval_arriba = integrando_difusa(beta_c_grid, alpha_c_grid, x, y, arriba = True)
        integrando_eval_abajo = integrando_difusa(beta_c_grid, alpha_c_grid, x, y, arriba = False)

        resultado_arriba = simpson(simpson(integrando_eval_arriba, angulos, axis = 0), angulos)
        resultado_abajo = simpson(simpson(integrando_eval_abajo, angulos, axis = 0), angulos)
        
        return resultado_arriba, resultado_abajo
    
    for j in range(num_celdas_y):
        for i in range(num_celdas_x):
            x = x_grid[j, i]
            y = y_grid[j, i]
            result_arriba, result_abajo = integrar_angulos(x, y)
            PAR_difusa_celda_arriba[j, i] = PARdf * result_arriba
            PAR_difusa_celda_abajo[j, i] = PARdf * albedo * result_abajo

    
    return PAR_difusa_celda_arriba + np.flipud(PAR_difusa_celda_abajo)

# --- CÁLCULOS DE TRNSMISIÓN DE LA LUZ (PRIMERA VERSIÓN, SECCIÓN RECTANGULAR, MENOS GENERAL) ---

def calcular_PAR_directa(param_arquitectura, param_modelo_luz, num_celdas_x, num_celdas_z, df_datos_doy, freq_min):
    """
    Calcula la PAR directa absorbida por celda en una hilera rectangular.
    
    Parámetros:
    -----------
    Wr : float
        Ancho de la hilera (m)
    Hr : float
        Altura de la hilera (m)
    Wp : float
        Ancho del camino entre hileras (m)
    num_celdas : int
        Número de celdas por lado para discretizar el rectangulo
    df_datos_doy : DataFrame
        Debe contener las columnas:
        - time: tiempo 
        - theta: elevación solar en grados
        - phi: azimut solar en grados (0=Norte, 90=Este)
        - PARdr: PAR directa incidente (W/m² o mol/m²/s)
        
    Retorna:
    --------
    PARdr_celda : ndarray
        Matriz num_celdas×num_celdas con la PAR directa integrada diaria por celda
    """
    Wr = param_arquitectura["Wr"]
    Hr = param_arquitectura["Hr"]
    Wp = param_arquitectura["Wp"]
    eta = np.radians(param_arquitectura["eta"])
    sigma = param_modelo_luz["O_av"]
    O_av = param_modelo_luz["O_av"]
    LAI = param_modelo_luz["LAI"]
    Ld = LAI / Hr
    
    # Crear malla de puntos (centros de celdas)
    x_edges = np.linspace(0, Wr, num_celdas_x + 1)
    z_edges = np.linspace(0, Hr, num_celdas_z + 1)
    x_centers = (x_edges[:-1] + x_edges[1:]) / 2
    z_centers = (z_edges[:-1] + z_edges[1:]) / 2
    
    
    x_grid, z_grid = np.meshgrid(x_centers, z_centers)
    
    # 4. Inicializar matriz de acumulación de PAR directa
    PAR_directa_acum = np.zeros((num_celdas_z, num_celdas_x))
    
    # 5. Iterar sobre cada paso de tiempo.
    alpha_c_list = []
    beta_c_list = []
    beta_list = []
    alpha_list = []

    eta = np.radians(eta)
    delta_t = int(pd.Timedelta(freq_min).total_seconds())  

    for idx, row in df_datos_doy.iterrows():
        # Saltar si no hay radiación directa
        if row['PARdr (mol PAR/m^2/s)'] <= 0:
            continue
            
        # Ángulos solares (convertir a radianes)
        beta = np.radians(row['theta'])  # Elevación
        phi_sun = np.radians(row['phi'])  # Azimut solar
        
        #alpha = np.abs((phi_sun - eta + np.pi) % (2*np.pi) - np.pi)  # Diferencia entre azimut solar y de hilera
        alpha = phi_sun-eta
        alpha_prime = alpha % (2*np.pi) # entre 0 y 2pi

        # Calcular ángulos transformados (ec. 1-2 de G&G)
        sin_beta_c = np.cos(alpha) * np.cos(beta)
        beta_c = np.arcsin(sin_beta_c) 
        
        cos_beta_c = np.cos(beta_c)
        if abs(cos_beta_c) < 1e-10:
            cos_beta_c = 1e-10  # Evitar división por cero

        cos_alpha_c = np.sin(beta)/cos_beta_c
        alpha_c = np.arccos(cos_alpha_c)
        if alpha_prime>np.pi:
            alpha_c = -alpha_c # criterio de signos: si eta=0,  alpha_c>0 en los cuadrantes NE y SE y alpha<0 en los cuadrnates NO y SO

            
        # Coeficiente de extinción para luz directa (ec. 13-14 de G&G)
        K_dir_bl = O_av / np.sin(beta) if np.sin(beta) > 0 else 100  # Evitar división por cero
        K_dir = K_dir_bl * np.sqrt(1 - sigma)
        
        PAR_dr_inst = row['PARdr (mol PAR/m^2/s)'] # PAR directa incidente en este instante
        
        # Calcular para cada celda
        alpha_c_list.append(alpha_c)
        beta_c_list.append(beta_c)
        beta_list.append(beta)
        alpha_list.append(alpha)

        
        for i in range(num_celdas_x):
            for j in range(num_celdas_z):
                x = x_grid[j, i]  # Distancia desde borde izquierdo de hilera
                z = z_grid[j, i]  # Distancia desde tope de hilera
                if x < 0 or x > Wr or z < 0 or z > Hr:  # Saltar si el punto está fuera de la hilera
                    continue
                
                L_t = calcular_Lt(Wr, Wp, Ld, alpha_c, beta_c, x, z)
                I_dir = PAR_dr_inst * np.exp(-K_dir * L_t) / np.sin(beta) # Usamos la forma simplificada equivalente a Beer-Lambert
                
                PAR_directa_acum[j, i] += I_dir*delta_t # Acumular (suponiendo que df_datos_doy tiene pasos de tiempo regulares)
                
    df_angulos = pd.DataFrame({
        "alpha": np.degrees(alpha_list),
        "beta": np.degrees(beta_list),
    'alpha_c_grados': np.degrees(alpha_c_list),
    'beta_c_grados': np.degrees(beta_c_list),
    "sin_beta": np.sin(beta_list),
    "cos_alpha_c_sin_alpha_c": np.cos(alpha_c_list)*np.cos(beta_c_list),
    "sin_beta_c": np.sin(beta_c_list),
    "cos_alpha_cos_beta": np.cos(alpha_list)*np.cos(beta_list),
    "tan_alpha_c": np.tan(alpha_c_list),
    "sin_alpha_c": np.sin(alpha_c_list),
    "cos_beta_c": np.cos(beta_c_list)
})  
    
    return PAR_directa_acum, df_angulos

def calcular_Lt(Wr, Wp, Ld,  alpha_c, beta_c, x, z):
     # Calcular longitud del camino en el plano XZ (ec. 3-7 de G&G)
    if abs(alpha_c) < 1e-6:
        alpha_c = np.sign(alpha_c) * 1e-6  # aseguras valor mínimo

    cos_beta_c = np.cos(beta_c)
    cos_beta_c = max(abs(cos_beta_c), 1e-6) * np.sign(cos_beta_c)

    x_hor = z * np.tan(alpha_c)
    total_dist = x_hor + x
    unit_size = Wp + Wr
    
    x_r = total_dist - np.floor(total_dist / unit_size) * unit_size
    N_u = np.floor(total_dist / unit_size)
    
    if x_r <= Wr:
        P_r_prime = (N_u * Wr - x + x_r) / np.sin(alpha_c)
    else:
        P_r_prime = ((N_u + 1) * Wr - x) / np.sin(alpha_c)
    
    # Longitud total del camino en 3D (ec. 8 de G&G)
    P_r = P_r_prime / cos_beta_c  # Es importante que -pi/2 < cos(beta_c) < pi/2
    L_t = P_r * Ld # Área foliar atravesada (ec. 9 de G&G)
    return L_t

def calcular_PAR_difusa(Wr, Hr, Wp, num_celdas_x, num_celdas_z, PARdf):
    # Crear malla de puntos (centros de celdas)
    x_edges = np.linspace(0, Wr, num_celdas_x + 1)
    z_edges = np.linspace(0, Hr, num_celdas_z + 1)
    x_centers = (x_edges[:-1] + x_edges[1:]) / 2
    z_centers = (z_edges[:-1] + z_edges[1:]) / 2
    

    x_grid, z_grid = np.meshgrid(x_centers, z_centers)
    
    sigma = 0.2  
    O_av = 0.5  
    rho = 0.1 # rho 0.05 - 0.10
    
    Ld = 2.0 / Hr # leaf area density 
    PAR_difusa_celda = np.zeros((num_celdas_z, num_celdas_x))
    
    def integrando_difusa(beta_c, alpha_c, x, z):
        Lt = calcular_Lt(Wr, Wp, Ld, alpha_c, beta_c, x, z)
        return 1/np.pi*(1-rho)*np.exp(-O_av*Lt*np.sqrt(1-sigma))*np.cos(beta_c)**2*np.cos(alpha_c)
    
    def integrar_simpson(x, z):
        num_divisiones = 20
        angulos = np.linspace(-np.pi/2, np.pi/2, num_divisiones)
        alpha_c_grid, beta_c_grid = np.meshgrid(angulos, angulos)
        integrando_eval = np.zeros((len(angulos), len(angulos)))
        for i in range(len(angulos)):
            for j in range(len(angulos)):
                integrando_eval[i,j] = integrando_difusa(beta_c_grid[i, j], alpha_c_grid[i, j], x, z)
        
        inner = simpson(integrando_eval, angulos, axis = 0)
        resultado = simpson(inner, angulos)
           
        return resultado
    
    for i in range(num_celdas_x):
            for j in range(num_celdas_z):
                x = x_grid[j, i]  # Distancia desde borde izquierdo de hilera
                z = z_grid[j, i]  # Distancia desde tope de hile
                result = integrar_simpson(x,z)
                PAR_difusa_celda[j, i] = PARdf*result
                                
    return PAR_difusa_celda






        