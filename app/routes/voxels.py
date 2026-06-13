import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm # Librería para barras de progreso, útil para procesos largos

# --- Funciones Base ---

def calcular_angulos_solares(latitud, dia_ano):
    """
    Calcula la hora de salida y puesta del sol para un día y una ubicación.
    """
    declinacion_rad = np.radians(23.45 * np.sin(np.radians(360.0 / 365.0 * (dia_ano - 81))))
    latitud_rad = np.radians(latitud)
    cos_h_salida = -np.tan(latitud_rad) * np.tan(declinacion_rad)
    if abs(cos_h_salida) > 1:
        return None, None
    h_salida_rad = np.arccos(cos_h_salida)
    duracion_dia_horas = 2 * np.degrees(h_salida_rad) / 15
    hora_salida = 12 - duracion_dia_horas / 2
    hora_puesta = 12 + duracion_dia_horas / 2
    return hora_salida, hora_puesta

def simular_irradiancia_hora(params, dia_ano, hora_solar):
    """Simula la irradiancia PAR (umol/m^2/s) incidente exterior y los ángulos solares."""
    h, w, t, d, eta, psi, latitud, k, sigma = params.values()
    declinacion_rad = np.radians(23.45 * np.sin(np.radians(360.0 / 365.0 * (dia_ano - 81))))
    latitud_rad = np.radians(latitud)
    hora_angular_rad = np.radians((hora_solar - 12) * 15)
    
    sin_elev = np.sin(latitud_rad) * np.sin(declinacion_rad) + np.cos(latitud_rad) * np.cos(declinacion_rad) * np.cos(hora_angular_rad)
    if sin_elev <= 0: return 0, 0, 0, 0
    elev_rad = np.arcsin(sin_elev)
    elev_deg = np.degrees(elev_rad)
    
    cos_phi_sol = (np.sin(declinacion_rad) - np.sin(latitud_rad) * np.sin(elev_rad)) / (np.cos(latitud_rad) * np.cos(elev_rad))
    phi_sol_rad = np.arccos(np.clip(cos_phi_sol, -1, 1))
    if hora_solar > 12: phi_sol_rad = 2 * np.pi - phi_sol_rad

    Idr_total = 1000 * (1 - np.exp(-0.06 * elev_deg)) # W/m^2 [cite: 651]
    Idf_total = 10 * (0.55 + 9.6 * (1 - np.exp(-0.05 * elev_deg))) # W/m^2 [cite: 652]
    
    FACTOR_CONVERSION_PAR = 0.5 * 4.57
    par_directa_incidente = Idr_total * FACTOR_CONVERSION_PAR
    par_difusa_incidente = Idf_total * FACTOR_CONVERSION_PAR
    
    return par_directa_incidente, par_difusa_incidente, phi_sol_rad, elev_rad

def simular_dia_completo(params, dia_ano, num_celdas_x, num_celdas_z):
    """Ejecuta la simulación para un día completo con cálculo de trayectoria combinado."""
    h, w, t, d, eta, psi, latitud, k, sigma = params.values()
    dosel_par_diario = np.zeros((num_celdas_z, num_celdas_x))
    
    hora_salida, hora_puesta = calcular_angulos_solares(latitud, dia_ano)
    if hora_salida is None: return dosel_par_diario

    for hora in range(int(np.ceil(hora_salida)), int(np.floor(hora_puesta)) + 1):
        par_dir, par_dif, phi_sol_rad, elev_rad = simular_irradiancia_hora(params, dia_ano, hora)
        if elev_rad <= 0: continue
            
        par_hora_actual = np.zeros((num_celdas_z, num_celdas_x))
        altura_dosel = h - t
        ancho_celda = w / num_celdas_x
        alto_celda = altura_dosel / num_celdas_z
        psi_rad, eta_rad = np.radians(psi), np.radians(eta)
        phi_rel_rad = np.abs(phi_sol_rad - eta_rad)
        
        par_dir_sup = par_dir * (np.sin(elev_rad) * np.sin(psi_rad) + np.cos(elev_rad) * np.cos(psi_rad) * np.cos(phi_rel_rad))
        par_dif_sup = par_dif * (1 - 0.5 * (psi / 90))
        par_total_superficie = par_dir_sup + par_dif_sup

        for iz in range(num_celdas_z):
            for ix in range(num_celdas_x):
                # --- INICIO DE LA CORRECCIÓN LÓGICA COMBINADA ---
                z_celda = (iz + 0.5) * alto_celda
                x_celda_centro = (ix + 0.5) * ancho_celda - w / 2
                
                # 1. Calcular trayectoria vertical (desde arriba)
                s_vertical = z_celda / np.sin(elev_rad) if np.sin(elev_rad) > 1e-6 else float('inf')

                # 2. Calcular trayectoria horizontal (desde el lateral)
                dist_horiz_borde = (w / 2) - np.abs(x_celda_centro)
                denominador_h = np.cos(elev_rad) * np.abs(np.sin(phi_rel_rad))
                s_horizontal = dist_horiz_borde / denominador_h if denominador_h > 1e-6 else float('inf')
                
                # 3. La trayectoria real 's' es la más corta de las dos
                s = min(s_vertical, s_horizontal)
                # --- FIN DE LA CORRECCIÓN LÓGICA COMBINADA ---

                tau = np.exp(-sigma * k * s)
                par_celda = tau * par_total_superficie
                par_hora_actual[iz, ix] = max(par_celda, 0)
        
        dosel_par_diario += par_hora_actual * 3600
        
    return dosel_par_diario / 1e6

def simular_anual_completo(params, num_celdas_x, num_celdas_z):
    """Ejecuta la simulación para los 365 días del año."""
    dosel_par_anual = np.zeros((num_celdas_z, num_celdas_x))
    print("Iniciando simulación anual...")
    for dia in tqdm(range(1, 366), desc="Simulando días del año"):
        dosel_diario = simular_dia_completo(params, dia, num_celdas_x, num_celdas_z)
        dosel_par_anual += dosel_diario
    return dosel_par_anual

def visualizar_heatmap_anual(dosel_par, params):
    """Muestra un mapa de calor con la distribución de PAR anual y ejes corregidos."""
    h, w, t, d, eta, psi, latitud, k, sigma = params.values()
    plt.figure(figsize=(8, 8))
    im = plt.imshow(dosel_par, cmap='inferno', extent=[0, w, t, h], aspect='auto')
    plt.colorbar(im, label='PAR Anual Acumulado (mol/m²/año)')
    plt.title('Distribución Anual de PAR en la Sección Transversal del Seto')
    plt.xlabel('Anchura del seto (m)')
    plt.ylabel('Altura desde el suelo (m)')
    plt.grid(True, color='white', linestyle='--', linewidth=0.5, alpha=0.3)
    plt.show()

# --- Parámetros de Simulación y Ejecución ---

if __name__ == '__main__':
    parametros_explotacion = {
        "h": 2.5, "w": 1.2, "t": 0.5, "d": 4.0, "eta": 90,
        "psi": 0, "latitud": 39.8, "k": -np.log((20 / 100) / 1.2),
        "sigma": 0.9
    }
    num_celdas_x = 20
    num_celdas_z = 20

    mapa_par_anual = simular_anual_completo(parametros_explotacion, num_celdas_x, num_celdas_z)
    visualizar_heatmap_anual(mapa_par_anual, parametros_explotacion)