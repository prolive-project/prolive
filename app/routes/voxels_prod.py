import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import time

# --- 1: CÁLCULOS SOLARES Y DE IRRADIANCIA ---

def calcular_angulos_solares(latitud, dia_ano):
    """Calcula la hora de salida y puesta del sol."""
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
    """Simula la irradiancia PAR incidente exterior y los ángulos solares."""
    latitud = params["latitud"]
    declinacion_rad = np.radians(23.45 * np.sin(np.radians(360.0 / 365.0 * (dia_ano - 81))))
    latitud_rad = np.radians(latitud)
    hora_angular_rad = np.radians((hora_solar - 12) * 15)
    
    sin_elev = np.sin(latitud_rad) * np.sin(declinacion_rad) + np.cos(latitud_rad) * np.cos(declinacion_rad) * np.cos(hora_angular_rad)
    if sin_elev <= 0: return 0, 0, 0, 0
    elev_rad = np.arcsin(sin_elev)
    
    cos_phi_sol = (np.sin(declinacion_rad) - np.sin(latitud_rad) * np.sin(elev_rad)) / (np.cos(latitud_rad) * np.cos(elev_rad))
    phi_sol_rad = np.arccos(np.clip(cos_phi_sol, -1, 1))
    if hora_solar > 12: phi_sol_rad = 2 * np.pi - phi_sol_rad
    
    elev_deg = np.degrees(elev_rad)
    Idr_total = 1000 * (1 - np.exp(-0.06 * elev_deg))
    Idf_total = 10 * (0.55 + 9.6 * (1 - np.exp(-0.05 * elev_deg)))
    
    FACTOR_CONVERSION_PAR = 0.5 * 4.57
    par_directa_incidente = Idr_total * FACTOR_CONVERSION_PAR
    par_difusa_incidente = Idf_total * FACTOR_CONVERSION_PAR
    
    return par_directa_incidente, par_difusa_incidente, phi_sol_rad, elev_rad

# --- 2: SIMULACIÓN DEL DOSEL ---

def simular_dia_completo_vectorizado(params, dia_ano, num_celdas_x, num_celdas_z):
    """Ejecuta la simulación para un día completo usando operaciones vectorizadas de NumPy."""
    h, w, t, d, eta, psi, latitud, k, sigma = params.values()
    hora_salida, hora_puesta = calcular_angulos_solares(latitud, dia_ano)
    if hora_salida is None: return np.zeros((num_celdas_z, num_celdas_x))

    altura_dosel = h - t
    ancho_celda = w / num_celdas_x
    alto_celda = altura_dosel / num_celdas_z
    x_coords = np.linspace(-w/2 + ancho_celda/2, w/2 - ancho_celda/2, num_celdas_x)
    z_coords = np.linspace(alto_celda/2, altura_dosel - alto_celda/2, num_celdas_z)
    xx, zz = np.meshgrid(x_coords, z_coords)

    psi_rad = np.radians(psi)
    ancho_superior_seto = w/2 - (zz * np.tan(psi_rad))
    mascara_trapezoidal = np.abs(xx) < ancho_superior_seto
    dosel_par_diario = np.zeros((num_celdas_z, num_celdas_x))

    for hora in range(int(np.ceil(hora_salida)), int(np.floor(hora_puesta)) + 1):
        par_dir, par_dif, phi_sol_rad, elev_rad = simular_irradiancia_hora(params, dia_ano, hora)
        if elev_rad <= 0: continue
        eta_rad, psi_rad = np.radians(eta), np.radians(psi)
        phi_rel_rad = np.abs(phi_sol_rad - eta_rad)

        par_dir_sup = par_dir * (np.sin(elev_rad) * np.sin(psi_rad) + np.cos(elev_rad) * np.cos(psi_rad) * np.cos(phi_rel_rad))
        par_dif_sup = par_dif * (1 - 0.5 * (psi / 90))
        par_total_superficie = par_dir_sup + par_dif_sup
        
        s_vertical = zz / np.sin(elev_rad) if np.sin(elev_rad) > 1e-6 else float('inf')
        dist_horiz_borde = ancho_superior_seto - np.abs(xx)
        denominador_h = np.cos(elev_rad) * np.abs(np.sin(phi_rel_rad))
        s_horizontal = dist_horiz_borde / denominador_h if denominador_h > 1e-6 else float('inf')
        s = np.minimum(s_vertical, s_horizontal)
        
        tau = np.exp(-sigma * k * s)
        par_hora_actual = tau * par_total_superficie
        par_hora_actual[~mascara_trapezoidal] = 0
        par_hora_actual = np.maximum(par_hora_actual, 0)
        dosel_par_diario += par_hora_actual * 3600

    return dosel_par_diario / 1e6

def simular_anual_completo(params, num_celdas_x, num_celdas_z):
    """Ejecuta la simulación para los 365 días del año usando la función vectorizada."""
    dosel_par_anual = np.zeros((num_celdas_z, num_celdas_x))
    # LA BARRA DE PROGRESO PRINCIPAL ESTÁ AQUÍ
    for dia in tqdm(range(1, 366), desc="[FASE 1/3] Simulando año completo"):
        dosel_diario = simular_dia_completo_vectorizado(params, dia, num_celdas_x, num_celdas_z)
        dosel_par_anual += dosel_diario
    return dosel_par_anual

# --- 3:  PRODUCTIVIDAD ---

def convertir_par_a_produccion(mapa_par_anual):
    """Convierte el mapa de PAR anual a mapas de potencial productivo."""
    par_mj = mapa_par_anual * 0.219 # 1 mol PAR ≈ 0.219 MJ
    
    mapa_tamano_fruto = 0.39 + 0.0051 * par_mj
    mapa_concentracion_aceite = 34.8 + 0.25 * par_mj
    mapa_densidad_frutos = -35 + 59 * par_mj
    
    mapas_produccion = {
        "Tamaño del Fruto (g)": np.maximum(mapa_tamano_fruto, 0),
        "Concentración de Aceite (%)": np.maximum(mapa_concentracion_aceite, 0),
        "Densidad de Frutos (frutos/m²)": np.maximum(mapa_densidad_frutos, 0)
    }
    return mapas_produccion

# --- MÓDULO 4: VISUALIZACIÓN (Sin cambios) ---

def visualizar_resultados(dosel_par, mapas_produccion, params):
    """Muestra el mapa de calor de PAR y los 3 mapas de producción."""
    h, w, t, _, _, _, _, _, _ = params.values()

    plt.figure(figsize=(7, 6))
    im_par = plt.imshow(dosel_par, cmap='inferno', extent=[0, w, t, h], aspect='auto', origin='lower')
    plt.colorbar(im_par, label='PAR Anual Acumulado (mol/m²/año)')
    plt.title('Distribución Anual de PAR en el Seto')
    plt.xlabel('Anchura del seto (m)')
    plt.ylabel('Altura desde el suelo (m)')
    plt.grid(True, color='white', linestyle='--', linewidth=0.5, alpha=0.3)
    plt.show()

    fig, axs = plt.subplots(1, 3, figsize=(20, 6), sharey=True)
    fig.suptitle('Mapas de Potencial Productivo del Seto', fontsize=16)
    for ax, (titulo, mapa) in zip(axs, mapas_produccion.items()):
        im = ax.imshow(mapa, cmap='viridis', extent=[0, w, t, h], aspect='auto', origin='lower')
        fig.colorbar(im, ax=ax, label=titulo)
        ax.set_title(titulo)
        ax.set_xlabel('Anchura del seto (m)')
    axs[0].set_ylabel('Altura desde el suelo (m)')
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()

# --- EJECUCIÓN  ---

if __name__ == '__main__':
    parametros_seto_trapezoidal = {
        "h": 3.0, "w": 1.5, "t": 0, "d": 5.0, "eta": 89.99,
        "psi": 15, "latitud": 38.88, "k": 0.8, "sigma": 0.9
    }
    num_celdas_x = 19
    num_celdas_z = 19

    # --- FASE 1: Simulación ---
    start_time = time.time()
    mapa_par_anual = simular_anual_completo(parametros_seto_trapezoidal, num_celdas_x, num_celdas_z)
    end_time = time.time()
    print(f" Simulación anual completada en {end_time - start_time:.2f} segundos.")

    # --- FASE 2: Conversión a producción ---
    print("\n[FASE 2/3] Convirtiendo mapa de luz a mapas de producción...")
    start_conv_time = time.time()
    mapas_productivos = convertir_par_a_produccion(mapa_par_anual)
    end_conv_time = time.time()
    print(f" Conversión finalizada en {end_conv_time - start_conv_time:.4f} segundos.")

    # --- FASE 3: Visualización ---
    print("\n[FASE 3/3] Generando visualizaciones...")
    visualizar_resultados(mapa_par_anual, mapas_productivos, parametros_seto_trapezoidal)
    print("Proceso finalizado.")