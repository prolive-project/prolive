import numpy as np
import time

# Se importa directamente desde skopt y sus submódulos
from skopt import gp_minimize
from skopt.space import Real, Integer
from skopt.utils import use_named_args

# --- MÓDULO 1: SIMULACIÓN ---

def calcular_angulos_solares(latitud, dia_ano):
    declinacion_rad = np.radians(23.45 * np.sin(np.radians(360.0 / 365.0 * (dia_ano - 81))))
    latitud_rad = np.radians(latitud)
    cos_h_salida = -np.tan(latitud_rad) * np.tan(declinacion_rad)
    if abs(cos_h_salida) > 1: return None, None
    h_salida_rad = np.arccos(cos_h_salida)
    duracion_dia_horas = 2 * np.degrees(h_salida_rad) / 15
    return 12 - duracion_dia_horas / 2, 12 + duracion_dia_horas / 2

def simular_irradiancia_hora(latitud, dia_ano, hora_solar):
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
    return Idr_total * FACTOR_CONVERSION_PAR, Idf_total * FACTOR_CONVERSION_PAR, phi_sol_rad, elev_rad

# ===== INICIO DE LA FUNCIÓN CORREGIDA =====
def simular_anual(params, num_celdas_x, num_celdas_z):
    """
    VERSIÓN DE DEPURACIÓN para capturar el estado de las variables antes del error.
    """
    h, w, t, d, eta, psi, latitud, k, sigma = params.values()
    dosel_par_anual = np.zeros((num_celdas_z, num_celdas_x))

    for dia in range(1, 366):
        dosel_par_diario = np.zeros((num_celdas_z, num_celdas_x))
        
        hora_salida, hora_puesta = calcular_angulos_solares(latitud, dia)
        if hora_salida is None: continue
        
        altura_dosel = h - t
        # Protección contra altura de dosel no positiva
        if altura_dosel <= 0: continue

        ancho_celda = w / num_celdas_x
        alto_celda = altura_dosel / num_celdas_z
        x_coords = np.linspace(-w/2 + ancho_celda/2, w/2 - ancho_celda/2, num_celdas_x)
        z_coords = np.linspace(alto_celda/2, altura_dosel - alto_celda/2, num_celdas_z)
        xx, zz = np.meshgrid(x_coords, z_coords)
        psi_rad = np.radians(psi)
        ancho_superior_seto = w/2 - (zz * np.tan(psi_rad))
        mascara_trapezoidal = np.abs(xx) < ancho_superior_seto

        for hora in range(int(np.ceil(hora_salida)), int(np.floor(hora_puesta)) + 1):
            par_dir, par_dif, phi_sol_rad, elev_rad = simular_irradiancia_hora(latitud, dia, hora)
            if elev_rad <= 0: continue
            
            eta_rad = np.radians(eta)
            phi_rel_rad = np.abs(phi_sol_rad - eta_rad)
            par_dir_sup = par_dir * (np.sin(elev_rad) * np.sin(psi_rad) + np.cos(elev_rad) * np.cos(psi_rad) * np.cos(phi_rel_rad))
            par_dif_sup = par_dif * (1 - 0.5 * (psi / 90))
            par_total_superficie = par_dir_sup + par_dif_sup
            
            s_vertical = zz / np.sin(elev_rad) if np.sin(elev_rad) > 1e-6 else float('inf')
            dist_horiz_borde = np.maximum(0, ancho_superior_seto - np.abs(xx))
            denominador_h = np.cos(elev_rad) * np.abs(np.sin(phi_rel_rad))
            s_horizontal = dist_horiz_borde / denominador_h if denominador_h > 1e-6 else float('inf')
            s = np.minimum(s_vertical, s_horizontal)
            
            tau = np.exp(-sigma * k * s)
            par_hora_actual = tau * par_total_superficie

            # --- INICIO DEL BLOQUE DE DEPURACIÓN ---
            try:
                # Esta es la línea que falla
                par_hora_actual[~mascara_trapezoidal] = 0
            except TypeError:
                # Si falla, imprimimos toda la información posible y detenemos el programa
                print("\n" + "!"*60)
                print("!!! ERROR FATAL DETECTADO. DUMP DE DEPURACIÓN: !!!")
                print(f"Parámetros actuales: h={h}, w={w}, d={d}, eta={eta}, psi={psi}")
                print(f"Día del año: {dia}, Hora solar: {hora}")
                print(f"Variable 'tau': tipo={type(tau)}, shape={getattr(tau, 'shape', 'N/A')}")
                print(f"Variable 'par_total_superficie': tipo={type(par_total_superficie)}, valor={par_total_superficie}")
                print(f"Variable 'par_hora_actual': tipo={type(par_hora_actual)}, shape={getattr(par_hora_actual, 'shape', 'N/A')}")
                print(f"Contenido de 'par_hora_actual':\n{par_hora_actual}")
                print("!"*60 + "\n")
                # Relanzamos el error para que el programa se detenga como antes
                raise
            # --- FIN DEL BLOQUE DE DEPURACIÓN ---

            dosel_par_diario += np.maximum(par_hora_actual, 0) * 3600

        dosel_par_anual += dosel_par_diario
    
    return dosel_par_anual / 1e6

# ===== FIN DE LA FUNCIÓN CORREGIDA =====


# --- MÓDULO 2: CONFIGURACIÓN DE LA OPTIMIZACIÓN ---

# Constantes Fisiológicas y de Simulación
UMBRAL_PAR_MJ = 1500
UMBRAL_PAR_MOL = UMBRAL_PAR_MJ / 0.219
NUM_CELDAS_X = 15
NUM_CELDAS_Z = 15
N_LLAMADAS_OPTIMIZADOR = 75
contador_llamadas = 0

# 1. Definir el espacio de búsqueda para los parámetros
espacio_busqueda = [
    Real(2.0, 4.0, name='h'),
    Real(0.8, 2.0, name='w'),
    Real(3.5, 7.0, name='d'),
    Integer(0, 90, name='eta'),
    Integer(0, 20, name='psi')
]

# 2. Definir la Función Objetivo a minimizar
@use_named_args(espacio_busqueda)
def funcion_objetivo(h, w, d, eta, psi):
    global contador_llamadas, PRIMERA_LLAMADA
    contador_llamadas += 1
    
    # Verificación para asegurarnos de que se usa la nueva función
    if contador_llamadas == 1:
        print("--- Verificación: Ejecutando la versión de código corregida. ---")

    parametros_simulacion = {
        "h": h, "w": w, "t": 0.5, "d": d, "eta": float(eta),
        "psi": float(psi), "latitud": LATITUD_FIJA, "k": 0.8, "sigma": 0.9
    }

    mapa_par_anual = simular_anual(parametros_simulacion, NUM_CELDAS_X, NUM_CELDAS_Z)
    
    volumen_productivo = np.sum(mapa_par_anual > UMBRAL_PAR_MOL)
    
    print(f"--> Llamada {contador_llamadas}/{N_LLAMADAS_OPTIMIZADOR}: [h={h:.2f}, w={w:.2f}, d={d:.2f}, eta={eta}, psi={psi}] - Volumen Productivo: {volumen_productivo}")
    
    return -volumen_productivo

# --- MÓDULO 3: EJECUCIÓN PRINCIPAL ---

if __name__ == '__main__':
    try:
        LATITUD_FIJA = float(input("Introduce la latitud de la finca (ej. para Badajoz: 38.88): "))
    except ValueError:
        print("Entrada no válida. Usando latitud por defecto de 38.88.")
        LATITUD_FIJA = 38.88

    print("\n" + "="*50)
    print("🚀 Iniciando Proceso de Optimización Bayesiana 🚀")
    print(f"Objetivo: Maximizar el 'Volumen Productivo' (Vóxeles > {UMBRAL_PAR_MOL:.0f} mol/m²/año)")
    print(f"Se realizarán {N_LLAMADAS_OPTIMIZADOR} simulaciones para encontrar el diseño óptimo.")
    print("="*50 + "\n")
    
    start_time = time.time()
    
    resultado_optimizacion = gp_minimize(
        func=funcion_objetivo,
        dimensions=espacio_busqueda,
        n_calls=N_LLAMADAS_OPTIMIZADOR,
        n_initial_points=10,
        acq_func='EI',
        random_state=42
    )
    
    end_time = time.time()
    
    mejor_volumen = -resultado_optimizacion.fun
    mejores_parametros = resultado_optimizacion.x
    
    print("\n" + "="*50)
    print("🎉 ¡Optimización Finalizada! 🎉")
    print(f"Tiempo total de ejecución: {((end_time - start_time) / 60):.2f} minutos.")
    print("-"*50)
    print(f"🏆 Mejor Volumen Productivo encontrado: {mejor_volumen:.0f} vóxeles")
    print("Parámetros del Diseño Óptimo:")
    print(f"  - Altura del seto (h):   {mejores_parametros[0]:.2f} m")
    print(f"  - Anchura de la base (w): {mejores_parametros[1]:.2f} m")
    print(f"  - Distancia entre hileras (d): {mejores_parametros[2]:.2f} m")
    print(f"  - Orientación (eta):     {mejores_parametros[3]}° (0=N-S, 90=E-O)")
    print(f"  - Ángulo de pared (psi): {mejores_parametros[4]}°")
    print("="*50)