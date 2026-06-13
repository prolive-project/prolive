from matplotlib import pyplot as plt
import numpy as np
from app.routes.modulo_2 import *
from app.routes.modulo_3_4 import *

# --- VISUALIZACIÓN ---

def plot_traps(matriz_par, param_arquitectura):
    w = param_arquitectura["w"]
    t = param_arquitectura["t"]
    h = param_arquitectura["h"]
    W = max(w,t)
    x_left_top  = (W - t) / 2
    x_right_top = (W + t) / 2
    x_left_bot  = (W - w) / 2
    x_right_bot = (W + w) / 2
    
    x_trap = [x_left_top, x_right_top, x_right_bot, x_left_bot, x_left_top]
    y_trap = [h,           h,           0,            0,          h]
    fig, ax = plt.subplots()
   
    im = ax.imshow(
        matriz_par,
        cmap='inferno',
        extent=[0, W, 0, h],
        #aspect='auto'
    )
    plt.colorbar(im, label='PAR diaria (mol PAR/m²)')
    # Dibujar el contorno del seto
    ax.plot(x_trap, y_trap, color='white', linewidth=2)
    ax.set_xlabel('Anchura del seto (m)')
    ax.set_ylabel('Altura del seto (m)')
    return ax

def plot_produccion(matriz_par, param_arquitectura):
    w = param_arquitectura["w"]
    t = param_arquitectura["t"]
    h = param_arquitectura["h"]
    W = max(w,t)
    x_left_top  = (W - t) / 2
    x_right_top = (W + t) / 2
    x_left_bot  = (W - w) / 2
    x_right_bot = (W + w) / 2
    
    x_trap = [x_left_top, x_right_top, x_right_bot, x_left_bot, x_left_top]
    y_trap = [h,           h,           0,            0,          h]
    fig, ax = plt.subplots()
    im = ax.imshow(
        matriz_par,
        cmap='viridis',
        extent=[0, W, 0, h],
        #aspect='auto'
    )
    ax.plot(x_trap, y_trap, color='white', linewidth=2)
    ax.set_xlabel('Anchura del seto (m)')
    ax.set_ylabel('Altura del seto (m)')
    return ax, im

# --- EJECUCIÓN EN TERMINAL ---

if __name__ == '__main__':
    DOY = 170
    lat = 38
    lon = 0 

    frec_min = "10min"
    num_celdas_x = 8
    num_celdas_y = 12

    arquitectura_olivar = {
        "w" : 2,   # base inferior (m)
        "t" : 2,   # base superior (m)
        "h" : 3,   # altura (m)
        "d" : 6,   # separacion entre setos (m)
        "eta" : 90 # orientacion (eta = 0 en NS, eta = 90 en EO)
    } 
    
    param_modelo_luz = {
        "sigma": 0.2, # Coeficiente de dispersión de la hoja
        "O_av": 0.5,  # Proyección promedio para distribución esférica de ángulos foliares
        "rho" : 0.1,  # coeficiente de reflexion de las hojas
        "LAI" : 8,    # Indice de area foliar
        "albedo" : 0.3 
    }
    
    arquitectura_GG = {
        "Wp" : 2,
        "Wr": 2,
        "Hr": 3,
        "eta": 90
    } # parámetros del modelo de seccion rectangular de Gijzen y Goudrian 
    
    # Datos de produccion
    param_produccion = { 
        "Wf" : 0.7,     # Peso medio del fruto (g)
        "aceite" : 45,  # Porcentaje medio de aceite (%)
        "Wt" : 4.9e3,   # Peso seco total (kg)
        "Ns" : 10,      # Número de hileras
        "Ls" : 100      # Longitud de la hilera
    }
    
      
    def simular_radiacion_diaria():
        PAR_medida_dia = 60 # mol PAR/m^2
        S_gd_medida = PAR_medida_dia/(0.5/0.217)# MJ/m^2
        datos_doy, PAR_dia, PARdr_dia, PARdf_dia = calculate_radiation_components(lat, lon, DOY, frec_min, None, clear_sky=True)
        print(f"Radiación global diaria: {PAR_dia/(0.5/0.217)} MJ/m²")
        print(f"PAR diaria: {PAR_dia} mol PAR/m^2")
        print(datos_doy[datos_doy["theta"]>0].head())
        print(
    f"PAR del día "
    f"{datos_doy['PAR (mol PAR/m^2/s)'].sum() * int(pd.Timedelta(frec_min).total_seconds())}"
)

        datos_doy.to_csv('datos_doy.csv', index=False)
    simular_radiacion_diaria()
    
    def simular_transmision():
        datos_doy = pd.read_csv('datos_doy.csv')
        PARdf_dia = datos_doy["PARdf (mol PAR/m^2/s)"].sum()*int(pd.Timedelta(frec_min).total_seconds())
        PARdr_dia = datos_doy["PARdr (mol PAR/m^2/s)"].sum()*int(pd.Timedelta(frec_min).total_seconds())
        print(f"PAR diaria directa: {PARdr_dia} mol PAR/m^2")
        print(f"PAR diaria difusa: {PARdf_dia} mol PAR/m^2")
        
        matriz_PAR_directa, df_angulos = calcular_PAR_directa(arquitectura_GG, param_modelo_luz, num_celdas_x, num_celdas_y, datos_doy[datos_doy["theta"]>0], frec_min)
        
        
        matriz_PAR_directa_traps = calcular_PAR_directa_traps(arquitectura_olivar, param_modelo_luz, num_celdas_x, num_celdas_y,
                                                              datos_doy[datos_doy["theta"]>0], frec_min)
        ax1 = plot_traps(matriz_PAR_directa_traps, arquitectura_olivar)
        ax1.set_title('Distribución diaria de PAR directa en el seto')
        
        matriz_PAR_difusa_traps = calcular_PAR_difusa_traps(arquitectura_olivar, param_modelo_luz, num_celdas_x, num_celdas_y, PARdf_dia)
        ax2 = plot_traps(matriz_PAR_difusa_traps, arquitectura_olivar)
        ax2.set_title('Distribución diaria de PAR difusa en el seto')

        matriz_PAR = matriz_PAR_directa_traps + matriz_PAR_difusa_traps
        np.save("matriz_PAR.npy", matriz_PAR)
        ax3 = plot_traps(matriz_PAR, arquitectura_olivar)
        ax3.set_title('Distribución diaria de PAR total en el seto')
        
        plt.show()
        print(f"Orientacion: {arquitectura_GG['eta']}")
        print(df_angulos.head(10))
        print(df_angulos.tail(10))
        print(f"PAR diaria directa: {PARdr_dia} mol PAR/m^2")
        print(f"PAR diaria difusa: {PARdf_dia} mol PAR/m^2")
        df_angulos.to_csv('datos_angulos.csv', index=False)
    simular_transmision()

    def simular_produccion_capas():
         matriz_PAR = np.load("matriz_PAR.npy")

         resultados = calcular_produccion_seto(
             matriz_PAR,
             arquitectura_olivar,
             longitud_seto=100.0
         )

         print("")
         print("=== PRODUCCIÓN DEL SETO (Connor 2016) ===")
         print("--- POR METRO DE SETO ---")
         print(f"Frutos por metro de fila: {resultados['frutos_por_metro']:.1f}")
         print(f"Peso seco por metro (kg): {resultados['peso_seco_por_metro_kg']:.3f}")
         print(f"Aceite por metro (kg): {resultados['aceite_por_metro_kg']:.3f}")
         print("--- POR SETO ---")
         print(f"Frutos totales en el seto: {resultados['frutos_totales']:.0f}")
         print(f"Peso seco total (kg): {resultados['peso_seco_total_kg']:.2f}")
         print(f"Peso aceite total (kg): {resultados['aceite_total_kg']:.2f}")
         print("--- POR HECTAREA ---")
         print(f"Frutos totales por hectárea: {resultados['frutos_hectarea']:.0f}")
         print(f"Peso seco total (kg/hectarea): {resultados['peso_seco_hectarea']:.2f}")
         print(f"Peso aceite total (kg/hectarea): {resultados['aceite_hectarea']:.2f}")


         matriz_PAR_capas = resultados["matriz_PAR_capas"]
         matriz_wfruto_capas = resultados["matriz_wfruto_capas"]
         matriz_densidad_capas = resultados["matriz_densidad_capas"]
         matriz_peso_capas = resultados["matriz_peso_capas"]
         matriz_concen_aceite_capas = resultados["matriz_concen_aceite_capas"]
         matriz_aceite_capas = resultados["matriz_aceite_capas"]
         ax1 = plot_traps(matriz_PAR_capas, arquitectura_olivar)
         ax1.set_title('Distribución diaria de PAR total en capas')
         
         
         ax2, im2 = plot_produccion(matriz_wfruto_capas, arquitectura_olivar)
         plt.colorbar(im2, label = "Peso del fruto seco (g)")
         ax2.set_title('Peso de fruto por capas')
        
         ax3, im3 = plot_produccion(matriz_peso_capas, arquitectura_olivar)
         plt.colorbar(im3, label = "Peso total de fruto por metro (g)")
         ax3.set_title('Peso seco por metro de hilera')

         ax4, im4 = plot_produccion(matriz_concen_aceite_capas, arquitectura_olivar)
         plt.colorbar(im4, label = "Concentracion de aceite (%)")
         ax4.set_title('Concentracion de aceite por capas')

         ax5, im5 = plot_produccion(matriz_aceite_capas, arquitectura_olivar)
         plt.colorbar(im5, label = "Produccion de aceite (g)")
         ax5.set_title('Peso de aceite por metro de hilera (g)')
         
         ax6, im6 = plot_produccion(matriz_densidad_capas, arquitectura_olivar)
         plt.colorbar(im6, label = "Densidad de frutos por m^2")
         ax6.set_title('Densidad de frutos por m^2')


         plt.show()   
    #simular_produccion_capas()

    def simular_produccion_voxeles():
        I_sat = 28 
        peso_medio_fruto_g = param_produccion["Wf"]
        aceite_medio = param_produccion["aceite"]
        Ns = param_produccion["Ns"]
        Ls = param_produccion["Ls"]
        peso_total_kg = param_produccion["Wt"]
        densidad_lineal_frutos = peso_total_kg/((peso_medio_fruto_g*1e-3)*Ns*Ls)

        matriz_PAR = np.load("matriz_PAR.npy")
        resultados = transformar_PAR_a_fisiologia_voxs(matriz_PAR, arquitectura_olivar, peso_medio_fruto_g, aceite_medio, densidad_lineal_frutos, I_sat)
        
        ax2, im2 = plot_produccion(resultados["matriz_wfruto_voxs"], arquitectura_olivar)
        plt.colorbar(im2, label = "Peso del fruto seco (g)")
        ax2.set_title('Peso de fruto por voxel')
        
        ax4, im4 = plot_produccion(resultados["matriz_concen_aceite_voxs"], arquitectura_olivar)
        plt.colorbar(im4, label = "Concentracion de aceite (%)")
        ax4.set_title('Concentracion de aceite por voxel')
     
        ax6, im6 = plot_produccion(resultados["matriz_densidad_voxs"], arquitectura_olivar)
        plt.colorbar(im6, label = "Densidad de frutos por metro de hilera")
        ax6.set_title('Densidad de frutos por m')

        plt.show()
    #simular_produccion_voxeles()

