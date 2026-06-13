from flask import Blueprint, render_template, request, jsonify
from flask_login import login_required
import traceback
import numpy as np
import pandas as pd
import pvlib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import io
import base64
from datetime import datetime
from scipy.integrate import trapezoid
from app.routes.rca import calculate_rca


def solar_angles_day(latitude, longitude, day_of_year, freq_min):
    """
    Calcula la posición solar (elevación y azimut) a lo largo de un día
    para una localización geográfica dada, usando pvlib.

    La función genera una serie temporal regular durante el día especificado
    (en UTC) y devuelve los ángulos solares correspondientes, adecuados
    para cálculos de radiación directa y difusa.

    Parámetros
    ----------
    latitude : float
        Latitud del emplazamiento en grados decimales (positiva Norte).
    longitude : float
        Longitud del emplazamiento en grados decimales (positiva Este).
    day_of_year : int
        Día del año (1-365 o 1-366 en año bisiesto).
    freq_min : str
        Frecuencia temporal en minutos, en formato compatible con pandas
        (por ejemplo '5min', '10min', '15min').

    Retorna
    -------
    df : pandas.DataFrame
        DataFrame con una fila por instante temporal y las columnas:
        - 'elevation' : elevación solar en grados sobre el horizonte.
        - 'azimuth'   : azimut solar en grados (0 = Norte, 90 = Este).
        - 'time'      : instante temporal en UTC (Timestamp).
    """
     
    date = pd.Timestamp('2024-01-01') + pd.Timedelta(days=day_of_year - 1)
    times = pd.date_range( 
        start=date, 
        end=date + pd.Timedelta(hours=23.9),
        freq= freq_min,
        tz='UTC'   # Para pvlib es mejor trabajar en UTC
    ) # Secuencia de tiempo cada 0.1 h (6 min)
    
    solpos = pvlib.solarposition.get_solarposition(times, latitude, longitude) # Cálculo de posiciones solares
    df = solpos[['elevation', 'azimuth']].copy() # Extraer elevación y azimut
    df['time'] = times
    
    return df.reset_index(drop=True)

def daily_radiation(angle_df, h, psi, d, eta):
    """
    Cálculo de la radiacion incidente y en la base de acuerdo con Connor 2006.

    Recibe: 
    Elevacion y azimut en funcion del tiempo a lo largo de un dia, DataFrame con columnas: ['time', 'elevation', 'azimuth']
    Parametros de arquitectura h, psi, d, eta

    Devuelve: 
    la radiación diaria incidente en MJ, energia_inc
    la radiacion diaria en la base en MJ, energia_base
    la ratio r = energia_base/energia_inc
    """
    def Idr(theta):
        return 1000 * (1 - np.exp(-0.06 * theta))
    def Idf(theta):
        return 10 * (0.55 + 9.6 * (1 - np.exp(-0.05 * theta)))

    theta = angle_df["elevation"].values  # en grados
    phi = angle_df["azimuth"].values 

    mask = theta > 0 # Solo cuando el sol está sobre el horizonte
    theta = theta[mask]
    phi = phi[mask]
    times = angle_df["time"].values[mask]
    times_sec = times.astype('datetime64[s]').astype(float)

    theta_c = np.degrees(np.arctan(h * np.cos(np.radians(phi - eta)) / (d + h * np.tan(np.radians(psi)))))
    theta_c_n = np.degrees(np.arctan(h / (d + h * np.tan(np.radians(psi)))))
    
    #============= Incidente:========================
    I_direct = Idr(theta) * np.sin(np.radians(theta))
    I_diffuse = Idf(theta)
    I_total = I_direct + I_diffuse  # en W/m²
    energia_inc = trapezoid(I_total, times_sec) / 1e6  # Integración: ∫ I(t) dt  (W/m² * s → J/m²) y convertir J/m² → MJ/m²
    #============== En la base:=====================
    Idr_base = Idr(theta) * (np.sin(np.radians(theta)) * np.sin(np.radians(psi)) +
                              np.cos(np.radians(theta)) * np.cos(np.radians(psi)) * np.cos(np.radians(phi - eta)))
    Idf_base = Idf(theta) * (90 - theta_c_n + psi) / 180
    
    energia_Idf_base = trapezoid(Idf_base, times_sec)
    #mask_direct = theta>theta_c 
    Idr_base = np.maximum(0,np.where(theta>theta_c, Idr_base, 0)) # disciminar si se supera el angulo critico (solo para la radiacion directa)
    energia_Idr_base = trapezoid(Idr_base, times_sec)
    energia_base = (energia_Idr_base + energia_Idf_base)/1e6
    r = energia_base/energia_inc # ratio a comparar con el umbral x
    
    return {"energia_diaria_incidente_MJ": energia_inc, "energia_diaria_base_MJ": energia_base, "ratio": r}


def calculate_critical_angle(solar_data, h, psi, d, eta):
    phi = solar_data['azimuth']
    theta_c = np.degrees(np.arctan(h * np.cos(np.radians(phi - eta)) / (d + h * np.tan(np.radians(psi)))))
    return theta_c


def simular_daily_radiation(range_DOY, lat, lon, freq_min, h, psi, d, eta):
    r = []
    for DOY in range_DOY:
        angulos_solares = solar_angles_day(lat, lon, DOY, freq_min)
        radiacion_diaria = daily_radiation(angulos_solares, h, psi, d, eta) 
        r.append(radiacion_diaria["ratio"])
    r = np.array(r)
    return {"DOY": range_DOY, "ratio": r}

def plot_daily_radiation(datos):
    """
    plot de ratio de radiacion en la base vs DOY
    La x es el umbral elegido para comparar con la ratio
    """
    plt.figure(figsize=(8, 5))
    
    DOY = datos["DOY"] 
    r = datos["ratio"]
       
    plt.plot(DOY, r*100, ".")
    
    plt.xlabel("Day of year")
    plt.ylabel("Ratio (%)")
    plt.title("Ratio de radiación diaria en la base")
    #plt.legend(loc = "best")
    plt.grid()

query_history = []
MAX_QUERIES = 3


def compute_fr_daily(
    lat, lon, rango_DOY, freq_min,
    h, psi, d, w, eta
):

    datos = simular_daily_radiation(
        rango_DOY, lat, lon, freq_min, h, psi, d, eta
    )

    RCA = calculate_rca(h, d, w, psi)

    return datos, RCA


def generate_fr_daily_plot_superpuesto(datos, label):
    global query_history

    query_history.append({
        "DOY": datos["DOY"],
        "ratio": datos["ratio"],
        "label": label
    })

    if len(query_history) > MAX_QUERIES:
        query_history.pop(0)

    img = io.BytesIO()
    fig, ax = plt.subplots(figsize=(8, 5))

    colores = ["tab:blue", "tab:orange", "tab:green"]

    for i, q in enumerate(query_history):
        ax.plot(
            q["DOY"],
            q["ratio"] * 100,
            ".",
            color=colores[i],
            label=q["label"]
        )

    ax.set_xlabel("Day of year")
    ax.set_ylabel("FR diaria (%)")
    ax.set_title("Fracción de radiación interceptada (Connor)")
    ax.set_ylim(0, 80)
    ax.grid(True)
    ax.legend()

    plt.savefig(img, format="png", bbox_inches="tight")
    plt.close(fig)

    img.seek(0)
    return base64.b64encode(img.read()).decode("utf-8")


fr_daily_bp = Blueprint("fr_daily", __name__)


@fr_daily_bp.route("/fr_daily")
@login_required
def fr_daily():
    return render_template("fr_daily.html")


@fr_daily_bp.route("/simulate_fr_daily", methods=["POST"])
def simulate_fr_daily():
    try:
        data = request.get_json()

        lat = float(data.get("lat", 38))
        lon = float(data.get("lon", 0))
        freq_min = data.get("freq_min", "6min")

        DOY_ini = int(data.get("DOY_ini", 150))
        DOY_fin = int(data.get("DOY_fin", 300))
        rango_DOY = np.arange(DOY_ini, DOY_fin)

        # Arquitectura
        h = float(data.get("h", 3))
        psi = float(data.get("psi", 15))
        d = float(data.get("d", 3))
        w = float(data.get("w", 3))
        eta = float(data.get("eta", 0))

        datos, RCA = compute_fr_daily(
    lat, lon, rango_DOY, freq_min,
    h, psi, d, w, 90-eta  # redefinir eta para que eta = 0 sea NS
)    
        label = f"h={h}, d={d}, w={w}, ψ={psi:.1f}°, η={eta}"
        
        img_data = generate_fr_daily_plot_superpuesto(datos, label)

        return jsonify({
            "RCA": RCA,
            "image": img_data
        })

    except Exception as e:
        print("🔥 ERROR EN simulate_fr_daily")
        traceback.print_exc()
        return jsonify({"error": str(e)}), 500
 