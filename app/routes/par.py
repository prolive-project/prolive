
from flask import Blueprint, request, jsonify, render_template
import traceback
from flask_login import login_required
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from app.routes.modulo_2 import *
from app.routes.modulo_5 import plot_traps
import io
import base64
from flask_login import login_required

par_bp = Blueprint("par", __name__)

def simular_PAR_diaria(DOY, lat, lon, arquitectura_olivar, param_modelo_luz, frec_min):
    """
    Devuelve PAR diaria total, directa y difusa (mol PAR/m^2)
    """

    datos_doy, PAR_dia, PARdr_dia, PARdf_dia = calculate_radiation_components(
        lat, lon, DOY, frec_min, None, clear_sky=True
    )

    # Integración diaria (solo cuando el sol está sobre el horizonte)
    mask = datos_doy["theta"] > 0
    segundos = int(pd.Timedelta(frec_min).total_seconds())

    PAR_total = datos_doy.loc[mask, "PAR (mol PAR/m^2/s)"].sum() * segundos
    PAR_directa = datos_doy.loc[mask, "PARdr (mol PAR/m^2/s)"].sum() * segundos
    PAR_difusa = datos_doy.loc[mask, "PARdf (mol PAR/m^2/s)"].sum() * segundos

    return PAR_total, PAR_directa, PAR_difusa, datos_doy.loc[mask]

def generate_par_map(datos_doy, arquitectura_olivar, param_modelo_luz,
                     num_celdas_x, num_celdas_y, frec_min, PARdf_dia):

    matriz_PAR_directa = calcular_PAR_directa_traps(
        arquitectura_olivar,
        param_modelo_luz,
        num_celdas_x,
        num_celdas_y,
        datos_doy,
        frec_min
    )

    matriz_PAR_difusa = calcular_PAR_difusa_traps(
        arquitectura_olivar,
        param_modelo_luz,
        num_celdas_x,
        num_celdas_y,
        PARdf_dia
    )

    matriz_PAR_total = matriz_PAR_directa + matriz_PAR_difusa

    # === Gráfico ===
    img = io.BytesIO()
    fig, ax = plt.subplots(figsize=(6, 6))
    plot_traps(matriz_PAR_total, arquitectura_olivar)
    ax.set_title("Distribución diaria de PAR total")

    plt.savefig(img, format="png", bbox_inches="tight")
    plt.close(fig)
    img.seek(0)

    return base64.b64encode(img.read()).decode("utf-8")


@par_bp.route("/par")
@login_required
def par():
    return render_template("par.html")

@par_bp.route("/simulate_radiation", methods=["POST"])
def simulate_radiation():
    try:
        data = request.get_json()

        DOY = int(data.get("DOY", 170))
        lat = float(data.get("lat", 38))
        lon = 0
        frec_min = data.get("frec_min", "10min")

        arquitectura_olivar = data.get("arquitectura_olivar", {
        "w": 2,
        "t": 2,
        "h": 3,
        "d": 6,
        "eta": 90})

        param_modelo_luz = data.get("param_modelo_luz", {})

        LAI = float(param_modelo_luz.get("LAI", 8))
        sigma = float(param_modelo_luz.get("sigma", 0.2))
        rho = float(param_modelo_luz.get("rho", 0.1))
        albedo = float(param_modelo_luz.get("albedo", 0.3))
        O_av = 0.5
        
        param_modelo_luz = {
            "LAI": LAI,
            "sigma": sigma,
            "rho": rho,
            "albedo": albedo,
            "O_av": O_av
        }


        PAR_total, PAR_directa, PAR_difusa, datos_doy = simular_PAR_diaria(
            DOY, lat, lon, arquitectura_olivar, param_modelo_luz, frec_min
        )

        img_data = generate_par_map(
            datos_doy,
            arquitectura_olivar,
            param_modelo_luz,
            num_celdas_x=20,
            num_celdas_y=20,
            frec_min=frec_min,
            PARdf_dia=PAR_difusa
        )

        return jsonify({
            "PAR_total": PAR_total,
            "PAR_directa": PAR_directa,
            "PAR_difusa": PAR_difusa,
            "image": img_data
        })

    except Exception as e:
        print("🔥 ERROR EN simulate_radiation")
        traceback.print_exc()   # 👈 ESTO ES CLAVE
        return jsonify({
            "error": str(e)
        }), 500
