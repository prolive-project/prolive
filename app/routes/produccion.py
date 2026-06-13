from flask import Blueprint, request, jsonify, render_template
import traceback
from flask_login import login_required

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import numpy as np
import io
import base64

from app.routes.par import simular_PAR_diaria
from app.routes.modulo_2 import (
    calcular_PAR_directa_traps,
    calcular_PAR_difusa_traps
)
from app.routes.modulo_3_4 import transformar_PAR_a_fisiologia_voxs
from app.routes.modulo_5 import plot_traps, plot_produccion

prod_bp = Blueprint("prod", __name__)
def calcular_PAR_media_rango(
    rango_DOY,
    lat,
    lon,
    arquitectura_olivar,
    param_modelo_luz,
    frec_min,
    num_celdas_x,
    num_celdas_y
):
    PAR_acumulada = None
    n_dias = 0

    for DOY in rango_DOY:
        PAR_total_dia, _, PAR_difusa, datos_doy = simular_PAR_diaria(
            DOY, lat, lon, arquitectura_olivar, param_modelo_luz, frec_min
        )

        PAR_dir = calcular_PAR_directa_traps(
            arquitectura_olivar,
            param_modelo_luz,
            num_celdas_x,
            num_celdas_y,
            datos_doy,
            frec_min
        )

        PAR_dif = calcular_PAR_difusa_traps(
            arquitectura_olivar,
            param_modelo_luz,
            num_celdas_x,
            num_celdas_y,
            PAR_difusa
        )

        PAR_total = PAR_dir + PAR_dif

        if PAR_acumulada is None:
            PAR_acumulada = PAR_total
        else:
            PAR_acumulada += PAR_total

        n_dias += 1

    PAR_media = PAR_acumulada / n_dias
    return PAR_media

def generate_prod_results(
    PAR_media,
    arquitectura_olivar,
    param_produccion
):
    densidad_lineal = (
        param_produccion["Wt"] /
        ((param_produccion["Wf"] * 1e-3)
         * param_produccion["Ns"]
         * param_produccion["Ls"])
    )

    resultados = transformar_PAR_a_fisiologia_voxs(
        PAR_media,
        arquitectura_olivar,
        param_produccion["Wf"],
        param_produccion["aceite"],
        densidad_lineal,
        param_produccion["I_sat"]
    )

    return resultados


def plot_prod_maps(PAR_total, resultados, arquitectura_olivar):

    images = {}

    def fig_to_base64(fig):
        buf = io.BytesIO()
        fig.savefig(buf, format="png", bbox_inches="tight")
        plt.close(fig)
        buf.seek(0)
        return base64.b64encode(buf.read()).decode("utf-8")

    # === Peso de fruto ===
    ax, im = plot_produccion(resultados["matriz_wfruto_voxs"], arquitectura_olivar)
    plt.colorbar(im, ax=ax, label="Peso del fruto seco (g)")
    ax.set_title("Peso de fruto por voxel")
    fig = ax.figure
    images["wfruto"] = fig_to_base64(fig)
    plt.close(fig)

    # === Concentración de aceite ===
    ax, im = plot_produccion(resultados["matriz_concen_aceite_voxs"], arquitectura_olivar)
    plt.colorbar(im, ax=ax, label="Concentración de aceite (%)")
    ax.set_title("Concentración de aceite por voxel")
    fig = ax.figure
    images["aceite"] = fig_to_base64(fig)
    plt.close(fig)

    # === Densidad de frutos ===
    ax, im = plot_produccion(resultados["matriz_densidad_voxs"], arquitectura_olivar)
    plt.colorbar(im, ax=ax, label="Densidad de frutos (frutos/m)")
    ax.set_title("Densidad de frutos por metro")
    fig = ax.figure
    images["densidad"] = fig_to_base64(fig)
    plt.close(fig)

    return images

@prod_bp.route("/prod")
@login_required
def prod():
    return render_template("prod.html")

@prod_bp.route("/simulate_prod", methods=["POST"])
@login_required
def simulate_prod():
    try:
        data = request.get_json()

        DOY_ini = int(data.get("DOY_ini", 150))
        DOY_fin = int(data.get("DOY_fin", 250))
     
        lat = float(data.get("lat", 38))
        lon = 0.0
        frec_min = data.get("frec_min", "10min")

        arquitectura_olivar = data.get("arquitectura_olivar", {
            "w": 2,
            "t": 2,
            "h": 3,
            "d": 6,
            "eta": 90
        })

        param_modelo_luz = data.get("param_modelo_luz", {})
        param_modelo_luz["O_av"] = 0.5  # 

        param_produccion = data.get("param_produccion", {
            "I_sat": 28,  # Radiación a la que satura la densidad de frutos (mol PAR/m^2)
            "Wf": 0.7,    # peso seco del fruto (g)
            "aceite": 45, # porcentaje medio de aceite respecto al peso seco (%)
            "Wt": 5e3,  # Peso seco total (kg)
            "Ns": 10,     # Numero de setos
            "Ls": 100     # Longitud de las hileras
        })

        PAR_media = calcular_PAR_media_rango(
            range(DOY_ini, DOY_fin + 1),
            lat,
            lon,
            arquitectura_olivar,
            param_modelo_luz,
            frec_min,
            num_celdas_x=8,
            num_celdas_y=12
        )

        resultados = generate_prod_results(
            PAR_media,
            arquitectura_olivar,
            param_produccion
        )

        images = plot_prod_maps(PAR_media, resultados, arquitectura_olivar)

        return jsonify({
            "PAR_media": float(np.mean(PAR_media)),
            "images": images
        })

    except Exception:
        print("🔥 ERROR EN simulate_prod")
        traceback.print_exc()
        return jsonify({"error": "Error interno en producción"}), 500
