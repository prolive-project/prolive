from flask import Blueprint, render_template, request, jsonify
import pandas as pd
import pvlib
import numpy as np
import matplotlib.pyplot as plt
import io
import base64
from datetime import datetime
from flask_login import login_required

hle_bp = Blueprint('hle', __name__)

# 🔹 Historial de las últimas 5 consultas
hle_history = []

def calculate_solar_angles(latitude, longitude):
    times = pd.date_range(start='2024-01-01', end='2024-12-31 23:00', freq='h', tz='Europe/Madrid')
    solar_position = pvlib.solarposition.get_solarposition(times, latitude, longitude)
    results = solar_position[['azimuth', 'elevation']]
    return results.reset_index()

def calculate_critical_angle(solar_data, h, psi, d, eta):
    phi = solar_data['azimuth']
    theta_c = np.degrees(np.arctan(h * np.cos(np.radians(phi - eta)) / (d + h * np.tan(np.radians(psi)))))
    return theta_c

def optimize_light_hours(solar_data, h, psi, d, eta, threshold=0.50):
    hours_above_threshold = 0
    for _, row in solar_data.iterrows():
        phi = row['azimuth']
        theta = row['elevation']
        theta_c = np.degrees(np.arctan(h * np.cos(np.radians(phi - eta)) / (d + h * np.tan(np.radians(psi)))))
        if theta > theta_c:
            Idr = 1000 * (1 - np.exp(-0.06 * theta))
            Idf = 10 * (0.55 + 9.6 * (1 - np.exp(-0.05 * theta)))
            R_direct = Idr * (np.sin(np.radians(theta)) * np.sin(np.radians(psi)) +
                              np.cos(np.radians(theta)) * np.cos(np.radians(psi)) * np.cos(np.radians(phi - eta)))
            R_diffuse = Idf * (90 - theta_c + psi) / 180
            R = (R_direct + R_diffuse) / (Idr + Idf)
            if R >= threshold:
                hours_above_threshold += 1
    return hours_above_threshold

def generate_hle_plot(h, psi, d, eta, lat, lon):
    solar_data = calculate_solar_angles(lat, lon)
    effective_hours = optimize_light_hours(solar_data, h, psi, d, eta)

    # 🔹 Guardar consulta en historial (máximo 5)
    hle_history.append({'lat': lat, 'lon': lon, 'h': h, 'psi': psi, 'd': d, 'eta': eta, 'hle_value': effective_hours})
    if len(hle_history) > 5:
        hle_history.pop(0)

    # 🔹 Crear gráfico de historial de consultas
    plt.figure(figsize=(8, 5))
    indices = list(range(1, len(hle_history) + 1))
    valores = [q['hle_value'] for q in hle_history]
    etiquetas = [f"Consulta {i}" for i in indices]

    plt.bar(indices, valores, tick_label=etiquetas, color='skyblue')
    plt.ylabel("HLE (horas)")
    plt.title("Historial de Consultas HLE")

    # Agregar anotaciones sobre cada barra
    for i, valor in enumerate(valores):
        plt.text(indices[i], valor, f"{valor} h", ha='center', va='bottom')

    img = io.BytesIO()
    plt.savefig(img, format='png')
    img.seek(0)
    img_base64 = base64.b64encode(img.getvalue()).decode('utf8')
    plt.close()

    return img_base64, effective_hours

@hle_bp.route('/HLE')
@login_required
def hle():
    return render_template('hle.html')

@hle_bp.route('/update_hle', methods=['POST'])
def update_hle():
    try:
        data = request.get_json()
        if not data:
            return jsonify({"error": "No se recibieron datos"}), 400

        lat = float(data.get('lat', 40))  
        lon = float(data.get('lon', -3))
        h = float(data.get('h', 1.5))
        psi = float(data.get('psi', 10))
        d = float(data.get('d', 3.0))
        eta = float(data.get('eta', 90))

        img_base64, hle_value = generate_hle_plot(h, psi, d, eta, lat, lon)

        return jsonify({'hle_value': hle_value, 'image': img_base64, 'history': hle_history})

    except Exception as e:
        return jsonify({"error": str(e)}), 500



