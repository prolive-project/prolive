from flask import Blueprint, render_template, request, jsonify
import io
import base64
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pvlib
from flask_login import login_required

# 🔹 Definir el Blueprint correctamente
rca_hle_bp = Blueprint('rca_hle', __name__)

# 🔹 Historial de consultas (últimas 5 consultas)
query_history = []

# Función para calcular ángulos solares
def calculate_solar_angles(latitude, longitude):
    times = pd.date_range(start='2024-01-01', end='2024-12-31 23:00', freq='h', tz='Europe/Madrid')
    solar_position = pvlib.solarposition.get_solarposition(times, latitude, longitude)
    return solar_position[['azimuth', 'elevation']].reset_index()

# Función para calcular RCA
def calculate_rca(h, d, w):
    return (2 * h / np.cos(np.radians(0)) + w - 2 * h * np.tan(np.radians(0))) / (d + w)

# Función para combinar HLE y RCA en una métrica optimizada
def objective_function(params, solar_data, threshold):
    h, d, w, eta, psi = params
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
    rca = calculate_rca(h, d, w)
    return hours_above_threshold * rca

# Generar una gráfica basada en consultas previas
def generate_hle_rca_plot(h, d, w, eta, psi, lat, lon):
    global query_history
    solar_data = calculate_solar_angles(lat, lon)
    hle_rca_value = objective_function([h, d, w, eta, psi], solar_data, 0.50)

    # Guardar consulta en historial
    query_history.append({
        "h": h, "d": d, "w": w, "eta": eta, "psi": psi, "lat": lat, "lon": lon, "hle_rca": hle_rca_value
    })
    if len(query_history) > 5:
        query_history.pop(0)

    # Generar gráfico con historial de consultas
    plt.figure(figsize=(8, 5))
    indices = list(range(1, len(query_history) + 1))
    values = [q['hle_rca'] for q in query_history]
    labels = [f"Consulta {i}" for i in indices]
    bars = plt.bar(indices, values, tick_label=labels, color='skyblue')
    plt.ylabel("Valor de HLE * RCA")
    plt.title("Historial de Optimización de HLE * RCA")

    for bar, query in zip(bars, query_history):
        plt.text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
                 f"h: {query['h']}, d: {query['d']}, w: {query['w']}, eta: {query['eta']}, psi: {query['psi']}",
                 ha='center', va='bottom', fontsize=8, rotation=45, color='black', bbox=dict(facecolor='white', alpha=0.7))

    img = io.BytesIO()
    plt.savefig(img, format='png')
    img.seek(0)
    img_base64 = base64.b64encode(img.getvalue()).decode('utf8')
    plt.close()
    return img_base64

# 🔹 Ruta para la vista de RCA * HLE
@rca_hle_bp.route('/rca_hle')
@login_required
def rca_hle():
    return render_template('rca_hle.html')

# 🔹 Ruta para actualizar la gráfica de RCA * HLE
@rca_hle_bp.route('/update_hle_rca', methods=['POST'])
def update_hle_rca():
    try:
        data = request.get_json()
        h = float(data['h'])
        d = float(data['d'])
        w = float(data['w'])
        psi = float(data['psi'])
        eta = float(data['eta'])
        lat = float(data['lat'])
        lon = float(data['lon'])

        img_base64 = generate_hle_rca_plot(h, d, w, eta, psi, lat, lon)
        return jsonify({'image': img_base64})

    except Exception as e:
        return jsonify({"error": str(e)}), 500
