import numpy as np
import matplotlib
matplotlib.use('Agg')  # Evita problemas con la GUI en Flask
import matplotlib.pyplot as plt
from flask import Blueprint, render_template, request, jsonify
import io
import base64
import scipy.optimize as opt
from flask_login import login_required

rca_bp = Blueprint('rca', __name__)

# Función para calcular RCA
def calculate_rca(h, d, w, psi):
    rca = (2 * h / np.cos(np.radians(psi)) + w - 2 * h * np.tan(np.radians(psi))) / (d + w)
    return rca

# Función para generar la gráfica en base64
def generate_plot(h, d, w, psi):
    rca_value = calculate_rca(h, d, w, psi)

    img = io.BytesIO()
    plt.figure(figsize=(8, 5))
    plt.plot([0, 1], [rca_value, rca_value], lw=2, label=f"RCA = {rca_value:.2f}")
    plt.ylim(0, 3.5)
    plt.title("Relative Canopy Area (RCA)")
    plt.xlabel("X")
    plt.ylabel("Valor de RCA")
    plt.legend()
    plt.grid(True)
    plt.savefig(img, format='png', bbox_inches='tight')
    plt.close()  # Asegura que Matplotlib no mantenga abierto el gráfico

    img.seek(0)
    return base64.b64encode(img.read()).decode('utf-8')

# Ruta principal para RCA
@rca_bp.route('/RCA')
@login_required
def rca():
    return render_template('rca.html')

# Endpoint para actualizar la gráfica dinámicamente
@rca_bp.route('/update_rca', methods=['POST'])
def update_rca():
    try:
        data = request.get_json()  # 🔹 Asegurar que Flask interpreta JSON
        if not data:
            return jsonify({"error": "No se recibieron datos"}), 400

        h = float(data.get('h', 1.5))
        d = float(data.get('d', 3.0))
        w = float(data.get('w', 1.0))
        psi = float(data.get('psi', 10))

        rca_value = calculate_rca(h, d, w, psi)
        img_data = generate_plot(h, d, w, psi)

        return jsonify({
            'rca_value': rca_value,
            'image': img_data
        })

    except Exception as e:
        return jsonify({"error": str(e)}), 500


# Función para maximizar el RCA
@rca_bp.route('/maximize_rca', methods=['POST'])
def maximize_rca():
    bounds = [(0.75, 2.7), (1.5, 6.0), (0.4, 2.0), (0, 45)]
    result = opt.differential_evolution(lambda x: -calculate_rca(*x), bounds)
    best_h, best_d, best_w, best_psi = result.x

    # Calcular el nuevo RCA y generar la gráfica con los valores óptimos
    rca_value = calculate_rca(best_h, best_d, best_w, best_psi)
    img_data = generate_plot(best_h, best_d, best_w, best_psi)

    return jsonify({
        'h': best_h,
        'd': best_d,
        'w': best_w,
        'psi': best_psi,
        'rca_value': rca_value,
        'image': img_data
    })

