from flask import Blueprint, render_template, request, jsonify
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import io
import base64
from flask_login import login_required

fr_bp = Blueprint('fr', __name__)

# Historial de las últimas 5 consultas
query_history = []

def calcular_angulos_solares(latitud, dia_ano, hora):
    declinacion = 23.45 * np.sin(np.radians((360/365) * (dia_ano - 81)))
    hora_angular = (hora - 12) * 15
    cos_zenith = (np.sin(np.radians(latitud)) * np.sin(np.radians(declinacion)) +
                  np.cos(np.radians(latitud)) * np.cos(np.radians(declinacion)) * np.cos(np.radians(hora_angular)))
    zenith = np.degrees(np.arccos(np.clip(cos_zenith, -1, 1)))
    sin_azimuth = (np.cos(np.radians(declinacion)) * np.sin(np.radians(hora_angular))) / np.cos(np.radians(zenith))
    azimuth = np.degrees(np.arcsin(np.clip(sin_azimuth, -1, 1)))
    return zenith, azimuth

def calcular_sombra(H, Wx, Wy, fR, latitud, dia_ano, hora):
    zenith, azimuth = calcular_angulos_solares(latitud, dia_ano, hora)
    fR_mod = fR % 180  
    fS_fR = ((azimuth - fR_mod + 180) % 360) - 180
    Lx = H * np.tan(np.radians(zenith)) * np.sin(np.radians(fS_fR))
    Ly = H * np.tan(np.radians(zenith)) * np.cos(np.radians(fS_fR))
    return max(Lx, 0), max(Ly, 0)

def fraccion_interceptada(H, Wx, Wy, Ex, Ey, fR, latitud, dia_ano, hora, Cp=0.3, alpha=0.85):
    Lx, Ly = calcular_sombra(H, Wx, Wy, fR, latitud, dia_ano, hora)
    fs_h = ((Lx + Wx) * Wy + (Ly + Wy) * Wx - Wx * Wy) / (Ex * Ey)
    fs_h = min(max(fs_h, 0), 1)
    Cp_eff = np.exp(np.log(Cp) * np.sqrt(alpha))
    fb_h = fs_h * (1 - Cp_eff)
    tb_h = 1 - fb_h
    return fb_h, tb_h

def simular_anual(H, Wx, Wy, Ex, Ey, fR, latitud, Cp=0.3, alpha=0.85):
    dias = np.arange(1, 366)
    horas = np.linspace(6, 18, 13)
    fracciones_par = np.zeros(len(dias))

    for i, dia in enumerate(dias):
        fb_h_vals = []
        tb_h_vals = []
        for h in horas:
            fb_h, tb_h = fraccion_interceptada(H, Wx, Wy, Ex, Ey, fR, latitud, dia, h, Cp, alpha)
            fb_h_vals.append(fb_h)
            tb_h_vals.append(tb_h)
        td_D = 2 * np.sum(np.array(tb_h_vals) * np.cos(np.radians(horas)) * np.sin(np.radians(horas))) / len(horas)
        fd_D = 1 - td_D
        f_IPAR_h = np.array(fb_h_vals) * 0.5 + fd_D * 0.5
        f_DIPAR = np.mean(f_IPAR_h)
        fracciones_par[i] = f_DIPAR
    return dias, fracciones_par

def generate_fr_plot(H, Wx, Wy, Ex, Ey, fR, latitud):
    global query_history
    dias, fracciones = simular_anual(H, Wx, Wy, Ex, Ey, fR, latitud)
    query_data = {"H": H, "Wx": Wx, "Wy": Wy, "Ex": Ex, "Ey": Ey, "fR": fR, "latitud": latitud, "fracciones": fracciones}
    query_history.append(query_data)
    if len(query_history) > 5:
        query_history.pop(0)
    plt.figure(figsize=(10, 5))
    colores = ['b', 'g', 'r', 'c', 'm']
    for i, q in enumerate(query_history):
        plt.plot(dias, q['fracciones'], label=f"Consulta {i+1}: H={q['H']}, Wx={q['Wx']}, Wy={q['Wy']}, Ex={q['Ex']}, Ey={q['Ey']}, fR={q['fR']}", color=colores[i])
    plt.xlabel("Día del año")
    plt.ylabel("Fracción de radiación interceptada")
    plt.title("Evolución de la intercepción de luz con cambios de parámetros")
    plt.ylim(0, 1)
    plt.legend()
    plt.grid()
    img = io.BytesIO()
    plt.savefig(img, format='png')
    img.seek(0)
    img_base64 = base64.b64encode(img.getvalue()).decode('utf8')
    plt.close()
    return img_base64

@fr_bp.route('/FR')
@login_required
def fr():
    return render_template('fr.html')

@fr_bp.route('/update_fr', methods=['POST'])
def update_fr():
    data = request.get_json()
    H = float(data['H'])
    Wx = float(data['Wx'])
    Wy = float(data['Wy'])
    Ex = float(data['Ex'])
    Ey = float(data['Ey'])
    fR = float(data['fR'])
    latitud = float(data['latitud'])
    img_base64 = generate_fr_plot(H, Wx, Wy, Ex, Ey, fR, latitud)
    return jsonify({'image': img_base64})
