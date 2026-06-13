from flask import Blueprint, render_template
from flask_login import login_required

home_bp = Blueprint('home', __name__)

@home_bp.route('/')
def landing():
    return render_template('landing.html')  # 🔹 Renderiza la nueva landing page

@home_bp.route('/home')
@login_required
def home():
    return render_template('home.html', title="Inicio")  # 🔹 Página protegida