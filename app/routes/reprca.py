from flask import Blueprint, render_template
from flask_login import login_required

reprca_bp = Blueprint('reprca', __name__)

@reprca_bp.route('/reprca')
@login_required  # 🔹 Bloquea acceso sin autenticación
def reprca():
    return render_template('reprca.html')