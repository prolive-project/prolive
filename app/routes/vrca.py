from flask import Blueprint, render_template

vrca_bp = Blueprint('vrca', __name__)

@vrca_bp.route('/VRCA')
def vrca():
    return render_template('vrca.html')
