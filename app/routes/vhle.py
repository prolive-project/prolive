from flask import Blueprint, render_template

vhle_bp = Blueprint('vhle', __name__)

@vhle_bp.route('/VHLE')
def vhle():
    return render_template('vhle.html')