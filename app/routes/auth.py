from flask import Blueprint, render_template, redirect, url_for, request, flash, jsonify
from flask_login import login_user, logout_user, login_required, current_user
from werkzeug.security import generate_password_hash, check_password_hash
from app.models.user import User, db

auth_bp = Blueprint('auth', __name__)

@auth_bp.route('/login', methods=['POST'])
def login():
    email = request.form.get('email')
    password = request.form.get('password')
    
    # Buscar usuario en la base de datos
    user = User.query.filter_by(email=email).first()

    if not user:
        return jsonify({"success": False, "message": "Este correo no está registrado"}), 401

    if not check_password_hash(user.password, password):
        return jsonify({"success": False, "message": "Contraseña incorrecta"}), 401

    # Si todo está correcto, iniciar sesión
    login_user(user)
    return jsonify({"success": True, "redirect_url": url_for('home.home')})


@auth_bp.route('/register', methods=['POST'])
def register():
    email = request.form.get('email')
    name = request.form.get('name')
    password = request.form.get('password')
    hashed_password = generate_password_hash(password, method='pbkdf2:sha256')

    existing_user = User.query.filter_by(email=email).first()
    if existing_user:
        return jsonify({"success": False, "message": "El correo ya está registrado"}), 400

    new_user = User(email=email, name=name, password=hashed_password)
    db.session.add(new_user)
    db.session.commit()

    login_user(new_user)  # 🔹 Iniciar sesión automáticamente tras el registro
    return jsonify({"success": True, "redirect_url": url_for('home.home')})

@auth_bp.route('/users')
def get_users():
    users = User.query.all()
    users_list = [{"id": user.id, "name": user.name, "email": user.email} for user in users]
    return jsonify(users_list)

@auth_bp.route('/logout')
@login_required
def logout():
    logout_user()
    flash('Sesión cerrada correctamente', 'info')
    return redirect(url_for('landing'))  # 🔹 Redirigir a landing al cerrar sesión
