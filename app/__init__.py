from flask import Flask
from flask_sqlalchemy import SQLAlchemy
from flask_login import LoginManager
from flask_migrate import Migrate
from app.config import Config
from app.models.user import User, db

def create_app():
    app = Flask(__name__)
    app.config.from_object(Config)

    # Inicializar base de datos
    db.init_app(app)

    # Flask-Migrate
    migrate = Migrate(app, db)

    # Flask-Login
    login_manager = LoginManager()
    login_manager.login_view = "auth.login"
    login_manager.init_app(app)

    @login_manager.user_loader
    def load_user(user_id):
        return User.query.get(int(user_id))

    # Crear tablas + usuario admin (IMPORTANTE)
    with app.app_context():
        db.create_all()

        # Crear usuario admin si no existe
        admin = User.query.filter_by(email="admin@test.com").first()
        if not admin:
            admin = User(
                name="admin",
                email="admin@test.com",
                password="admin123"
            )
            db.session.add(admin)
            db.session.commit()

    # Blueprints
    from app.routes.home import home_bp
    from app.routes.rca import rca_bp
    from app.routes.hle import hle_bp
    from app.routes.hlerca import rca_hle_bp
    from app.routes.fr import fr_bp
    from app.routes.reprca import reprca_bp
    from app.routes.auth import auth_bp
    from app.routes.vrca import vrca_bp
    from app.routes.vhle import vhle_bp 
    from app.routes.par import par_bp
    from app.routes.fr_daily_rad import fr_daily_bp
    from app.routes.produccion import prod_bp

    from app.routes.dash_reprca import create_dash_app

    app.register_blueprint(home_bp)
    app.register_blueprint(rca_bp)
    app.register_blueprint(hle_bp)
    app.register_blueprint(rca_hle_bp)
    app.register_blueprint(fr_bp)
    app.register_blueprint(reprca_bp)
    app.register_blueprint(auth_bp)
    app.register_blueprint(vrca_bp)
    app.register_blueprint(vhle_bp)
    app.register_blueprint(par_bp)
    app.register_blueprint(fr_daily_bp)
    app.register_blueprint(prod_bp)

    # Dash
    create_dash_app(app)

    return app
