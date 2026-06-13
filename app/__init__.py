from flask import Flask
from flask_sqlalchemy import SQLAlchemy
from flask_login import LoginManager
from flask_migrate import Migrate
from app.config import Config
from app.models.user import User, db

def create_app():
    app = Flask(__name__)
    app.config.from_object(Config)

    # Inicializar la base de datos
    db.init_app(app)
    
    # Inicializar Flask-Migrate
    migrate = Migrate(app, db)

    # Configurar Flask-Login
    login_manager = LoginManager()
    login_manager.login_view = "auth.login"
    login_manager.init_app(app)

    @login_manager.user_loader
    def load_user(user_id):
        return User.query.get(int(user_id))

    # Importar y registrar rutas
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

    from app.routes.dash_reprca import create_dash_app # 🔹 Importar Dash

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

    # 🔹 Iniciar Dash dentro de Flask
    create_dash_app(app)
    return app
