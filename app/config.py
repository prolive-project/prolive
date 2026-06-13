import os

BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

uri = os.getenv("DATABASE_URL")

if uri and uri.startswith("postgres://"):
    uri = uri.replace("postgres://", "postgresql://", 1)

SQLALCHEMY_DATABASE_URI = uri or 'sqlite:///' + os.path.join(BASE_DIR, 'instance', 'database.db')

SECRET_KEY = os.getenv("SECRET_KEY", "dev_key")
SQLALCHEMY_TRACK_MODIFICATIONS = False
