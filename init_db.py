from app import create_app, db
from app.user import User

app = create_app()

with app.app_context():
    db.create_all()

    if not User.query.filter_by(email="admin@test.com").first():
        user = User(
            name="admin",
            email="admin@test.com",
            password="admin123"
        )
        db.session.add(user)
        db.session.commit()

    print("DB inicializada")
