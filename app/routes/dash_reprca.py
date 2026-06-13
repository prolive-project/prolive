from dash import Dash, dcc, html, Input, Output
import plotly.graph_objects as go

# Parámetros iniciales
default_h = 1.5  # Altura del dosel en metros
default_d = 6  # Separación entre hileras
default_w = 4  # Ancho del dosel
spacing_within_row = 1.3  # Separación entre árboles en la misma fila
trunk_height = 1.5  # Altura del tronco en metros
plot_size_x = 100  # 🔹 Mantenemos X en el rango [0,100]
plot_size_y = 100  # 🔹 Mantenemos Y en el rango [0,100]

def create_dash_app(flask_app):
    """Crea la aplicación Dash integrada en Flask."""
    dash_app = Dash(
        server=flask_app,
        routes_pathname_prefix='/co2_dash/',
    )

    # Layout de la app Dash con inputs dinámicos
    dash_app.layout = html.Div([
        html.H3("Representación de hectárea"),
        
        html.Div([
            html.Label("Altura del dosel (h en metros):"),
            dcc.Input(id='h', type='number', value=default_h, step=0.1, min=0.1),

            html.Label("Espaciado entre hileras (d en metros):"),
            dcc.Input(id='d', type='number', value=default_d, step=0.1, min=0.1),

            html.Label("Ancho del dosel (w en metros):"),
            dcc.Input(id='w', type='number', value=default_w, step=0.1, min=0.1),
        ], style={'display': 'flex', 'gap': '10px', 'margin-bottom': '20px'}),

        dcc.Graph(id='co2-graph'),
    ])

    @dash_app.callback(
        Output('co2-graph', 'figure'),
        [Input('h', 'value'),
         Input('d', 'value'),
         Input('w', 'value')]
    )
    def update_graph(h, d, w):
        """Actualiza la representación gráfica en función de los valores ingresados."""
        if not h or not d or not w:
            return go.Figure()

        # 🔹 PASAMOS 'h' A generate_hedgerows()
        canopies, trunks, row_count, trees_per_row = generate_hedgerows(plot_size_x, plot_size_y, h, d, w, spacing_within_row)

        # 🔹 Cálculo de absorción de CO₂
        total_trees = row_count * trees_per_row

        fig = go.Figure()

        # Dibujar los troncos
        for trunk in trunks:
            x = [trunk[0][0], trunk[1][0]]
            y = [trunk[0][1], trunk[1][1]]
            z = [trunk[0][2], trunk[1][2]]
            fig.add_trace(go.Scatter3d(x=x, y=y, z=z, mode='lines', line=dict(color='brown', width=5), showlegend=False))

        # Dibujar los doseles (pirámides)
        for canopy in canopies:
            x = [v[0] for v in canopy] + [canopy[0][0]]  # Cerrar el polígono
            y = [v[1] for v in canopy] + [canopy[0][1]]
            z = [v[2] for v in canopy] + [canopy[0][2]]
            fig.add_trace(go.Mesh3d(x=x, y=y, z=z, color='green', opacity=0.5, showlegend=False))

        # Representar el suelo
        fig.add_trace(go.Surface(
            z=[[0, 0], [0, 0]],
            x=[[0, plot_size_x], [0, plot_size_x]],
            y=[[0, plot_size_y], [0, plot_size_y]],
            colorscale=[(0, 'burlywood'), (1, 'burlywood')],
            showscale=False
        ))

        # 🔹 Ajustamos la proporción de los ejes asegurando que Z se escale proporcionalmente
        max_z = h + trunk_height
        aspect_ratio_z = max_z / max(plot_size_x, plot_size_y)  # Relación correcta con X e Y

        fig.update_layout(
            title=f"Olivar en Seto - {row_count} Hileras x {trees_per_row} Árboles ({total_trees} Árboles en Total)",
            scene=dict(
                xaxis=dict(range=[0, plot_size_x], backgroundcolor="lightblue"),
                yaxis=dict(range=[0, plot_size_y], backgroundcolor="lightblue"),
                zaxis=dict(range=[0, max_z], backgroundcolor="burlywood"),
                aspectmode='manual',  # 🔹 Define la proporción entre los ejes
                aspectratio=dict(x=1, y=1, z=aspect_ratio_z)  # 🔹 Relación correcta
            )
        )
        return fig

    return dash_app

# Funciones auxiliares para la representación
def create_canopy(x, y, z=0, h=1.5, w=4):
    """Crea la geometría de un dosel de árbol en forma de pirámide de base cuadrada."""
    base_bottom_left = [x - w / 2, y - w / 2, z]
    base_bottom_right = [x + w / 2, y - w / 2, z]
    base_top_left = [x - w / 2, y + w / 2, z]
    base_top_right = [x + w / 2, y + w / 2, z]
    top = [x, y, z + h]
    return [base_bottom_left, base_bottom_right, base_top_right, base_top_left, top]

def create_trunk(x, y, z=0, height=trunk_height):
    """Crea la geometría del tronco del árbol como una línea vertical."""
    return [[x, y, z], [x, y, z + height]]

def generate_hedgerows(plot_size_x, plot_size_y, h, d, w, spacing_within_row):
    """Genera la disposición de los árboles en el seto."""
    canopies = []
    trunks = []
    row_count = plot_size_x // (w + d)  # Número máximo de filas considerando separación entre hileras (w + d)
    trees_per_row = int(plot_size_y // spacing_within_row)  # Número máximo de árboles por fila

    for row in range(row_count):
        for tree in range(trees_per_row):
            x = row * (w + d)  # Separación entre hileras considerando w + d
            y = tree * spacing_within_row  # Separación entre árboles dentro de una fila
            canopies.append(create_canopy(x, y, z=trunk_height, h=h, w=w))  # 🔹 PASAMOS 'h'
            trunks.append(create_trunk(x, y))

    return canopies, trunks, row_count, trees_per_row
