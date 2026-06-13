import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Polygon, LineString
from matplotlib.animation import FuncAnimation

# --- Cálculo con Shapely ---

def calcular_Lt_traps(w, t, h, d, Ld, alpha_c, beta_c, x0, y0):
        vecinos = 2
        theta = np.pi-alpha_c
        P = w+d     # periodo horizontal

        dx = np.sin(theta)
        dy = -np.cos(theta)
        
        # segmento largo de rayo
        if dy == 0:
            Tmax = 10
        else:
            Tmax = (h - y0) / dy   # dy negativo → Tmax>0
        
        ray = LineString([(x0, y0), (x0 + dx*Tmax, y0 + dy*Tmax)])
        
        if abs(alpha_c) < 1e-6:
            alpha_c = np.sign(alpha_c) * 1e-6  # aseguras valor mínimo
        cos_beta_c = np.cos(beta_c)
        cos_beta_c = max(abs(cos_beta_c), 1e-6) * np.sign(cos_beta_c)
        
        # -------------------------------
        # Construir trapezoides periódicos
        # -------------------------------
        def trapezoid_at(k):
            """Devuelve el trapezoide número k (entero)."""
            x_shift = k * P
            return Polygon([
                (x_shift - w/2, 0),
                (x_shift + w/2, 0),
                (x_shift + t/2, h),
                (x_shift - t/2, h)
            ])
        
        K = range(-vecinos, vecinos+1)
        traps = [trapezoid_at(k) for k in K]
        # -------------------------------
        # Cálculo de distancia dentro de trapezoides
        # -------------------------------
        total_inside = 0
        segments_inside = []
        
        for tr in traps:
            inter = ray.intersection(tr)
            if inter.is_empty:
                continue
            if isinstance(inter, LineString):
                seg_len = inter.length
                total_inside += seg_len
                segments_inside.append(inter)
            elif hasattr(inter, "geoms"):
                for part in inter.geoms:
                    total_inside += part.length
                    segments_inside.append(part)

        Pr = total_inside/cos_beta_c
        Lt = Pr*Ld
        
        #print("Distancia total dentro de los trapezoides (proyectada en el plano transversal):", total_inside)
        #print("Distancia total dentro del seto:", Pr)

        return {"Lt": Lt , "distancia_atravesada": Pr, "trapezoides": traps, "rayo": ray, "segmentos": segments_inside}

def init_plot(ax, datos_traps):
    traps = datos_traps["trapezoides"]
    ray = datos_traps["rayo"]
    segments_inside = datos_traps["segmentos"]

    trap_patches = []
    for tr in traps:
        xs, ys = tr.exterior.xy
        patch = ax.fill(xs, ys, alpha=0.3, color="gray")[0]
        trap_patches.append(patch)

    # rayo
    rx, ry = ray.xy
    ray_line, = ax.plot(rx, ry, '--', linewidth=2)

    # segmentos dentro
    seg_lines = []
    for seg in segments_inside:
        sx, sy = seg.xy
        line, = ax.plot(sx, sy, linewidth=3)
        seg_lines.append(line)

    ax.set_aspect('equal')
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.grid(True)

    return trap_patches, ray_line, seg_lines

def update_plot(alpha_c_deg, ax, artists, params):
    w, t, h, d, beta_c, x0, y0 = params
    trap_patches, ray_line, seg_lines = artists

    datos_traps = calcular_Lt_traps(
        w, t, h, d, 2/h,
        np.radians(alpha_c_deg),
        np.radians(beta_c),
        x0, y0
    )

    # actualizar rayo
    rx, ry = datos_traps["rayo"].xy
    ray_line.set_data(rx, ry)

    # borrar segmentos anteriores
    for line in seg_lines:
        line.remove()
    seg_lines.clear()

    # dibujar nuevos segmentos
    for seg in datos_traps["segmentos"]:
        sx, sy = seg.xy
        line, = ax.plot(sx, sy, linewidth=3)
        seg_lines.append(line)

    Pr = datos_traps["distancia_atravesada"]
    ax.set_title(f"α_c = {alpha_c_deg}°, L = {Pr:.2f} m")

    return ray_line, *seg_lines




if __name__ == "__main__":

    # -------------------------------
    # Parámetros del problema
    # -------------------------------
    w = 2   # base inferior
    t = 0.5 # base superior
    h = 3       # altura
    d = 2  # separacion entre setos
    #alpha_c = 60  # -90 < alpha_c < 90
    
    beta_c = 0
    x0, y0 = 0, 0
    
    fig, ax = plt.subplots(figsize=(8, 6))
    
    xlim = 3*(d+w)
    ax.set_xlim(-xlim, xlim)
    ax.set_ylim(-0.2, h + 0.2)
    
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")

    # frame inicial
    alpha_init = -80
    datos_init = calcular_Lt_traps(
        w, t, h, d, 2/h,
        np.radians(alpha_init),
        np.radians(beta_c),
        x0, y0
    )
    
    artists = init_plot(ax, datos_init)
    
    params = (w, t, h, d, beta_c, x0, y0)
    alphas = list(range(-80, 90, 10))
    
    ani = FuncAnimation(
        fig,
        update_plot,
        frames=alphas,
        fargs=(ax, artists, params),
        interval=500,   # ms entre frames
        blit=False
    )
    
    plt.show()
         
         
    
    
    

    