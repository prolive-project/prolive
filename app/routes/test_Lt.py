import numpy as np

def longitud_rayo_en_trapecio(
    x0, y0, dx, dy,
    w, t, h,
    eps=1e-9
):
    """
    Longitud geométrica del rayo dentro de un trapecio centrado en x=0.
    dx, dy pueden ser arrays.
    """

    # Inicializar intervalo [s_min, s_max]
    s_min = np.zeros_like(dx)
    s_max = np.full_like(dx, np.inf)

    # ---------
    # y >= 0
    # ---------
    mask = np.abs(dy) > eps
    s1 = (0 - y0) / dy
    s2 = (h - y0) / dy
    s_low = np.minimum(s1, s2)
    s_high = np.maximum(s1, s2)

    s_min = np.where(mask, np.maximum(s_min, s_low), s_min)
    s_max = np.where(mask, np.minimum(s_max, s_high), s_max)

    # Pendiente de los lados
    m = (w - t) / (2 * h)

    # -----------------
    # Lado izquierdo
    # x >= -w/2 + m y
    # -----------------
    A = dx - m * dy
    B = (-w/2 + m * y0) - x0

    mask = np.abs(A) > eps
    s = B / A
    s_min = np.where(mask & (A > 0), np.maximum(s_min, s), s_min)
    s_max = np.where(mask & (A < 0), np.minimum(s_max, s), s_max)

    # -----------------
    # Lado derecho
    # x <= w/2 - m y
    # -----------------
    A = dx + m * dy
    B = (w/2 - m * y0) - x0

    mask = np.abs(A) > eps
    s = B / A
    s_min = np.where(mask & (A < 0), np.maximum(s_min, s), s_min)
    s_max = np.where(mask & (A > 0), np.minimum(s_max, s), s_max)

    # Longitud geométrica
    length = np.maximum(0.0, s_max - s_min)

    return length

def calcular_Pr_traps_vec(
    w, t, h, d,
    alpha_c, beta_c,
    x0, y0,
    vecinos=2,
    eps=1e-9
):
    """
    Calcula la distancia efectiva recorrida por un rayo al atravesar un conjunto
    periódico de setos de sección trapezoidal.
    
    Para cada dirección definida por (alpha_c, beta_c), se determina de forma
    analítica la longitud total del rayo contenida en el follaje, sumando las
    contribuciones de varios setos vecinos. Para generalizar fuera dl plano 
    transversal dividimos por cos(beta_c).
    
    La función es vectorizable y no utiliza geometría discreta.
    """


    alpha_c = np.asarray(alpha_c)
    beta_c  = np.asarray(beta_c)
    x0      = np.asarray(x0)
    y0      = np.asarray(y0)

    theta = np.pi - alpha_c
    dx = np.sin(theta)
    dy = -np.cos(theta)
    dy = np.where(np.abs(dy) < eps, -eps, dy)

    cos_beta = np.cos(beta_c)
    cos_beta = np.where(np.abs(cos_beta) < eps, eps, cos_beta)

    P = w + d

    out_shape = np.broadcast(x0, dx).shape
    total_length = np.zeros(out_shape)

    for k in range(-vecinos, vecinos + 1):
        x_shift = x0 - k * P

        Lk = longitud_rayo_en_trapecio(
            x_shift, y0,
            dx, dy,
            w, t, h,
            eps
        )

        total_length += Lk

    Pr = total_length / cos_beta
 

    return Pr

if __name__ == "__main__":
    w = 2   # base inferior
    t = 0.5 # base superior
    h = 3       # altura
    d = 2  # separacion entre setos
    #alpha_c = 60  # -90 < alpha_c < 90
    Ld = 2/h
    num_divisiones = 20
    angulos = np.linspace(-np.pi/2, np.pi/2, num_divisiones)
    alpha_c_grid, beta_c_grid = np.meshgrid(angulos, angulos)
    x0, y0 = 0, 0
    Lt = calcular_Pr_traps_vec(w, t, h, d, alpha_c_grid, beta_c_grid, x0, y0, vecinos=2)
    print(Lt)