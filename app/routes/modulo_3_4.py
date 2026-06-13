import numpy as np

# --- PRODUCCIÓN ---

def densidad_frutos_connor(PAR):
    """
    Relación bilineal simplificada tipo Connor 2016
    PAR en mol PAR m-2 d-1
    """
    PAR_sat = 24  # frutos m-2 por mol PAR
    a = 59
    b = -35
    dens = a * PAR +b 
    dens = np.minimum(dens, a * PAR_sat +b )

    return dens

def peso_fruto_seco_connor(PAR):
    """
    Relación lineal Connor 2016
    """
    return 0.0051 * PAR + 0.39

def concentracion_aceite_connor(PAR):
    """
    Relación lineal Connor 2016
    """
    return 0.25 * PAR + 34.8

def calcular_produccion_seto(matriz_PAR, arquitectura, longitud_seto=100.0):
    """
    Calcula número total de frutos y peso total producido
    siguiendo Connor et al. (2016)
    """

    h = arquitectura["h"]
    w = arquitectura["w"]  # = t (rectangular)
    d = arquitectura["d"]
    nz, nx = matriz_PAR.shape

    dz = h / nz
    dx = w / nx

    # Área de pared por capa y por metro de fila
    area_capa = dz * 1.0  # m2 (1 m de longitud de fila)

    frutos_por_metro = 0.0
    peso_seco_por_metro = 0.0
    peso_aceite_por_metro = 0.0

    # Definición de lados
    mitad_x = nx // 2
    lados = {
        "izquierdo": range(0, mitad_x),
        "derecho": range(mitad_x, nx)
    }
    matriz_PAR_capas = np.zeros((nz, nx))
    for lado, idx_x in lados.items():
        for iz in range(nz):
            # PAR media de la capa y lado
            PAR_capa = matriz_PAR[iz, idx_x].mean()
            matriz_PAR_capas[iz, idx_x] = PAR_capa
            # Densidad de frutos (frutos m-2 de seto)
            dens_f = densidad_frutos_connor(PAR_capa)
            concen_aceite = concentracion_aceite_connor(PAR_capa)
            # Número de frutos por metro de fila
            frutos_capa = dens_f * area_capa

            # Peso seco medio
            w_seco = peso_fruto_seco_connor(PAR_capa)
            
            # Peso seco total capa
            peso_seco_capa = frutos_capa * w_seco

            aceite_capa = peso_seco_capa*concen_aceite/100

            frutos_por_metro += frutos_capa
            peso_seco_por_metro += peso_seco_capa
            peso_aceite_por_metro += aceite_capa

    # Escalado a longitud real del seto
    frutos_totales = frutos_por_metro * longitud_seto
    peso_seco_total = peso_seco_por_metro * longitud_seto
    aceite_total = peso_aceite_por_metro * longitud_seto
    
    matriz_wfruto_capas = peso_fruto_seco_connor(matriz_PAR_capas)
    matriz_densidad_capas = densidad_frutos_connor(matriz_PAR_capas)
    matriz_peso_capas = matriz_densidad_capas*matriz_wfruto_capas*area_capa
    matriz_concen_aceite_capas = concentracion_aceite_connor(matriz_PAR_capas)
    matriz_aceite_capas = matriz_concen_aceite_capas*matriz_peso_capas/100

    return {
        "frutos_por_metro": frutos_por_metro,
        "peso_seco_por_metro_kg": peso_seco_por_metro / 1000,
        "aceite_por_metro_kg": peso_aceite_por_metro/1000,
        "frutos_totales": frutos_totales,
        "peso_seco_total_kg": peso_seco_total / 1000,
        "aceite_total_kg": aceite_total/ 1000,
        "frutos_hectarea": frutos_por_metro/(w+d)*1e4,
        "peso_seco_hectarea": peso_seco_por_metro/(w+d)*1e4,
        "aceite_hectarea": peso_aceite_por_metro/(w+d)*1e4,
        "matriz_PAR_capas": matriz_PAR_capas,
        "matriz_wfruto_capas": matriz_wfruto_capas,
        "matriz_densidad_capas": matriz_densidad_capas,
        "matriz_peso_capas": matriz_peso_capas,
        "matriz_concen_aceite_capas": matriz_concen_aceite_capas,
        "matriz_aceite_capas": matriz_aceite_capas
    }

# --- REPARTIR PRODUCCION ---

def mascara_trapecio(matriz_PAR, arquitectura):
    """
    Devuelve una máscara booleana con True en las celdas
    contenidas dentro del seto trapezoidal.
    """

    w = arquitectura["w"]   # base inferior
    t = arquitectura["t"]   # base superior
    h = arquitectura["h"]
    nz, nx = matriz_PAR.shape
    W = max(w, t)

    dx = W / nx
    dz = h / nz

    mask = np.zeros_like(matriz_PAR, dtype=bool)
    
    for iz in range(nz):
        z = h - (iz + 0.5) * dz
        # anchura del seto a esa altura
        w_z = w + (t - w) * z / h
        x_left = (W - w_z) / 2
        x_right = x_left + w_z

        for ix in range(nx):
            x = (ix + 0.5) * dx
            if x_left <= x <= x_right:
                mask[iz, ix] = True
    return mask


def normalizar_pesos(w, mask):
    w_eff = np.where(mask, w, 0.0)
    s = np.sum(w_eff)

    if s == 0:
        raise ValueError("La suma de pesos dentro del seto es cero")
    return w_eff / s

def reparto_lineal_por_PAR(matriz_PAR, valor_medio, arquitectura):
    """
    Reparte un valor medio según PAR de forma proporcional.
    Devuelve una matriz del mismo tamaño.
    """
    mask = mascara_trapecio(matriz_PAR, arquitectura)
    pesos = normalizar_pesos(matriz_PAR, mask)
    n_voxels_seto = mask.sum()

    return valor_medio * pesos * n_voxels_seto

def peso_bilineal_PAR(PAR, I_sat=28.0):
    """
    Peso relativo bilineal con saturación.
    """
    return np.minimum(PAR, I_sat)

def reparto_densidad_por_PAR(
    matriz_PAR,
    densidad_lineal_total,
    arquitectura,
    I_sat=28.0
):
    """
    Reparte frutos/m según PAR con saturación.
    Conserva exactamente la densidad total.
    Devuelve una matriz con la densidad lineal de fruto en cada voxel
    """
    mask = mascara_trapecio(matriz_PAR, arquitectura)
    pesos_brutos = peso_bilineal_PAR(matriz_PAR, I_sat)
    pesos = normalizar_pesos(pesos_brutos, mask)
    n_voxels_seto = mask.sum()
    densidad_voxel = densidad_lineal_total * pesos * n_voxels_seto
    return densidad_voxel

def transformar_PAR_a_fisiologia_voxs(
    matriz_PAR,
    arquitectura,
    peso_medio_fruto_g,
    aceite_medio_pct,
    densidad_lineal_frutos,
    I_sat=28.0
):
    """
    Devuelve:
    - peso de fruto por voxel (g/fruto)
    - concentración de aceite por voxel (% d.b.)
    - densidad de frutos por voxel (frutos/m)
    """

    peso_fruto_voxel = reparto_lineal_por_PAR(
        matriz_PAR,
        peso_medio_fruto_g,
        arquitectura
    )
    
    concen_aceite_voxel = reparto_lineal_por_PAR(
        matriz_PAR,
        aceite_medio_pct,
        arquitectura
    )

    densidad_voxel = reparto_densidad_por_PAR(
        matriz_PAR,
        densidad_lineal_frutos,
        arquitectura,
        I_sat
    )

    return {
        "matriz_wfruto_voxs": peso_fruto_voxel,
        "matriz_densidad_voxs": densidad_voxel,
        "matriz_concen_aceite_voxs": concen_aceite_voxel,
    }
        
    
