"""
formas_nuevas.py — Figuras propuestas, construidas con piezas.py.

Cada forma ilustra una simetría distinta y, por tanto, un patrón distinto de
coeficientes no nulos:

    dado        8 máximos en los vértices de un cubo   -> simetría cúbica: solo C_40, C_44, C_60, C_64
    estrella    5 lóbulos ecuatoriales                  -> solo m múltiplo de 5 (y zonales)
    cacahuete   dos polos y una cintura ecuatorial      -> solo m = 0 y n par
    corazon     dos lóbulos, punta y hendidura          -> solo simetría φ -> -φ: todos los S_nm = 0

Uso:
    python formas_nuevas.py            # todas
    python formas_nuevas.py dado       # una sola
"""

import sys

import numpy as np

from piezas import banda, bulto, construir, hoyo, lobulos

LAT_VERTICE_CUBO = np.degrees(np.arctan(1 / np.sqrt(2)))  # 35.26°

FORMAS = {
    "dado": (
        sum(bulto(s * LAT_VERTICE_CUBO, 45 + 90 * k, ancho=30) for s in (1, -1) for k in range(4)),
        2.0,
        [(25, -30, "Oblicua"), (90, -90, "Vista polar"), (0, -45, "Vista ecuatorial")],
    ),
    "estrella": (
        lobulos(5, lat=0, ancho=25) + 0.6 * banda(0, ancho=25),
        2.0,
        [(35, -60, "Oblicua"), (90, -90, "Vista polar"), (0, -90, "Vista ecuatorial")],
    ),
    "cacahuete": (
        bulto(90, 0, ancho=30) + bulto(-90, 0, ancho=30) - 1.2 * banda(0, ancho=18),
        3.0,
        [(30, -60, "Oblicua"), (90, -90, "Vista polar"), (0, -90, "Vista ecuatorial")],
    ),
    "corazon": (
        bulto(25, -45, ancho=30) + bulto(25, 45, ancho=30)   # lóbulos superiores
        + 1.2 * bulto(-55, 0, ancho=25, punta=True)          # punta inferior
        + 0.8 * hoyo(55, 0, ancho=20)                        # hendidura entre lóbulos
        + 0.5 * hoyo(0, 180, ancho=60),                      # cara posterior más plana
        1.5,
        [(0, 0, "Frontal"), (90, -90, "Vista polar"), (0, -90, "Lateral")],
    ),
}

if __name__ == "__main__":
    for nombre in sys.argv[1:] or FORMAS:
        forma, q, vistas = FORMAS[nombre]
        construir(nombre, forma, q, vistas)
