"""
explorar.py — Generación de figuras nuevas por muestreo aleatorio estructurado.

Sortear los 49 coeficientes de forma independiente produce un campo aleatorio isótropo
dominado por los grados altos (el grado n aporta 2n + 1 términos): todas las muestras
son irregulares y de aspecto similar, y la probabilidad de obtener una simetría
reconocible es despreciable. Aquí se sortean en cambio unos pocos parámetros con
significado geométrico, organizados en familias:

    zonal           solo C_n0: cuerpos de revolución
    simetria_m      solo m múltiplo de k: invariancia bajo giros de 360°/k
    espiral         m lóbulos con fase dependiente de la latitud
    bultos          2-5 máximos o mínimos gaussianos en posiciones aleatorias
    dispersa        2-3 armónicos aislados, para estudiar su interferencia
    espectral       varianza por grado ∝ n^(-β), espectro decreciente como el del geoide
    aleatorio_puro  referencia: 49 coeficientes U(-10, 10) independientes

Uso:
    python explorar.py                                   # todas las familias, 8 candidatos
    python explorar.py --familia espiral --n 12 --semilla 3 --guardar-dat

Salida: resultados/exploracion/<familia>_s<semilla>.png (y los .dat con --guardar-dat).
"""

import argparse

import numpy as np

import armonicos as A
from piezas import espiral as pieza_espiral

TAU = 2 * np.pi
DIR_EXPLORACION = A.RESULTADOS / "exploracion"


def _distancia(th, ph, t0, p0):
    c = np.sin(th) * np.sin(t0) * np.cos(ph - p0) + np.cos(th) * np.cos(t0)
    return np.arccos(np.clip(c, -1, 1))


def _vacios():
    return np.zeros((7, 7)), np.zeros((7, 7))


def zonal(rng):
    C, S = _vacios()
    n = np.arange(1, 7)
    C[1:, 0] = rng.normal(size=6) / n ** rng.uniform(0, 1.5)
    return "zonal", C, S


def simetria_m(rng):
    k = int(rng.integers(2, 6))
    C, S = _vacios()
    for n in range(7):
        for m in range(0, n + 1, k):
            if rng.random() < 0.6:
                C[n, m] = rng.normal()
    return f"simetría k={k}", C, S


def espiral(rng):
    m = int(rng.integers(2, 6))
    k = rng.uniform(0.5, 1.4) * rng.choice([-1, 1])
    C, S = A.proyectar(pieza_espiral(m, k, lat=rng.uniform(-20, 20), ancho=rng.uniform(30, 50)))
    return f"espiral m={m} k={k:.2f}", C, S


def bultos(rng):
    nb = int(rng.integers(2, 6))
    params = [(np.arccos(rng.uniform(-1, 1)), rng.uniform(0, TAU),          # punto uniforme en la esfera
               np.radians(rng.uniform(15, 35)), rng.choice([-1, 1]) * rng.uniform(0.4, 1))
              for _ in range(nb)]

    def f(th, ph):
        return sum(a * np.exp(-_distancia(th, ph, t0, p0) ** 2 / (2 * s**2)) for t0, p0, s, a in params)

    C, S = A.proyectar(f)
    return f"bultos ×{nb}", C, S


def dispersa(rng):
    C, S = _vacios()
    usados = []
    for _ in range(int(rng.integers(2, 4))):
        n = int(rng.integers(1, 7))
        m = int(rng.integers(0, n + 1))
        C[n, m] = rng.normal()
        S[n, m] = rng.normal() if m else 0.0
        usados.append(f"Y{n}{m}")
    return " + ".join(usados), C, S


def espectral(rng):
    beta = rng.uniform(2, 4)
    C, S = _vacios()
    for n in range(1, 7):
        sigma = n ** (-beta / 2)
        C[n, : n + 1] = rng.normal(scale=sigma, size=n + 1)
        S[n, 1 : n + 1] = rng.normal(scale=sigma, size=n)
    return f"espectral β={beta:.1f}", C, S


def aleatorio_puro(rng):
    C = np.tril(rng.uniform(-10, 10, (7, 7)))
    S = np.tril(rng.uniform(-10, 10, (7, 7)))
    S[:, 0] = 0
    return "aleatorio puro", C, S


FAMILIAS = {f.__name__: f for f in (zonal, simetria_m, espiral, bultos, dispersa, espectral, aleatorio_puro)}


def hoja(candidatos, archivo, vistas=((30, -60), (90, -90)), columnas=4, dpi=80):
    """Hoja de miniaturas: por candidato, una vista oblicua y una polar."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    nv = len(vistas)
    filas = int(np.ceil(len(candidatos) / columnas))
    fig = plt.figure(figsize=(2.3 * nv * columnas, 2.5 * filas))
    for i, (etiqueta, C, S) in enumerate(candidatos):
        r, c = divmod(i, columnas)
        for j, (elev, azim) in enumerate(vistas):
            ax = fig.add_subplot(filas, nv * columnas, r * nv * columnas + c * nv + j + 1, projection="3d")
            A.dibujar(ax, C, S, elev, azim)
            if j == 0:
                ax.text2D(0.02, 0.97, f"[{i}] {etiqueta}", transform=ax.transAxes, fontsize=7, va="top")
    fig.subplots_adjust(0, 0, 1, 1, 0, 0.05)
    archivo.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(archivo, dpi=dpi)
    plt.close(fig)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--familia", default="todas", choices=["todas", *FAMILIAS])
    p.add_argument("--n", type=int, default=8, help="candidatos por familia")
    p.add_argument("--semilla", type=int, default=0)
    p.add_argument("--q", type=float, default=1.5, help="relieve: mínimo en radio 1 - 0.1·q")
    p.add_argument("--guardar-dat", action="store_true", help="guardar el .dat de cada candidato")
    a = p.parse_args(argv)

    for familia in (FAMILIAS if a.familia == "todas" else [a.familia]):
        rng = np.random.default_rng(a.semilla)
        candidatos = []
        for k in range(a.n):
            etiqueta, C, S = FAMILIAS[familia](rng)
            if familia != "aleatorio_puro":  # la referencia se deja sin modificar
                C, S = A.fijar_relieve(C, S, a.q)
            C, S, _ = A.limpiar(C, S, umbral=0.005)
            if a.guardar_dat:
                A.escribir_dat(DIR_EXPLORACION / f"{familia}_s{a.semilla}_{k}.dat", C, S)
            candidatos.append((etiqueta, C, S))
        salida = DIR_EXPLORACION / f"{familia}_s{a.semilla}.png"
        hoja(candidatos, salida)
        print("->", salida.relative_to(A.RAIZ))


if __name__ == "__main__":
    main()
