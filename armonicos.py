"""
armonicos.py — Núcleo de cálculo para la suma de armónicos esféricos.

Reproduce exactamente la convención del script de la práctica (original/spheriharm1.py):

    f(θ, φ) = Σ_{n=0}^{6} Σ_{m=0}^{n} K_nm · P_nm(cos θ) · [C_nm cos(mφ) + S_nm sin(mφ)]

    K_nm = sqrt( (2n+1)/(4π) · (n-m)!/(n+m)! )
    P_nm = scipy.special.lpmv(m, n, ·)        (incluye la fase de Condon-Shortley (-1)^m)
    ρ    = 0.1 · (10 + f / max f)             (radio de la superficie dibujada)

con θ la colatitud y φ la longitud. Sobre esa base añade el problema inverso
(obtener C y S a partir de una función objetivo) y utilidades de E/S y dibujo.

Funciones principales
    sintetizar(C, S)         coeficientes -> f en una malla
    proyectar(func)          función objetivo -> (C, S), mejor aproximación con n <= 6
    bulto_axial(a, θ0, φ0)   patrón axisimétrico en torno a un eje -> (C, S) (teorema de adición)
    fijar_relieve(C, S, q)   ajusta C_00 para fijar la profundidad de las zonas mínimas
    limpiar(C, S)            normaliza, anula términos despreciables y redondea
    leer_dat / escribir_dat  E/S en el formato de spheriharm1.dat
    dibujar / render_vistas  representación 3D con la misma geometría y colores que el script
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import scipy.special as ss

DMAX = 6

RAIZ = Path(__file__).resolve().parent
RESULTADOS = RAIZ / "resultados"
DIR_DATS = RESULTADOS / "dats"
DIR_FIGURAS = RESULTADOS / "figuras"


# ---------------------------------------------------------------------------
# Base de armónicos
# ---------------------------------------------------------------------------
def knorm(n: int, m: int) -> float:
    """Factor de normalización K_nm del script."""
    f = ss.factorial(n - m) / ss.factorial(n + m)
    return float(np.sqrt((2 * n + 1) / (4 * np.pi) * f))


def malla(nphi: int = 90, ntheta: int = 45):
    """Malla del script: φ en [0, 2π] y θ en [0, π]. Las matrices tienen forma (nphi, ntheta)."""
    vphi = np.linspace(0, 2 * np.pi, nphi)
    vtheta = np.linspace(0, np.pi, ntheta)
    theta, phi = np.meshgrid(vtheta, vphi)
    return vtheta, vphi, theta, phi


def sintetizar(C, S, vtheta=None, vphi=None, dmax: int = DMAX):
    """Evalúa f en la malla producto (vphi × vtheta). Por defecto usa la malla del script."""
    if vtheta is None or vphi is None:
        vtheta, vphi, _, _ = malla()
    f = np.zeros((len(vphi), len(vtheta)))
    x = np.cos(vtheta)
    for n in range(dmax + 1):
        for m in range(n + 1):
            if C[n, m] == 0 and S[n, m] == 0:
                continue
            parte_angular = C[n, m] * np.cos(m * vphi) + S[n, m] * np.sin(m * vphi)
            f += knorm(n, m) * np.outer(parte_angular, ss.lpmv(m, n, x))
    return f


def proyectar(func, dmax: int = DMAX, nquad: int = 64):
    """
    Análisis armónico de func(θ, φ): devuelve (C, S) en la convención del script.

    Por la ortogonalidad de la base, la proyección es la mejor aproximación por
    mínimos cuadrados sobre la esfera con grado <= dmax:

        C_nm = (2 - δ_m0) ∫ f · K_nm P_nm(cos θ) cos(mφ) dΩ
        S_nm =      2     ∫ f · K_nm P_nm(cos θ) sin(mφ) dΩ      (m > 0)

    Las integrales se calculan con cuadratura de Gauss-Legendre en x = cos θ y regla
    del trapecio periódica en φ; ambas son exactas para funciones de grado < 2·nquad.
    """
    x, w = np.polynomial.legendre.leggauss(nquad)
    th = np.arccos(x)
    nph = 2 * nquad
    ph = np.arange(nph) * 2 * np.pi / nph
    dphi = 2 * np.pi / nph
    TH, PH = np.meshgrid(th, ph, indexing="ij")
    F = np.asarray(func(TH, PH), dtype=float)

    C = np.zeros((dmax + 1, dmax + 1))
    S = np.zeros((dmax + 1, dmax + 1))
    for n in range(dmax + 1):
        for m in range(n + 1):
            Y = knorm(n, m) * ss.lpmv(m, n, x)
            norma2 = 1.0 if m == 0 else 0.5  # ||K_nm P_nm cos(mφ)||² con dΩ
            C[n, m] = np.einsum("i,i,ij,j->", w, Y, F, np.cos(m * ph)) * dphi / norma2
            if m > 0:
                S[n, m] = np.einsum("i,i,ij,j->", w, Y, F, np.sin(m * ph)) * dphi / norma2
    return C, S


def bulto_axial(a, theta0: float, phi0: float, dmax: int = DMAX):
    """
    Coeficientes de f = Σ_n a_n P_n(cos γ), con γ el ángulo al punto (θ0, φ0).

    Por el teorema de adición,
        P_n(cos γ) = Σ_m ε_m (n-m)!/(n+m)! P_nm(cos θ) P_nm(cos θ0) cos(m(φ - φ0)),
    con ε_0 = 1 y ε_m = 2, de donde
        C_nm = a_n ε_m (n-m)!/(n+m)! P_nm(cos θ0) cos(mφ0) / K_nm   (y S_nm con sin).
    """
    C = np.zeros((dmax + 1, dmax + 1))
    S = np.zeros((dmax + 1, dmax + 1))
    for n in range(min(dmax, len(a) - 1) + 1):
        for m in range(n + 1):
            eps = 1.0 if m == 0 else 2.0
            fac = ss.factorial(n - m) / ss.factorial(n + m)
            base = a[n] * eps * fac * ss.lpmv(m, n, np.cos(theta0)) / knorm(n, m)
            C[n, m] = base * np.cos(m * phi0)
            if m > 0:
                S[n, m] = base * np.sin(m * phi0)
    return C, S


# ---------------------------------------------------------------------------
# Relieve, normalización y limpieza de coeficientes
# ---------------------------------------------------------------------------
def fijar_relieve(C, S, q: float = 1.0):
    """
    Ajusta solo C_00 para que, con ρ = 0.1·(10 + f/max f), el radio vaya de 1 - 0.1·q
    (mínimo de f) a 1.1 (máximo de f).

    C_00 multiplica a Y_00 = 1/sqrt(4π), una constante: modificarlo desplaza f entera,
    lo que no altera los colores (escala min-max) pero sí el radio, porque max f está
    en el denominador. Restando c = (min f + q·max f)/(1 + q) se cumple la condición, y
    max f - c = (max f - min f)/(1 + q) > 0, de modo que la figura nunca se invierte.
    """
    C = C.copy()
    f = sintetizar(C, S)
    fmax, fmin = f.max(), f.min()
    c = (fmin + q * fmax) / (1 + q)
    C[0, 0] -= c / knorm(0, 0)
    return C, S


def radios_script(C, S):
    """Radio mínimo y máximo de la superficie tal como la dibuja spheriharm1.py."""
    f = sintetizar(C, S)
    rho = 0.1 * (10 + f / f.max())
    return float(rho.min()), float(rho.max())


def normalizar(C, S):
    """Escala los coeficientes a max|c| = 10. La figura no cambia: el script divide por max f."""
    s = max(np.abs(C).max(), np.abs(S).max())
    return C * 10 / s, S * 10 / s


def limpiar(C, S, umbral: float = 0.01, cifras: int = 3):
    """
    Normaliza a max|c| = 10, anula los coeficientes menores que umbral·10, redondea a
    'cifras' cifras significativas y fija S_n0 = 0 (sin(0·φ) = 0).

    Devuelve (C, S, error), con error = ||f_limpia - f|| / ||f|| en la malla del script.
    """
    C, S = normalizar(C, S)
    C2, S2 = C.copy(), S.copy()
    for M in (C2, S2):
        M[np.abs(M) < umbral * 10] = 0.0
        nz = M != 0
        M[nz] = np.array([float(f"{v:.{cifras}g}") for v in M[nz]])
    S2[:, 0] = 0.0
    f0, f1 = sintetizar(C, S), sintetizar(C2, S2)
    error = float(np.linalg.norm(f1 - f0) / np.linalg.norm(f0))
    return C2, S2, error


def coeficientes_no_nulos(C, S, dmax: int = DMAX):
    """Lista de (n, m, C_nm, S_nm) con algún coeficiente distinto de cero."""
    filas = []
    for n in range(dmax + 1):
        for m in range(n + 1):
            if C[n, m] != 0 or (m > 0 and S[n, m] != 0):
                filas.append((n, m, float(C[n, m]), float(S[n, m]) if m > 0 else 0.0))
    return filas


# ---------------------------------------------------------------------------
# Entrada / salida en el formato de spheriharm1.dat
# ---------------------------------------------------------------------------
def leer_dat(ruta, dmax: int = DMAX):
    """Lee un .dat con el mismo criterio que spheriharm1.py."""
    C = np.zeros((dmax + 1, dmax + 1))
    S = np.zeros((dmax + 1, dmax + 1))
    with open(ruta) as fichero:
        for nlinea, linea in enumerate(fichero):
            campos = linea.split()
            if nlinea > 1:
                for k in range(1, len(campos)):
                    C[nlinea - 2, k - 1] = float(campos[k].split("(")[1].split(",")[0])
                    S[nlinea - 2, k - 1] = float(campos[k].split(",")[1].split(")")[0])
    return C, S


def escribir_dat(ruta, C, S, dmax: int = DMAX, cifras: int = 5):
    """Escribe la tabla triangular de tuplas (C_nm,S_nm) en el formato de spheriharm1.dat."""

    def fmt(v):
        v = float(v)
        return "1e-6" if abs(v) < 1e-9 else f"{v:.{cifras}g}"  # 1e-6 es el cero del fichero original

    lineas = [
        "Please fill in the list with the C and S coefficients of your choice",
        "n\\m \t" + "\t".join(str(m) for m in range(dmax + 1)),
    ]
    for n in range(dmax + 1):
        lineas.append(f"{n}\t" + "\t".join(f"({fmt(C[n, m])},{fmt(S[n, m])})" for m in range(n + 1)))
    Path(ruta).parent.mkdir(parents=True, exist_ok=True)
    Path(ruta).write_text("\n".join(lineas) + "\n")


# ---------------------------------------------------------------------------
# Representación
# ---------------------------------------------------------------------------
def geometria(C, S, nphi: int = 120, ntheta: int = 60):
    """Superficie ρ(θ, φ) = 0.1·(10 + f/max f), idéntica a la del script."""
    vtheta, vphi, theta, phi = malla(nphi, ntheta)
    f = sintetizar(C, S, vtheta, vphi)
    rho = 0.1 * (10 + f / np.amax(f))
    x = rho * np.sin(theta) * np.cos(phi)
    y = rho * np.sin(theta) * np.sin(phi)
    z = rho * np.cos(theta)
    return x, y, z, f


def dibujar(ax, C, S, elev, azim, titulo=None, luz=True, ejes=False):
    """
    Dibuja la suma en un eje 3D con la geometría y la escala de colores del script.

    Diferencias de presentación respecto al script, que no cambian la figura:
      - aspecto 1:1:1 (matplotlib usa por defecto una caja 4:4:3 que comprime z un 25 %);
      - iluminación direccional (luz=True), para distinguir zonas hundidas de salientes.
    """
    import matplotlib.colors as mcolors
    from matplotlib import cm

    x, y, z, f = geometria(C, S)
    colores = cm.jet(mcolors.Normalize(vmin=f.min(), vmax=f.max())(f))
    opciones = dict(rstride=1, cstride=1, linewidth=0, antialiased=False, facecolors=colores)
    if luz:
        fuente = mcolors.LightSource(azdeg=azim + 60, altdeg=min(elev + 35, 89))
        opciones.update(shade=True, lightsource=fuente)
    else:
        opciones.update(shade=False)
    ax.plot_surface(x, y, z, **opciones)

    r = max(np.abs(x).max(), np.abs(y).max(), np.abs(z).max())
    ax.set_xlim(-r, r)
    ax.set_ylim(-r, r)
    ax.set_zlim(-r, r)
    ax.set_box_aspect((1, 1, 1), zoom=1.35)
    ax.view_init(elev, azim)
    if not ejes:
        ax.set_axis_off()
    if titulo:
        ax.set_title(titulo, fontsize=10)


VISTAS_POR_DEFECTO = [(30, -60, "Oblicua"), (90, -90, "Vista polar"), (0, -90, "Vista ecuatorial")]


def render_vistas(C, S, archivo, titulo="", vistas=VISTAS_POR_DEFECTO, dpi=110):
    """Guarda una figura con varias vistas (INC, AZ, etiqueta) de la misma suma."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(3.3 * len(vistas), 3.9))
    for i, (elev, azim, nombre) in enumerate(vistas):
        ax = fig.add_subplot(1, len(vistas), i + 1, projection="3d")
        dibujar(ax, C, S, elev, azim, f"{nombre}\nINC={elev}  AZ={azim}")
    if titulo:
        fig.suptitle(titulo, fontsize=13)
    fig.subplots_adjust(left=0, right=1, bottom=0, top=0.80 if titulo else 0.88, wspace=0)
    Path(archivo).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(archivo, dpi=dpi)
    plt.close(fig)


def imprimir_coeficientes(C, S):
    print(f"{'n':>2} {'m':>2} {'C_nm':>8} {'S_nm':>8}")
    for n, m, c, s in coeficientes_no_nulos(C, S):
        print(f"{n:>2} {m:>2} {c:>8.3g} {s:>8.3g}")
