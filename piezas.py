"""
piezas.py — Construcción de funciones objetivo f(θ, φ) a partir de piezas simples.

Cada pieza es una función sobre la esfera que admite suma, resta y producto por un
escalar, de modo que una forma se describe como combinación lineal de piezas. Como la
proyección es lineal, los coeficientes de la suma son la suma de los coeficientes.

Posiciones en grados geográficos (lat ∈ [-90, 90], lon ∈ [0, 360)); anchuras en grados.

    bulto(lat, lon, ancho, punta=False)        gaussiana centrada en (lat, lon); con punta=True, exponencial
    hoyo(lat, lon, ancho, punta=False)         -bulto
    crater(lat, lon, ancho, borde, radio_borde) hoyo rodeado de un anillo elevado
    banda(lat, ancho)                          anillo zonal a una latitud
    lobulos(m, lat, ancho, fase)               m lóbulos en longitud alrededor de una latitud
    espiral(m, k, lat, ancho)                  m lóbulos cuya fase varía con la latitud: cos(m(φ - k·lat))

construir(nombre, forma, q) proyecta la forma, fija el relieve, guarda el .dat y la figura
e imprime un diagnóstico (fracción no representable con n <= 6, simetrías, coeficientes).
"""

from __future__ import annotations

from pathlib import Path

import armonicos as A
import numpy as np

R = np.radians


class Pieza:
    """Función f(θ, φ) con aritmética: suma, resta, producto por escalar y producto de piezas."""

    def __init__(self, f, nombre: str = ""):
        self.f, self.nombre = f, nombre

    def __call__(self, th, ph):
        return self.f(th, ph)

    def __add__(self, otra):
        otra = otra if isinstance(otra, Pieza) else constante(otra)
        return Pieza(
            lambda th, ph: self(th, ph) + otra(th, ph), f"{self.nombre} + {otra.nombre}"
        )

    __radd__ = __add__

    def __neg__(self):
        return Pieza(lambda th, ph: -self(th, ph), f"-({self.nombre})")

    def __sub__(self, otra):
        return self + (-otra if isinstance(otra, Pieza) else -otra)

    def __rsub__(self, otra):
        return (-self) + otra

    def __mul__(self, k):
        if isinstance(k, Pieza):
            return Pieza(
                lambda th, ph: self(th, ph) * k(th, ph), f"({self.nombre})·({k.nombre})"
            )
        return Pieza(lambda th, ph: k * self(th, ph), f"{k}·{self.nombre}")

    __rmul__ = __mul__


def constante(c):
    return Pieza(lambda th, ph: c + 0 * th, str(c))


def distancia_angular(th, ph, lat, lon):
    """Ángulo γ entre (θ, φ) y el punto (lat, lon): cos γ = sin θ sin θ0 cos(φ - φ0) + cos θ cos θ0."""
    t0, p0 = R(90 - lat), R(lon)
    c = np.sin(th) * np.sin(t0) * np.cos(ph - p0) + np.cos(th) * np.cos(t0)
    return np.arccos(np.clip(c, -1, 1))


def bulto(lat, lon=0.0, ancho=30.0, punta=False):
    """exp(-γ²/2σ²); con punta=True, exp(-γ/σ), que tiene un máximo anguloso."""
    s = R(ancho)
    if punta:
        return Pieza(
            lambda th, ph: np.exp(-distancia_angular(th, ph, lat, lon) / s),
            f"bulto({lat},{lon})",
        )
    return Pieza(
        lambda th, ph: np.exp(-(distancia_angular(th, ph, lat, lon) ** 2) / (2 * s**2)),
        f"bulto({lat},{lon})",
    )


def hoyo(lat, lon=0.0, ancho=30.0, punta=False):
    return -bulto(lat, lon, ancho, punta)


def crater(lat, lon=0.0, ancho=18.0, borde=0.35, radio_borde=None):
    """
    Hoyo gaussiano de anchura 'ancho' rodeado de un anillo de altura relativa 'borde'
    a una distancia angular 'radio_borde' (por defecto 2.6·ancho).

    Un mínimo de f por sí solo no produce una depresión visible: cerca de un punto, la
    esfera ya desciende como z ≈ ρ(1 - γ²/2), de modo que el fondo debe descender más
    rápido que eso. El anillo elevado aporta esa variación y hace visible la depresión.
    """
    s = R(ancho)
    r0 = R(radio_borde) if radio_borde is not None else 2.6 * s
    w = 0.7 * s

    def f(th, ph):
        g = distancia_angular(th, ph, lat, lon)
        return -np.exp(-(g**2) / (2 * s**2)) + borde * np.exp(
            -((g - r0) ** 2) / (2 * w**2)
        )

    return Pieza(f, f"crater({lat},{lon})")


def banda(lat, ancho=20.0):
    s = R(ancho)
    return Pieza(
        lambda th, ph: np.exp(-(((np.pi / 2 - th) - R(lat)) ** 2) / (2 * s**2)),
        f"banda({lat})",
    )


def lobulos(m, lat=0.0, ancho=40.0, fase=0.0):
    s = R(ancho)
    return Pieza(
        lambda th, ph: (
            np.exp(-(((np.pi / 2 - th) - R(lat)) ** 2) / (2 * s**2))
            * np.cos(m * (ph - R(fase)))
        ),
        f"lobulos({m})",
    )


def espiral(m, k=1.0, lat=0.0, ancho=40.0):
    """Envolvente gaussiana en latitud por cos(m(φ - k·lat)); k en radianes de giro por radián de latitud."""
    s = R(ancho)

    def f(th, ph):
        la = np.pi / 2 - th
        return np.exp(-((la - R(lat)) ** 2) / (2 * s**2)) * np.cos(m * (ph - k * la))

    return Pieza(f, f"espiral({m},{k})")


# ---------------------------------------------------------------------------
# Diagnóstico
# ---------------------------------------------------------------------------
def error_truncamiento(func, C, S, nquad: int = 64) -> float:
    """||f - f_6|| / ||f - <f>||: fracción de la variación de f que no cabe en grado <= 6."""
    x, w = np.polynomial.legendre.leggauss(nquad)
    th = np.arccos(x)
    ph = np.arange(2 * nquad) * np.pi / nquad
    F = func(*np.meshgrid(th, ph, indexing="xy"))
    F6 = A.sintetizar(C, S, th, ph)
    media = np.sum(w * F) / np.sum(w * np.ones_like(F))
    return float(np.sqrt(np.sum(w * (F - F6) ** 2) / np.sum(w * (F - media) ** 2)))


def simetrias(C, S) -> list[str]:
    """Simetrías de f que se deducen de qué coeficientes son distintos de cero."""
    filas = A.coeficientes_no_nulos(C, S)
    terminos = [(n, m) for n, m, _, _ in filas if (n, m) != (0, 0)]
    ms = sorted({m for _, m in terminos})
    salida = [
        f"órdenes m presentes: {ms}",
        f"grados n presentes: {sorted({n for n, _ in terminos})}",
    ]
    if ms == [0]:
        salida.append("solo m = 0: simetría de revolución alrededor del eje z")
    else:
        positivos = [m for m in ms if m > 0]
        k = min(positivos) if positivos else 0
        if k > 1 and all(m % k == 0 for m in ms):
            salida.append(
                f"m múltiplos de {k}: invariante bajo giros de 360°/{k} = {360 // k}° en torno a z"
            )
    if all(abs(s) < 1e-12 for _, m, _, s in filas if m > 0):
        salida.append("todos los S_nm = 0: simétrica respecto al plano φ = 0 (φ -> -φ)")
    if terminos and all((n + m) % 2 == 0 for n, m in terminos):
        salida.append(
            "solo n + m par: simétrica respecto al ecuador (z -> -z), pues P_nm(-x) = (-1)^(n+m) P_nm(x)"
        )
    return salida


def construir(
    nombre, forma, q=1.5, vistas=None, umbral=0.01, carpeta: Path | None = None
):
    """Forma -> coeficientes -> .dat y figura. Devuelve (C, S)."""
    vistas = vistas or A.VISTAS_POR_DEFECTO
    carpeta = Path(carpeta) if carpeta else A.RESULTADOS

    C, S = A.proyectar(forma)
    err_trunc = error_truncamiento(forma, C, S)
    C, S = A.fijar_relieve(C, S, q)
    C, S, err_red = A.limpiar(C, S, umbral=umbral, cifras=3)

    ruta_dat = carpeta / "dats" / f"{nombre}.dat"
    ruta_png = carpeta / "figuras" / f"{nombre}.png"
    A.escribir_dat(ruta_dat, C, S)
    A.render_vistas(C, S, ruta_png, nombre, vistas=vistas)

    r_min, r_max = A.radios_script(C, S)
    print(f"\n=== {nombre} ===")
    print(
        f"fracción no representable con n <= 6: {100 * err_trunc:.0f} %"
        + (
            "  (aviso: piezas demasiado estrechas para grado 6)"
            if err_trunc > 0.35
            else ""
        )
    )
    print(f"error introducido al redondear: {100 * err_red:.2f} %")
    print(
        f"radio en spheriharm1.py: mínimo {r_min:.2f}, máximo {r_max:.2f}"
        + ("  (aviso: relieve insuficiente, aumentar q)" if r_min > 0.95 else "")
    )
    for linea in simetrias(C, S):
        print("  ·", linea)
    A.imprimir_coeficientes(C, S)
    print(
        "->",
        ruta_dat.relative_to(A.RAIZ) if ruta_dat.is_relative_to(A.RAIZ) else ruta_dat,
    )
    return C, S
