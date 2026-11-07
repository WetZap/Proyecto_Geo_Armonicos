"""
Comprobaciones del método. Ejecutar desde la raíz del repositorio:

    python -m pytest -q          (o, sin pytest: python tests/test_armonicos.py)
"""

import sys
import tempfile
from pathlib import Path

import numpy as np
import scipy.special as ss

RAIZ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RAIZ))

import armonicos as A  # noqa: E402
import figuras_guion as G  # noqa: E402
import formas_nuevas as F  # noqa: E402
from piezas import bulto, crater  # noqa: E402

DAT_ORIGINAL = RAIZ / "original" / "spheriharm1.dat"


# ---------------------------------------------------------------------------
# Referencia: lectura y suma copiadas literalmente de original/spheriharm1.py
# ---------------------------------------------------------------------------
def _leer_como_el_original(ruta):
    C = np.zeros((7, 7))
    S = np.zeros((7, 7))
    nline = 0
    with open(ruta) as f:
        for line in f:
            data = line.split()
            if nline > 1:
                for nfield in range(1, len(data)):
                    C[nline - 2, nfield - 1] = data[nfield].split("(")[1].split(",")[0]
                    S[nline - 2, nfield - 1] = data[nfield].split(",")[1].split(")")[0]
            nline += 1
    return C, S


def _sumar_como_el_original(C, S):
    vphi = np.linspace(0, 2 * np.pi, 90)
    vtheta = np.linspace(0, np.pi, 45)
    RR2 = np.zeros((len(vphi), len(vtheta)))
    for degree in range(7):
        for order in range(degree + 1):
            Plm = ss.lpmv(order, degree, np.cos(vtheta))
            SC = C[degree, order] * np.cos(vphi * order) + S[degree, order] * np.sin(vphi * order)
            f = ss.factorial(degree - order) / ss.factorial(degree + order)
            K_norm = np.sqrt((2 * degree + 1) / (4 * np.pi) * f)
            RR2 += K_norm * np.outer(SC, Plm)
    return RR2


def _coeficientes_aleatorios(semilla=0):
    rng = np.random.default_rng(semilla)
    C = np.tril(rng.normal(size=(7, 7)))
    S = np.tril(rng.normal(size=(7, 7)))
    S[:, 0] = 0
    return C, S


def _terminos(C, S):
    return {(n, m) for n, m, _, _ in A.coeficientes_no_nulos(C, S) if (n, m) != (0, 0)}


# ---------------------------------------------------------------------------
# Equivalencia con el script de la práctica
# ---------------------------------------------------------------------------
def test_sintesis_identica_al_script():
    C, S = _leer_como_el_original(DAT_ORIGINAL)
    assert np.allclose(A.sintetizar(C, S), _sumar_como_el_original(C, S), atol=1e-12)


def test_lectura_identica_al_script():
    for ruta in [DAT_ORIGINAL, *sorted((RAIZ / "resultados" / "dats").glob("*.dat"))]:
        C1, S1 = _leer_como_el_original(ruta)
        C2, S2 = A.leer_dat(ruta)
        assert np.array_equal(C1, C2) and np.array_equal(S1, S2), ruta.name


def test_escritura_y_lectura(tmp_path=None):
    ruta = Path(tmp_path or tempfile.gettempdir()) / "prueba.dat"
    C, S = _coeficientes_aleatorios(1)
    A.escribir_dat(ruta, C, S, cifras=12)
    C2, S2 = A.leer_dat(ruta)
    # los ceros se escriben como 1e-6, igual que en el fichero original
    assert np.allclose(C, C2, rtol=1e-10, atol=2e-6) and np.allclose(S, S2, rtol=1e-10, atol=2e-6)


# ---------------------------------------------------------------------------
# Proyección y teorema de adición
# ---------------------------------------------------------------------------
def test_proyeccion_recupera_los_coeficientes():
    """Si f es ya una suma de armónicos con n <= 6, la proyección es exacta."""
    C0, S0 = _coeficientes_aleatorios(0)

    def f(th, ph):
        total = np.zeros_like(th)
        for n in range(7):
            for m in range(n + 1):
                total += A.knorm(n, m) * ss.lpmv(m, n, np.cos(th)) * (C0[n, m] * np.cos(m * ph) + S0[n, m] * np.sin(m * ph))
        return total

    C1, S1 = A.proyectar(f)
    assert np.abs(C1 - C0).max() < 1e-10 and np.abs(S1 - S0).max() < 1e-10


def test_teorema_de_adicion():
    a = [1.0, 0.8, 0.6, 0.4, 0.3, 0.2, 0.1]
    theta0, phi0 = 1.1, 0.7
    C, S = A.bulto_axial(a, theta0, phi0)
    vtheta, vphi, T, P = A.malla()
    cos_g = np.sin(T) * np.sin(theta0) * np.cos(P - phi0) + np.cos(T) * np.cos(theta0)
    directa = sum(a[n] * ss.eval_legendre(n, cos_g) for n in range(7))
    assert np.abs(A.sintetizar(C, S) - directa).max() < 1e-10


# ---------------------------------------------------------------------------
# Relieve
# ---------------------------------------------------------------------------
def test_fijar_relieve():
    C, S = _coeficientes_aleatorios(2)
    for q in (0.5, 1.0, 2.0, 3.0):
        C2, S2 = A.fijar_relieve(C, S, q)
        r_min, r_max = A.radios_script(C2, S2)
        assert abs(r_max - 1.1) < 1e-12 and abs(r_min - (1 - 0.1 * q)) < 1e-12
        diferencia = A.sintetizar(C2, S2) - A.sintetizar(C, S)
        assert np.ptp(diferencia) < 1e-12  # solo cambia una constante: los colores no varían
        assert np.array_equal(np.delete(C2.ravel(), 0), np.delete(C.ravel(), 0)) and np.array_equal(S2, S)


def test_crater_produce_depresion_en_la_silueta():
    """Con el borde, z = ρ cos θ crece al alejarse del polo: hay una depresión visible."""
    def muesca(forma, q=3.0):
        C, S = A.proyectar(forma)
        C, S = A.fijar_relieve(C, S, q)
        th = np.linspace(0, np.pi / 2, 200)
        fmax = A.sintetizar(C, S).max()
        z = 0.1 * (10 + A.sintetizar(C, S, th, np.array([np.pi]))[0] / fmax) * np.cos(th)
        return z.max() - z[0]

    lateral = 0.9 * bulto(0, 0, ancho=60)
    sin_borde = -(bulto(90, 0, 18) + bulto(-90, 0, 18)) + lateral
    con_borde = crater(90, 0, 18, 0.35, 40) + crater(-90, 0, 18, 0.35, 40) + lateral
    assert muesca(sin_borde) < 0.03 < 0.08 < muesca(con_borde)


# ---------------------------------------------------------------------------
# Simetrías -> coeficientes no nulos
# ---------------------------------------------------------------------------
def test_simetria_cubica_del_dado():
    C, S = A.proyectar(F.FORMAS["dado"][0])
    C[0, 0] = 0  # se descarta la media para que no domine la normalización
    C, S, _ = A.limpiar(C, S, umbral=0.01)
    assert _terminos(C, S) == {(4, 0), (4, 4), (6, 0), (6, 4)}
    # la simetría cúbica fija el cociente, sea cual sea la forma de los máximos
    C, S = A.proyectar(F.FORMAS["dado"][0])
    assert abs(C[4, 4] / C[4, 0] - 2 * np.sqrt(5 / 14)) < 1e-8


def test_cacahuete_zonal_y_ecuatorial():
    C, S = A.proyectar(F.FORMAS["cacahuete"][0])
    C, S, _ = A.limpiar(C, S)
    assert all(m == 0 and n % 2 == 0 for n, m in _terminos(C, S))


def test_estrella_orden_5():
    C, S = A.proyectar(F.FORMAS["estrella"][0])
    C, S, _ = A.limpiar(C, S)
    assert all(m % 5 == 0 for _, m in _terminos(C, S)) and np.all(S == 0)


def test_figuras_del_guion():
    C, S, _ = G.calcular("bombilla")
    assert all(m == 0 for _, m in _terminos(C, S))                       # revolución
    C, S, _ = G.calcular("ojo")
    assert all((n + m) % 2 == 0 for n, m in _terminos(C, S)) and np.all(S == 0)
    C, S, _ = G.calcular("manzana")
    assert all((n + m) % 2 == 0 for n, m in _terminos(C, S)) and np.all(S == 0)
    C, S, _ = G.calcular("espiral")
    assert {m for _, m in _terminos(C, S)} == {4}


def test_dats_publicados_coinciden_con_el_codigo():
    """Los .dat de resultados/ son los que genera el código actual."""
    for nombre in G.FIGURAS:
        C, S, _ = G.calcular(nombre)
        C2, S2 = A.leer_dat(RAIZ / "resultados" / "dats" / f"{nombre}.dat")
        assert np.allclose(C, C2, atol=2e-6) and np.allclose(S, S2, atol=2e-6), nombre


if __name__ == "__main__":
    pruebas = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    fallos = 0
    for prueba in pruebas:
        try:
            prueba()
            print(f"OK     {prueba.__name__}")
        except AssertionError as e:
            fallos += 1
            print(f"FALLO  {prueba.__name__} {e}")
    print(f"\n{len(pruebas) - fallos}/{len(pruebas)} pruebas superadas")
    sys.exit(1 if fallos else 0)
