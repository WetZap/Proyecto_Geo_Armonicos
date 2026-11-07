"""
generar_todo.py — Regenera todos los resultados de resultados/ a partir del código.

    python generar_todo.py

Pasos:
  1. Figuras del guion (figuras_guion.py): .dat, figuras, comparativa y tabla LaTeX.
  2. Formas nuevas (formas_nuevas.py).
  3. Hojas de exploración de todas las familias (explorar.py, semilla 0).
  4. Figuras explicativas del método: efecto de q y efecto del borde de un cráter.
  5. Ejecución del script original y de la variante corregida sobre dos .dat.
"""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402

import armonicos as A  # noqa: E402
import explorar  # noqa: E402
import figuras_guion  # noqa: E402
import formas_nuevas  # noqa: E402
from ejecutar_original import ejecutar  # noqa: E402
from piezas import bulto, construir, crater  # noqa: E402


def panel(casos, archivo, vistas):
    """casos: lista de (etiqueta, C, S). Una fila por caso y una columna por vista."""
    fig = plt.figure(figsize=(3.2 * len(vistas), 3.2 * len(casos)))
    for i, (etiqueta, C, S) in enumerate(casos):
        r_min, _ = A.radios_script(C, S)
        for j, (elev, azim) in enumerate(vistas):
            ax = fig.add_subplot(len(casos), len(vistas), i * len(vistas) + j + 1, projection="3d")
            A.dibujar(ax, C, S, elev, azim, f"{etiqueta}  (radio mínimo {r_min:.2f})" if j == 0 else None)
    fig.subplots_adjust(0, 0, 1, 0.96, 0, 0.08)
    fig.savefig(archivo, dpi=90)
    plt.close(fig)


def figuras_metodo():
    vistas = [(25, -60), (5, -60)]
    manzana = figuras_guion.MANZANA

    C0, S0 = A.proyectar(manzana)
    casos = []
    for q in (0.5, 1.0, 2.5):
        C, S = A.fijar_relieve(C0, S0, q)
        casos.append((f"q = {q}", C, S))
    panel(casos, A.DIR_FIGURAS / "metodo_efecto_q.png", vistas)

    lateral = 0.9 * bulto(0, 0, ancho=60)
    sin_borde = -(bulto(90, 0, ancho=18) + bulto(-90, 0, ancho=18)) + lateral
    casos = []
    for etiqueta, forma in (("hoyo sin borde", sin_borde), ("cráter con borde", manzana)):
        C, S = A.fijar_relieve(*A.proyectar(forma), 2.5)
        casos.append((etiqueta, C, S))
    panel(casos, A.DIR_FIGURAS / "metodo_efecto_borde.png", vistas)


def main():
    print("1. Figuras del guion")
    figuras_guion.main(list(figuras_guion.FIGURAS))

    print("\n2. Formas nuevas")
    for nombre, (forma, q, vistas) in formas_nuevas.FORMAS.items():
        construir(nombre, forma, q, vistas)

    print("\n3. Exploración")
    explorar.main(["--n", "8", "--semilla", "0"])

    print("\n4. Figuras del método")
    figuras_metodo()

    print("\n5. Script original y variante corregida")
    for nombre, elev, azim in (("manzana", 25, -60), ("espiral", 30, -60)):
        for corregido in (False, True):
            ruta = ejecutar(A.DIR_DATS / f"{nombre}.dat", elev, azim, corregido)
            print("->", ruta.relative_to(A.RAIZ))


if __name__ == "__main__":
    main()
