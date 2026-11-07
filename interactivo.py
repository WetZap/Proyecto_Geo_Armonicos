"""
interactivo.py — Ajuste manual de coeficientes con deslizadores.

Uso:
    python interactivo.py resultados/dats/espiral.dat C44 S54 C64 C00
    python interactivo.py original/spheriharm1.dat C20 C32 S32

Cada deslizador controla un coeficiente (Cnm o Snm, con 0 <= m <= n <= 6; S_n0 no existe
porque sin(0·φ) = 0). Las dos vistas se redibujan al cambiar un valor y pueden girarse
con el ratón. «Guardar» escribe <fichero>_ajustado.dat en el mismo formato.
"""

import re
import sys

import matplotlib.pyplot as plt
from matplotlib.widgets import Button, Slider

import armonicos as A


def interpretar(nombre):
    coincidencia = re.fullmatch(r"([CS])(\d)(\d)", nombre.upper())
    if not coincidencia:
        raise SystemExit(f"Nombre no válido: {nombre} (formato Cnm o Snm, p. ej. C20, S54)")
    tipo, n, m = coincidencia.group(1), int(coincidencia.group(2)), int(coincidencia.group(3))
    if not (0 <= m <= n <= A.DMAX) or (tipo == "S" and m == 0):
        raise SystemExit(f"{nombre}: se requiere 0 <= m <= n <= 6 y, para S, m > 0")
    return tipo, n, m


def main(argv):
    if len(argv) < 2:
        raise SystemExit(__doc__)
    ruta, nombres = argv[0], argv[1:]
    C, S = A.leer_dat(ruta)
    parametros = [interpretar(x) for x in nombres]

    alto = 5.5 + 0.35 * len(parametros)
    fig = plt.figure(figsize=(10, alto))
    base = 0.35 * len(parametros) / alto + 0.08
    ax1 = fig.add_axes([0.0, base, 0.5, 0.95 - base], projection="3d")
    ax2 = fig.add_axes([0.5, base, 0.5, 0.95 - base], projection="3d")
    A.dibujar(ax1, C, S, 30, -60)
    A.dibujar(ax2, C, S, 90, -90)

    def redibujar():
        for ax in (ax1, ax2):
            elev, azim = ax.elev, ax.azim  # conserva la orientación elegida con el ratón
            ax.clear()
            A.dibujar(ax, C, S, elev, azim)
        fig.canvas.draw_idle()

    deslizadores = []
    for i, (tipo, n, m) in enumerate(parametros):
        matriz = C if tipo == "C" else S
        limite = max(10.0, 2 * abs(matriz[n, m]))
        eje = fig.add_axes([0.12, 0.04 + 0.35 * i / alto, 0.6, 0.25 / alto])
        deslizador = Slider(eje, f"{tipo}{n}{m}", -limite, limite, valinit=matriz[n, m])

        def al_cambiar(valor, matriz=matriz, n=n, m=m):
            matriz[n, m] = valor
            redibujar()

        deslizador.on_changed(al_cambiar)
        deslizadores.append(deslizador)

    boton = Button(fig.add_axes([0.8, 0.04, 0.12, 0.05]), "Guardar")

    def guardar(_):
        salida = re.sub(r"\.dat$", "", ruta) + "_ajustado.dat"
        A.escribir_dat(salida, C, S)
        print("Guardado en", salida)

    boton.on_clicked(guardar)
    fig._widgets = (deslizadores, boton)  # mantiene las referencias vivas
    plt.show()
    return fig, deslizadores, boton


if __name__ == "__main__":
    main(sys.argv[1:])
