"""
figuras_guion.py — Coeficientes (C_nm, S_nm) para las figuras de la Figura 4 del guion:
Bombilla, Ojo, Manzana y Espiral.

Método: para cada figura se define una función objetivo f(θ, φ) cuyas simetrías y
rasgos coinciden con los observados en el guion, y se proyecta sobre los armónicos de
grado <= 6. Después se fija el relieve con C_00 y se normalizan los coeficientes.
La justificación de cada función objetivo está en docs/figuras.md.

Uso:
    python figuras_guion.py              # las cuatro figuras
    python figuras_guion.py espiral      # una sola

Salida (en resultados/):
    dats/<figura>.dat                    coeficientes en el formato de spheriharm1.dat
    figuras/<figura>.png                 tres vistas de la superficie
    figuras/comparativa_guion.png        Figura 4 del guion junto a las figuras obtenidas
    coeficientes_figuras.tex             tablas de coeficientes no nulos para LaTeX
"""

import sys

import armonicos as A
from piezas import bulto, crater, espiral

# ---------------------------------------------------------------------------
# Funciones objetivo
# ---------------------------------------------------------------------------
# Bombilla: simetría de revolución (solo m = 0). Máximo estrecho en el polo sur y
# pequeña depresión en el polo norte.
BOMBILLA = bulto(-90, 0, ancho=20) + 0.3 * crater(90, 0, ancho=15, borde=0.33, radio_borde=30)

# Ojo: máximo con punta centrado en el ecuador (lat 0°, lon 0°). Las curvas de nivel son
# planos perpendiculares a ese eje, de modo que f depende solo de la distancia a él.
OJO = bulto(0, 0, ancho=40, punta=True)

# Manzana: depresiones con borde en ambos polos (simetría ecuatorial) y un lado más alto
# en lon 0°, que rompe la simetría de revolución.
MANZANA = (crater(90, 0, ancho=18, borde=0.35, radio_borde=40)
           + crater(-90, 0, ancho=18, borde=0.35, radio_borde=40)
           + 0.9 * bulto(0, 0, ancho=60))

# Espiral: cuatro lóbulos (m = 4) cuya fase gira con la latitud.
ESPIRAL = espiral(4, k=1.0, lat=10, ancho=40)

# nombre: (función objetivo, q, vistas (INC, AZ, etiqueta))
# q fija el relieve: el máximo de f queda en radio 1.1 y el mínimo en 1 - 0.1·q.
FIGURAS = {
    "bombilla": (BOMBILLA, 1.0, [(20, 30, "Oblicua"), (90, 0, "Vista polar"), (0, 0, "Vista ecuatorial")]),
    "ojo": (OJO, 2.0, [(20, -60, "Oblicua"), (90, -90, "Vista polar"), (10, -60, "Vista ecuatorial")]),
    "manzana": (MANZANA, 2.5, [(25, -60, "Oblicua"), (90, -90, "Vista polar"), (5, -60, "Vista ecuatorial")]),
    "espiral": (ESPIRAL, 3.0, [(30, -60, "Oblicua"), (90, -90, "Vista polar"), (5, -60, "Vista ecuatorial")]),
}

DIR_REFERENCIA = A.RAIZ / "original" / "figura4"


def calcular(nombre):
    """Devuelve (C, S, error de redondeo) para una figura."""
    forma, q, _ = FIGURAS[nombre]
    C, S = A.proyectar(forma)
    C, S = A.fijar_relieve(C, S, q)
    return A.limpiar(C, S, umbral=0.01, cifras=3)


def tabla_latex(resultados):
    lineas = []
    for nombre, (C, S) in resultados.items():
        lineas += [
            r"\begin{table}[htbp]", r"\centering", r"\small",
            rf"\caption{{Coeficientes no nulos de la figura \emph{{{nombre.capitalize()}}} "
            r"(normalizados a $\max|c| = 10$).}",
            rf"\label{{tab:coef-{nombre}}}",
            r"\begin{tabular}{rrrr}", r"\hline", r"$n$ & $m$ & $C_{nm}$ & $S_{nm}$ \\", r"\hline",
        ]
        lineas += [f"{n} & {m} & {c:.3g} & {s:.3g} \\\\" for n, m, c, s in A.coeficientes_no_nulos(C, S)]
        lineas += [r"\hline", r"\end{tabular}", r"\end{table}", ""]
    return "\n".join(lineas)


def comparativa_con_guion(nombres, archivo):
    """Coloca cada figura de referencia del guion junto a la obtenida."""
    from PIL import Image

    filas = []
    for nombre in nombres:
        ref = DIR_REFERENCIA / f"{nombre}.png"
        obt = A.DIR_FIGURAS / f"{nombre}.png"
        if not ref.exists():
            continue
        a, b = Image.open(ref).convert("RGB"), Image.open(obt).convert("RGB")
        h = 380
        a = a.resize((int(a.width * h / a.height), h))
        b = b.resize((int(b.width * h / b.height), h))
        fila = Image.new("RGB", (a.width + b.width + 20, h), "white")
        fila.paste(a, (0, 0))
        fila.paste(b, (a.width + 20, 0))
        filas.append(fila)
    if not filas:
        return
    lienzo = Image.new("RGB", (max(f.width for f in filas), sum(f.height for f in filas)), "white")
    y = 0
    for fila in filas:
        lienzo.paste(fila, (0, y))
        y += fila.height
    lienzo.save(archivo)


def main(nombres):
    resultados = {}
    for nombre in nombres:
        C, S, err = calcular(nombre)
        resultados[nombre] = (C, S)
        A.escribir_dat(A.DIR_DATS / f"{nombre}.dat", C, S)
        A.render_vistas(C, S, A.DIR_FIGURAS / f"{nombre}.png", nombre.capitalize(),
                        vistas=FIGURAS[nombre][2], dpi=120)
        r_min, r_max = A.radios_script(C, S)
        print(f"\n=== {nombre.upper()} ===")
        print(f"error de redondeo: {100 * err:.2f} %  |  radio en el script: {r_min:.2f} – {r_max:.2f}")
        A.imprimir_coeficientes(C, S)
        elev, azim, _ = FIGURAS[nombre][2][0]
        print(f"con el script original: python ejecutar_original.py resultados/dats/{nombre}.dat {elev} {azim}")

    if set(nombres) == set(FIGURAS):
        (A.RESULTADOS / "coeficientes_figuras.tex").write_text(tabla_latex(resultados))
        comparativa_con_guion(nombres, A.DIR_FIGURAS / "comparativa_guion.png")


if __name__ == "__main__":
    main(sys.argv[1:] or list(FIGURAS))
