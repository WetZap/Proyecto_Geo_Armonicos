"""
mi_forma.py — Plantilla para construir una forma nueva.

Procedimiento:
  1. Describir la forma como combinación de piezas (ver piezas.py).
  2. Ajustar posiciones, anchuras, pesos y q comparando la figura generada con la deseada.
     Anchuras por debajo de ~25° no son representables con grado <= 6.
  3. Para el ajuste fino de coeficientes concretos:
         python interactivo.py resultados/dats/<NOMBRE>.dat C20 C40

Uso:
    python mi_forma.py
"""

from piezas import banda, bulto, construir, crater, espiral, hoyo, lobulos  # noqa: F401

NOMBRE = "mi_forma"

# Ejemplo: dos máximos opuestos en el ecuador y un cráter en el polo norte.
FORMA = bulto(0, 0, ancho=30) + bulto(0, 180, ancho=30) + 0.8 * crater(90, 0, ancho=20)

# Relieve: el máximo de f queda en radio 1.1 y el mínimo en 1 - 0.1·q (solo se modifica C_00).
Q = 2.0

VISTAS = [(30, -60, "Oblicua"), (90, -90, "Vista polar"), (0, -90, "Vista ecuatorial")]

if __name__ == "__main__":
    construir(NOMBRE, FORMA, Q, VISTAS)
