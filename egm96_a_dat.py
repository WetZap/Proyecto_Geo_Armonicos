#!/usr/bin/env python3
"""
egm96_a_dat.py — rellena spheriharm1.dat con los coeficientes de EGM96 hasta grado y orden 5.

    python egm96_a_dat.py EGM96.gfc spheriharm1_EGM96.dat [--directo]

Pasos:
 1. Lee C̄nm, S̄nm (n <= 5) del .gfc (totalmente normalizados 4π, sin fase de Condon-Shortley).
 2. Quita el elipsoide de referencia WGS84 (C̄20 y C̄40 del campo normal) -> ondulaciones.
    C00 = 1 también se quita (la masa del elipsoide es la de la Tierra). Grado 1 ya es 0 en EGM96.
 3. Convierte a la convención de spheriharm1.py (K_nm con 1/(4π) y lpmv con (-1)^m):
        C_script = (-1)^m · sqrt(2 - δ_m0) · C̄   (el factor global sqrt(4π) no cambia la figura)
    Con --directo NO se convierte (valores tal cual, solo para comparar).
 4. Multiplica por 1e6 (unidades de 1e-6) para que el '1e-6' de las casillas vacías sea despreciable.
"""
import sys
import numpy as np
import scipy.special as ss

NMAX, DMAX = 5, 6
# Campo normal del elipsoide WGS84 (coeficientes totalmente normalizados)
C20_ELIP = -0.484166774985e-3
C40_ELIP = 0.790303733511e-6
ESCALA = 1e6


def leer_gfc(ruta, nmax=NMAX):
    C = np.zeros((DMAX + 1, DMAX + 1)); S = np.zeros_like(C)
    with open(ruta) as f:
        for line in f:
            if line.startswith("end_of_head"):
                break
        for line in f:
            p = line.split()
            if not p or p[0] != "gfc":
                continue
            n, m = int(p[1]), int(p[2])
            if n > nmax:
                break
            C[n, m], S[n, m] = float(p[3]), float(p[4])
    return C, S


def a_script(C, S):
    C, S = C.copy(), S.copy()
    for m in range(DMAX + 1):
        f = (-1) ** m * np.sqrt(2.0 if m else 1.0)
        C[:, m] *= f; S[:, m] *= f
    return C, S


def escribir_dat(ruta, C, S, cifras=6):
    def fmt(v):
        return "1e-6" if abs(v) < 1e-12 else f"{v:.{cifras}g}"
    lin = ["Please fill in the list with the C and S coefficients of your choice",
           "n\\m \t" + "\t".join(str(m) for m in range(DMAX + 1))]
    for n in range(DMAX + 1):
        lin.append(f"{n}\t" + "\t".join(f"({fmt(C[n, m])},{fmt(S[n, m])})" for m in range(n + 1)))
    with open(ruta, "w") as f:
        f.write("\n".join(lin) + "\n")


def main(argv):
    gfc, salida = argv[0], argv[1]
    C, S = leer_gfc(gfc)
    C[0, 0] = 0.0
    C[2, 0] -= C20_ELIP
    C[4, 0] -= C40_ELIP
    if "--directo" not in argv:
        C, S = a_script(C, S)
    escribir_dat(salida, C * ESCALA, S * ESCALA)
    print("Escrito", salida)


if __name__ == "__main__":
    main(sys.argv[1:])
