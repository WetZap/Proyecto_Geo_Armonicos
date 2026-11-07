"""
ejecutar_original.py — Ejecuta el script de la práctica con un fichero de coeficientes.

spheriharm1.py lee siempre ./spheriharm1.dat y guarda en ./Imagenes/. Para no modificar
original/spheriharm1.dat, este lanzador copia el .dat indicado a un directorio temporal,
ejecuta allí el script y mueve la imagen a resultados/script_original/.

Uso:
    python ejecutar_original.py resultados/dats/manzana.dat 25 -60
    python ejecutar_original.py resultados/dats/manzana.dat 25 -60 --corregido

--corregido usa spheriharm1_corregido.py (misma suma y misma superficie; solo añade
iluminación y aspecto 1:1:1).
"""

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

RAIZ = Path(__file__).resolve().parent
SALIDA = RAIZ / "resultados" / "script_original"


def ejecutar(dat, elev, azim, corregido=False, nombre=None):
    script = RAIZ / ("spheriharm1_corregido.py" if corregido else "original/spheriharm1.py")
    dat = Path(dat).resolve()
    nombre = nombre or dat.stem + ("_corregido" if corregido else "")
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        shutil.copy(dat, tmp / "spheriharm1.dat")
        (tmp / "Imagenes").mkdir()
        subprocess.run([sys.executable, str(script), f"{elev:g}", f"{azim:g}", nombre],
                       cwd=tmp, check=True, env={**os.environ, "MPLBACKEND": "Agg"})
        imagen = next((tmp / "Imagenes").glob("*.png"))
        SALIDA.mkdir(parents=True, exist_ok=True)
        destino = SALIDA / imagen.name
        shutil.move(imagen, destino)
    return destino


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("dat", help="fichero de coeficientes")
    p.add_argument("inc", type=float, help="inclinación de la vista (90 = polar)")
    p.add_argument("az", type=float, help="azimut de la vista")
    p.add_argument("--corregido", action="store_true", help="usar spheriharm1_corregido.py")
    a = p.parse_args()
    print("->", ejecutar(a.dat, a.inc, a.az, a.corregido).relative_to(RAIZ))
