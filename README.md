# GROMACS MD Pipeline Automatizado (Proteina + Ligando)

Pipeline para simulaciones de dinamica molecular de complejos proteina-ligando con GROMACS.

## Requisitos

Dependencias validadas al inicio de `md_pipeline_mejorado.sh`:

- `gmx` o `gmx_mpi`
- `awk`
- `sed`
- `grep`
- `bc`

Dependencias recomendadas para analisis y reportes:

- `python3`
- `numpy`
- `matplotlib`

## Estructura esperada

```text
Automatizacion_proteina_ligando/
|- md_pipeline_mejorado.sh
|- plot_analysis.py
|- README.md
|- proteinas/
|  `- mi_proteina/
|     |- proteina.gro
|     |- topol.top
|     `- posre.itp
|- ligandos/
|  `- mi_ligando/
|     |- ligando.gro
|     |- ligando.itp
|     `- ligando.prm   (opcional)
|- mdp/
|  |- ions.mdp
|  |- em.mdp
|  |- nvt.mdp
|  |- npt.mdp
|  `- md_prod.mdp
`- MD_RUN/             (salida automatica)
```

## Ejecucion interactiva

```bash
chmod +x md_pipeline_mejorado.sh
./md_pipeline_mejorado.sh
```

## Modo no interactivo (flags CLI)

El pipeline entra automaticamente a modo no interactivo cuando detecta los obligatorios:

- `--prot`
- `--lig`
- `--ff`

Si faltan obligatorios, el script hace fallback a modo interactivo (salvo que uses `--non-interactive`, donde falla de forma explicita).

### Ejemplo minimo no interactivo

```bash
./md_pipeline_mejorado.sh \
  --prot caspasa9 \
  --lig M4-A \
  --ff charmm36-jul2022.ff
```

### Ejemplo no interactivo completo

```bash
./md_pipeline_mejorado.sh \
  --prot caspasa9 \
  --lig M4-A \
  --ff charmm36-jul2022.ff \
  --box dodecahedron \
  --box-dist 1.2 \
  --water tip3p \
  --ion 0.15 \
  --prod-ns 50
```

Tambien se soportan alias:

- `--box` (alias de `--box-type`)
- `--ion` (alias de `--ion-conc`)

## Archivo de configuracion (`--config`)

Puedes cargar parametros desde archivo:

```bash
./md_pipeline_mejorado.sh --config pipeline.conf
```

Formato:

```bash
PROT=caspasa9
LIG=M4-A
FF=charmm36-jul2022.ff
BOX_TYPE=dodecahedron
BOX_DIST=1.2
WATER_MODEL=tip3p
ION_CONC=0.15
PROD_NS=50
```

## Modo `--dry-run`

`--dry-run` ejecuta validaciones, armado del sistema y verificaciones de `grompp`, pero no ejecuta etapas largas de `mdrun`.

```bash
./md_pipeline_mejorado.sh \
  --prot caspasa9 \
  --lig M4-A \
  --ff charmm27.ff \
  --dry-run
```

En dry-run, revisa especialmente:

- `logs/dryrun_grompp_em.log`
- `logs/dryrun_grompp_nvt.log`
- `logs/dryrun_grompp_npt.log`
- `logs/dryrun_grompp_md.log`

## Reanudar corridas

```bash
./md_pipeline_mejorado.sh --resume
```

El checkpoint se guarda en cada paso exitoso y se valida con parser seguro antes de cargarse.

## Seguridad y robustez implementadas

- Se elimino parseo fragil con `grep -P` / `grep -oP`; se usan alternativas con `awk`/`sed`.
- Se valida `bc` explicitamente al inicio para comparaciones de punto flotante.
- La carga de checkpoints evita `source` directo de contenido no confiable.
- La edicion de `topol.top` ahora es encapsulada, idempotente y con verificacion estricta post-edicion.
- Manejo de errores reforzado con `set -E` + `trap ERR` + limpieza de temporales.

## Reproducibilidad

Al final de cada corrida se genera `SUMMARY.txt` con:

- Fecha/hora exacta
- Hostname
- Binario y version de GROMACS
- Info basica de CPU/GPU
- Parametros clave de simulacion
- Modo de ejecucion (`normal` o `dry-run`)

## Ejemplos para SLURM/SGE

### SLURM

```bash
#!/bin/bash
#SBATCH -J md_prot_lig
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=48:00:00

./md_pipeline_mejorado.sh \
  --prot caspasa9 \
  --lig M4-A \
  --ff charmm36-jul2022.ff \
  --box dodecahedron \
  --water tip3p \
  --ion 0.15 \
  --prod-ns 100
```

### SGE

```bash
#!/bin/bash
#$ -N md_prot_lig
#$ -pe smp 16
#$ -l h_rt=48:00:00

./md_pipeline_mejorado.sh \
  --prot caspasa9 \
  --lig M4-A \
  --ff charmm36-jul2022.ff \
  --prod-ns 100
```

## Licencia

Uso libre para investigacion academica.
