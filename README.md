# GROMACS MD Pipeline Automatizado v4.5 (Proteína + Ligando)

Pipeline para simulaciones de dinámica molecular de complejos proteína-ligando con GROMACS.

## Novedades v4.5

### `md_pipeline_mejorado.sh`
- **Soporte GPU (`--gpu-id`)**: Selecciona GPU específica para `mdrun` (esencial para Colab/HPC)
- **Grupos de energía automáticos**: Inyecta `energygrps = Protein LIG` en el MDP de producción para garantizar análisis de energía de interacción
- **Extender simulaciones (`--extend`)**: Continúa una simulación existente sin re-preparar el sistema
- **Progreso visual con ETA**: Muestra porcentaje y progreso durante `mdrun`
- **Validación de ligandos (v4.5)**: Verifica consistencia entre `.gro` y `.itp` del ligando antes de simular
- **Modo `--analysis-only`**: Re-ejecuta solo el análisis sobre una corrida existente terminada
- **14 pasos** (antes 13): Se añade validación del ligando como paso 2

### `plot_analysis.py`
- **Mapa de contactos (heatmap)**: Visualización residuo vs tiempo de contactos proteína-ligando
- **Running average**: Overlay de media móvil en todas las gráficas de series de tiempo
- **Histograma H-bonds**: Distribución y fracción acumulada de puentes de hidrógeno proteína-ligando
- **Reporte HTML (`--format html`)**: Genera reporte interactivo con gráficas embebidas
- **Exportación CSV**: `statistics_summary.csv` con todas las métricas para comparaciones batch

### `run_mmpbsa.sh`
- **Auto-detección de grupos**: Ya no requiere `--receptor`/`--ligand` si los nombres estándar existen
- **Ventana temporal (`--start-time`/`--end-time`)**: Selecciona qué porción de la trayectoria analizar
- **Convergencia energética**: Evalúa automáticamente si el cálculo MM-PBSA convergió (σ < 2 kcal/mol)
- **`--skip-recenter`**: Omite re-centrado PBC si la trayectoria ya está centrada

## Novedades v4.0

- **Equilibración mejorada**: NVT y NPT ahora usan 1 ns (antes 100 ps)
- **Fix tau_t**: Eliminada sobreescritura de `tau_t=0.1` que degradaba el muestreo canónico
- **Auto-detección de threads**: `NT` se detecta automáticamente con `nproc`, o se configura con `--nthreads`
- **`--maxwarn` configurable**: Controla el máximo de warnings de `grompp` desde CLI
- **`--prod-ns` acepta decimales**: Ej: `--prod-ns 0.5` para pruebas rápidas
- **Timestamps en logs**: Cada mensaje incluye hora para diagnóstico de rendimiento
- **Warnings de grompp visibles**: Se muestran automáticamente las warnings al compilar
- **Análisis extendido v4.0**:
  - Energía de interacción proteína-ligando (Coulomb + LJ)
  - Contactos nativos (d < 0.45 nm)
  - Evaluación de convergencia (RMSD block averaging)
- **`plot_analysis.py` v4.0**:
  - CLI con `argparse`: `--output-dir`, `--no-pca`, `--no-dccm`
  - B-factors calculados desde RMSF
  - Tabla resumen de estadísticas en PDF
  - Evaluación de convergencia visual
  - Stride automático para DCCM en trayectorias largas
  - Gráficas para energía de interacción y contactos nativos
- **`run_mmpbsa.sh` v4.0**:
  - Modo no interactivo completo con flags CLI
  - `set -E` para trap ERR correcto
  - `grep -oP` reemplazado por `awk` (portable)
  - Archivos temporales seguros con `mktemp`
  - Soporte de entropía con `--entropy` (Normal Mode Analysis)

## Requisitos

Dependencias validadas al inicio de `md_pipeline_mejorado.sh`:

- `gmx` o `gmx_mpi`
- `awk`
- `sed`
- `grep`
- `bc`

Dependencias recomendadas para análisis y reportes:

- `python3`
- `numpy`
- `matplotlib`

## Estructura esperada

```text
Automatizacion_proteina_ligando/
|- md_pipeline_mejorado.sh
|- plot_analysis.py
|- run_mmpbsa.sh
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
|  |- nvt.mdp          (1 ns equilibración)
|  |- npt.mdp          (1 ns equilibración)
|  `- md_prod.mdp
`- MD_RUN/             (salida automática)
```

## Ejecución interactiva

```bash
chmod +x md_pipeline_mejorado.sh
./md_pipeline_mejorado.sh
```

## Modo no interactivo (flags CLI)

El pipeline entra automáticamente a modo no interactivo cuando detecta los obligatorios:

- `--prot`
- `--lig`
- `--ff`

Si faltan obligatorios, el script hace fallback a modo interactivo (salvo que uses `--non-interactive`, donde falla de forma explícita).

### Ejemplo mínimo no interactivo

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
  --prod-ns 50 \
  --nthreads 8 \
  --maxwarn 2
```

También se soportan alias:

- `--box` (alias de `--box-type`)
- `--ion` (alias de `--ion-conc`)

### Parámetros nuevos v4.0

| Parámetro | Default | Descripción |
|-----------|---------|-------------|
| `--nthreads` | auto (`nproc`) | Número de threads para mdrun |
| `--maxwarn` | 1 | Máximo de warnings aceptadas por grompp |
| `--prod-ns` | 10 | Tiempo de producción (acepta decimales: 0.5, 1, 50) |

### Parámetros nuevos v4.5

| Parámetro | Default | Descripción |
|-----------|---------|-------------|
| `--gpu-id N` | auto | ID de GPU para mdrun |
| `--extend N` | — | Extender simulación por N ns (con `--resume`) |
| `--analysis-only DIR` | — | Re-ejecutar solo análisis sobre corrida existente |

### Ejemplos v4.5

```bash
# Usar GPU específica (HPC / Colab):
./md_pipeline_mejorado.sh --prot cyp51 --lig eugenol --ff charmm36.ff --gpu-id 0

# Extender simulación existente 50 ns más:
./md_pipeline_mejorado.sh --extend 50 --resume

# Re-ejecutar análisis sobre corrida anterior:
./md_pipeline_mejorado.sh --analysis-only MD_RUN/cyp51_eugenol_20260326_120000

# MM-PBSA con ventana temporal (solo últimos 50 ns):
./run_mmpbsa.sh --rundir MD_RUN/mi_corrida --calc gb_decomp \
  --start-time 50 --end-time 100

# MM-PBSA con auto-detección de grupos:
./run_mmpbsa.sh --rundir MD_RUN/mi_corrida --calc gb_decomp
```

## Archivo de configuración (`--config`)

Puedes cargar parámetros desde archivo:

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
NT=16
MAXWARN=1
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

## Batch mode (múltiples ligandos)

Para screening de múltiples ligandos:

```bash
for lig in M4-A M4-B M4-C; do
  ./md_pipeline_mejorado.sh \
    --prot caspasa9 \
    --lig "$lig" \
    --ff charmm36-jul2022.ff \
    --prod-ns 50
done
```

## MM-PB(GB)SA (`run_mmpbsa.sh`)

### Modo interactivo

```bash
./run_mmpbsa.sh
```

### Modo no interactivo (nuevo v4.0)

```bash
./run_mmpbsa.sh \
  --rundir MD_RUN/caspasa9_M4-A_20260224_120000 \
  --calc gb_decomp \
  --interval 5 \
  --salt 0.15 \
  --igb 5 \
  --receptor 1 \
  --ligand 13
```

#### Opciones de cálculo (`--calc`)

| Tipo | Descripción |
|------|-------------|
| `gb_only` | GB solamente (rápido, ~5-30 min) |
| `pb_only` | PB solamente (lento, ~1-4 hrs) |
| `gb_pb` | GB + PB |
| `gb_decomp` | GB + Descomposición por residuo (recomendado) |
| `gb_pb_decomp` | GB + PB + Descomposición (análisis completo) |

#### Entropía (nuevo v4.0)

```bash
./run_mmpbsa.sh \
  --rundir MD_RUN/mi_corrida \
  --calc gb_decomp \
  --receptor 1 --ligand 13 \
  --entropy
```

## Generación de reportes (`plot_analysis.py`)

### Uso básico

```bash
python3 plot_analysis.py MD_RUN/mi_corrida/
```

### Opciones CLI (nuevo v4.0)

```bash
python3 plot_analysis.py MD_RUN/mi_corrida/ --output-dir ./reportes/
python3 plot_analysis.py MD_RUN/mi_corrida/ --no-pca --no-dccm
python3 plot_analysis.py --help
```

### Gráficas incluidas

| Categoría | Métricas |
|-----------|----------|
| Minimización | Energía potencial |
| Equilibración | Temperatura, Presión, Densidad |
| Estabilidad | RMSD backbone, RMSD proteína, Radio de giro |
| Ligando | RMSD ligando, Distancia mínima, Contactos nativos, Energía de interacción |
| Flexibilidad | RMSF por residuo, **B-factors** (nuevo v4.0) |
| Superficie | SASA |
| Interacciones | H-bonds proteína, H-bonds prot-lig, Contactos por residuo |
| PCA | Eigenvalores, PC1, PC2, FEL |
| DCCM | Mapa de correlación dinámica |
| **Convergencia** | Evaluación visual de convergencia RMSD (nuevo v4.0) |
| **Resumen** | Tabla estadística de todas las métricas (nuevo v4.0) |

## Seguridad y robustez implementadas

- Se eliminó parseo frágil con `grep -P` / `grep -oP`; se usan alternativas con `awk`/`sed`.
- Se valida `bc` explícitamente al inicio para comparaciones de punto flotante.
- La carga de checkpoints evita `source` directo de contenido no confiable.
- La edición de `topol.top` ahora es encapsulada, idempotente y con verificación estricta post-edición.
- Manejo de errores reforzado con `set -E` + `trap ERR` + limpieza de temporales.
- Archivos temporales usan `mktemp` en lugar de rutas fijas en `/tmp/`.

## Reproducibilidad

Al final de cada corrida se genera `SUMMARY.txt` con:

- Fecha/hora exacta
- Hostname
- Binario y versión de GROMACS
- Info básica de CPU/GPU
- Parámetros clave de simulación
- Modo de ejecución (`normal` o `dry-run`)

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
  --prod-ns 100 \
  --nthreads 16

# MM-PBSA post-simulación (v4.0: modo no interactivo)
./run_mmpbsa.sh \
  --rundir MD_RUN/caspasa9_M4-A_* \
  --calc gb_decomp \
  --receptor 1 --ligand 13
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
  --prod-ns 100 \
  --nthreads 16
```

## Licencia

Uso libre para investigación académica.
