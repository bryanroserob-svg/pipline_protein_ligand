# GROMACS MD Pipeline Automatizado v4.6 (Proteína + Ligando)

Pipeline para simulaciones de dinámica molecular de complejos proteína-ligando con GROMACS.

## Novedades v4.6 — Robustez HPC

### `md_pipeline_mejorado.sh`
- **`--resume-dir <ruta>`**: Reanuda o extiende una simulación **específica** por ruta, sin depender de la búsqueda automática. Permite elegir exactamente qué corrida continuar aunque haya varias incompletas.
- **Auto-selección en HPC**: En modo `--non-interactive`, si hay múltiples checkpoints, selecciona automáticamente el más reciente (evita que `read` bloquee en nodos SLURM donde `stdin = /dev/null`).
- **Validación de RUNDIR**: Al cargar un checkpoint, verifica que el directorio de la corrida es accesible en el nodo actual (crítico en sistemas NFS/Lustre).
- **Checkpoint atómico**: El archivo `.checkpoint` se escribe primero en un `.tmp` y luego se reemplaza con `mv` (evita checkpoint corrupto si el job recibe SIGKILL durante la escritura).
- **Checkpoint completo**: Guarda `MAXWARN`, `NT` y `GPU_ID` para que el resume use exactamente los mismos parámetros.
- **SLURM-aware threads**: Usa `$SLURM_CPUS_PER_TASK` para auto-detectar threads en lugar de `nproc` (que devuelve los CPUs del nodo completo, no los del job).
- **Logs HPC limpios**: Los colores ANSI se desactivan automáticamente cuando el output se redirige a un archivo de log (sin caracteres basura).
- **Cleanup de procesos monitor**: El PID del proceso de monitoreo se registra para limpieza correcta si el job es cancelado.

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
| `--extend N` | — | Extender simulación por N ns adicionales |
| `--analysis-only DIR` | — | Re-ejecutar solo análisis sobre corrida existente |

### Parámetros nuevos v4.6 (HPC)

| Parámetro | Default | Descripción |
|-----------|---------|-------------|
| `--resume-dir DIR` | — | Ruta explícita a la corrida a reanudar o extender. Implica `--resume` automáticamente. Permite elegir exactamente qué simulación continuar, incluso si está marcada como COMPLETADA (para `--extend`). |

### Ejemplos v4.5 y v4.6

```bash
# Usar GPU específica (HPC / Colab):
./md_pipeline_mejorado.sh --prot cyp51 --lig eugenol --ff charmm36.ff --gpu-id 0

# Reanudar la simulación incompleta más reciente (auto-detecta):
./md_pipeline_mejorado.sh --resume --non-interactive

# Reanudar una simulación específica que se cayó (v4.6):
./md_pipeline_mejorado.sh \
  --resume-dir MD_RUN/caspasa9_M4-A_20260607_153000 \
  --nthreads 16 --gpu-id 0 --non-interactive

# Extender una simulación específica 100 ns más (v4.6):
./md_pipeline_mejorado.sh \
  --extend 100 \
  --resume-dir MD_RUN/caspasa9_M4-A_20260607_153000 \
  --nthreads 16 --gpu-id 0 --non-interactive

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

### Reanudación automática (busca la más reciente incompleta)

```bash
./md_pipeline_mejorado.sh --resume --non-interactive
```

### Reanudación explícita de una corrida específica (v4.6)

Útil cuando hay múltiples simulaciones incompletas y quieres elegir cuál continuar:

```bash
./md_pipeline_mejorado.sh \
  --resume-dir MD_RUN/caspasa9_M4-A_20260607_153000 \
  --nthreads 16 \
  --gpu-id 0 \
  --non-interactive
```

El checkpoint se guarda en cada paso exitoso y se valida con parser seguro antes de cargarse.
En HPC, si el job muere por OOM KILL o walltime, el checkpoint queda en el último paso
exitoso y la simulación puede reanudarse exactamente desde ahí.

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
- **(v4.6)** Checkpoint escrito atómicamente con `mv` — no puede quedar corrupto por SIGKILL.
- **(v4.6)** `RUNDIR` del checkpoint se valida antes de reanudar — falla rápido si el filesystem no está montado.
- **(v4.6)** Threads auto-detectados desde `$SLURM_CPUS_PER_TASK` para respetar el job allocation.

## Reproducibilidad

Al final de cada corrida se genera `SUMMARY.txt` con:

- Fecha/hora exacta
- Hostname
- Binario y versión de GROMACS
- Info básica de CPU/GPU
- Parámetros clave de simulación
- Modo de ejecución (`normal` o `dry-run`)

## Ejemplos para SLURM/SGE

### SLURM — Simulación nueva

```bash
#!/bin/bash
#SBATCH --job-name=md_caspasa9_m4a
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --gres=gpu:a100-sxm4-40gb:1
#SBATCH --time=48:00:00
#SBATCH --output=md_out_%j.log
#SBATCH --error=md_err_%j.log

enroot start --root --rw --mount "$HOME:$HOME" gromacs2026-cuda13 bash -c "
    source /opt/gromacs2026.2/bin/GMXRC
    cd \"\$HOME/pipline_protein_ligand\"
    ./md_pipeline_mejorado.sh \\
      --prot caspasa9 --lig M4-A --ff charmm36-jul2022.ff \\
      --box dodecahedron --water tip3p --ion 0 \\
      --prod-ns 100 --nthreads 16 --gpu-id 0 --non-interactive
"
```

### SLURM — Reanudar simulación que se cayó (auto-detecta la más reciente)

```bash
#!/bin/bash
#SBATCH --job-name=resume_#NOMBRE_DINAMICA
#SBATCH --partition=gpu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=16 --mem=32G
#SBATCH --gres=gpu:a100-sxm4-40gb:1
#SBATCH --time=48:00:00
#SBATCH --output=md_out_%j.log
#SBATCH --error=md_err_%j.log

enroot start --root --rw --mount "$HOME:$HOME" gromacs2026-cuda13 bash -c "
    source /opt/gromacs2026.2/bin/GMXRC
    cd \"\$HOME/pipline_protein_ligand\"
    ./md_pipeline_mejorado.sh \\
      --resume-dir MD_RUN/#NOMBRE_DINAMICA \\
      --nthreads 16 --gpu-id 0 --non-interactive
"
```

### SLURM — Extender simulación específica (v4.6)

```bash
#!/bin/bash
#SBATCH --job-name=extend_#NOMBRE_DINAMICA
#SBATCH --partition=gpu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=16 --mem=32G
#SBATCH --gres=gpu:a100-sxm4-40gb:1
#SBATCH --time=48:00:00
#SBATCH --output=md_out_%j.log
#SBATCH --error=md_err_%j.log

enroot start --root --rw --mount "$HOME:$HOME" gromacs2026-cuda13 bash -c "
    source /opt/gromacs2026.2/bin/GMXRC
    cd \"\$HOME/pipline_protein_ligand\"
    ./md_pipeline_mejorado.sh \\
      --extend 100 \\
      --resume-dir MD_RUN/#NOMBRE_DINAMICA \\
      --nthreads 16 --gpu-id 0 --non-interactive
"
```

### SLURM — Re-generar análisis y resumen

```bash
#!/bin/bash
#SBATCH --job-name=analysis_#NOMBRE_DINAMICA
#SBATCH --partition=gpu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=16 --mem=32G
#SBATCH --gres=gpu:a100-sxm4-40gb:1
#SBATCH --time=4:00:00
#SBATCH --output=md_out_%j.log
#SBATCH --error=md_err_%j.log

enroot start --root --rw --mount "$HOME:$HOME" gromacs2026-cuda13 bash -c "
    source /opt/gromacs2026.2/bin/GMXRC
    cd \"\$HOME/pipline_protein_ligand\"
    ./md_pipeline_mejorado.sh \\
      --analysis-only MD_RUN/#NOMBRE_DINAMICA
"
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
