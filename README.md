# 🧬 GROMACS MD Pipeline Automatizado

Pipeline automatizado para simulaciones de Dinámica Molecular de complejos proteína-ligando usando **GROMACS 2025.x**.

## 📋 Requisitos

| Software | Versión | Uso |
|---|---|---|
| **GROMACS** | 2025.x | Simulación MD |
| **Python 3** | 3.8+ | Gráficas de análisis |
| **matplotlib + numpy** | - | `pip install matplotlib numpy` |
| **gmx_MMPBSA** *(opcional)* | 1.6+ | Energía libre de unión |
| **AmberTools** *(opcional)* | 23+ | Requerido por gmx_MMPBSA |

## 📁 Estructura del Proyecto

```
gromacs_automatizado/
├── md_pipeline_mejorado.sh    # Pipeline principal de MD
├── run_mmpbsa.sh              # Análisis MM-PB(GB)SA
├── plot_analysis.py           # Generador de gráficas PDF
├── README.md
├── proteinas/                 # Carpeta de proteínas
│   └── mi_proteina/
│       ├── proteina.gro
│       ├── topol.top
│       └── posre.itp
├── ligandos/                  # Carpeta de ligandos
│   └── mi_ligando/
│       ├── ligando.gro
│       ├── ligando.itp
│       └── ligando.prm        # (opcional, para CHARMM)
├── mdp/                       # Archivos de parámetros
│   ├── ions.mdp
│   ├── em.mdp
│   ├── nvt.mdp
│   ├── npt.mdp
│   └── md_prod.mdp
├── charmm36.ff/               # (opcional) Force field
└── MD_RUN/                    # Resultados (generado automáticamente)
    └── proteina_ligando_YYYYMMDD_HHMMSS/
        ├── 00_setup/
        ├── 01_minimization/
        ├── 02_equilibration/
        ├── 03_production/
        ├── 04_analysis/
        ├── 05_mmpbsa_*/       # (si se ejecuta run_mmpbsa.sh)
        ├── logs/
        ├── mdp_used/
        ├── .checkpoint         # Estado para resume
        ├── SUMMARY.txt
        └── REPORT_MD.pdf      # (si se ejecuta plot_analysis.py)
```

## 🚀 Uso Rápido

### 1. Pipeline de MD (simulación completa)

```bash
chmod +x md_pipeline_mejorado.sh
./md_pipeline_mejorado.sh
```

El script preguntará interactivamente:
- Nombre de la proteína y ligando
- Tipo de caja (cúbica, triclínica, dodecaedro, octaedro)
- Distancia a bordes (nm)
- Modelo de agua (TIP3P, SPC, SPC/E, TIP4P, TIP5P)

**Flujo del pipeline (13 pasos con checkpoint automático):**

```
Paso  1: Copiar archivos de entrada
Paso  2: Generar restraints del ligando
Paso  3: Actualizar topología
Paso  4: Ensamblar complejo proteína-ligando
Paso  5: Solvatar sistema
Paso  6: Neutralizar con iones
Paso  7: Crear grupos de índice
Paso  8: Minimización de energía
Paso  9: Equilibración NVT (calentamiento)
Paso 10: Equilibración NPT (presión)
Paso 11: Producción (sin restricciones)
Paso 12: Análisis post-producción
Paso 13: Generar resumen
```

### 2. Reanudar una corrida interrumpida

Si el pipeline falla o se interrumpe, puedes reanudarlo **sin empezar de cero**:

```bash
./md_pipeline_mejorado.sh --resume
```

El sistema de checkpoints:
- **Guarda progreso** automáticamente después de cada paso
- **No crea una carpeta nueva** — reanuda en la misma
- **No pide datos otra vez** — carga proteína, ligando, etc. del checkpoint
- **Si hay varias corridas incompletas**, muestra un menú para elegir cuál reanudar

### 3. Gráficas automáticas

```bash
python3 plot_analysis.py                          # Menú interactivo
python3 plot_analysis.py MD_RUN/mi_corrida/       # Directo
```

Genera `REPORT_MD.pdf` con ~17 gráficas organizadas por secciones:
- Minimización, Equilibración, Estabilidad, Ligando, Flexibilidad, Interacciones, PCA

### 4. Energía libre de unión (MM-PB/GBSA)

```bash
conda activate gmx_mmpbsa_env   # (se activa automáticamente si existe)
chmod +x run_mmpbsa.sh
./run_mmpbsa.sh
```

Opciones de cálculo:
1. **GB** — Rápido (~5-30 min)
2. **PB** — Lento (~1-4 hrs), más preciso
3. **GB + PB** — Ambos métodos
4. **GB + Descomposición** — Identifica residuos clave
5. **GB + PB + Descomposición** — Análisis completo

## 📊 Archivos de Análisis Generados

### En `04_analysis/`

| Archivo | Descripción |
|---|---|
| `rmsd_backbone.xvg` | Estabilidad del backbone vs tiempo |
| `rmsd_protein.xvg` | RMSD de toda la proteína |
| `rmsd_ligand.xvg` | RMSD del ligando (¿se mueve del sitio?) |
| `gyrate.xvg` | Radio de giro (compactación) |
| `rmsf_residue.xvg` | Flexibilidad por residuo |
| `mindist_prot_lig.xvg` | Distancia mínima prot-lig (disociación) |
| `sasa_protein.xvg` | Superficie accesible al solvente |
| `hbond_protein_ligand.xvg` | Puentes de H proteína↔ligando |
| `eigenval.xvg` | Eigenvalores del PCA |
| `proj_pc1.xvg` / `proj_pc2.xvg` | Componentes principales |
| `clusters.pdb` | Conformaciones representativas |
| `rmsd_matrix.xpm` | Mapa de calor RMSD vs RMSD |
| `contacts_prot_lig.xvg` | Contactos por residuo |
| `dssp.dat` | Estructura secundaria vs tiempo |

### En `05_mmpbsa_*/`

| Archivo | Descripción |
|---|---|
| `FINAL_RESULTS_MMPBSA.dat` | ΔG_bind (energía libre de unión) |
| `FINAL_DECOMP_MMPBSA.dat` | Contribución por residuo |
| `_GMXMMPBSA_info` | Para `gmx_MMPBSA_ana` (visualización interactiva) |

## 🔧 Preparación de Archivos de Entrada

### Proteína
```bash
gmx pdb2gmx -f proteina.pdb -o proteina.gro -water tip3p
# Copiar proteina.gro, topol.top, posre.itp a proteinas/nombre/
```

### Ligando
1. Obtener la estructura del ligando (`.mol2`)
2. Generar topología en [CGenFF](https://cgenff.umaryland.edu/) (ParamChem):
   - Subir el `.mol2` a CGenFF
   - Descargar el `.str` (stream file)
   - Convertir a formato GROMACS con `cgenff_charmm2gmx_py3.py` o herramientas similares
3. **Verificar manualmente** la topología generada (penaltys, cargas, tipos de átomos)
4. Copiar `ligando.gro`, `ligando.itp` (y `ligando.prm`) a `ligandos/nombre/`

### Archivos MDP
Los parámetros `tc-grps`, `tau_t`, `ref_t` y `define` se ajustan automáticamente por el pipeline.

## ⚙️ Instalación de gmx_MMPBSA

```bash
conda create -n gmx_mmpbsa_env python=3.11 -y
conda activate gmx_mmpbsa_env
conda install -c conda-forge ambertools -y
sudo apt install libopenmpi-dev openmpi-bin   # Solo OpenMPI (NO instalar MPICH también)
pip install gmx_MMPBSA
gmx_MMPBSA --help                            # Verificar
```

## 📝 Notas Importantes

- **GROMACS 2025.x**: Usa `gmx dssp` (integrado), no requiere `mkdssp` externo
- **Grupos de termostato**: Se crean automáticamente como `Protein_Ligand` y `Solvent`
- **Restricciones**: NVT y NPT usan `-DPOSRES -DPOSRES_LIG`; producción sin restricciones
- **`-maxwarn 1`**: Se permite máximo 1 warning en grompp
- **Checkpoint**: Se guarda en `$RUNDIR/.checkpoint` después de cada paso exitoso
- **Conda**: `run_mmpbsa.sh` detecta y activa automáticamente entornos con "mmpbsa" en el nombre

## 📄 Licencia

Uso libre para investigación académica.
