#!/bin/bash
set -Eeuo pipefail

#==========================================
# MM-PB(GB)SA ANÁLISIS DE ENERGÍA LIBRE v4.5
# Requiere: gmx_MMPBSA, AmberTools, GROMACS
# Uso: ./run_mmpbsa.sh [directorio_corrida_MD]
#      ./run_mmpbsa.sh --rundir DIR --calc TYPE [flags]
#==========================================

#==========================================
# COLORES PARA OUTPUT
#==========================================
readonly RED='\033[0;31m'
readonly GREEN='\033[0;32m'
readonly YELLOW='\033[1;33m'
readonly BLUE='\033[0;34m'
readonly CYAN='\033[0;36m'
readonly NC='\033[0m'

CURRENT_STAGE="inicialización"

# CLI inputs for non-interactive mode
INPUT_RUNDIR=""
INPUT_CALC_TYPE=""
INPUT_INTERVAL=""
INPUT_SALT=""
INPUT_IGB=""
INPUT_RECEPTOR=""
INPUT_LIGAND=""
INPUT_ENTROPY=false
INPUT_SKIP_RECENTER=false
INPUT_START_TIME=""
INPUT_END_TIME=""
NON_INTERACTIVE=false

#==========================================
# FUNCIONES AUXILIARES
#==========================================
log_step() {
    CURRENT_STAGE="$1"
    echo -e "\n${BLUE}=========================================${NC}"
    echo -e "${BLUE}>>> $1${NC}"
    echo -e "${BLUE}=========================================${NC}\n"
}

log_success() { echo -e "${GREEN}✓${NC} $1"; }
log_error()   { echo -e "${RED}❌ ERROR:${NC} $1" >&2; }
log_warning() { echo -e "${YELLOW}⚠${NC} $1"; }
log_info()    { echo -e "${CYAN}ℹ${NC} $1"; }

cleanup_on_error() {
    local exit_code=$?
    trap - ERR
    if [ $exit_code -ne 0 ]; then
        echo ""
        log_error "Script interrumpido durante: ${CURRENT_STAGE}"
        log_warning "Revisa los logs en: ${MMPBSA_DIR:-desconocido}/"
    fi
}
trap 'cleanup_on_error' ERR
trap 'cleanup_on_error' EXIT

print_mmpbsa_usage() {
    cat <<EOF
Uso:
  $0 [directorio_corrida_MD]
  $0 --rundir DIR --calc TYPE [flags]

Opciones (modo no interactivo):
  --rundir DIR          Directorio de la corrida MD
  --calc TYPE           Tipo: gb_only|pb_only|gb_pb|gb_decomp|gb_pb_decomp
  --interval N          Intervalo de frames (default: 5)
  --salt M              Concentración salina en M (default: 0.150)
  --igb N               Modelo GB: 1|2|5|7|8 (default: 5)
  --receptor IDX        Índice del grupo receptor (auto-detecta si omitido)
  --ligand IDX          Índice del grupo ligando (auto-detecta si omitido)
  --entropy             Incluir cálculo de entropía (Normal Mode Analysis)
  --start-time NS       Tiempo inicial en ns para análisis (v4.5)
  --end-time NS         Tiempo final en ns para análisis (v4.5)
  --skip-recenter       Omitir re-centrado PBC (si ya está centrada) (v4.5)
  --help, -h            Mostrar esta ayuda

Ejemplos v4.5:
  # Análisis de los últimos 50 ns solamente:
  $0 --rundir MD_RUN/mi_corrida --calc gb_decomp \
    --receptor 1 --ligand 13 --start-time 50 --end-time 100

  # Auto-detección de grupos (sin --receptor/--ligand):
  $0 --rundir MD_RUN/mi_corrida --calc gb_decomp

  # Saltar re-centrado:
  $0 --rundir MD_RUN/mi_corrida --calc gb_only --skip-recenter
EOF
}

parse_mmpbsa_args() {
    while [ $# -gt 0 ]; do
        case "$1" in
            --rundir)
                [ -n "${2:-}" ] || { log_error "Falta valor para --rundir"; exit 1; }
                INPUT_RUNDIR="$2"; shift 2 ;;
            --calc)
                [ -n "${2:-}" ] || { log_error "Falta valor para --calc"; exit 1; }
                INPUT_CALC_TYPE="$2"; shift 2 ;;
            --interval)
                [ -n "${2:-}" ] || { log_error "Falta valor para --interval"; exit 1; }
                INPUT_INTERVAL="$2"; shift 2 ;;
            --salt)
                [ -n "${2:-}" ] || { log_error "Falta valor para --salt"; exit 1; }
                INPUT_SALT="$2"; shift 2 ;;
            --igb)
                [ -n "${2:-}" ] || { log_error "Falta valor para --igb"; exit 1; }
                INPUT_IGB="$2"; shift 2 ;;
            --receptor)
                [ -n "${2:-}" ] || { log_error "Falta valor para --receptor"; exit 1; }
                INPUT_RECEPTOR="$2"; shift 2 ;;
            --ligand)
                [ -n "${2:-}" ] || { log_error "Falta valor para --ligand"; exit 1; }
                INPUT_LIGAND="$2"; shift 2 ;;
            --entropy)
                INPUT_ENTROPY=true; shift ;;
            --skip-recenter)
                INPUT_SKIP_RECENTER=true; shift ;;
            --start-time)
                [ -n "${2:-}" ] || { log_error "Falta valor para --start-time"; exit 1; }
                INPUT_START_TIME="$2"; shift 2 ;;
            --end-time)
                [ -n "${2:-}" ] || { log_error "Falta valor para --end-time"; exit 1; }
                INPUT_END_TIME="$2"; shift 2 ;;
            --help|-h)
                print_mmpbsa_usage; exit 0 ;;
            *)
                # Legacy positional argument
                if [ -z "$INPUT_RUNDIR" ] && [ -d "$1" ]; then
                    INPUT_RUNDIR="$1"
                fi
                shift ;;
        esac
    done

    # Check if enough args for non-interactive (v4.5: receptor/ligand optional)
    if [ -n "$INPUT_RUNDIR" ] && [ -n "$INPUT_CALC_TYPE" ]; then
        NON_INTERACTIVE=true
    fi
}

#==========================================
# ACTIVACIÓN DE ENTORNO CONDA
#==========================================
activate_conda_env() {
    # Si gmx_MMPBSA ya está disponible, no hacer nada
    if command -v gmx_MMPBSA &> /dev/null; then
        return 0
    fi

    log_info "gmx_MMPBSA no encontrado en PATH, buscando entornos conda..."

    # Detectar la instalación de conda/mamba
    local conda_exe=""
    if command -v conda &> /dev/null; then
        conda_exe="conda"
    elif command -v mamba &> /dev/null; then
        conda_exe="mamba"
    elif [ -f "$HOME/miniconda3/bin/conda" ]; then
        conda_exe="$HOME/miniconda3/bin/conda"
    elif [ -f "$HOME/anaconda3/bin/conda" ]; then
        conda_exe="$HOME/anaconda3/bin/conda"
    else
        return 1
    fi

    # Inicializar conda para este shell si no está inicializado
    local conda_base
    conda_base=$("$conda_exe" info --base 2>/dev/null)
    if [ -n "$conda_base" ] && [ -f "$conda_base/etc/profile.d/conda.sh" ]; then
        source "$conda_base/etc/profile.d/conda.sh"
    fi

    # Buscar entornos que contengan "mmpbsa" en el nombre
    local envs
    envs=$(conda env list 2>/dev/null | grep -i "mmpbsa" | awk '{print $1}' || true)

    if [ -z "$envs" ]; then
        # Buscar en todos los entornos si alguno tiene gmx_MMPBSA instalado
        envs=$(conda env list 2>/dev/null | grep -v "^#" | grep -v "^$" | awk '{print $1}' || true)
        for env_name in $envs; do
            local env_bin
            env_bin=$(conda run -n "$env_name" which gmx_MMPBSA 2>/dev/null || true)
            if [ -n "$env_bin" ]; then
                log_success "gmx_MMPBSA encontrado en entorno: $env_name"
                conda activate "$env_name"
                return 0
            fi
        done
        return 1
    fi

    # Si encontramos entornos con "mmpbsa" en el nombre
    local env_name
    env_name=$(echo "$envs" | head -1)

    log_success "Entorno conda encontrado: $env_name"
    log_info "Activando entorno: $env_name"
    conda activate "$env_name"

    if command -v gmx_MMPBSA &> /dev/null; then
        log_success "gmx_MMPBSA activado correctamente desde entorno: $env_name"
        return 0
    else
        log_warning "El entorno $env_name no contiene gmx_MMPBSA"
        return 1
    fi
}

#==========================================
# VALIDACIÓN DE DEPENDENCIAS
#==========================================
validate_dependencies() {
    log_step "Validando dependencias para MM-PB(GB)SA"

    # Intentar activar entorno conda si gmx_MMPBSA no está disponible
    if ! command -v gmx_MMPBSA &> /dev/null; then
        if ! activate_conda_env; then
            log_error "gmx_MMPBSA no encontrado y no se pudo activar automáticamente"
            echo ""
            echo "Opciones para resolver:"
            echo ""
            echo "  Opción 1: Activar manualmente antes de ejecutar el script:"
            echo "    conda activate gmx_mmpbsa_env"
            echo "    ./run_mmpbsa.sh"
            echo ""
            echo "  Opción 2: Instalar gmx_MMPBSA:"
            echo "    conda create -n gmx_mmpbsa_env python=3.11 -y"
            echo "    conda activate gmx_mmpbsa_env"
            echo "    conda install -c conda-forge ambertools -y"
            echo "    pip install gmx_MMPBSA"
            echo ""
            exit 1
        fi
    fi
    log_success "gmx_MMPBSA encontrado: $(which gmx_MMPBSA)"

    if ! command -v ante-MMPBSA.py &> /dev/null; then
        log_error "AmberTools no encontrado (ante-MMPBSA.py)"
        echo "  Instalar con: conda install -c conda-forge ambertools -y"
        exit 1
    fi
    log_success "AmberTools encontrado"

    # Detectar GROMACS
    if command -v gmx &> /dev/null; then
        GMX=gmx
    elif command -v gmx_mpi &> /dev/null; then
        GMX=gmx_mpi
    else
        log_error "GROMACS no encontrado"
        exit 1
    fi
    log_success "GROMACS encontrado: $(which $GMX)"
}

#==========================================
# SELECCIÓN DEL DIRECTORIO DE CORRIDA
#==========================================
select_run_directory() {
    log_step "Selección del directorio de corrida MD"

    if [ -n "${1:-}" ] && [ -d "$1" ]; then
        RUNDIR="$(cd "$1" && pwd)"
        log_success "Directorio proporcionado: $RUNDIR"
    else
        # Buscar corridas disponibles en MD_RUN/
        if [ -d "MD_RUN" ]; then
            echo "Corridas MD disponibles:"
            echo ""
            local i=1
            local runs=()
            for dir in MD_RUN/*/; do
                if [ -f "${dir}03_production/md.tpr" ]; then
                    runs+=("$(cd "$dir" && pwd)")
                    echo "  $i) $(basename "$dir")"
                    ((i++))
                fi
            done

            if [ ${#runs[@]} -eq 0 ]; then
                log_error "No se encontraron corridas MD completadas en MD_RUN/"
                exit 1
            fi

            echo ""
            echo "Seleccione la corrida [1-${#runs[@]}]:"
            read -r RUN_CHOICE

            if [[ "$RUN_CHOICE" =~ ^[0-9]+$ ]] && [ "$RUN_CHOICE" -ge 1 ] && [ "$RUN_CHOICE" -le ${#runs[@]} ]; then
                RUNDIR="${runs[$((RUN_CHOICE - 1))]}"
            else
                log_error "Selección inválida"
                exit 1
            fi
        else
            echo "Ruta al directorio de la corrida MD:"
            read -r RUNDIR
            RUNDIR="$(cd "$RUNDIR" && pwd)"
        fi
    fi

    # Validar que existen los archivos necesarios
    local required_files=(
        "$RUNDIR/03_production/md.tpr"
        "$RUNDIR/03_production/md.xtc"
        "$RUNDIR/00_setup/topol.top"
        "$RUNDIR/00_setup/index.ndx"
    )

    for f in "${required_files[@]}"; do
        if [ ! -f "$f" ]; then
            log_error "Archivo no encontrado: $f"
            exit 1
        fi
    done

    log_success "Archivos de producción validados"
}

#==========================================
# MENÚ DE TIPO DE CÁLCULO
#==========================================
select_calculation_type() {
    log_step "Configuración del cálculo MM-PB(GB)SA"

    echo -e "${CYAN}Tipo de cálculo:${NC}"
    echo ""
    echo "  1) GB solamente"
    echo "     └─ Rápido (~5-30 min). Ideal para screening de ligandos"
    echo ""
    echo "  2) PB solamente"
    echo "     └─ Lento (~1-4 hrs). Más preciso, para resultados finales"
    echo ""
    echo "  3) GB + PB"
    echo "     └─ Ambos métodos para comparar"
    echo ""
    echo "  4) GB + Descomposición por residuo"
    echo "     └─ Identifica qué residuos contribuyen más a la unión"
    echo ""
    echo "  5) GB + PB + Descomposición por residuo"
    echo "     └─ Análisis completo (más lento, más información)"
    echo ""
    echo -e "Seleccione el tipo de cálculo [1-5] (default: 4):"
    read -r CALC_CHOICE
    CALC_CHOICE=${CALC_CHOICE:-4}

    case $CALC_CHOICE in
        1) CALC_TYPE="gb_only"     ; CALC_DESC="GB solamente" ;;
        2) CALC_TYPE="pb_only"     ; CALC_DESC="PB solamente" ;;
        3) CALC_TYPE="gb_pb"       ; CALC_DESC="GB + PB" ;;
        5) CALC_TYPE="gb_pb_decomp"; CALC_DESC="GB + PB + Descomposición" ;;
        *) CALC_TYPE="gb_decomp"   ; CALC_DESC="GB + Descomposición por residuo" ;;
    esac
    log_success "Tipo de cálculo: $CALC_DESC"

    # Intervalo de frames
    echo ""
    echo -e "${CYAN}Intervalo de frames:${NC}"
    echo "  Un intervalo mayor = más rápido, menos precisión estadística"
    echo "  Recomendado: 1 (todos), 5 (rápido), 10 (muy rápido)"
    echo ""
    echo "Intervalo de frames [default: 5]:"
    read -r FRAME_INTERVAL
    FRAME_INTERVAL=${FRAME_INTERVAL:-5}
    log_success "Intervalo: cada $FRAME_INTERVAL frames"

    # Concentración salina
    echo ""
    echo -e "${CYAN}Concentración salina (M):${NC}"
    echo "  0.150 = fisiológica (recomendado)"
    echo ""
    echo "Concentración salina [default: 0.150]:"
    read -r SALT_CONC
    SALT_CONC=${SALT_CONC:-0.150}
    log_success "Concentración salina: $SALT_CONC M"

    # Modelo GB
    if [[ "$CALC_TYPE" == *"gb"* ]]; then
        echo ""
        echo -e "${CYAN}Modelo GB (Generalized Born):${NC}"
        echo "  1) igb=1 - HCT (Hawkins, Cramer, Truhlar)"
        echo "  2) igb=2 - OBC1 (Onufriev, Bashford, Case)"
        echo "  5) igb=5 - OBC2/GBneck2 (más preciso, recomendado)"
        echo "  7) igb=7 - GBneck (bueno para ácidos nucleicos)"
        echo "  8) igb=8 - GBneck2 alternativo"
        echo ""
        echo "Modelo GB [default: 5]:"
        read -r IGB_CHOICE
        IGB_CHOICE=${IGB_CHOICE:-5}

        case $IGB_CHOICE in
            1) IGB=1 ;;
            2) IGB=2 ;;
            7) IGB=7 ;;
            8) IGB=8 ;;
            *) IGB=5 ;;
        esac
        log_success "Modelo GB: igb=$IGB"
    fi
}

#==========================================
# DETECTAR GRUPOS DE ÍNDICE
#==========================================
detect_index_groups() {
    log_step "Detectando grupos en index.ndx"

    local index_file="$RUNDIR/00_setup/index.ndx"

    # Listar grupos disponibles
    local tmp_ndx
    tmp_ndx=$(mktemp /tmp/gmxmmpbsa_ndx_XXXXXX.ndx)
    groups_output=$(echo q | "$GMX" make_ndx -f "$RUNDIR/03_production/md.tpr" \
        -n "$index_file" -o "$tmp_ndx" 2>&1 || true)
    rm -f "$tmp_ndx"

    echo "$groups_output" | grep "^ *[0-9]"
    echo ""

    # Intentar detectar automáticamente
    local protein_idx=""
    local ligand_idx=""

    # Buscar grupo "Protein" (normalmente índice 1)
    protein_idx=$(echo "$groups_output" | awk '/Protein[[:space:]]/ && /^[[:space:]]*[0-9]/ {print $1; exit}')
    # Buscar grupo con "LIG" en el nombre
    ligand_idx=$(echo "$groups_output" | awk '/LIG[[:space:]]/ && /^[[:space:]]*[0-9]/ {print $1; exit}')

    if [ -n "$protein_idx" ] && [ -n "$ligand_idx" ]; then
        log_success "Auto-detectados: Proteína = grupo $protein_idx, Ligando = grupo $ligand_idx"
        echo ""
        echo "¿Usar estos grupos? [S/n]:"
        read -r USE_AUTO
        USE_AUTO=${USE_AUTO:-S}

        if [[ "$USE_AUTO" =~ ^[Ss]$ ]]; then
            GRP_RECEPTOR=$protein_idx
            GRP_LIGAND=$ligand_idx
            return
        fi
    else
        log_warning "No se pudieron auto-detectar los grupos"
    fi

    echo ""
    echo "Ingrese el índice del grupo RECEPTOR (proteína):"
    read -r GRP_RECEPTOR
    echo "Ingrese el índice del grupo LIGANDO:"
    read -r GRP_LIGAND

    log_success "Receptor: grupo $GRP_RECEPTOR, Ligando: grupo $GRP_LIGAND"
}

#==========================================
# RE-CENTRADO DE TRAYECTORIA (PBC FIX)
#==========================================
recenter_trajectory() {
    # v4.5: Skip recenter if requested
    if [ "$INPUT_SKIP_RECENTER" = true ]; then
        log_info "Re-centrado omitido (--skip-recenter)"
        TRAJ_CENTERED="$RUNDIR/03_production/md.xtc"
        return 0
    fi

    log_step "Re-centrando trayectoria (corrección PBC)"

    local traj_in="$RUNDIR/03_production/md.xtc"
    local tpr_in="$RUNDIR/03_production/md.tpr"
    local index_in="$RUNDIR/00_setup/index.ndx"
    TRAJ_CENTERED="$MMPBSA_DIR/md_center.xtc"

    log_info "Esto corrige artefactos de condiciones periódicas de contorno"
    log_info "para que el complejo proteína-ligando permanezca intacto."

    # Paso 1: Hacer whole
    log_info "Paso 1/3: Reconstruyendo moléculas (whole)..."
    echo "System" | "$GMX" trjconv \
        -f "$traj_in" \
        -s "$tpr_in" \
        -n "$index_in" \
        -o "$MMPBSA_DIR/_tmp_whole.xtc" \
        -pbc whole \
        2>&1 | tail -5

    # Paso 2: Centrar en el complejo
    log_info "Paso 2/3: Centrando complejo proteína-ligando..."
    echo -e "Protein_Ligand\nSystem" | "$GMX" trjconv \
        -f "$MMPBSA_DIR/_tmp_whole.xtc" \
        -s "$tpr_in" \
        -n "$index_in" \
        -o "$MMPBSA_DIR/_tmp_center.xtc" \
        -center \
        -pbc mol \
        -ur compact \
        2>&1 | tail -5

    # Paso 3: Cluster final
    log_info "Paso 3/3: Agrupando moléculas (cluster)..."
    echo -e "Protein_Ligand\nSystem" | "$GMX" trjconv \
        -f "$MMPBSA_DIR/_tmp_center.xtc" \
        -s "$tpr_in" \
        -n "$index_in" \
        -o "$TRAJ_CENTERED" \
        -pbc cluster \
        2>&1 | tail -5

    # Limpiar temporales
    rm -f "$MMPBSA_DIR/_tmp_whole.xtc" "$MMPBSA_DIR/_tmp_center.xtc"

    if [ -f "$TRAJ_CENTERED" ]; then
        local size_mb
        size_mb=$(du -m "$TRAJ_CENTERED" | cut -f1)
        log_success "Trayectoria re-centrada: md_center.xtc (${size_mb} MB)"
    else
        log_error "Falló el re-centrado de la trayectoria"
        log_warning "Continuando con trayectoria original (puede causar errores)"
        TRAJ_CENTERED="$traj_in"
    fi

    # v4.5: Recortar por ventana temporal si se especificó
    trim_trajectory_by_time
}

#==========================================
# GENERAR ARCHIVO mmpbsa.in
#==========================================
generate_mmpbsa_input() {
    log_step "Generando archivo de configuración mmpbsa.in"

    local mmpbsa_in="$MMPBSA_DIR/mmpbsa.in"

    # Sección general (siempre presente)
    cat > "$mmpbsa_in" <<EOF
# Archivo generado automáticamente por run_mmpbsa.sh
# Tipo de cálculo: $CALC_DESC
# Fecha: $(date)

&general
sys_name="MMPBSA_$(basename "$RUNDIR")",
startframe=1,
endframe=9999999,
interval=$FRAME_INTERVAL,
forcefields="leaprc.gaff2",
temperature=298.15,
PBRadii=3,
keep_files=2,
/
EOF

    # Sección GB
    if [[ "$CALC_TYPE" == *"gb"* ]]; then
        cat >> "$mmpbsa_in" <<EOF
&gb
igb=$IGB, saltcon=$SALT_CONC,
/
EOF
        log_success "Sección GB añadida (igb=$IGB, salt=$SALT_CONC)"
    fi

    # Sección PB
    if [[ "$CALC_TYPE" == *"pb"* ]]; then
        cat >> "$mmpbsa_in" <<EOF
&pb
istrng=$SALT_CONC,
indi=2.0,
exdi=80.0,
ipb=2,
inp=1,
radiopt=0,
fillratio=4.0,
scale=2.0,
linit=1000,
prbrad=1.4,
/
EOF
        log_success "Sección PB añadida (istrng=$SALT_CONC, ipb=2, inp=1)"
    fi

    # Sección descomposición por residuo
    if [[ "$CALC_TYPE" == *"decomp"* ]]; then
        cat >> "$mmpbsa_in" <<EOF
&decomp
idecomp=2, dec_verbose=3,
print_res="within 10"
/
EOF
        log_success "Sección de descomposición añadida (residuos dentro de 10 Å)"
    fi

    # Sección entropía (Normal Mode Analysis)
    if [ "$INPUT_ENTROPY" = true ]; then
        cat >> "$mmpbsa_in" <<EOF
&nmode
nmstartframe=1, nmendframe=100, nminterval=10,
maxcyc=10000, drms=0.001,
/
EOF
        log_success "Sección de entropía añadida (Normal Mode Analysis)"
    fi

    echo ""
    log_info "Archivo generado: $mmpbsa_in"
    echo "--- Contenido ---"
    cat "$mmpbsa_in"
    echo "-----------------"
}

#==========================================
# DETECTAR MODO DE EJECUCIÓN (MPI o serie)
#==========================================
detect_mpi_mode() {
    # Prioridad de CPUs: SLURM > nproc > 1
    local ncpus
    if [ -n "${SLURM_CPUS_PER_TASK:-}" ]; then
        ncpus="$SLURM_CPUS_PER_TASK"
        log_info "CPUs detectados desde SLURM: $ncpus"
    else
        ncpus=$(nproc 2>/dev/null || echo 1)
        log_info "CPUs detectados desde nproc: $ncpus"
    fi

    # Verificar disponibilidad de mpirun y mpi4py
    if command -v mpirun &>/dev/null && \
       python3 -c "from mpi4py import MPI" &>/dev/null 2>&1; then
        log_success "MPI disponible (mpirun + mpi4py) — usando $ncpus procesos"
        MMPBSA_USE_MPI=true
        MMPBSA_NCPUS="$ncpus"
    else
        log_warning "MPI no disponible (mpirun o mpi4py ausentes) — modo serie"
        if ! command -v mpirun &>/dev/null; then
            log_info "  Motivo: mpirun no encontrado en PATH"
        else
            log_info "  Motivo: mpi4py no instalado (conda install -c conda-forge mpi4py)"
        fi
        MMPBSA_USE_MPI=false
        MMPBSA_NCPUS=1
    fi
}

#==========================================
# EJECUTAR gmx_MMPBSA
#==========================================
run_mmpbsa() {
    log_step "Ejecutando gmx_MMPBSA ($CALC_DESC)"

    cd "$MMPBSA_DIR" || exit 1

    # Detectar modo MPI
    detect_mpi_mode

    local start_time
    start_time=$(date +%s)

    echo -e "${YELLOW}Esto puede tardar desde minutos hasta horas según el cálculo...${NC}"
    echo -e "${YELLOW}Progreso detallado en: $MMPBSA_DIR/gmx_mmpbsa.log${NC}\n"

    # Argumentos comunes
    local mmpbsa_args=(
        -O
        -i mmpbsa.in
        -cs "$RUNDIR/03_production/md.tpr"
        -ci "$RUNDIR/00_setup/index.ndx"
        -cg "$GRP_RECEPTOR" "$GRP_LIGAND"
        -ct "$TRAJ_CENTERED"
        -cp "$RUNDIR/00_setup/topol.top"
        -o FINAL_RESULTS_MMPBSA.dat
        -eo FINAL_RESULTS_MMPBSA.csv
        -do FINAL_DECOMP_MMPBSA.dat
        -deo FINAL_DECOMP_MMPBSA.csv
    )

    if [ "$MMPBSA_USE_MPI" = true ]; then
        log_info "Lanzando con MPI: mpirun -np $MMPBSA_NCPUS gmx_MMPBSA MPI"
        mpirun -np "$MMPBSA_NCPUS" gmx_MMPBSA MPI "${mmpbsa_args[@]}" \
            > gmx_mmpbsa.log 2>&1
    else
        log_info "Lanzando en modo serie (1 proceso)"
        gmx_MMPBSA "${mmpbsa_args[@]}" \
            > gmx_mmpbsa.log 2>&1
    fi

    local end_time
    end_time=$(date +%s)
    local elapsed=$(( end_time - start_time ))
    local minutes=$(( elapsed / 60 ))
    local seconds=$(( elapsed % 60 ))

    if [ "$MMPBSA_USE_MPI" = true ]; then
        log_success "gmx_MMPBSA (MPI x${MMPBSA_NCPUS}) completado en ${minutes}m ${seconds}s"
    else
        log_success "gmx_MMPBSA (serie) completado en ${minutes}m ${seconds}s"
    fi
}

#==========================================
# GENERAR RESUMEN DE RESULTADOS
#==========================================
generate_results_summary() {
    log_step "Procesando resultados"

    cd "$MMPBSA_DIR" || exit 1

    # Mostrar resultados principales
    if [ -f "FINAL_RESULTS_MMPBSA.dat" ]; then
        echo -e "\n${CYAN}=========================================${NC}"
        echo -e "${CYAN}  RESULTADOS DE ENERGÍA LIBRE DE UNIÓN${NC}"
        echo -e "${CYAN}=========================================${NC}\n"

        cat FINAL_RESULTS_MMPBSA.dat

        log_success "Resultados completos: FINAL_RESULTS_MMPBSA.dat"
    fi

    if [ -f "FINAL_RESULTS_MMPBSA.csv" ]; then
        log_success "Resultados CSV: FINAL_RESULTS_MMPBSA.csv"
    fi

    if [ -f "FINAL_DECOMP_MMPBSA.dat" ]; then
        log_success "Descomposición por residuo: FINAL_DECOMP_MMPBSA.dat"
    fi

    if [ -f "FINAL_DECOMP_MMPBSA.csv" ]; then
        log_success "Descomposición CSV: FINAL_DECOMP_MMPBSA.csv"
    fi

    # Generar resumen
    cat > "$MMPBSA_DIR/SUMMARY_MMPBSA.txt" <<EOF
====================================
RESUMEN DE ANÁLISIS MM-PB(GB)SA
====================================
Corrida MD:           $(basename "$RUNDIR")
Tipo de cálculo:      $CALC_DESC
Fecha:                $(date)
Directorio:           $MMPBSA_DIR

PARÁMETROS:
- Intervalo de frames: $FRAME_INTERVAL
- Concentración salina: $SALT_CONC M
$([ "${IGB:-}" ] && echo "- Modelo GB:           igb=$IGB")
- Grupo receptor:      $GRP_RECEPTOR
- Grupo ligando:       $GRP_LIGAND

ARCHIVOS GENERADOS:
- FINAL_RESULTS_MMPBSA.dat   : Resultados principales (ΔG_bind)
- FINAL_RESULTS_MMPBSA.csv   : Resultados en formato CSV
$([ -f "$MMPBSA_DIR/FINAL_DECOMP_MMPBSA.dat" ] && echo "- FINAL_DECOMP_MMPBSA.dat   : Descomposición por residuo")
$([ -f "$MMPBSA_DIR/FINAL_DECOMP_MMPBSA.csv" ] && echo "- FINAL_DECOMP_MMPBSA.csv   : Descomposición CSV")
- mmpbsa.in                  : Archivo de configuración usado
- gmx_mmpbsa.log             : Log completo de ejecución

INTERPRETACIÓN:
- ΔG_bind negativo = unión favorable
- ΔG_bind < -10 kcal/mol = unión muy fuerte
- ΔG_bind entre -5 y -10 = unión moderada
- ΔG_bind > -5 = unión débil

VISUALIZACIÓN INTERACTIVA:
  gmx_MMPBSA_ana -p _GMXMMPBSA_info

PRÓXIMOS PASOS:
1. Revisar FINAL_RESULTS_MMPBSA.dat para ΔG_bind global
2. Si hay descomposición, identificar residuos clave en FINAL_DECOMP_MMPBSA.dat
3. Usar gmx_MMPBSA_ana para visualización gráfica interactiva
4. Comparar con datos experimentales si están disponibles
EOF

    cat "$MMPBSA_DIR/SUMMARY_MMPBSA.txt"
    log_success "Resumen guardado: SUMMARY_MMPBSA.txt"

    echo ""
    echo -e "${GREEN}=========================================${NC}"
    echo -e "${GREEN}  Para visualización interactiva ejecuta:${NC}"
    echo -e "${GREEN}  gmx_MMPBSA_ana -p $MMPBSA_DIR/_GMXMMPBSA_info${NC}"
    echo -e "${GREEN}=========================================${NC}"
}

#==========================================
# RECORTE TEMPORAL DE TRAYECTORIA (v4.5)
#==========================================
trim_trajectory_by_time() {
    local start_time="${INPUT_START_TIME:-}"
    local end_time="${INPUT_END_TIME:-}"

    # Si no se especificó ninguna ventana, no hacer nada
    if [ -z "$start_time" ] && [ -z "$end_time" ]; then
        return 0
    fi

    log_step "Recortando trayectoria por ventana temporal (v4.5)"

    local traj_input="$TRAJ_CENTERED"
    local traj_trimmed="$MMPBSA_DIR/md_trimmed.xtc"
    local trjconv_args=""

    # Convertir ns a ps para gmx trjconv
    if [ -n "$start_time" ]; then
        local start_ps
        start_ps=$(awk "BEGIN{printf \"%d\", $start_time * 1000}")
        trjconv_args="$trjconv_args -b $start_ps"
        log_info "Tiempo inicial: $start_time ns ($start_ps ps)"
    fi

    if [ -n "$end_time" ]; then
        local end_ps
        end_ps=$(awk "BEGIN{printf \"%d\", $end_time * 1000}")
        trjconv_args="$trjconv_args -e $end_ps"
        log_info "Tiempo final: $end_time ns ($end_ps ps)"
    fi

    echo "System" | "$GMX" trjconv \
        -f "$traj_input" \
        -s "$RUNDIR/03_production/md.tpr" \
        -n "$RUNDIR/00_setup/index.ndx" \
        -o "$traj_trimmed" \
        $trjconv_args \
        2>&1 | tail -5

    if [ -f "$traj_trimmed" ]; then
        local size_mb
        size_mb=$(du -m "$traj_trimmed" | cut -f1)
        TRAJ_CENTERED="$traj_trimmed"
        log_success "Trayectoria recortada: md_trimmed.xtc (${size_mb} MB)"
    else
        log_warning "No se pudo recortar; usando trayectoria completa"
    fi
}

#==========================================
# AUTO-DETECCIÓN DE GRUPOS (v4.5)
#==========================================
auto_detect_groups() {
    local index_file="$RUNDIR/00_setup/index.ndx"

    local tmp_ndx
    tmp_ndx=$(mktemp /tmp/gmxmmpbsa_ndx_XXXXXX.ndx)
    local groups_output
    groups_output=$(echo q | "$GMX" make_ndx -f "$RUNDIR/03_production/md.tpr" \
        -n "$index_file" -o "$tmp_ndx" 2>&1 || true)
    rm -f "$tmp_ndx"

    # Buscar grupo "Protein"
    local protein_idx
    protein_idx=$(echo "$groups_output" | awk '/Protein[[:space:]]/ && /^[[:space:]]*[0-9]/ {print $1; exit}')
    # Buscar grupo "LIG"
    local ligand_idx
    ligand_idx=$(echo "$groups_output" | awk '/LIG[[:space:]]/ && /^[[:space:]]*[0-9]/ {print $1; exit}')

    if [ -n "$protein_idx" ] && [ -n "$ligand_idx" ]; then
        GRP_RECEPTOR=$protein_idx
        GRP_LIGAND=$ligand_idx
        log_success "Auto-detectados: Proteína = grupo $GRP_RECEPTOR, Ligando = grupo $GRP_LIGAND"
        return 0
    fi

    return 1
}

#==========================================
# REPORTE DE CONVERGENCIA ENERGÉTICA (v4.5)
#==========================================
analyze_mmpbsa_convergence() {
    log_step "Evaluando convergencia del cálculo MM-PB(GB)SA (v4.5)"

    local csv_file="$MMPBSA_DIR/FINAL_RESULTS_MMPBSA.csv"
    if [ ! -f "$csv_file" ]; then
        log_warning "No se encontró CSV de resultados para evaluar convergencia"
        return
    fi

    # Extraer ΔG values from CSV (buscar columna TOTAL)
    local dg_values
    dg_values=$(awk -F',' '
        NR==1 {
            for(i=1;i<=NF;i++) {
                if($i ~ /TOTAL/) { col=i; break }
            }
            next
        }
        col && NF>=col { print $col }
    ' "$csv_file" 2>/dev/null | head -200)

    if [ -z "$dg_values" ]; then
        log_warning "No se pudieron extraer valores de \u0394G del CSV"
        return
    fi

    # Calcular media y desviación estándar con awk
    local stats
    stats=$(echo "$dg_values" | awk '
        { sum+=$1; sumsq+=$1*$1; n++ }
        END {
            mean=sum/n
            std=sqrt(sumsq/n - mean*mean)
            printf "%.4f %.4f %d", mean, std, n
        }
    ')

    local mean std nframes
    mean=$(echo "$stats" | awk '{print $1}')
    std=$(echo "$stats" | awk '{print $2}')
    nframes=$(echo "$stats" | awk '{print $3}')

    echo ""
    echo -e "${CYAN}Evaluación de convergencia (v4.5):${NC}"
    echo -e "  Frames analizados: $nframes"
    echo -e "  \u0394G medio:          $mean kcal/mol"
    echo -e "  Desviación est.:  $std kcal/mol"

    # Evaluar convergencia (\u03c3 < 2 kcal/mol = convergido)
    local converged
    converged=$(awk "BEGIN{print ($std < 2.0) ? 1 : 0}")
    if [ "$converged" -eq 1 ]; then
        echo -e "  Estado:            ${GREEN}\u2713 CONVERGIDO (\u03c3 < 2 kcal/mol)${NC}"
    else
        echo -e "  Estado:            ${YELLOW}\u26a0 CONVERGENCIA DUDOSA (\u03c3 \u2265 2 kcal/mol)${NC}"
        echo -e "  ${YELLOW}Considere: usar m\u00e1s frames (--interval 1) o simulaci\u00f3n m\u00e1s larga${NC}"
    fi
    echo ""
}

#==========================================
# FUNCIÓN PRINCIPAL
#==========================================
main() {
    echo -e "${BLUE}=========================================${NC}"
    echo -e "${BLUE}  ANÁLISIS MM-PB(GB)SA${NC}"
    echo -e "${BLUE}  Energía libre de unión proteína-ligando${NC}"
    echo -e "${BLUE}=========================================${NC}\n"

    parse_mmpbsa_args "$@"
    validate_dependencies

    if [ "$NON_INTERACTIVE" = true ]; then
        # Non-interactive mode
        RUNDIR="$(cd "$INPUT_RUNDIR" && pwd)"
        log_success "Directorio: $RUNDIR"

        # Validate files
        for f in "$RUNDIR/03_production/md.tpr" "$RUNDIR/03_production/md.xtc" \
                 "$RUNDIR/00_setup/topol.top" "$RUNDIR/00_setup/index.ndx"; do
            [ -f "$f" ] || { log_error "Archivo no encontrado: $f"; exit 1; }
        done

        CALC_TYPE="$INPUT_CALC_TYPE"
        case "$CALC_TYPE" in
            gb_only) CALC_DESC="GB solamente" ;;
            pb_only) CALC_DESC="PB solamente" ;;
            gb_pb) CALC_DESC="GB + PB" ;;
            gb_decomp) CALC_DESC="GB + Descomposición por residuo" ;;
            gb_pb_decomp) CALC_DESC="GB + PB + Descomposición" ;;
            *) log_error "Tipo de cálculo inválido: $CALC_TYPE"; exit 1 ;;
        esac

        FRAME_INTERVAL="${INPUT_INTERVAL:-5}"
        SALT_CONC="${INPUT_SALT:-0.150}"
        IGB="${INPUT_IGB:-5}"

        # v4.5: Auto-detección de grupos si no se proporcionaron
        if [ -n "$INPUT_RECEPTOR" ] && [ -n "$INPUT_LIGAND" ]; then
            GRP_RECEPTOR="$INPUT_RECEPTOR"
            GRP_LIGAND="$INPUT_LIGAND"
        else
            log_info "Grupos no especificados, intentando auto-detección (v4.5)..."
            if ! auto_detect_groups; then
                log_error "No se pudieron auto-detectar grupos. Use --receptor y --ligand"
                exit 1
            fi
        fi

        log_success "Modo no interactivo: $CALC_DESC"
    else
        # v4.5: Preguntar ventana temporal en modo interactivo
        select_run_directory "${INPUT_RUNDIR:-${1:-}}"
        select_calculation_type
        detect_index_groups

        if [ -z "${INPUT_START_TIME:-}" ] && [ -z "${INPUT_END_TIME:-}" ]; then
            echo ""
            echo -e "${CYAN}Ventana temporal de análisis (v4.5):${NC}"
            echo "  1) Analizar TODA la trayectoria (default)"
            echo "  2) Seleccionar ventana temporal (ej: solo los últimos 50 ns)"
            echo ""
            echo "Selección [1-2] (default: 1):"
            read -r TIME_CHOICE
            TIME_CHOICE=${TIME_CHOICE:-1}

            if [ "$TIME_CHOICE" = "2" ]; then
                echo "Tiempo inicial en ns (dejar vacío = desde el inicio):"
                read -r INPUT_START_TIME
                echo "Tiempo final en ns (dejar vacío = hasta el final):"
                read -r INPUT_END_TIME
                [ -n "$INPUT_START_TIME" ] && log_info "Desde: $INPUT_START_TIME ns"
                [ -n "$INPUT_END_TIME" ] && log_info "Hasta: $INPUT_END_TIME ns"
            fi
        fi
    fi

    # Crear directorio de trabajo
    MMPBSA_DIR="$RUNDIR/05_mmpbsa_$(date +%Y%m%d_%H%M%S)"
    mkdir -p "$MMPBSA_DIR"
    log_success "Directorio de trabajo: $MMPBSA_DIR"

    generate_mmpbsa_input
    recenter_trajectory
    run_mmpbsa
    generate_results_summary

    # v4.5: Análisis de convergencia
    analyze_mmpbsa_convergence

    log_step "¡ANÁLISIS MM-PB(GB)SA FINALIZADO CON ÉXITO!"
    echo -e "${GREEN}Resultados en: $MMPBSA_DIR${NC}\n"
}

# Ejecutar
main "$@"
