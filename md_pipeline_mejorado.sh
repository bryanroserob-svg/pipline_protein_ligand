#!/bin/bash
set -Eeuo pipefail

# GROMACS MD Pipeline Automatizado v4.5 (Proteína + Ligando)
#==========================================
# CONFIGURACIÓN
#==========================================
NT=""
GMX=""
USE_MPI=false
MDRUN=""
GPU_ID=""

readonly BASE_PROT=proteinas
readonly BASE_LIG=ligandos
readonly BASE_MDP=mdp
readonly WORKDIR=MD_RUN
readonly INITIAL_DIR="$(pwd)"
TOTAL_STEPS=13

# Flag para indicar si el force field es local (copiado) o nativo de GROMACS
FF_IS_LOCAL=false

#==========================================
# COLORES PARA OUTPUT
#==========================================
readonly RED='\033[0;31m'
readonly GREEN='\033[0;32m'
readonly YELLOW='\033[1;33m'
readonly BLUE='\033[0;34m'
readonly CYAN='\033[0;36m'
readonly NC='\033[0m' # No Color

# Variable para rastrear etapa actual (para trap)
CURRENT_STAGE="inicialización"
NON_INTERACTIVE=false
FORCE_NON_INTERACTIVE=false
RESUME_MODE=false
RESUME_FROM=1
DRY_RUN=false
SHOW_HELP=false
ANALYSIS_ONLY=false
EXTEND_NS=""
TEMP_FILES=()

# Parámetros opcionales para modo no interactivo
INPUT_PROT=""
INPUT_LIG=""
INPUT_FF=""
INPUT_BOX_TYPE=""
INPUT_BOX_DIST=""
INPUT_WATER_MODEL=""
INPUT_ION_CONC=""
INPUT_PROD_NS=""
INPUT_NT=""
INPUT_MAXWARN=""
INPUT_GPU_ID=""

#==========================================
# FUNCIONES AUXILIARES
#==========================================
log_step() {
    CURRENT_STAGE="$1"
    local ts
    ts=$(date '+%H:%M:%S')
    echo -e "\n${BLUE}[$ts] =========================================${NC}"
    echo -e "${BLUE}[$ts] >>> $1${NC}"
    echo -e "${BLUE}[$ts] =========================================${NC}\n"
}

log_success() {
    echo -e "${GREEN}✓${NC} [$(date '+%H:%M:%S')] $1"
}

log_error() {
    echo -e "${RED}❌ ERROR:${NC} [$(date '+%H:%M:%S')] $1" >&2
}

log_warning() {
    echo -e "${YELLOW}⚠${NC} [$(date '+%H:%M:%S')] $1"
}

log_info() {
    echo -e "${CYAN}ℹ${NC} [$(date '+%H:%M:%S')] $1"
}

create_dir() {
    if [ ! -d "$1" ]; then
        mkdir -p "$1"
        log_success "Carpeta creada: $1"
    fi
}

validate_command() {
    if ! command -v "$1" &> /dev/null; then
        log_error "$1 no está instalado o no está en PATH"
        exit 1
    fi
}

init_gmx() {
    if command -v gmx &> /dev/null; then
        GMX=gmx
        USE_MPI=false
    elif command -v gmx_mpi &> /dev/null; then
        GMX=gmx_mpi
        USE_MPI=true
    else
        log_error "No se encontró GROMACS (gmx o gmx_mpi)"
        exit 1
    fi

    # Auto-detect threads if not set
    if [ -z "$NT" ]; then
        NT=$(nproc 2>/dev/null || echo 4)
    fi

    if [ "$USE_MPI" = true ]; then
        # 1 rank MPI + N threads OpenMP: configuración óptima para 1 GPU / 1 nodo
        # Evita que gmx_mpi lance N ranks automáticamente (cada uno reserva su propia RAM)
        MDRUN="$GMX mdrun -ntmpi 1 -ntomp $NT"
    else
        MDRUN="$GMX mdrun -nt $NT"
    fi

    # GPU support (v4.5)
    if [ -n "$GPU_ID" ]; then
        MDRUN="$MDRUN -gpu_id $GPU_ID"
        log_info "GPU seleccionada: gpu_id=$GPU_ID"
    fi
}

# Wrapper simplificado para ejecutar comandos GMX
run_gmx() {
    "$GMX" "$@"
}

#==========================================
# ADAPTAR MDP SEGÚN FORCE FIELD SELECCIONADO
# CHARMM36 usa force-switch, otros FF no
#==========================================
adapt_mdp_for_forcefield() {
    local mdp_dir="$1"

    # Detectar si es CHARMM desde FF_DIR (ej: charmm36-jul2022.ff, charmm27.ff)
    local is_charmm=false
    case "${FF_DIR:-}" in
        charmm*) is_charmm=true ;;
    esac

    if [ "$is_charmm" = true ]; then
        log_info "Force field CHARMM detectado: MDP sin modificaciones (force-switch correcto)"
        return 0
    fi

    log_step "Adaptando archivos MDP para force field ${FF_DIR:-desconocido} (no-CHARMM)"

    for mdp_file in "$mdp_dir"/*.mdp; do
        [ -f "$mdp_file" ] || continue

        # Reemplazar vdw-modifier: force-switch → Potential-shift
        if grep -q 'force-switch' "$mdp_file"; then
            sed -i 's/vdw-modifier.*=.*force-switch/vdw-modifier             = Potential-shift/' "$mdp_file"
            # Eliminar rvdw-switch (no se usa con Potential-shift)
            sed -i '/^rvdw-switch/d' "$mdp_file"
            # Para AMBER/GROMOS/OPLS: DispCorr = EnerPres
            sed -i 's/^DispCorr.*=.*no/DispCorr                 = EnerPres/' "$mdp_file"
            # Actualizar comentario
            sed -i 's/; VAN DER WAALS (CHARMM36: force-switch)/; VAN DER WAALS (adaptado para '"${FF_DIR:-}"')/' "$mdp_file"
            log_success "$(basename "$mdp_file") adaptado: Potential-shift + DispCorr=EnerPres"
        fi
    done
}

#==========================================
# TRAP PARA LIMPIEZA EN CASO DE ERROR
#==========================================
register_temp_file() {
    TEMP_FILES+=("$1")
}

cleanup_temp_files() {
    local tmp
    for tmp in "${TEMP_FILES[@]}"; do
        if [ -n "$tmp" ] && [ -f "$tmp" ]; then
            rm -f "$tmp"
        fi
    done
}

cleanup_on_error() {
    local line_no=$1
    local cmd=$2
    local exit_code=$3

    trap - ERR

    echo ""
    log_error "Script interrumpido durante: ${CURRENT_STAGE}"
    log_error "Comando fallido (línea ${line_no}): ${cmd}"
    log_warning "Código de salida: ${exit_code}"
    cleanup_temp_files
    log_warning "Archivos parciales en: ${RUNDIR:-desconocido}"
    log_warning "Revisa los logs en: ${RUNDIR:-desconocido}/logs/"
}

cleanup_on_exit() {
    local exit_code=$1
    if [ "$exit_code" -eq 0 ]; then
        log_success "Finalización limpia del pipeline"
        cleanup_temp_files
    fi
}

trap 'cleanup_on_error ${LINENO} "${BASH_COMMAND}" "$?"' ERR
trap 'cleanup_on_exit "$?"' EXIT

#==========================================
# VALIDACIÓN DE DEPENDENCIAS
#==========================================
validate_dependencies() {
    log_step "Validando dependencias del sistema"
    local deps=(awk sed grep bc)
    local dep
    for dep in "${deps[@]}"; do
        validate_command "$dep"
    done
    log_success "Dependencias validadas: ${deps[*]}"
    log_success "GROMACS encontrado: $(command -v "$GMX") $([ "$USE_MPI" = true ] && echo "(MPI)" || echo "(serial)")"
}

is_int() {
    [[ "$1" =~ ^[0-9]+$ ]]
}

is_float() {
    [[ "$1" =~ ^[-+]?[0-9]+([.][0-9]+)?$ || "$1" =~ ^[-+]?[.][0-9]+$ ]]
}

is_positive_float() {
    is_float "$1" && (( $(echo "$1 > 0" | bc -l) ))
}

is_non_negative_float() {
    is_float "$1" && (( $(echo "$1 >= 0" | bc -l) ))
}

extract_checkpoint_value() {
    local ckpt_file="$1"
    local key="$2"
    awk -F'=' -v key="$key" '
        $1==key {
            value=$2
            sub(/^"/, "", value)
            sub(/"$/, "", value)
            print value
            exit
        }
    ' "$ckpt_file"
}

parse_checkpoint_file() {
    local ckpt_file="$1"
    local line key value

    while IFS= read -r line || [ -n "$line" ]; do
        [[ "$line" =~ ^[[:space:]]*$ ]] && continue
        [[ "$line" =~ ^[[:space:]]*# ]] && continue

        if [[ "$line" =~ ^([A-Z_]+)=\"([^\"]*)\"$ ]]; then
            key="${BASH_REMATCH[1]}"
            value="${BASH_REMATCH[2]}"
        else
            log_error "Checkpoint con formato inválido: $line"
            return 1
        fi

        case "$key" in
            LAST_STEP)
                if ! is_int "$value"; then
                    log_error "LAST_STEP inválido en checkpoint: $value"
                    return 1
                fi
                LAST_STEP="$value"
                ;;
            LAST_STEP_NAME) LAST_STEP_NAME="$value" ;;
            PROT) PROT="$value" ;;
            LIG) LIG="$value" ;;
            BOX_TYPE) BOX_TYPE="$value" ;;
            BOX_DIST) BOX_DIST="$value" ;;
            ION_CONC) ION_CONC="$value" ;;
            WATER_MODEL) WATER_MODEL="$value" ;;
            WATER_FILE) WATER_FILE="$value" ;;
            PROD_NS)
                if ! is_int "$value"; then
                    log_error "PROD_NS inválido en checkpoint: $value"
                    return 1
                fi
                PROD_NS="$value"
                ;;
            PROD_NSTEPS)
                if ! is_int "$value"; then
                    log_error "PROD_NSTEPS inválido en checkpoint: $value"
                    return 1
                fi
                PROD_NSTEPS="$value"
                ;;
            RUNDIR) RUNDIR="$value" ;;
            FF_DIR) FF_DIR="$value" ;;
            FF_IS_LOCAL)
                case "$value" in
                    true|false) FF_IS_LOCAL="$value" ;;
                    *) log_error "FF_IS_LOCAL inválido en checkpoint: $value"; return 1 ;;
                esac
                ;;
            DRY_RUN)
                case "$value" in
                    true|false) DRY_RUN="$value" ;;
                    *) log_error "DRY_RUN inválido en checkpoint: $value"; return 1 ;;
                esac
                ;;
            TOTAL_STEPS)
                if ! is_int "$value"; then
                    log_error "TOTAL_STEPS inválido en checkpoint: $value"
                    return 1
                fi
                TOTAL_STEPS="$value"
                ;;
            *)
                log_error "Clave no permitida en checkpoint: $key"
                return 1
                ;;
        esac
    done < "$ckpt_file"

    [ -n "${LAST_STEP:-}" ] || { log_error "Checkpoint sin LAST_STEP"; return 1; }
    [ -n "${RUNDIR:-}" ] || { log_error "Checkpoint sin RUNDIR"; return 1; }
    return 0
}

print_usage() {
    cat <<EOF
Uso:
  $0 [--resume|-r]
  $0 [--config <archivo.conf>] [flags]

Opciones:
  --prot <nombre>         Carpeta dentro de proteinas/
  --lig <nombre>          Carpeta dentro de ligandos/
  --ff <ff_dir>           Force field (ej: charmm36-jul2022.ff)
  --box <tipo>            Alias de --box-type
  --box-type <tipo>       cubic|triclinic|dodecahedron|octahedron
  --box-dist <nm>         Distancia de caja (default: 1.2)
  --water <modelo>        tip3p|spc|spce|tip4p|tip5p
  --ion <M>               Alias de --ion-conc
  --ion-conc <M>          Concentración NaCl (default: 0)
  --prod-ns <ns>          Tiempo de producción en ns (acepta decimales, default: 10)
  --nthreads <N>          Número de threads para mdrun (default: auto)
  --maxwarn <N>           Máximo de warnings para grompp (default: 1)
  --gpu-id <N>            ID de GPU para mdrun (default: auto-detección)
  --extend <ns>           Extender simulación existente por N ns adicionales
  --analysis-only <dir>   Re-ejecutar solo el análisis sobre una corrida existente
  --config <archivo>      Archivo KEY=VALUE (sin espacios)
  --non-interactive       Fuerza ejecución sin prompts
  --dry-run               Ejecuta preparación + grompp, sin mdrun
  --resume, -r            Reanuda una corrida incompleta
  --help, -h              Mostrar esta ayuda

Modo no interactivo automático:
  Si pasas --prot --lig --ff el pipeline entra en modo no interactivo.
  Si faltan obligatorios y NO usas --non-interactive, hace fallback a prompts.

Ejemplos v4.5:
  # Usar GPU específica:
  $0 --prot caspasa9 --lig M4-A --ff charmm36-jul2022.ff --gpu-id 0

  # Extender simulación existente 50 ns:
  $0 --extend 50 --resume

  # Re-ejecutar análisis sobre corrida anterior:
  $0 --analysis-only MD_RUN/caspasa9_M4-A_20260224_120000

Formato de --config:
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
  GPU_ID=0
EOF
}

load_config_file() {
    local cfg_file="$1"
    if [ ! -f "$cfg_file" ]; then
        log_error "Archivo de configuración no encontrado: $cfg_file"
        exit 1
    fi

    while IFS='=' read -r key value; do
        [ -z "${key:-}" ] && continue
        [[ "$key" =~ ^[[:space:]]*# ]] && continue
        key=$(echo "$key" | awk '{$1=$1; print}')
        value=$(echo "$value" | awk '{$1=$1; print}')
        case "$key" in
            PROT) INPUT_PROT="$value" ;;
            LIG) INPUT_LIG="$value" ;;
            FF) INPUT_FF="$value" ;;
            BOX_TYPE) INPUT_BOX_TYPE="$value" ;;
            BOX_DIST) INPUT_BOX_DIST="$value" ;;
            WATER_MODEL) INPUT_WATER_MODEL="$value" ;;
            ION_CONC) INPUT_ION_CONC="$value" ;;
            PROD_NS) INPUT_PROD_NS="$value" ;;
            NT) INPUT_NT="$value" ;;
            MAXWARN) INPUT_MAXWARN="$value" ;;
            GPU_ID) INPUT_GPU_ID="$value" ;;
            "") ;;
            *) log_warning "Clave desconocida en config: $key (ignorada)" ;;
        esac
    done < "$cfg_file"
}

parse_args() {
    while [ $# -gt 0 ]; do
        case "$1" in
            --resume|-r)
                RESUME_MODE=true
                shift
                ;;
            --config)
                [ -n "${2:-}" ] || { log_error "Falta archivo para --config"; exit 1; }
                load_config_file "$2"
                shift 2
                ;;
            --non-interactive)
                FORCE_NON_INTERACTIVE=true
                NON_INTERACTIVE=true
                shift
                ;;
            --dry-run)
                DRY_RUN=true
                shift
                ;;
            --prot)
                [ -n "${2:-}" ] || { log_error "Falta valor para --prot"; exit 1; }
                INPUT_PROT="$2"
                shift 2
                ;;
            --lig)
                [ -n "${2:-}" ] || { log_error "Falta valor para --lig"; exit 1; }
                INPUT_LIG="$2"
                shift 2
                ;;
            --ff)
                [ -n "${2:-}" ] || { log_error "Falta valor para --ff"; exit 1; }
                INPUT_FF="$2"
                shift 2
                ;;
            --box|--box-type)
                [ -n "${2:-}" ] || { log_error "Falta valor para --box-type"; exit 1; }
                INPUT_BOX_TYPE="$2"
                shift 2
                ;;
            --box-dist)
                [ -n "${2:-}" ] || { log_error "Falta valor para --box-dist"; exit 1; }
                INPUT_BOX_DIST="$2"
                shift 2
                ;;
            --water)
                [ -n "${2:-}" ] || { log_error "Falta valor para --water"; exit 1; }
                INPUT_WATER_MODEL="$2"
                shift 2
                ;;
            --ion|--ion-conc)
                [ -n "${2:-}" ] || { log_error "Falta valor para --ion"; exit 1; }
                INPUT_ION_CONC="$2"
                shift 2
                ;;
            --prod-ns)
                [ -n "${2:-}" ] || { log_error "Falta valor para --prod-ns"; exit 1; }
                INPUT_PROD_NS="$2"
                shift 2
                ;;
            --nthreads)
                [ -n "${2:-}" ] || { log_error "Falta valor para --nthreads"; exit 1; }
                INPUT_NT="$2"
                shift 2
                ;;
            --maxwarn)
                [ -n "${2:-}" ] || { log_error "Falta valor para --maxwarn"; exit 1; }
                INPUT_MAXWARN="$2"
                shift 2
                ;;
            --gpu-id)
                [ -n "${2:-}" ] || { log_error "Falta valor para --gpu-id"; exit 1; }
                INPUT_GPU_ID="$2"
                shift 2
                ;;
            --extend)
                [ -n "${2:-}" ] || { log_error "Falta valor para --extend"; exit 1; }
                EXTEND_NS="$2"
                shift 2
                ;;
            --analysis-only)
                [ -n "${2:-}" ] || { log_error "Falta directorio para --analysis-only"; exit 1; }
                ANALYSIS_ONLY=true
                ANALYSIS_ONLY_DIR="$2"
                shift 2
                ;;
            --help|-h)
                SHOW_HELP=true
                shift
                ;;
            *)
                log_error "Opción desconocida: $1"
                print_usage
                exit 1
                ;;
        esac
    done
}

has_any_non_interactive_input() {
    [ -n "$INPUT_PROT" ] || [ -n "$INPUT_LIG" ] || [ -n "$INPUT_FF" ] || [ -n "$INPUT_BOX_TYPE" ] || [ -n "$INPUT_BOX_DIST" ] || [ -n "$INPUT_WATER_MODEL" ] || [ -n "$INPUT_ION_CONC" ] || [ -n "$INPUT_PROD_NS" ] || [ -n "$INPUT_NT" ] || [ -n "$INPUT_MAXWARN" ]
}

has_required_non_interactive_input() {
    [ -n "$INPUT_PROT" ] && [ -n "$INPUT_LIG" ] && [ -n "$INPUT_FF" ]
}

set_water_file_from_model() {
    case "$WATER_MODEL" in
        tip3p|spc|spce) WATER_FILE="spc216.gro" ;;
        tip4p) WATER_FILE="tip4p.gro" ;;
        tip5p) WATER_FILE="tip5p.gro" ;;
        *) log_error "WATER_MODEL inválido: $WATER_MODEL"; exit 1 ;;
    esac
}

validate_runtime_parameters() {
    case "$BOX_TYPE" in
        cubic|triclinic|dodecahedron|octahedron) ;;
        *) log_error "BOX_TYPE inválido: $BOX_TYPE (use cubic|triclinic|dodecahedron|octahedron)"; exit 1 ;;
    esac

    if ! is_positive_float "$BOX_DIST"; then
        log_error "BOX_DIST inválido: '$BOX_DIST'. Debe ser un float > 0"
        exit 1
    fi

    if ! is_non_negative_float "$ION_CONC"; then
        log_error "ION_CONC inválido: '$ION_CONC'. Debe ser un float >= 0"
        exit 1
    fi

    if ! is_positive_float "$PROD_NS"; then
        log_error "PROD_NS inválido: '$PROD_NS'. Debe ser un número positivo"
        exit 1
    fi

    set_water_file_from_model
    PROD_NSTEPS=$(awk "BEGIN{printf \"%d\", $PROD_NS * 500000}")

    # Aplicar NT si fue proporcionado
    if [ -n "${INPUT_NT:-}" ]; then
        if is_int "$INPUT_NT" && [ "$INPUT_NT" -gt 0 ]; then
            NT="$INPUT_NT"
        else
            log_error "NT inválido: '$INPUT_NT'. Debe ser entero positivo"
            exit 1
        fi
    fi

    # Aplicar MAXWARN si fue proporcionado
    MAXWARN="${INPUT_MAXWARN:-2}"
    if ! is_int "$MAXWARN" || [ "$MAXWARN" -lt 0 ]; then
        log_error "MAXWARN inválido: '$MAXWARN'. Debe ser entero >= 0"
        exit 1
    fi

    # Aplicar GPU_ID si fue proporcionado (v4.5)
    if [ -n "${INPUT_GPU_ID:-}" ]; then
        if is_int "$INPUT_GPU_ID" && [ "$INPUT_GPU_ID" -ge 0 ]; then
            GPU_ID="$INPUT_GPU_ID"
        else
            log_error "GPU_ID inválido: '$INPUT_GPU_ID'. Debe ser entero >= 0"
            exit 1
        fi
    fi
}

resolve_execution_mode() {
    if has_required_non_interactive_input; then
        NON_INTERACTIVE=true
        apply_non_interactive_inputs
        return
    fi

    if [ "$FORCE_NON_INTERACTIVE" = true ]; then
        log_error "Modo no interactivo forzado, pero faltan obligatorios: --prot --lig --ff"
        exit 1
    fi

    if has_any_non_interactive_input || [ "$NON_INTERACTIVE" = true ]; then
        log_warning "Parámetros CLI/config incompletos; cambiando a modo interactivo"
    fi

    NON_INTERACTIVE=false
    get_user_input
}

#==========================================
# INPUT DEL USUARIO CON VALIDACIÓN
#==========================================
get_user_input() {
    echo -e "${BLUE}=========================================${NC}"
    echo -e "${BLUE}  SIMULACIÓN DE DINÁMICA MOLECULAR${NC}"
    echo -e "${BLUE}  (Proteína + Ligando)${NC}"
    echo -e "${BLUE}=========================================${NC}\n"

    # ========================================
    # Selección de proteína (auto-detección)
    # ========================================
    echo -e "\n${CYAN}═══════════════════════════════════════════${NC}"
    echo -e "${CYAN}  SELECCIÓN DE PROTEÍNA${NC}"
    echo -e "${CYAN}═══════════════════════════════════════════${NC}\n"

    local PROT_DIRS=()
    local PROT_PDBS=()

    if [ ! -d "$INITIAL_DIR/$BASE_PROT" ]; then
        log_error "No se encontró la carpeta '$BASE_PROT/' en el directorio actual"
        log_info "Crea la carpeta y coloca subcarpetas con archivos .pdb dentro"
        exit 1
    fi

    for prot_dir in "$INITIAL_DIR/$BASE_PROT"/*/; do
        if [ -d "$prot_dir" ]; then
            local pdb_found
            pdb_found=$(find "$prot_dir" -maxdepth 1 -name "*.pdb" -type f | head -1)
            if [ -n "$pdb_found" ]; then
                local dir_name
                dir_name=$(basename "$prot_dir")
                local pdb_name
                pdb_name=$(basename "$pdb_found")
                PROT_DIRS+=("$dir_name")
                PROT_PDBS+=("$pdb_name")
            fi
        fi
    done

    if [ ${#PROT_DIRS[@]} -eq 0 ]; then
        log_error "No se encontraron proteínas (carpetas con .pdb) en '$BASE_PROT/'"
        exit 1
    fi

    echo "  Proteínas disponibles en $BASE_PROT/:"
    echo ""
    for i in "${!PROT_DIRS[@]}"; do
        local idx=$((i + 1))
        echo "   ${idx}) ${PROT_DIRS[$i]}  (${PROT_PDBS[$i]})"
    done
    echo ""
    echo -e "  ${YELLOW}Seleccione la proteína [1-${#PROT_DIRS[@]}]:${NC}"
    read -r PROT_CHOICE

    if [[ "$PROT_CHOICE" =~ ^[0-9]+$ ]] && [ "$PROT_CHOICE" -ge 1 ] && [ "$PROT_CHOICE" -le ${#PROT_DIRS[@]} ]; then
        local sel=$((PROT_CHOICE - 1))
        PROT="${PROT_DIRS[$sel]}"
    else
        log_error "Selección inválida"
        exit 1
    fi

    log_success "Proteína seleccionada: $PROT"

    # ========================================
    # Selección de ligando (auto-detección)
    # ========================================
    echo -e "\n${CYAN}═══════════════════════════════════════════${NC}"
    echo -e "${CYAN}  SELECCIÓN DE LIGANDO${NC}"
    echo -e "${CYAN}═══════════════════════════════════════════${NC}\n"

    local LIG_DIRS=()
    local LIG_FILES=()

    if [ ! -d "$INITIAL_DIR/$BASE_LIG" ]; then
        log_error "No se encontró la carpeta '$BASE_LIG/' en el directorio actual"
        log_info "Crea la carpeta y coloca subcarpetas con archivos ligando.itp y ligando.gro dentro"
        exit 1
    fi

    for lig_dir in "$INITIAL_DIR/$BASE_LIG"/*/; do
        if [ -d "$lig_dir" ]; then
            # Buscar ligando.itp y ligando.gro
            if [ -f "$lig_dir/ligando.itp" ] && [ -f "$lig_dir/ligando.gro" ]; then
                local dir_name
                dir_name=$(basename "$lig_dir")
                local extras=""
                [ -f "$lig_dir/ligando.prm" ] && extras=" +prm"
                [ -f "$lig_dir/ligando.pdb" ] && extras="${extras} +pdb"
                LIG_DIRS+=("$dir_name")
                LIG_FILES+=("ligando.itp, ligando.gro${extras}")
            fi
        fi
    done

    if [ ${#LIG_DIRS[@]} -eq 0 ]; then
        log_error "No se encontraron ligandos (carpetas con ligando.itp + ligando.gro) en '$BASE_LIG/'"
        exit 1
    fi

    echo "  Ligandos disponibles en $BASE_LIG/:"
    echo ""
    for i in "${!LIG_DIRS[@]}"; do
        local idx=$((i + 1))
        echo "   ${idx}) ${LIG_DIRS[$i]}  (${LIG_FILES[$i]})"
    done
    echo ""
    echo -e "  ${YELLOW}Seleccione el ligando [1-${#LIG_DIRS[@]}]:${NC}"
    read -r LIG_CHOICE

    if [[ "$LIG_CHOICE" =~ ^[0-9]+$ ]] && [ "$LIG_CHOICE" -ge 1 ] && [ "$LIG_CHOICE" -le ${#LIG_DIRS[@]} ]; then
        local sel=$((LIG_CHOICE - 1))
        LIG="${LIG_DIRS[$sel]}"
    else
        log_error "Selección inválida"
        exit 1
    fi

    log_success "Ligando seleccionado: $LIG"

    # ========================================
    # Selección de Force Field (auto-detección)
    # ========================================
    echo -e "\n${CYAN}═══════════════════════════════════════════${NC}"
    echo -e "${CYAN}  SELECCIÓN DE CAMPO DE FUERZA${NC}"
    echo -e "${CYAN}═══════════════════════════════════════════${NC}\n"

    # Detectar force fields locales (*.ff en el directorio raíz)
    local LOCAL_FF_DIRS=()
    local LOCAL_FF_NAMES=()

    for ff_dir in "$INITIAL_DIR"/*.ff; do
        if [ -d "$ff_dir" ]; then
            local dir_basename
            dir_basename=$(basename "$ff_dir")
            local ff_desc="$dir_basename"
            if [ -f "$ff_dir/forcefield.doc" ]; then
                ff_desc=$(grep -m1 '[^ ]' "$ff_dir/forcefield.doc" 2>/dev/null || echo "$dir_basename")
            fi
            LOCAL_FF_DIRS+=("$dir_basename")
            LOCAL_FF_NAMES+=("$ff_desc")
        fi
    done

    # Detectar force fields nativos de GROMACS
    local NATIVE_FF_DIRS=()
    local NATIVE_FF_NAMES=()
    local gmx_topdir=""
    local gmx_datadir=""
    gmx_datadir=$($GMX --version 2>/dev/null | awk -F': *' '/^Data prefix/ {print $2; exit}' || true)
    if [ -n "$gmx_datadir" ]; then
        gmx_topdir="$gmx_datadir/share/gromacs/top"
    elif [ -d "/usr/local/gromacs/share/gromacs/top" ]; then
        gmx_topdir="/usr/local/gromacs/share/gromacs/top"
    elif [ -d "/usr/share/gromacs/top" ]; then
        gmx_topdir="/usr/share/gromacs/top"
    fi

    if [ -n "$gmx_topdir" ] && [ -d "$gmx_topdir" ]; then
        for ff_dir in "$gmx_topdir"/*.ff; do
            if [ -d "$ff_dir" ]; then
                local dir_basename
                dir_basename=$(basename "$ff_dir")
                # Evitar duplicados con los locales
                local is_dup=false
                for local_ff in "${LOCAL_FF_DIRS[@]+"${LOCAL_FF_DIRS[@]}"}"; do
                    [ "$local_ff" = "$dir_basename" ] && is_dup=true && break
                done
                [ "$is_dup" = true ] && continue

                local ff_desc="$dir_basename"
                if [ -f "$ff_dir/forcefield.doc" ]; then
                    ff_desc=$(grep -m1 '[^ ]' "$ff_dir/forcefield.doc" 2>/dev/null || echo "$dir_basename")
                fi
                NATIVE_FF_DIRS+=("$dir_basename")
                NATIVE_FF_NAMES+=("$ff_desc")
            fi
        done
    fi

    local total_ff=$(( ${#LOCAL_FF_DIRS[@]} + ${#NATIVE_FF_DIRS[@]} ))

    if [ "$total_ff" -eq 0 ]; then
        log_error "No se encontraron force fields (ni locales ni de GROMACS)"
        exit 1
    fi

    # Mostrar force fields locales
    local ff_idx=1
    if [ ${#LOCAL_FF_DIRS[@]} -gt 0 ]; then
        echo -e "  ${GREEN}── Force fields locales (descargados) ──${NC}"
        for i in "${!LOCAL_FF_DIRS[@]}"; do
            echo "   ${ff_idx}) ${LOCAL_FF_NAMES[$i]}"
            echo -e "      ${CYAN}[${LOCAL_FF_DIRS[$i]}]${NC}"
            ff_idx=$((ff_idx + 1))
        done
        echo ""
    fi

    # Mostrar force fields nativos
    if [ ${#NATIVE_FF_DIRS[@]} -gt 0 ]; then
        echo -e "  ${YELLOW}── Force fields nativos de GROMACS ──${NC}"
        for i in "${!NATIVE_FF_DIRS[@]}"; do
            echo "   ${ff_idx}) ${NATIVE_FF_NAMES[$i]}"
            echo -e "      ${CYAN}[${NATIVE_FF_DIRS[$i]}]${NC}"
            ff_idx=$((ff_idx + 1))
        done
        echo ""
    fi

    echo -e "  ${YELLOW}Seleccione el force field [1-${total_ff}]:${NC}"
    read -r FF_CHOICE

    if [[ "$FF_CHOICE" =~ ^[0-9]+$ ]] && [ "$FF_CHOICE" -ge 1 ] && [ "$FF_CHOICE" -le "$total_ff" ]; then
        if [ "$FF_CHOICE" -le ${#LOCAL_FF_DIRS[@]} ]; then
            local sel=$((FF_CHOICE - 1))
            SELECTED_FF_DIR="${LOCAL_FF_DIRS[$sel]}"
            SELECTED_FF_LOCAL=true
        else
            local sel=$((FF_CHOICE - ${#LOCAL_FF_DIRS[@]} - 1))
            SELECTED_FF_DIR="${NATIVE_FF_DIRS[$sel]}"
            SELECTED_FF_LOCAL=false
        fi
    else
        log_error "Selección inválida"
        exit 1
    fi

    log_success "Force field seleccionado: $SELECTED_FF_DIR (local=$SELECTED_FF_LOCAL)"

    # ========================================
    # Tipo de caja
    # ========================================
    echo -e "\n${CYAN}═══════════════════════════════════════════${NC}"
    echo -e "${CYAN}  PARÁMETROS DE SIMULACIÓN${NC}"
    echo -e "${CYAN}═══════════════════════════════════════════${NC}\n"

    echo "  Tipo de caja de simulación:"
    echo "   1) cubic       - Caja cúbica (más volumen, más moléculas de agua)"
    echo "   2) triclinic   - Romboedro truncado (~71% del volumen de cubic)"
    echo "   3) dodecahedron - Dodecaedro rómbico (~71%, recomendado)"
    echo "   4) octahedron  - Octaedro truncado (~77% del volumen de cubic)"
    echo -e "\n  ${YELLOW}Seleccione el tipo de caja [1-4] (default: 3 dodecahedron):${NC}"
    read -r BOX_CHOICE

    case $BOX_CHOICE in
        1) BOX_TYPE="cubic" ;;
        2) BOX_TYPE="triclinic" ;;
        4) BOX_TYPE="octahedron" ;;
        *) BOX_TYPE="dodecahedron" ;;
    esac

    # Distancia a bordes
    while true; do
        echo -e "\n  Distancia mínima a los bordes de la caja (nm) [default: 1.2]:"
        read -r BOX_DIST
        BOX_DIST=${BOX_DIST:-1.2}
        if is_positive_float "$BOX_DIST"; then
            break
        fi
        log_error "BOX_DIST inválido: '$BOX_DIST'. Debe ser un float > 0"
    done

    # Modelo de agua
    echo -e "\n  Modelo de agua a utilizar:"
    echo "   1) tip3p  - TIP3P (rápido, recomendado para la mayoría)"
    echo "   2) spc    - SPC (simple point charge)"
    echo "   3) spce   - SPC/E (extended simple point charge)"
    echo "   4) tip4p  - TIP4P (4 sitios)"
    echo "   5) tip5p  - TIP5P (5 sitios, más preciso pero más lento)"
    echo -e "\n  ${YELLOW}Seleccione el modelo de agua [1-5] (default: 1):${NC}"
    read -r WATER_CHOICE

    case $WATER_CHOICE in
        2) WATER_MODEL="spc" ;;
        3) WATER_MODEL="spce" ;;
        4) WATER_MODEL="tip4p" ;;
        5) WATER_MODEL="tip5p" ;;
        *) WATER_MODEL="tip3p" ;;
    esac
    set_water_file_from_model

    # Concentración iónica
    while true; do
        echo -e "\n  Concentración de iones NaCl (mol/L) [default: 0]:"
        echo -e "  ${CYAN}(0 = solo neutralizar, 0.15 = fisiológica)${NC}"
        read -r ION_CONC
        ION_CONC=${ION_CONC:-0}
        if is_non_negative_float "$ION_CONC"; then
            break
        fi
        log_error "ION_CONC inválido: '$ION_CONC'. Debe ser un float >= 0"
    done

    # Tiempo de producción
    while true; do
        echo -e "\n  Tiempo de simulación de producción en nanosegundos [default: 10]:"
        echo -e "  ${CYAN}(acepta decimales: 0.5, 1, 10, 50, 100)${NC}"
        read -r PROD_NS
        PROD_NS=${PROD_NS:-10}
        if is_positive_float "$PROD_NS"; then
            break
        fi
        log_error "PROD_NS inválido: '$PROD_NS'. Debe ser un número positivo"
    done

    # Threads
    echo -e "\n  Número de threads para mdrun [default: auto ($(nproc 2>/dev/null || echo 4))]:"
    read -r USER_NT
    if [ -n "$USER_NT" ]; then
        INPUT_NT="$USER_NT"
    fi

    validate_runtime_parameters
    log_info "Producción: $PROD_NS ns ($PROD_NSTEPS steps)"
    log_info "Threads: $NT"
}

#==========================================
# VALIDACIÓN DE ARCHIVOS MDP
#==========================================
validate_mdp_files() {
    log_step "Verificando archivos de entrada"

    if [ ! -d "$BASE_MDP" ]; then
        log_error "Carpeta '$BASE_MDP/' no encontrada"
        echo "   Debe contener los archivos .mdp necesarios:" >&2
        echo "   - ions.mdp, em.mdp, nvt.mdp, npt.mdp, md_prod.mdp" >&2
        exit 1
    fi

    local mdp_files=("ions.mdp" "em.mdp" "nvt.mdp" "npt.mdp" "md_prod.mdp")
    for mdp in "${mdp_files[@]}"; do
        if [ ! -f "$BASE_MDP/$mdp" ]; then
            log_error "Archivo '$BASE_MDP/$mdp' no encontrado"
            exit 1
        fi
    done

    log_success "Proteína: $PROT"
    log_success "Ligando: $LIG"
    log_success "Tipo de caja: $BOX_TYPE"
    log_success "Distancia a bordes: $BOX_DIST nm"
    log_success "Modelo de agua: $WATER_MODEL"
    log_success "Concentración iónica: $ION_CONC M"
    log_success "Archivos MDP validados"
}

apply_non_interactive_inputs() {
    log_step "Aplicando parámetros no interactivos"

    PROT="$INPUT_PROT"
    LIG="$INPUT_LIG"
    SELECTED_FF_DIR="$INPUT_FF"

    if [ -d "$INITIAL_DIR/$INPUT_FF" ]; then
        SELECTED_FF_LOCAL=true
    else
        SELECTED_FF_LOCAL=false
    fi

    BOX_TYPE="${INPUT_BOX_TYPE:-dodecahedron}"
    BOX_DIST="${INPUT_BOX_DIST:-1.2}"
    WATER_MODEL="${INPUT_WATER_MODEL:-tip3p}"
    ION_CONC="${INPUT_ION_CONC:-0}"
    PROD_NS="${INPUT_PROD_NS:-10}"

    validate_runtime_parameters

    [ -d "$INITIAL_DIR/$BASE_PROT/$PROT" ] || { log_error "Proteína no encontrada: $BASE_PROT/$PROT"; exit 1; }
    [ -f "$INITIAL_DIR/$BASE_PROT/$PROT/proteina.gro" ] || { log_error "Falta proteina.gro en $BASE_PROT/$PROT"; exit 1; }
    [ -f "$INITIAL_DIR/$BASE_PROT/$PROT/topol.top" ] || { log_error "Falta topol.top en $BASE_PROT/$PROT"; exit 1; }
    [ -f "$INITIAL_DIR/$BASE_PROT/$PROT/posre.itp" ] || { log_error "Falta posre.itp en $BASE_PROT/$PROT"; exit 1; }

    [ -d "$INITIAL_DIR/$BASE_LIG/$LIG" ] || { log_error "Ligando no encontrado: $BASE_LIG/$LIG"; exit 1; }
    [ -f "$INITIAL_DIR/$BASE_LIG/$LIG/ligando.gro" ] || { log_error "Falta ligando.gro en $BASE_LIG/$LIG"; exit 1; }
    [ -f "$INITIAL_DIR/$BASE_LIG/$LIG/ligando.itp" ] || { log_error "Falta ligando.itp en $BASE_LIG/$LIG"; exit 1; }

    if [ "$SELECTED_FF_LOCAL" = false ]; then
        if ! [ -d "$INITIAL_DIR/$INPUT_FF" ]; then
            log_info "Force field asumido como nativo de GROMACS: $INPUT_FF"
        fi
    fi

    log_success "Modo no interactivo configurado"
    log_info "  PROT=$PROT"
    log_info "  LIG=$LIG"
    log_info "  FF=$SELECTED_FF_DIR (local=$SELECTED_FF_LOCAL)"
    log_info "  BOX_TYPE=$BOX_TYPE"
    log_info "  BOX_DIST=$BOX_DIST"
    log_info "  WATER_MODEL=$WATER_MODEL"
    log_info "  ION_CONC=$ION_CONC"
    log_info "  PROD_NS=$PROD_NS"
    if [ -n "$GPU_ID" ]; then
        log_info "  GPU_ID=$GPU_ID"
    fi
    if [ "$DRY_RUN" = true ]; then
        log_info "  MODO=DRY-RUN (sin mdrun)"
    fi
}

#==========================================
# VALIDACIÓN DE INTEGRIDAD DEL LIGANDO (v4.5)
#==========================================
validate_ligand_files() {
    log_step "Validando integridad de archivos del ligando (v4.5)"

    local lig_gro="$INITIAL_DIR/$BASE_LIG/$LIG/ligando.gro"
    local lig_itp="$INITIAL_DIR/$BASE_LIG/$LIG/ligando.itp"

    # Verificar que el .gro contiene residuo LIG
    if ! grep -q 'LIG' "$lig_gro"; then
        log_error "ligando.gro no contiene residuo 'LIG'. Verifica la parametrización."
        log_info "Tip: el residuo del ligando debe llamarse 'LIG' en el archivo .gro"
        exit 1
    fi
    log_success "ligando.gro contiene residuo LIG"

    # Contar átomos pesados en .gro (línea 2 = número total de átomos)
    local natoms_gro
    natoms_gro=$(sed -n '2p' "$lig_gro" | awk '{print $1}')
    if [ -z "$natoms_gro" ] || ! is_int "$natoms_gro" || [ "$natoms_gro" -lt 1 ]; then
        log_error "ligando.gro tiene formato inválido: no se pudo leer número de átomos"
        exit 1
    fi

    # Contar átomos en .itp (sección [ atoms ])
    local natoms_itp
    natoms_itp=$(awk '
        /^\[/ { in_atoms=0 }
        /^\[[[:space:]]*atoms[[:space:]]*\]/ { in_atoms=1; next }
        in_atoms && /^[[:space:]]*[0-9]/ { count++ }
        END { print count+0 }
    ' "$lig_itp")

    if [ "$natoms_gro" -ne "$natoms_itp" ]; then
        log_error "Inconsistencia: ligando.gro tiene $natoms_gro átomos, ligando.itp tiene $natoms_itp"
        log_info "Regenera los archivos del ligando con CGenFF/GAFF2"
        exit 1
    fi
    log_success "Consistencia .gro/.itp verificada: $natoms_gro átomos"
}

#==========================================
# ESTRUCTURA DE CARPETAS
#==========================================
setup_directory_structure() {
    log_step "Creando estructura de carpetas"

    RUNDIR="$INITIAL_DIR/$WORKDIR/${PROT}_${LIG}_$(date +%Y%m%d_%H%M%S)"

    local dirs=(
        "$RUNDIR"
        "$RUNDIR/00_setup"
        "$RUNDIR/01_minimization"
        "$RUNDIR/02_equilibration"
        "$RUNDIR/03_production"
        "$RUNDIR/04_analysis"
        "$RUNDIR/logs"
        "$RUNDIR/mdp_used"
    )

    for dir in "${dirs[@]}"; do
        create_dir "$dir"
    done

    log_success "Directorio de trabajo: $RUNDIR"

    # Copiar archivos MDP
    log_step "Copiando archivos MDP para referencia"
    cp "$INITIAL_DIR/$BASE_MDP"/*.mdp "$RUNDIR/mdp_used/"
    log_success "Archivos MDP guardados en mdp_used/"
}

#==========================================
# SETUP INICIAL
#==========================================
setup_initial_files() {
    log_step "Copiando archivos de proteína y ligando"

    cd "$RUNDIR/00_setup" || exit 1

    # Proteína
    cp "$INITIAL_DIR/$BASE_PROT/$PROT/proteina.gro" .
    cp "$INITIAL_DIR/$BASE_PROT/$PROT/topol.top" .
    cp "$INITIAL_DIR/$BASE_PROT/$PROT/posre.itp" .

    # Ligando - copiar todos los archivos disponibles
    for f in ligando.itp ligando.gro ligando.pdb; do
        if [ -f "$INITIAL_DIR/$BASE_LIG/$LIG/$f" ]; then
            cp "$INITIAL_DIR/$BASE_LIG/$LIG/$f" .
        fi
    done

    if [ -f "$INITIAL_DIR/$BASE_LIG/$LIG/ligando.prm" ]; then
        cp "$INITIAL_DIR/$BASE_LIG/$LIG/ligando.prm" .
        log_success "Archivo ligando.prm encontrado (parámetros CHARMM)"
    fi

    # Usar el force field seleccionado por el usuario
    FF_DIR="${SELECTED_FF_DIR:-}"
    FF_IS_LOCAL="${SELECTED_FF_LOCAL:-false}"

    if [ -n "$FF_DIR" ]; then
        if [ "$FF_IS_LOCAL" = true ] && [ -d "$INITIAL_DIR/$FF_DIR" ]; then
            log_step "Copiando force field local: $FF_DIR"
            cp -r "$INITIAL_DIR/$FF_DIR" .
            log_success "$FF_DIR copiado a 00_setup/"
        else
            log_success "Force field nativo de GROMACS: $FF_DIR"
        fi
    else
        log_warning "No se seleccionó force field, usando el referenciado en topol.top"
        # Fallback: detectar desde topol.top
        FF_DIR=$(sed -n 's#.*"\./\([^/]*\.ff\)/.*#\1#p' topol.top | head -1)
        if [ -z "$FF_DIR" ]; then
            FF_DIR=$(sed -n 's|^#include[[:space:]]\+"\([^/]*\.ff\)/.*|\1|p' topol.top | head -1)
        fi
        if [ -n "$FF_DIR" ] && [ -d "$INITIAL_DIR/$FF_DIR" ]; then
            cp -r "$INITIAL_DIR/$FF_DIR" .
            FF_IS_LOCAL=true
        fi
    fi

    log_success "Archivos copiados a 00_setup/"

    # Adaptar MDPs si el force field no es CHARMM
    adapt_mdp_for_forcefield "$RUNDIR/mdp_used"
}

#==========================================
# RESTRAINTS DEL LIGANDO
#==========================================
generate_ligand_restraints() {
    log_step "Generando restraints del ligando (sin H)"

    cd "$RUNDIR/00_setup" || exit 1

    # Determinar índice del nuevo grupo dinámicamente
    local last_idx
    last_idx=$(echo q | "$GMX" make_ndx -f ligando.gro 2>&1 | awk '/^[[:space:]]*[0-9]+[[:space:]]/ {idx=$1} END{print idx+0}')
    local new_grp=$((last_idx + 1))

    # Crear índice excluyendo hidrógenos
    run_gmx make_ndx -f ligando.gro -o lig_noh.ndx &> "$RUNDIR/logs/make_ndx_lig.log" <<EOF
r LIG & !a H*
name ${new_grp} LIG-H
q
EOF

    if [ ! -f lig_noh.ndx ]; then
        log_error "lig_noh.ndx no se creó correctamente"
        exit 1
    fi

    # Generar posre
    echo "LIG-H" | run_gmx genrestr -f ligando.gro -n lig_noh.ndx -o posre_ligando.itp \
        -fc 1000 1000 1000 &> "$RUNDIR/logs/genrestr.log"

    log_success "posre_ligando.itp generado"
}

#==========================================
# MODIFICAR TOPOLOGÍA
#==========================================
has_include_line() {
    local top_file="$1"
    local include_file="$2"
    awk -v inc="$include_file" '
        $0 ~ "^[[:space:]]*#include[[:space:]]+\"" inc "\"[[:space:]]*$" {found=1}
        END {exit !found}
    ' "$top_file"
}

count_include_line() {
    local top_file="$1"
    local include_file="$2"
    awk -v inc="$include_file" '
        $0 ~ "^[[:space:]]*#include[[:space:]]+\"" inc "\"[[:space:]]*$" {count++}
        END {print count+0}
    ' "$top_file"
}

ensure_unique_include_after_anchor() {
    local top_file="$1"
    local include_file="$2"
    local anchor_regex="$3"
    local include_line="#include \"${include_file}\""
    local tmp_file="${top_file}.tmp.include.$$"
    register_temp_file "$tmp_file"

    awk -v inc="$include_file" '
        $0 ~ "^[[:space:]]*#include[[:space:]]+\"" inc "\"[[:space:]]*$" {
            seen++
            if (seen > 1) next
        }
        {print}
    ' "$top_file" > "$tmp_file"
    mv "$tmp_file" "$top_file"

    if ! has_include_line "$top_file" "$include_file"; then
        local anchor_line
        anchor_line=$(awk -v re="$anchor_regex" '$0 ~ re {print NR; exit}' "$top_file")
        [ -n "$anchor_line" ] || {
            log_error "No se encontró línea ancla para incluir ${include_file}"
            return 1
        }

        tmp_file="${top_file}.tmp.insert.$$"
        register_temp_file "$tmp_file"
        awk -v n="$anchor_line" -v ins="$include_line" '
            NR == n {print; print ins; next}
            {print}
        ' "$top_file" > "$tmp_file"
        mv "$tmp_file" "$top_file"
    fi

    local include_count
    include_count=$(count_include_line "$top_file" "$include_file")
    if [ "$include_count" -ne 1 ]; then
        log_error "El include ${include_file} quedó duplicado o ausente (count=${include_count})"
        return 1
    fi
}

ensure_single_ligand_molecule_entry() {
    local top_file="$1"
    local tmp_file="${top_file}.tmp.molecules.$$"
    register_temp_file "$tmp_file"

    awk '
        BEGIN {in_mol=0; has_lig=0; seen_section=0}
        /^\[[[:space:]]*molecules[[:space:]]*\]/ {
            seen_section=1
            in_mol=1
            print
            next
        }
        /^\[/ {
            if (in_mol && !has_lig) {
                print "LIG                 1"
                has_lig=1
            }
            in_mol=0
            print
            next
        }
        {
            if (in_mol) {
                if ($0 ~ /^[[:space:]]*;/ || $0 ~ /^[[:space:]]*$/) {
                    print
                    next
                }

                if ($1 == "LIG") {
                    if (has_lig) next
                    has_lig=1
                    print "LIG                 1"
                    next
                }
            }
            print
        }
        END {
            if (!seen_section) {
                exit 2
            }
            if (in_mol && !has_lig) {
                print "LIG                 1"
            }
        }
    ' "$top_file" > "$tmp_file" || return 1

    mv "$tmp_file" "$top_file"
}

verify_topology_integrity() {
    local top_file="$1"
    local include_count lig_count

    awk '/^\[[[:space:]]*molecules[[:space:]]*\]/{found=1} END{exit !found}' "$top_file" || {
        log_error "topol.top corrupto: no se encontró sección [ molecules ]"
        return 1
    }

    include_count=$(count_include_line "$top_file" "ligando.itp")
    [ "$include_count" -eq 1 ] || {
        log_error "topol.top corrupto: include de ligando.itp esperado=1, actual=${include_count}"
        return 1
    }

    if [ -f "ligando.prm" ]; then
        include_count=$(count_include_line "$top_file" "ligando.prm")
        [ "$include_count" -eq 1 ] || {
            log_error "topol.top corrupto: include de ligando.prm esperado=1, actual=${include_count}"
            return 1
        }
    fi

    lig_count=$(awk '
        /^\[[[:space:]]*molecules[[:space:]]*\]/{in_mol=1; next}
        /^\[/{in_mol=0}
        in_mol && $0 !~ /^[[:space:]]*;/ && $1=="LIG" {count++}
        END{print count+0}
    ' "$top_file")
    [ "$lig_count" -eq 1 ] || {
        log_error "topol.top corrupto: entrada LIG esperada=1, actual=${lig_count}"
        return 1
    }

    return 0
}

update_topology() {
    log_step "Actualizando archivos de topología"

    cd "$RUNDIR/00_setup" || exit 1

    if ! grep -q "POSRES_LIG" ligando.itp; then
        cat <<'EOF' >> ligando.itp

; Ligand position restraints
#ifdef POSRES_LIG
#include "posre_ligando.itp"
#endif
EOF
        log_success "Restraints añadidos a ligando.itp"
    fi

    log_step "Modificando topol.top (inserción segura e idempotente)"
    cp topol.top topol.top.backup

    local forcefield_anchor='^[[:space:]]*#include[[:space:]]+".*forcefield\.itp"[[:space:]]*$'
    local prm_anchor='^[[:space:]]*#include[[:space:]]+"ligando\.prm"[[:space:]]*$'

    if [ -f "ligando.prm" ]; then
        ensure_unique_include_after_anchor topol.top "ligando.prm" "$forcefield_anchor"
        log_success "Include de ligando.prm validado"
    fi

    if [ -f "ligando.prm" ]; then
        ensure_unique_include_after_anchor topol.top "ligando.itp" "$prm_anchor"
    else
        ensure_unique_include_after_anchor topol.top "ligando.itp" "$forcefield_anchor"
    fi
    log_success "Include de ligando.itp validado"

    ensure_single_ligand_molecule_entry topol.top || {
        log_error "No se pudo insertar entrada LIG en [ molecules ]"
        cp topol.top.backup topol.top
        exit 1
    }

    if ! verify_topology_integrity topol.top; then
        log_error "Falló verificación de topología; restaurando backup"
        cp topol.top.backup topol.top
        exit 1
    fi

    log_success "topol.top actualizado correctamente (idempotente)"
}

#==========================================
# CONSTRUCCIÓN DEL SISTEMA
#==========================================
build_complex() {
    log_step "Ensamblando complejo proteína-ligando"

    cd "$RUNDIR/00_setup" || exit 1

    run_gmx editconf -f proteina.gro -o prot.gro &> "$RUNDIR/logs/editconf_prot.log"

    run_gmx insert-molecules -f prot.gro -ci ligando.gro -o complex.gro -nmol 1 \
        &> "$RUNDIR/logs/insert_molecules.log"

    log_success "Complejo ensamblado: complex.gro"
}

#==========================================
# SOLVACIÓN
#==========================================
solvate_system() {
    log_step "Definiendo caja ($BOX_TYPE) y solvatando con $WATER_MODEL"

    cd "$RUNDIR/00_setup" || exit 1

    run_gmx editconf -f complex.gro -o boxed.gro -d "$BOX_DIST" -bt "$BOX_TYPE" \
        &> "$RUNDIR/logs/editconf_box.log"
    log_success "Caja $BOX_TYPE definida (d=$BOX_DIST nm)"

    run_gmx solvate -cp boxed.gro -cs "$WATER_FILE" -o solv.gro -p topol.top \
        &> "$RUNDIR/logs/solvate.log"
    log_success "Sistema solvatado con $WATER_MODEL"
}

#==========================================
# NEUTRALIZACIÓN
#==========================================
neutralize_system() {
    log_step "Añadiendo iones (concentración: $ION_CONC M NaCl)"

    cd "$RUNDIR/00_setup" || exit 1

    if ! run_gmx grompp -f "$RUNDIR/mdp_used/ions.mdp" -c solv.gro -p topol.top -o ions.tpr \
        -maxwarn "$MAXWARN" > "$RUNDIR/logs/grompp_ions.log" 2>&1; then
        log_error "grompp falló al generar ions.tpr"
        cat "$RUNDIR/logs/grompp_ions.log"
        exit 1
    fi

    if is_positive_float "$ION_CONC"; then
        echo "SOL" | run_gmx genion -s ions.tpr -o system.gro -p topol.top \
            -pname NA -nname CL -neutral -conc "$ION_CONC" \
            &> "$RUNDIR/logs/genion.log"
        log_success "Sistema neutralizado con $ION_CONC M NaCl"
    else
        echo "SOL" | run_gmx genion -s ions.tpr -o system.gro -p topol.top -neutral \
            &> "$RUNDIR/logs/genion.log"
        log_success "Sistema neutralizado (solo contraiones)"
    fi
}

#==========================================
# GRUPOS DE ÍNDICE
#==========================================
create_index_groups() {
    log_step "Creando grupos para termostatos"

    cd "$RUNDIR/00_setup" || exit 1

    # Determinar dinámicamente los índices de los nuevos grupos
    local last_idx
    last_idx=$(echo q | "$GMX" make_ndx -f system.gro 2>&1 | awk '/^[[:space:]]*[0-9]+[[:space:]]/ {idx=$1} END{print idx+0}')
    local new_grp1=$((last_idx + 1))
    local new_grp2=$((last_idx + 2))

    run_gmx make_ndx -f system.gro -o index.ndx &> "$RUNDIR/logs/make_ndx_system.log" <<EOF
1 | r LIG
r SOL | r NA | r CL
name ${new_grp1} Protein_Ligand
name ${new_grp2} Solvent
q
EOF

    log_success "index.ndx creado (Protein_Ligand y Solvent)"
}

#==========================================
# HELPER: Copiar archivos comunes desde 00_setup
# (cp en vez de ln -sf: symlinks no funcionan en /mnt/c/ NTFS/WSL)
#==========================================
link_setup_files() {
    local target_dir="$1"
    cd "$target_dir" || exit 1

    cp -f ../00_setup/topol.top .
    cp -f ../00_setup/index.ndx .

    # Copiar archivos .itp explícitamente
    for itp_file in posre.itp ligando.itp posre_ligando.itp; do
        [ -f "../00_setup/$itp_file" ] && cp -f "../00_setup/$itp_file" .
    done

    # Copiar .prm si existe
    [ -f "../00_setup/ligando.prm" ] && cp -f "../00_setup/ligando.prm" .

    # Copiar el force field solo si es local (no nativo de GROMACS)
    if [ "$FF_IS_LOCAL" = true ]; then
        for ff_candidate in ../00_setup/*.ff; do
            [ -d "$ff_candidate" ] && cp -r "$ff_candidate" .
        done
    fi
}

#==========================================
# MINIMIZACIÓN
#==========================================
run_minimization() {
    log_step "Ejecutando minimización de energía"

    link_setup_files "$RUNDIR/01_minimization"
    cp -f ../00_setup/system.gro .

    run_gmx grompp -f "$RUNDIR/mdp_used/em.mdp" -c system.gro -p topol.top -n index.ndx \
        -o em.tpr -maxwarn "$MAXWARN" &> "$RUNDIR/logs/grompp_em.log"
    # Log grompp warnings
    grep -i 'WARNING' "$RUNDIR/logs/grompp_em.log" | head -5 | while read -r w; do log_warning "grompp-em: $w"; done || true

    $MDRUN -deffnm em &> "$RUNDIR/logs/mdrun_em.log"
    log_success "Minimización completada"

    echo "Potential" | run_gmx energy -f em.edr -o energy_em.xvg &> "$RUNDIR/logs/energy_em.log"
    log_success "Energía potencial guardada: energy_em.xvg"
}

#==========================================
# EQUILIBRACIÓN NVT
#==========================================
run_nvt_equilibration() {
    log_step "Equilibración NVT (calentamiento) - CON RESTRICCIONES"

    link_setup_files "$RUNDIR/02_equilibration"
    cp -f ../01_minimization/em.gro .

    # Modificar MDP para incluir restricciones Y grupos correctos
    local nvt_mdp="nvt_temp.mdp"
    cp "$RUNDIR/mdp_used/nvt.mdp" "$nvt_mdp"

    # Añadir define
    if ! grep -q "^define" "$nvt_mdp"; then
        sed -i '1i define = -DPOSRES -DPOSRES_LIG' "$nvt_mdp"
    else
        sed -i 's/^define.*/define = -DPOSRES -DPOSRES_LIG/' "$nvt_mdp"
    fi

    # Reemplazar grupos de temperatura (tau_t se deja como está en el MDP original)
    sed -i 's/^tc-grps.*/tc-grps                  = Protein_Ligand Solvent/' "$nvt_mdp"
    sed -i 's/^ref_t.*/ref_t                    = 300     300/' "$nvt_mdp"

    # === CORRECCIÓN reproducibilidad: inyectar semilla determinista ===
    # Derivada del timestamp del RUNDIR → única por corrida, reproducible si se conoce RUNDIR
    local nvt_seed
    nvt_seed=$(echo "${RUNDIR##*_}" | tr -dc '0-9' | head -c 8)
    nvt_seed="${nvt_seed:-91337}"
    sed -i "s/__SEED_PLACEHOLDER__/${nvt_seed}/" "$nvt_mdp"
    log_info "Semilla NVT inyectada: gen_seed = ${nvt_seed} (reproducible con este RUNDIR)"

    log_success "Grupos de termostato actualizados en nvt_temp.mdp"

    run_gmx grompp -f "$nvt_mdp" -c em.gro -r em.gro -p topol.top -n index.ndx \
        -o nvt.tpr -maxwarn "$MAXWARN" &> "$RUNDIR/logs/grompp_nvt.log"
    grep -i 'WARNING' "$RUNDIR/logs/grompp_nvt.log" | head -5 | while read -r w; do log_warning "grompp-nvt: $w"; done || true

    $MDRUN -deffnm nvt &> "$RUNDIR/logs/mdrun_nvt.log"
    log_success "NVT completado (proteína y ligando restringidos)"

    echo "Temperature" | run_gmx energy -f nvt.edr -o temperature_nvt.xvg \
        &> "$RUNDIR/logs/energy_nvt.log"
    log_success "Temperatura guardada: temperature_nvt.xvg"
}

#==========================================
# EQUILIBRACIÓN NPT
#==========================================
run_npt_equilibration() {
    log_step "Equilibración NPT (presión) - CON RESTRICCIONES"

    cd "$RUNDIR/02_equilibration" || exit 1

    local npt_mdp="npt_temp.mdp"
    cp "$RUNDIR/mdp_used/npt.mdp" "$npt_mdp"

    # Añadir define
    if ! grep -q "^define" "$npt_mdp"; then
        sed -i '1i define = -DPOSRES -DPOSRES_LIG' "$npt_mdp"
    else
        sed -i 's/^define.*/define = -DPOSRES -DPOSRES_LIG/' "$npt_mdp"
    fi

    # Reemplazar grupos de temperatura (tau_t se deja como está en el MDP original)
    sed -i 's/^tc-grps.*/tc-grps                  = Protein_Ligand Solvent/' "$npt_mdp"
    sed -i 's/^ref_t.*/ref_t                    = 300     300/' "$npt_mdp"

    log_success "Grupos de termostato actualizados en npt_temp.mdp"

    run_gmx grompp -f "$npt_mdp" -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top \
        -n index.ndx -o npt.tpr -maxwarn "$MAXWARN" &> "$RUNDIR/logs/grompp_npt.log"
    grep -i 'WARNING' "$RUNDIR/logs/grompp_npt.log" | head -5 | while read -r w; do log_warning "grompp-npt: $w"; done || true

    $MDRUN -deffnm npt &> "$RUNDIR/logs/mdrun_npt.log"
    log_success "NPT completado (proteína y ligando restringidos)"

    echo "Pressure" | run_gmx energy -f npt.edr -o pressure_npt.xvg \
        &> "$RUNDIR/logs/energy_npt_press.log"
    echo "Density" | run_gmx energy -f npt.edr -o density_npt.xvg \
        &> "$RUNDIR/logs/energy_npt_dens.log"
    log_success "Presión y densidad guardadas"
}

#==========================================
# PRODUCCIÓN
#==========================================
run_production() {
    log_step "Iniciando simulación de producción - SIN RESTRICCIONES"

    link_setup_files "$RUNDIR/03_production"
    cp -f ../02_equilibration/npt.gro .
    cp -f ../02_equilibration/npt.cpt .

    # Modificar MDP para grupos correctos (sin restricciones)
    local md_mdp="md_prod_temp.mdp"
    cp "$RUNDIR/mdp_used/md_prod.mdp" "$md_mdp"

    # Actualizar nsteps según el tiempo elegido por el usuario
    sed -i "s/^nsteps.*/nsteps                   = ${PROD_NSTEPS}      ; ${PROD_NS} ns/" "$md_mdp"

    # Reemplazar grupos de temperatura (tau_t se deja como está en el MDP original)
    sed -i 's/^tc-grps.*/tc-grps                  = Protein_Ligand Solvent/' "$md_mdp"
    sed -i 's/^ref_t.*/ref_t                    = 300     300/' "$md_mdp"

    # v4.6: Producción SIN energygrps → GPU full speed
    # energygrps es incompatible con cálculos no-bonded en GPU (GROMACS lo fuerza a CPU)
    # Solución: producción en GPU limpia + rerun posterior para extraer energías Protein↔LIG
    if grep -q '^energygrps' "$md_mdp"; then
        sed -i '/^energygrps/d' "$md_mdp"
        sed -i '/^; GRUPOS DE ENERGÍA/d' "$md_mdp"
        log_info "energygrps eliminado del MDP de producción (incompatible con GPU)"
    fi

    log_success "Producción configurada: $PROD_NS ns ($PROD_NSTEPS steps)"

    # --- TPR PRINCIPAL: sin energygrps, GPU-compatible ---
    run_gmx grompp -f "$md_mdp" -c npt.gro -t npt.cpt -p topol.top \
        -n index.ndx -o md.tpr -maxwarn "$MAXWARN" &> "$RUNDIR/logs/grompp_md.log"
    grep -i 'WARNING' "$RUNDIR/logs/grompp_md.log" | head -5 | while read -r w; do log_warning "grompp-md: $w"; done || true

    # --- TPR RERUN: con energygrps = Protein LIG (para post-producción) ---
    local md_rerun_mdp="md_rerun_temp.mdp"
    cp "$md_mdp" "$md_rerun_mdp"
    printf '\n; GRUPOS DE ENERGÍA (rerun post-producción, no afecta velocidad GPU)\nenerygrps_placeholder\n' >> "$md_rerun_mdp"
    sed -i 's/^enerygrps_placeholder/energygrps               = Protein LIG/' "$md_rerun_mdp"
    run_gmx grompp -f "$md_rerun_mdp" -c npt.gro -t npt.cpt -p topol.top \
        -n index.ndx -o md_rerun.tpr -maxwarn "$MAXWARN" &> "$RUNDIR/logs/grompp_md_rerun.log"
    grep -i 'WARNING' "$RUNDIR/logs/grompp_md_rerun.log" | head -5 | while read -r w; do log_warning "grompp-rerun: $w"; done || true
    rm -f "$md_rerun_mdp"
    log_success "TPR para rerun preparado: md_rerun.tpr (energygrps = Protein LIG)"

    log_success "Sistema preparado para producción (GPU)"
    echo -e "\n${YELLOW}Ejecutando producción en GPU (esto puede tardar)...${NC}\n"

    # Monitor de progreso con ETA
    monitor_mdrun "$RUNDIR/logs/mdrun_md.log" "$PROD_NS" &
    local monitor_pid=$!

    $MDRUN -deffnm md &> "$RUNDIR/logs/mdrun_md.log"

    # Terminar monitor
    kill $monitor_pid 2>/dev/null || true
    wait $monitor_pid 2>/dev/null || true
    echo ""
    log_success "Simulación de producción completada"

    # --- RERUN: recalcular energías Protein↔LIG sobre la trayectoria ---
    log_step "Rerun: extrayendo energías Protein↔LIG de la trayectoria (sin re-simular)"
    if [ -f "md.xtc" ] && [ -f "md_rerun.tpr" ]; then
        $MDRUN -s md_rerun.tpr -rerun md.xtc -deffnm md_rerun \
            &> "$RUNDIR/logs/mdrun_rerun.log"
        log_success "Rerun completado → md_rerun.edr contiene energías Protein↔LIG"
    else
        log_warning "Rerun omitido: no se encontró md.xtc o md_rerun.tpr"
        log_warning "  Puedes ejecutarlo manualmente:"
        log_warning "  gmx_mpi mdrun -s md_rerun.tpr -rerun md.xtc -deffnm md_rerun"
    fi
}

#==========================================
# MONITOR DE PROGRESO CON ETA (v4.5)
#==========================================
monitor_mdrun() {
    local log_file="$1"
    local total_ns="$2"
    local total_ps
    total_ps=$(awk "BEGIN{printf \"%d\", $total_ns * 1000}")

    sleep 10  # Esperar a que mdrun genere output

    while true; do
        if [ ! -f "$log_file" ]; then
            sleep 15
            continue
        fi

        # Extraer el último paso completado del log
        local current_step
        current_step=$(tail -50 "$log_file" 2>/dev/null | awk '/^[[:space:]]*step[[:space:]]+[0-9]/{step=$2} /Step[[:space:]]+Time/{getline; if($1 ~ /^[0-9]/) step=$1} END{print step+0}' 2>/dev/null || echo 0)

        if [ "$current_step" -gt 0 ] && [ "$total_ps" -gt 0 ]; then
            local current_ps
            current_ps=$(awk "BEGIN{printf \"%.0f\", $current_step * 0.002}")
            local pct
            pct=$(awk "BEGIN{printf \"%.1f\", ($current_ps / $total_ps) * 100}")

            # ETA basado en tiempo transcurrido
            local elapsed
            elapsed=$(awk '/Performance:/{t=$2} END{if(t) print t+0; else print 0}' "$log_file" 2>/dev/null || echo 0)

            printf "\r  ${CYAN}⏳ Progreso: %.1f%% (%s / %s ps)${NC}" "$pct" "$current_ps" "$total_ps"
        fi

        sleep 30
    done
}

#==========================================
# EXTENDER SIMULACIÓN (v4.5)
#==========================================
extend_simulation() {
    log_step "Extendiendo simulación por $EXTEND_NS ns adicionales"

    cd "$RUNDIR/03_production" || exit 1

    if [ ! -f "md.tpr" ]; then
        log_error "No se encontró md.tpr en 03_production/. ¿Se completó la producción original?"
        exit 1
    fi

    if [ ! -f "md.cpt" ]; then
        log_error "No se encontró md.cpt en 03_production/. No se puede extender sin checkpoint."
        exit 1
    fi

    local extend_ps
    extend_ps=$(awk "BEGIN{printf \"%d\", $EXTEND_NS * 1000}")

    # Extender el TPR
    run_gmx convert-tpr -s md.tpr -extend "$extend_ps" -o md_extended.tpr \
        &> "$RUNDIR/logs/extend_tpr.log"
    mv md.tpr md_original.tpr
    mv md_extended.tpr md.tpr
    log_success "TPR extendido por $EXTEND_NS ns ($extend_ps ps)"

    # Continuar la simulación
    echo -e "\n${YELLOW}Continuando simulación extendida...${NC}\n"

    $MDRUN -deffnm md -cpi md.cpt -noappend \
        &> "$RUNDIR/logs/mdrun_extend.log"
    log_success "Simulación extendida completada"

    # Concatenar trayectorias si se usó -noappend
    if ls md.part*.xtc &>/dev/null; then
        run_gmx trjcat -f md.part*.xtc -o md_full.xtc -cat \
            &> "$RUNDIR/logs/trjcat_extend.log"
        mv md.xtc md_original.xtc 2>/dev/null || true
        mv md_full.xtc md.xtc
        log_success "Trayectorias concatenadas: md.xtc (completa)"
    fi
}

#==========================================
# MODO ANÁLISIS-ONLY (v4.5)
#==========================================
run_analysis_only() {
    log_step "Modo analysis-only: re-ejecutando análisis sobre corrida existente"

    RUNDIR=$(cd "$ANALYSIS_ONLY_DIR" && pwd)

    if [ ! -d "$RUNDIR/03_production" ]; then
        log_error "Directorio no es una corrida MD válida: falta 03_production/"
        exit 1
    fi

    if [ ! -f "$RUNDIR/03_production/md.tpr" ] || [ ! -f "$RUNDIR/03_production/md.xtc" ]; then
        log_error "Faltan archivos de producción (md.tpr, md.xtc)"
        exit 1
    fi

    # Crear/recrear directorio de análisis
    mkdir -p "$RUNDIR/04_analysis"
    mkdir -p "$RUNDIR/logs"

    # Detectar GROMACS y FF
    FF_DIR=$(sed -n 's#.*"\./\([^/]*\.ff\)/.*#\1#p' "$RUNDIR/00_setup/topol.top" 2>/dev/null | head -1 || echo "")
    FF_IS_LOCAL=false
    if [ -n "$FF_DIR" ] && [ -d "$RUNDIR/00_setup/$FF_DIR" ]; then
        FF_IS_LOCAL=true
    fi

    # Leer info del checkpoint si existe
    if [ -f "$RUNDIR/.checkpoint" ]; then
        PROT=$(extract_checkpoint_value "$RUNDIR/.checkpoint" "PROT")
        LIG=$(extract_checkpoint_value "$RUNDIR/.checkpoint" "LIG")
    else
        PROT=$(basename "$RUNDIR" | cut -d'_' -f1)
        LIG=$(basename "$RUNDIR" | cut -d'_' -f2)
    fi

    log_success "Corrida detectada: $PROT + $LIG"
    log_success "Directorio: $RUNDIR"

    run_analysis
    log_step "Análisis re-ejecutado con éxito"
    echo -e "${GREEN}Resultados en: $RUNDIR/04_analysis/${NC}\n"
}

run_dry_run_grompp_checks() {
    log_step "Dry-run: validando topologías con grompp (sin mdrun)"

    # EM
    link_setup_files "$RUNDIR/01_minimization"
    cp -f ../00_setup/system.gro dryrun_system.gro
    run_gmx grompp -f "$RUNDIR/mdp_used/em.mdp" -c dryrun_system.gro -p topol.top -n index.ndx \
        -o em_dry.tpr -maxwarn "$MAXWARN" &> "$RUNDIR/logs/dryrun_grompp_em.log"
    log_success "Dry-run grompp EM OK"

    # NVT
    link_setup_files "$RUNDIR/02_equilibration"
    cp -f ../00_setup/system.gro dryrun_system.gro

    local nvt_mdp="nvt_temp_dry.mdp"
    cp "$RUNDIR/mdp_used/nvt.mdp" "$nvt_mdp"
    if ! grep -q "^define" "$nvt_mdp"; then
        sed -i '1i define = -DPOSRES -DPOSRES_LIG' "$nvt_mdp"
    else
        sed -i 's/^define.*/define = -DPOSRES -DPOSRES_LIG/' "$nvt_mdp"
    fi
    sed -i 's/^tc-grps.*/tc-grps                  = Protein_Ligand Solvent/' "$nvt_mdp"
    sed -i 's/^ref_t.*/ref_t                    = 300     300/' "$nvt_mdp"

    run_gmx grompp -f "$nvt_mdp" -c dryrun_system.gro -r dryrun_system.gro -p topol.top -n index.ndx \
        -o nvt_dry.tpr -maxwarn "$MAXWARN" &> "$RUNDIR/logs/dryrun_grompp_nvt.log"
    log_success "Dry-run grompp NVT OK"

    # NPT
    local npt_mdp="npt_temp_dry.mdp"
    cp "$RUNDIR/mdp_used/npt.mdp" "$npt_mdp"
    if ! grep -q "^define" "$npt_mdp"; then
        sed -i '1i define = -DPOSRES -DPOSRES_LIG' "$npt_mdp"
    else
        sed -i 's/^define.*/define = -DPOSRES -DPOSRES_LIG/' "$npt_mdp"
    fi
    sed -i 's/^tc-grps.*/tc-grps                  = Protein_Ligand Solvent/' "$npt_mdp"
    sed -i 's/^ref_t.*/ref_t                    = 300     300/' "$npt_mdp"

    run_gmx grompp -f "$npt_mdp" -c dryrun_system.gro -r dryrun_system.gro -p topol.top -n index.ndx \
        -o npt_dry.tpr -maxwarn "$MAXWARN" &> "$RUNDIR/logs/dryrun_grompp_npt.log"
    log_success "Dry-run grompp NPT OK"

    # MD producción
    link_setup_files "$RUNDIR/03_production"
    cp -f ../00_setup/system.gro dryrun_system.gro

    local md_mdp="md_prod_temp_dry.mdp"
    cp "$RUNDIR/mdp_used/md_prod.mdp" "$md_mdp"
    sed -i "s/^nsteps.*/nsteps                   = ${PROD_NSTEPS}      ; ${PROD_NS} ns/" "$md_mdp"
    sed -i 's/^tc-grps.*/tc-grps                  = Protein_Ligand Solvent/' "$md_mdp"
    sed -i 's/^ref_t.*/ref_t                    = 300     300/' "$md_mdp"

    run_gmx grompp -f "$md_mdp" -c dryrun_system.gro -p topol.top -n index.ndx \
        -o md_dry.tpr -maxwarn "$MAXWARN" &> "$RUNDIR/logs/dryrun_grompp_md.log"
    log_success "Dry-run grompp producción OK"
    log_warning "Dry-run activo: no se ejecutó ningún mdrun"
}

#==========================================
# ANÁLISIS POST-PRODUCCIÓN
#==========================================
run_analysis() {
    log_step "Iniciando análisis post-producción"

    cd "$RUNDIR/04_analysis" || exit 1

    # Copiar archivos de producción (cp en vez de ln: symlinks no funcionan en NTFS/WSL)
    for f in md.tpr md.xtc md.edr; do
        [ -f "../03_production/$f" ] && cp -f "../03_production/$f" .
    done
    cp -f ../00_setup/index.ndx .

    # RMSD inicial
    echo "Backbone Backbone" | run_gmx rms -s md.tpr -f md.xtc -o rmsd_backbone_raw.xvg \
        -tu ns &> "$RUNDIR/logs/analysis_rmsd_raw.log"
    log_success "RMSD inicial: rmsd_backbone_raw.xvg"

    # Procesamiento de trayectorias
    echo -e "\n${YELLOW}Procesando trayectorias...${NC}"

    echo "Protein System" | run_gmx trjconv -s md.tpr -f md.xtc -o md_clean_temp.xtc \
        -pbc nojump -ur compact -center &> "$RUNDIR/logs/analysis_trjconv1.log"
    log_success "Trayectoria centrada: md_clean_temp.xtc"

    echo "Protein System" | run_gmx trjconv -s md.tpr -f md_clean_temp.xtc \
        -o md_clean_full.xtc -fit rot+trans &> "$RUNDIR/logs/analysis_trjconv2.log"
    log_success "Trayectoria completa alineada: md_clean_full.xtc"

    echo "Protein non-Water" | run_gmx trjconv -s md.tpr -f md_clean_temp.xtc \
        -o md_clean_nowat.xtc -fit rot+trans &> "$RUNDIR/logs/analysis_trjconv3.log"
    log_success "Trayectoria sin agua: md_clean_nowat.xtc"

    echo "Protein non-Water" | run_gmx trjconv -s md.tpr -f md_clean_temp.xtc \
        -o md_frame0_nowat.pdb -pbc nojump -ur compact -center -b 0 -e 0 \
        &> "$RUNDIR/logs/analysis_frame0.log"
    log_success "Primer frame: md_frame0_nowat.pdb"

    # Último frame: obtener tiempo final dinámicamente
    local end_time
    end_time=$("$GMX" check -f md_clean_temp.xtc 2>&1 \
        | awk '/Last frame/ {for (i=1;i<=NF;i++) if ($i=="time") {print $(i+1); exit}}' || echo "")
    if [ -z "$end_time" ]; then
        end_time=1000000000
        log_warning "No se pudo determinar tiempo final, usando fallback para último frame"
    fi

    echo "Protein Protein" | run_gmx trjconv -s md.tpr -f md_clean_temp.xtc \
        -o md_protein_lastframe.pdb -pbc nojump -ur compact -center \
        -dump "$end_time" &> "$RUNDIR/logs/analysis_lastframe.log"
    log_success "Último frame: md_protein_lastframe.pdb"

    # TPR sin agua
    echo "non-Water" | run_gmx convert-tpr -s md.tpr -o tpr_nowat.tpr \
        &> "$RUNDIR/logs/analysis_tpr_nowat.log"
    log_success "TPR sin agua: tpr_nowat.tpr"

    # Filtrado
    echo "Protein" | run_gmx filter -s tpr_nowat.tpr -f md_clean_nowat.xtc \
        -ol md_clean_nowat_filtered.xtc -all -fit &> "$RUNDIR/logs/analysis_filter.log"
    log_success "Trayectoria filtrada: md_clean_nowat_filtered.xtc"

    # RMSD sin agua
    echo -e "\n${YELLOW}Calculando métricas estructurales...${NC}"

    echo "Backbone Backbone" | run_gmx rms -s tpr_nowat.tpr -f md_clean_nowat.xtc \
        -o rmsd_backbone.xvg -tu ns &> "$RUNDIR/logs/analysis_rmsd_backbone.log"
    log_success "RMSD backbone: rmsd_backbone.xvg"

    echo "Protein Protein" | run_gmx rms -s tpr_nowat.tpr -f md_clean_nowat.xtc \
        -o rmsd_protein.xvg -tu ns &> "$RUNDIR/logs/analysis_rmsd_protein.log"
    log_success "RMSD proteína: rmsd_protein.xvg"

    # Radio de giro
    echo "Backbone" | run_gmx gyrate -s tpr_nowat.tpr -f md_clean_nowat.xtc \
        -o gyrate.xvg &> "$RUNDIR/logs/analysis_gyrate.log"
    log_success "Radio de giro: gyrate.xvg"

    # RMSF
    echo "Backbone" | run_gmx rmsf -s tpr_nowat.tpr -f md_clean_nowat.xtc \
        -o rmsf_residue.xvg -res &> "$RUNDIR/logs/analysis_rmsf.log"
    log_success "RMSF: rmsf_residue.xvg"

    # Puentes de hidrógeno (GROMACS 2025: nueva sintaxis con -r y -t)
    echo -e "\n${YELLOW}Analizando interacciones...${NC}"

    run_gmx hbond -s tpr_nowat.tpr -f md_clean_nowat.xtc \
        -r 'protein' -t 'protein' \
        -num hbond_protein.xvg -tu ns &> "$RUNDIR/logs/analysis_hbond_protein.log"
    log_success "H-bonds proteína: hbond_protein.xvg"

    run_gmx hbond -s tpr_nowat.tpr -f md_clean_nowat.xtc \
        -r 'protein' -t 'resname LIG' \
        -num hbond_protein_ligand.xvg -tu ns &> "$RUNDIR/logs/analysis_hbond_prot_lig.log"
    log_success "H-bonds proteína-ligando: hbond_protein_ligand.xvg"

    run_gmx hbond -s tpr_nowat.tpr -f md_clean_nowat.xtc \
        -r 'resname LIG' -t 'resname LIG' \
        -num hbond_ligand.xvg -tu ns &> "$RUNDIR/logs/analysis_hbond_ligand.log" 2>/dev/null || true
    log_success "H-bonds ligando: hbond_ligand.xvg"

    # Estructura secundaria (GROMACS 2025: gmx dssp integrado, no requiere mkdssp externo)
    echo -e "\n${YELLOW}Analizando estructura secundaria...${NC}"

    if run_gmx dssp -s tpr_nowat.tpr -f md_clean_nowat.xtc \
        -o dssp.dat -tu ns &> "$RUNDIR/logs/analysis_dssp.log"; then
        log_success "Estructura secundaria: dssp.dat (gmx dssp integrado)"
    else
        log_warning "gmx dssp falló o no disponible - análisis de estructura secundaria omitido"
    fi

    # ==========================================
    # ANÁLISIS EXTENDIDO: LIGANDO Y COMPLEJO
    # ==========================================
    echo -e "\n${YELLOW}Análisis extendido del complejo proteína-ligando...${NC}"

    # 1. RMSD del ligando (¿se mantiene en el sitio de unión?)
    echo "LIG LIG" | run_gmx rms -s tpr_nowat.tpr -f md_clean_nowat.xtc \
        -o rmsd_ligand.xvg -tu ns &> "$RUNDIR/logs/analysis_rmsd_ligand.log"
    log_success "RMSD ligando: rmsd_ligand.xvg"

    # 2. Distancia mínima proteína-ligando (detecta disociación)
    echo "Protein LIG" | run_gmx mindist -s tpr_nowat.tpr -f md_clean_nowat.xtc \
        -od mindist_prot_lig.xvg -tu ns &> "$RUNDIR/logs/analysis_mindist.log"
    log_success "Distancia mínima prot-lig: mindist_prot_lig.xvg"

    # 3. SASA - Superficie accesible al solvente
    echo "Protein" | run_gmx sasa -s tpr_nowat.tpr -f md_clean_nowat.xtc \
        -o sasa_protein.xvg -tu ns &> "$RUNDIR/logs/analysis_sasa.log"
    log_success "SASA proteína: sasa_protein.xvg"

    # 4. PCA - Análisis de componentes principales (movimientos dominantes)
    echo -e "\n${YELLOW}Análisis de componentes principales (PCA)...${NC}"

    echo "Backbone Backbone" | run_gmx covar -s tpr_nowat.tpr -f md_clean_nowat.xtc \
        -o eigenval.xvg -v eigenvec.trr &> "$RUNDIR/logs/analysis_covar.log"
    log_success "Covarianza calculada: eigenval.xvg, eigenvec.trr"

    echo "Backbone Backbone" | run_gmx anaeig -s tpr_nowat.tpr -f md_clean_nowat.xtc \
        -v eigenvec.trr -eig eigenval.xvg -proj proj_pc1.xvg \
        -first 1 -last 1 &> "$RUNDIR/logs/analysis_pc1.log"
    log_success "Proyección PC1: proj_pc1.xvg"

    echo "Backbone Backbone" | run_gmx anaeig -s tpr_nowat.tpr -f md_clean_nowat.xtc \
        -v eigenvec.trr -eig eigenval.xvg -proj proj_pc2.xvg \
        -first 2 -last 2 &> "$RUNDIR/logs/analysis_pc2.log"
    log_success "Proyección PC2: proj_pc2.xvg"

    # 4b. Proyección 2D (PC1 vs PC2) para Free Energy Landscape
    echo "Backbone Backbone" | run_gmx anaeig -s tpr_nowat.tpr -f md_clean_nowat.xtc \
        -v eigenvec.trr -eig eigenval.xvg -2d proj_2d.xvg \
        -first 1 -last 2 &> "$RUNDIR/logs/analysis_2d.log"
    log_success "Proyección 2D (PC1 vs PC2): proj_2d.xvg"

    # 4c. DCCM: Extraer coordenadas C-alpha para mapa de correlación dinámica
    echo "C-alpha" | run_gmx traj -s tpr_nowat.tpr -f md_clean_nowat.xtc \
        -ox ca_positions.xvg &> "$RUNDIR/logs/analysis_ca_coords.log"
    log_success "Coordenadas C-alpha: ca_positions.xvg"

    # 5. Cluster analysis (conformaciones representativas)
    echo -e "\n${YELLOW}Análisis de clusters...${NC}"

    echo "Backbone System" | run_gmx cluster -s tpr_nowat.tpr -f md_clean_nowat.xtc \
        -cl clusters.pdb -cutoff 0.2 -method gromos \
        -g cluster.log &> "$RUNDIR/logs/analysis_cluster.log"
    log_success "Clusters: clusters.pdb, cluster.log"

    # 6. RMSD matrix (mapa de cambios conformacionales)
    echo "Backbone Backbone" | run_gmx rms -s tpr_nowat.tpr -f md_clean_nowat.xtc \
        -m rmsd_matrix.xpm -tu ns &> "$RUNDIR/logs/analysis_rmsd_matrix.log"
    log_success "RMSD matrix: rmsd_matrix.xpm"

    if run_gmx xpm2ps -f rmsd_matrix.xpm -o rmsd_matrix.eps \
        &> "$RUNDIR/logs/analysis_rmsd_matrix_ps.log"; then
        log_success "RMSD matrix gráfico: rmsd_matrix.eps"
    fi

    # 7. Contactos proteína-ligando por residuo (qué residuos interactúan más)
    echo "Protein LIG" | run_gmx mindist -s tpr_nowat.tpr -f md_clean_nowat.xtc \
        -on contacts_prot_lig.xvg -d 0.4 -group -tu ns \
        &> "$RUNDIR/logs/analysis_contacts.log"
    log_success "Contactos por residuo: contacts_prot_lig.xvg"

    # ==========================================
    # ANÁLISIS EXTENDIDO v4.0
    # ==========================================
    echo -e "\n${YELLOW}Análisis extendido v4.0...${NC}"

    # 8. Energía de interacción proteína-ligando (Coulomb + LJ)
    echo -e "\n${YELLOW}Calculando energía de interacción proteína-ligando...${NC}"
    echo "Coul-SR:Protein-LIG LJ-SR:Protein-LIG" | \
        run_gmx energy -f md.edr -o interaction_energy.xvg \
        &> "$RUNDIR/logs/analysis_interaction_energy.log" 2>&1 || true
    if [ -f "interaction_energy.xvg" ]; then
        log_success "Energía de interacción: interaction_energy.xvg"
    else
        log_warning "No se pudo extraer energía de interacción (grupos de energía no disponibles en MDP)"
    fi

    # 9. Contactos nativos (fraction of native contacts - Q)
    echo "Protein LIG" | run_gmx mindist -s tpr_nowat.tpr -f md_clean_nowat.xtc \
        -od native_contacts_dist.xvg -d 0.45 -tu ns \
        &> "$RUNDIR/logs/analysis_native_contacts.log" 2>&1 || true
    log_success "Contactos nativos (d<0.45nm): native_contacts_dist.xvg"

    # 10. Convergencia: RMSD block averaging (evaluación de equilibrio)
    echo -e "\n${YELLOW}Evaluando convergencia...${NC}"
    if [ -f "rmsd_backbone.xvg" ]; then
        # Generar RMSD de segunda mitad vs primera mitad para evaluar convergencia
        local total_frames
        total_frames=$(grep -cv '^[#@]' rmsd_backbone.xvg || echo 0)
        local half_frame=$(( total_frames / 2 ))
        if [ "$half_frame" -gt 10 ]; then
            local half_time
            half_time=$(grep -v '^[#@]' rmsd_backbone.xvg | awk -v hf="$half_frame" 'NR==hf {print $1; exit}')
            if [ -n "$half_time" ]; then
                echo "Backbone Backbone" | run_gmx rms -s tpr_nowat.tpr -f md_clean_nowat.xtc \
                    -o rmsd_backbone_2ndhalf.xvg -tu ns -b "$half_time" \
                    &> "$RUNDIR/logs/analysis_rmsd_convergence.log" 2>&1 || true
                log_success "RMSD segunda mitad: rmsd_backbone_2ndhalf.xvg (evaluación convergencia)"
            fi
        fi
    fi
}

#==========================================
# GENERAR RESUMEN
#==========================================
generate_summary() {
    log_step "Generando resumen final"

    local run_date host_name gmx_version cpu_model cpu_cores gpu_info exec_mode
    run_date=$(date '+%Y-%m-%d %H:%M:%S')
    host_name=$(hostname 2>/dev/null || echo "desconocido")
    gmx_version=$($GMX --version 2>/dev/null | awk -F': *' '/^GROMACS version/ {print $2; exit}')
    gmx_version=${gmx_version:-"no-detectada"}
    cpu_model=$(lscpu 2>/dev/null | awk -F': *' '/Model name/ {print $2; exit}')
    cpu_model=${cpu_model:-"no-detectado"}
    cpu_cores=$(nproc 2>/dev/null || echo "no-detectado")
    if command -v nvidia-smi &> /dev/null; then
        gpu_info=$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | awk 'NR==1{print; exit}')
        gpu_info=${gpu_info:-"detectada pero no disponible"}
    else
        gpu_info="no-detectada"
    fi
    exec_mode=$([ "$DRY_RUN" = true ] && echo "dry-run" || echo "normal")

    cat > "$RUNDIR/SUMMARY.txt" <<EOF
====================================
RESUMEN DE SIMULACIÓN MD
====================================
Proteína:     $PROT
Ligando:      $LIG
Fecha:        $run_date
Directorio:   $RUNDIR

TRAZABILIDAD DEL ENTORNO:
- Modo ejecución:       $exec_mode
- Hostname:             $host_name
- GROMACS binario:      $(command -v "$GMX")
- GROMACS versión:      $gmx_version
- CPU modelo:           $cpu_model
- CPU cores (nproc):    $cpu_cores
- GPU:                  $gpu_info

PARÁMETROS DE SIMULACIÓN:
- GROMACS:             $GMX $([ "$USE_MPI" = true ] && echo "(MPI)" || echo "(serial)")
- Force field:         $FF_DIR (local=$FF_IS_LOCAL)
- Tipo de caja:        $BOX_TYPE
- Distancia a bordes:  $BOX_DIST nm
- Modelo de agua:      $WATER_MODEL
- Concentración iónica:$ION_CONC M NaCl
- Producción objetivo:  $PROD_NS ns ($PROD_NSTEPS steps)
- Número de threads:   $NT

RESTRICCIONES POSICIONALES:
- Minimización (em):      Sin restricciones
- NVT:                    -DPOSRES -DPOSRES_LIG
- NPT:                    -DPOSRES -DPOSRES_LIG
- Producción (md_prod):   Sin restricciones

GRUPOS DE TERMOSTATO:
- Protein_Ligand (proteína + ligando)
- Solvent (agua + iones)

ESTRUCTURA DE CARPETAS:
00_setup/          : Configuración inicial
01_minimization/   : Minimización de energía
02_equilibration/  : NVT y NPT
03_production/     : MD de producción
04_analysis/       : Análisis post-producción
logs/              : Logs completos
mdp_used/          : Archivos .mdp originales

ARCHIVOS CLAVE:
Producción:
- 03_production/md.tpr : Sistema de producción
- 03_production/md.xtc : Trayectoria comprimida
- 03_production/md.edr : Energías

Análisis:
- 04_analysis/md_clean_nowat.xtc         : Trayectoria sin agua (análisis)
- 04_analysis/md_clean_nowat_filtered    : Trayectoria filtrada
- 04_analysis/rmsd_backbone.xvg          : RMSD backbone
- 04_analysis/rmsd_ligand.xvg            : RMSD del ligando
- 04_analysis/gyrate.xvg                 : Radio de giro
- 04_analysis/rmsf_residue.xvg           : Fluctuaciones por residuo
- 04_analysis/hbond_protein_ligand.xvg   : H-bonds proteína-ligando
- 04_analysis/mindist_prot_lig.xvg       : Distancia mínima prot-lig
- 04_analysis/sasa_protein.xvg           : Superficie accesible al solvente
- 04_analysis/eigenval.xvg               : Eigenvalores PCA
- 04_analysis/proj_pc1.xvg               : Proyección componente principal 1
- 04_analysis/proj_pc2.xvg               : Proyección componente principal 2
- 04_analysis/clusters.pdb               : Conformaciones representativas
- 04_analysis/rmsd_matrix.xpm            : Matriz RMSD
- 04_analysis/contacts_prot_lig.xvg      : Contactos por residuo

PRÓXIMOS PASOS:
1. Verificar convergencia: rmsd_backbone.xvg
2. Identificar regiones flexibles: rmsf_residue.xvg
3. Analizar estabilidad: gyrate.xvg
4. Evaluar H-bonds: hbond_protein_ligand.xvg
5. Visualizar: PyMOL/VMD con md_clean_nowat.xtc

NOTA DRY-RUN:
- Si Modo ejecución = dry-run, no se ejecutaron etapas largas de mdrun ni análisis.
- Revise logs dryrun_grompp_* para confirmar topología y preparación.
EOF

    cat "$RUNDIR/SUMMARY.txt"
    log_success "Resumen guardado: SUMMARY.txt"
}

#==========================================
# SISTEMA DE CHECKPOINTS
#==========================================
save_checkpoint() {
    local step_num=$1
    local step_name=$2
    cat > "$RUNDIR/.checkpoint" <<CKPT
# Checkpoint generado automáticamente - NO EDITAR
LAST_STEP="$step_num"
LAST_STEP_NAME="$step_name"
PROT="$PROT"
LIG="$LIG"
BOX_TYPE="$BOX_TYPE"
BOX_DIST="$BOX_DIST"
ION_CONC="$ION_CONC"
WATER_MODEL="$WATER_MODEL"
WATER_FILE="$WATER_FILE"
PROD_NS="$PROD_NS"
PROD_NSTEPS="$PROD_NSTEPS"
FF_DIR="$FF_DIR"
FF_IS_LOCAL="$FF_IS_LOCAL"
DRY_RUN="$DRY_RUN"
TOTAL_STEPS="$TOTAL_STEPS"
RUNDIR="$RUNDIR"
CKPT
}

load_checkpoint() {
    local ckpt_file="$1"
    if [ ! -f "$ckpt_file" ]; then
        log_error "Checkpoint no encontrado: $ckpt_file"
        exit 1
    fi

    if ! parse_checkpoint_file "$ckpt_file"; then
        log_error "Checkpoint inválido o inseguro: $ckpt_file"
        exit 1
    fi

    RESUME_FROM=$((LAST_STEP + 1))

    log_success "Checkpoint cargado:"
    log_info "  Proteína: $PROT"
    log_info "  Ligando: $LIG"
    log_info "  Último paso completado: $LAST_STEP ($LAST_STEP_NAME)"
    log_info "  Reanudando desde paso: $RESUME_FROM"
    log_info "  Directorio: $RUNDIR"
}

find_checkpoints() {
    local runs=()
    local infos=()

    if [ -d "$INITIAL_DIR/$WORKDIR" ]; then
        for dir in "$INITIAL_DIR/$WORKDIR"/*/; do
            if [ -f "${dir}.checkpoint" ]; then
                local ckpt_file="${dir}.checkpoint"
                local ckpt_last_step
                local ckpt_last_name
                local ckpt_total_steps
                ckpt_last_step=$(extract_checkpoint_value "$ckpt_file" "LAST_STEP")
                ckpt_last_name=$(extract_checkpoint_value "$ckpt_file" "LAST_STEP_NAME")
                ckpt_total_steps=$(extract_checkpoint_value "$ckpt_file" "TOTAL_STEPS")
                ckpt_total_steps=${ckpt_total_steps:-$TOTAL_STEPS}
                [ -z "$ckpt_last_step" ] && continue
                # Solo mostrar si NO está completado
                if [ "${ckpt_last_name:-}" != "COMPLETADO" ] && [ "$ckpt_last_step" -lt "$ckpt_total_steps" ]; then
                    runs+=("$ckpt_file")
                    infos+=("$(basename "$dir") (detenido en: ${ckpt_last_name:-desconocido}, paso $ckpt_last_step/$ckpt_total_steps)")
                fi
            fi
        done
    fi

    if [ ${#runs[@]} -eq 0 ]; then
        log_error "No se encontraron corridas incompletas con checkpoint"
        exit 1
    fi

    if [ ${#runs[@]} -eq 1 ]; then
        load_checkpoint "${runs[0]}"
        return
    fi

    echo -e "\n${BLUE}Corridas con checkpoint disponible:${NC}\n"
    for i in "${!infos[@]}"; do
        echo "  $((i + 1))) ${infos[$i]}"
    done
    echo ""
    echo "¿Cuál reanudar? [1-${#runs[@]}]:"
    read -r RESUME_CHOICE

    if [[ "$RESUME_CHOICE" =~ ^[0-9]+$ ]] && [ "$RESUME_CHOICE" -ge 1 ] && [ "$RESUME_CHOICE" -le ${#runs[@]} ]; then
        load_checkpoint "${runs[$((RESUME_CHOICE - 1))]}"
    else
        log_error "Selección inválida"
        exit 1
    fi
}

run_step() {
    local step_num=$1
    local step_name=$2
    local step_func=$3

    # En modo resume, saltar pasos ya completados
    if [ "${RESUME_MODE:-false}" = true ] && [ "$step_num" -lt "$RESUME_FROM" ]; then
        log_warning "Paso $step_num ($step_name) ya completado - saltando"
        return 0
    fi

    log_step "Paso $step_num/$TOTAL_STEPS: $step_name"
    "$step_func"

    # Guardar checkpoint después de cada paso exitoso
    if [ -n "${RUNDIR:-}" ] && [ -d "${RUNDIR:-}" ]; then
        save_checkpoint "$step_num" "$step_name"
    fi
}

#==========================================
# FUNCIÓN PRINCIPAL
#==========================================
main() {
    RESUME_MODE=false
    RESUME_FROM=1
    TOTAL_STEPS=14

    parse_args "$@"

    if [ "$SHOW_HELP" = true ]; then
        print_usage
        exit 0
    fi

    init_gmx
    validate_dependencies

    # v4.5: Modo analysis-only (re-ejecutar análisis sobre corrida existente)
    if [ "$ANALYSIS_ONLY" = true ]; then
        run_analysis_only
        exit 0
    fi

    # Modo resume
    if [ "$RESUME_MODE" = true ]; then
        log_step "Modo RESUME activado"
        find_checkpoints

        # v4.5: Si se pidió extend junto con resume
        if [ -n "$EXTEND_NS" ]; then
            extend_simulation
            # Re-ejecutar análisis después de extender
            run_step 12 "Análisis post-producción"     run_analysis
            run_step 13 "Generar resumen"              generate_summary
            save_checkpoint "$TOTAL_STEPS" "COMPLETADO"
            log_step "¡EXTENSIÓN + ANÁLISIS FINALIZADOS CON ÉXITO!"
            echo -e "${GREEN}Todos los resultados en: $RUNDIR${NC}\n"
            exit 0
        fi
    fi

    if [ "$RESUME_MODE" = false ]; then
        resolve_execution_mode
        validate_mdp_files
        setup_directory_structure
    fi

    if [ "$DRY_RUN" = true ]; then
        TOTAL_STEPS=10
        run_step 1  "Copiar archivos"                setup_initial_files
        run_step 2  "Validar ligando (v4.5)"         validate_ligand_files
        run_step 3  "Restraints del ligando"         generate_ligand_restraints
        run_step 4  "Actualizar topología"           update_topology
        run_step 5  "Ensamblar complejo"             build_complex
        run_step 6  "Solvatar sistema"               solvate_system
        run_step 7  "Neutralizar sistema"            neutralize_system
        run_step 8  "Crear grupos de índice"         create_index_groups
        run_step 9  "Validar grompp (dry-run)"       run_dry_run_grompp_checks
        run_step 10 "Generar resumen"                generate_summary
    else
        run_step 1  "Copiar archivos"              setup_initial_files
        run_step 2  "Validar ligando (v4.5)"       validate_ligand_files
        run_step 3  "Restraints del ligando"       generate_ligand_restraints
        run_step 4  "Actualizar topología"         update_topology
        run_step 5  "Ensamblar complejo"           build_complex
        run_step 6  "Solvatar sistema"             solvate_system
        run_step 7  "Neutralizar sistema"          neutralize_system
        run_step 8  "Crear grupos de índice"       create_index_groups
        run_step 9  "Minimización de energía"      run_minimization
        run_step 10 "Equilibración NVT"            run_nvt_equilibration
        run_step 11 "Equilibración NPT"            run_npt_equilibration
        run_step 12 "Producción"                   run_production
        run_step 13 "Análisis post-producción"     run_analysis
        run_step 14 "Generar resumen"              generate_summary
    fi

    # Marcar como completado
    save_checkpoint "$TOTAL_STEPS" "COMPLETADO"

    log_step "¡SIMULACIÓN FINALIZADA CON ÉXITO!"
    echo -e "${GREEN}Todos los resultados en: $RUNDIR${NC}"
    echo -e "${GREEN}Para generar gráficas: python3 plot_analysis.py $RUNDIR${NC}\n"
}

# Ejecutar
main "$@"
