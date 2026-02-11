#!/bin/bash
set -euo pipefail

#==========================================
# CONFIGURACIÓN
#==========================================
# Auto-detectar gmx o gmx_mpi
if command -v gmx &> /dev/null; then
    readonly GMX=gmx
    readonly USE_MPI=false
elif command -v gmx_mpi &> /dev/null; then
    readonly GMX=gmx_mpi
    readonly USE_MPI=true
else
    echo "ERROR: No se encontró GROMACS (gmx o gmx_mpi)" >&2
    exit 1
fi

readonly NT=16
# Para gmx_mpi, mdrun se ejecuta diferente
if [ "$USE_MPI" = true ]; then
    readonly MDRUN="$GMX mdrun"
else
    readonly MDRUN="$GMX mdrun -nt $NT"
fi

readonly BASE_PROT=proteinas
readonly BASE_LIG=ligandos
readonly BASE_MDP=mdp
readonly WORKDIR=MD_RUN
readonly INITIAL_DIR="$(pwd)"

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

#==========================================
# FUNCIONES AUXILIARES
#==========================================
log_step() {
    CURRENT_STAGE="$1"
    echo -e "\n${BLUE}=========================================${NC}"
    echo -e "${BLUE}>>> $1${NC}"
    echo -e "${BLUE}=========================================${NC}\n"
}

log_success() {
    echo -e "${GREEN}✓${NC} $1"
}

log_error() {
    echo -e "${RED}❌ ERROR:${NC} $1" >&2
}

log_warning() {
    echo -e "${YELLOW}⚠${NC} $1"
}

log_info() {
    echo -e "${CYAN}ℹ${NC} $1"
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

# Wrapper simplificado para ejecutar comandos GMX
run_gmx() {
    "$GMX" "$@"
}

#==========================================
# TRAP PARA LIMPIEZA EN CASO DE ERROR
#==========================================
cleanup_on_error() {
    local exit_code=$?
    if [ $exit_code -ne 0 ]; then
        echo ""
        log_error "Script interrumpido durante: ${CURRENT_STAGE}"
        log_warning "Archivos parciales en: ${RUNDIR:-desconocido}"
        log_warning "Revisa los logs en: ${RUNDIR:-desconocido}/logs/"
    fi
}
trap cleanup_on_error EXIT

#==========================================
# VALIDACIÓN DE DEPENDENCIAS
#==========================================
validate_dependencies() {
    log_step "Validando dependencias del sistema"
    log_success "GROMACS encontrado: $(which $GMX) $([ "$USE_MPI" = true ] && echo "(MPI)" || echo "(serial)")"
}

#==========================================
# INPUT DEL USUARIO CON VALIDACIÓN
#==========================================
get_user_input() {
    echo -e "${BLUE}=========================================${NC}"
    echo -e "${BLUE}  SIMULACIÓN DE DINÁMICA MOLECULAR${NC}"
    echo -e "${BLUE}=========================================${NC}\n"

    # Proteína
    while true; do
        echo "Nombre de la proteína (carpeta en proteinas/):"
        read -r PROT
        if [ -d "$BASE_PROT/$PROT" ]; then
            break
        else
            log_error "Proteína '$PROT' no encontrada en $BASE_PROT/"
        fi
    done

    # Ligando
    while true; do
        echo -e "\nNombre del ligando (carpeta en ligandos/):"
        read -r LIG
        if [ -d "$BASE_LIG/$LIG" ]; then
            break
        else
            log_error "Ligando '$LIG' no encontrado en $BASE_LIG/"
        fi
    done

    # Tipo de caja
    echo -e "\nTipo de caja de simulación:"
    echo "  1) cubic       - Caja cúbica (más volumen, más moléculas de agua)"
    echo "  2) triclinic   - Romboedro truncado (~71% del volumen de cubic)"
    echo "  3) dodecahedron - Dodecaedro rómbico (~71%, recomendado)"
    echo "  4) octahedron  - Octaedro truncado (~77% del volumen de cubic)"
    echo -e "\nSeleccione el tipo de caja [1-4] (default: 3 dodecahedron):"
    read -r BOX_CHOICE

    case $BOX_CHOICE in
        1) BOX_TYPE="cubic" ;;
        2) BOX_TYPE="triclinic" ;;
        4) BOX_TYPE="octahedron" ;;
        *) BOX_TYPE="dodecahedron" ;;
    esac

    # Distancia a bordes
    echo -e "\nDistancia mínima a los bordes de la caja (nm) [default: 1.2]:"
    read -r BOX_DIST
    BOX_DIST=${BOX_DIST:-1.2}

    # Modelo de agua
    echo -e "\nModelo de agua a utilizar:"
    echo "  1) tip3p  - TIP3P (rápido, recomendado para la mayoría)"
    echo "  2) spc    - SPC (simple point charge)"
    echo "  3) spce   - SPC/E (extended simple point charge)"
    echo "  4) tip4p  - TIP4P (4 sitios)"
    echo "  5) tip5p  - TIP5P (5 sitios, más preciso pero más lento)"
    echo -e "\nSeleccione el modelo de agua [1-5] (default: 1):"
    read -r WATER_CHOICE

    case $WATER_CHOICE in
        2) WATER_MODEL="spc"; WATER_FILE="spc216.gro" ;;
        3) WATER_MODEL="spce"; WATER_FILE="spc216.gro" ;;
        4) WATER_MODEL="tip4p"; WATER_FILE="tip4p.gro" ;;
        5) WATER_MODEL="tip5p"; WATER_FILE="tip5p.gro" ;;
        *) WATER_MODEL="tip3p"; WATER_FILE="spc216.gro" ;;
    esac

    # Tiempo de producción
    echo -e "\nTiempo de simulación de producción en nanosegundos [default: 10]:"
    read -r PROD_NS
    PROD_NS=${PROD_NS:-10}
    # dt = 0.002 ps → 500000 steps por ns
    PROD_NSTEPS=$((PROD_NS * 500000))
    log_info "Producción: $PROD_NS ns ($PROD_NSTEPS steps)"
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
    log_success "Archivos MDP validados"
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

    # Copiar CHARMM36.ff si existe
    if [ -d "$INITIAL_DIR/charmm36.ff" ]; then
        log_step "Copiando CHARMM36.ff a 00_setup"
        cp -r "$INITIAL_DIR/charmm36.ff" .
        log_success "CHARMM36.ff copiado a 00_setup/"
    fi

    log_success "Archivos copiados a 00_setup/"
}

#==========================================
# RESTRAINTS DEL LIGANDO
#==========================================
generate_ligand_restraints() {
    log_step "Generando restraints del ligando (sin H)"

    cd "$RUNDIR/00_setup" || exit 1

    # Determinar índice del nuevo grupo dinámicamente
    local last_idx
    last_idx=$(echo q | "$GMX" make_ndx -f ligando.gro 2>&1 | grep -oP '^\s*\K\d+(?=\s)' | tail -1)
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
update_topology() {
    log_step "Actualizando archivos de topología"

    cd "$RUNDIR/00_setup" || exit 1

    # Añadir restraints a ligando.itp si no existen
    if ! grep -q "POSRES_LIG" ligando.itp; then
        cat <<'EOF' >> ligando.itp

; Ligand position restraints
#ifdef POSRES_LIG
#include "posre_ligando.itp"
#endif
EOF
        log_success "Restraints añadidos a ligando.itp"
    fi

    # CRÍTICO: Modificar topol.top para incluir parámetros del ligando
    log_step "Modificando topol.top para incluir ligando y sus parámetros"

    # Crear backup
    cp topol.top topol.top.backup

    # Si existe ligando.prm, incluirlo PRIMERO (antes del ligando.itp)
    if [ -f "ligando.prm" ]; then
        log_success "Incluyendo parámetros del ligando (ligando.prm)"

        # Buscar la línea del forcefield para insertar DESPUÉS
        local ff_line=$(grep -n '#include.*forcefield\.itp"' topol.top | cut -d: -f1)

        if [ -n "$ff_line" ]; then
            sed -i "${ff_line}a\\
; Include ligand parameters (CHARMM format)\\
#include \"ligando.prm\"" topol.top
            log_success "ligando.prm incluido después del forcefield"
        else
            log_error "No se encontró línea de forcefield.itp en topol.top"
            exit 1
        fi
    fi

    # Ahora incluir ligando.itp
    if ! grep -q 'ligando.itp' topol.top; then
        if [ -f "ligando.prm" ]; then
            local insert_line=$(grep -n 'ligando\.prm' topol.top | cut -d: -f1)
        else
            local insert_line=$(grep -n '#include.*forcefield\.itp"' topol.top | cut -d: -f1)
        fi

        if [ -n "$insert_line" ]; then
            sed -i "${insert_line}a\\
; Include ligand topology\\
#include \"ligando.itp\"" topol.top
            log_success "ligando.itp incluido en topol.top"
        else
            log_error "No se pudo determinar dónde insertar ligando.itp"
            exit 1
        fi
    else
        log_warning "ligando.itp ya estaba incluido en topol.top"
    fi

    # Añadir molécula LIG en la sección [ molecules ]
    if ! grep -q '^LIG' topol.top; then
        if grep -q '^\[ molecules \]' topol.top; then
            local prot_name=$(awk '/^\[ molecules \]/{flag=1; next} flag && NF>0 && !/^;/{print $1; exit}' topol.top)
            sed -i "/^\[ molecules \]/,/^${prot_name}/s/^\(${prot_name}.*\)/\1\nLIG                 1/" topol.top
            log_success "Molécula LIG añadida a topol.top"
        else
            log_error "No se encontró sección [ molecules ] en topol.top"
            exit 1
        fi
    else
        log_warning "Molécula LIG ya estaba en topol.top"
    fi

    # Verificar el resultado
    log_step "Verificando topol.top modificado"
    echo "--- Includes encontrados ---"
    grep '#include' topol.top | head -10
    echo -e "\n--- Sección molecules ---"
    sed -n '/\[ molecules \]/,/^$/p' topol.top
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
    log_step "Añadiendo iones para neutralizar el sistema"

    cd "$RUNDIR/00_setup" || exit 1

    if ! run_gmx grompp -f "$RUNDIR/mdp_used/ions.mdp" -c solv.gro -p topol.top -o ions.tpr \
        -maxwarn 1 > "$RUNDIR/logs/grompp_ions.log" 2>&1; then
        log_error "grompp falló al generar ions.tpr"
        cat "$RUNDIR/logs/grompp_ions.log"
        exit 1
    fi

    echo "SOL" | run_gmx genion -s ions.tpr -o system.gro -p topol.top -neutral \
        &> "$RUNDIR/logs/genion.log"

    log_success "Sistema neutralizado"
}

#==========================================
# GRUPOS DE ÍNDICE
#==========================================
create_index_groups() {
    log_step "Creando grupos para termostatos"

    cd "$RUNDIR/00_setup" || exit 1

    # Determinar dinámicamente los índices de los nuevos grupos
    local last_idx
    last_idx=$(echo q | "$GMX" make_ndx -f system.gro 2>&1 | grep -oP '^\s*\K\d+(?=\s)' | tail -1)
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

    # Copiar charmm36.ff si existe
    if [ -d "../00_setup/charmm36.ff" ]; then
        cp -r ../00_setup/charmm36.ff .
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
        -o em.tpr -maxwarn 1 &> "$RUNDIR/logs/grompp_em.log"

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

    # CRÍTICO: Reemplazar grupos de temperatura
    sed -i 's/^tc-grps.*/tc-grps                  = Protein_Ligand Solvent/' "$nvt_mdp"
    sed -i 's/^tau_t.*/tau_t                    = 0.1     0.1/' "$nvt_mdp"
    sed -i 's/^ref_t.*/ref_t                    = 300     300/' "$nvt_mdp"

    log_success "Grupos de termostato actualizados en nvt_temp.mdp"

    run_gmx grompp -f "$nvt_mdp" -c em.gro -r em.gro -p topol.top -n index.ndx \
        -o nvt.tpr -maxwarn 1 &> "$RUNDIR/logs/grompp_nvt.log"

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

    # CRÍTICO: Reemplazar grupos de temperatura
    sed -i 's/^tc-grps.*/tc-grps                  = Protein_Ligand Solvent/' "$npt_mdp"
    sed -i 's/^tau_t.*/tau_t                    = 0.1     0.1/' "$npt_mdp"
    sed -i 's/^ref_t.*/ref_t                    = 300     300/' "$npt_mdp"

    log_success "Grupos de termostato actualizados en npt_temp.mdp"

    run_gmx grompp -f "$npt_mdp" -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top \
        -n index.ndx -o npt.tpr -maxwarn 1 &> "$RUNDIR/logs/grompp_npt.log"

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

    # Reemplazar grupos de temperatura
    sed -i 's/^tc-grps.*/tc-grps                  = Protein_Ligand Solvent/' "$md_mdp"
    sed -i 's/^tau_t.*/tau_t                    = 0.1     0.1/' "$md_mdp"
    sed -i 's/^ref_t.*/ref_t                    = 300     300/' "$md_mdp"

    log_success "Producción configurada: $PROD_NS ns ($PROD_NSTEPS steps)"

    run_gmx grompp -f "$md_mdp" -c npt.gro -t npt.cpt -p topol.top \
        -n index.ndx -o md.tpr -maxwarn 1 &> "$RUNDIR/logs/grompp_md.log"

    log_success "Sistema preparado para producción"
    echo -e "\n${YELLOW}Ejecutando producción (esto puede tardar)...${NC}\n"

    $MDRUN -deffnm md &> "$RUNDIR/logs/mdrun_md.log"
    log_success "Simulación de producción completada"
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
        | grep -oP 'Last frame\s+\d+\s+time\s+\K[\d.]+' || echo "")
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
}

#==========================================
# GENERAR RESUMEN
#==========================================
generate_summary() {
    log_step "Generando resumen final"

    cat > "$RUNDIR/SUMMARY.txt" <<EOF
====================================
RESUMEN DE SIMULACIÓN MD
====================================
Proteína:     $PROT
Ligando:      $LIG
Fecha:        $(date)
Directorio:   $RUNDIR

PARÁMETROS DE SIMULACIÓN:
- GROMACS:             $GMX $([ "$USE_MPI" = true ] && echo "(MPI)" || echo "(serial)")
- Tipo de caja:        $BOX_TYPE
- Distancia a bordes:  $BOX_DIST nm
- Modelo de agua:      $WATER_MODEL
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
WATER_MODEL="$WATER_MODEL"
WATER_FILE="$WATER_FILE"
PROD_NS="$PROD_NS"
PROD_NSTEPS="$PROD_NSTEPS"
RUNDIR="$RUNDIR"
CKPT
}

load_checkpoint() {
    local ckpt_file="$1"
    if [ ! -f "$ckpt_file" ]; then
        log_error "Checkpoint no encontrado: $ckpt_file"
        exit 1
    fi

    source "$ckpt_file"
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
                source "${dir}.checkpoint"
                # Solo mostrar si NO está completado (paso 14 = generate_summary)
                if [ "$LAST_STEP" -lt 14 ]; then
                    runs+=("${dir}.checkpoint")
                    infos+=("$(basename "$dir") (detenido en: $LAST_STEP_NAME, paso $LAST_STEP/14)")
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

    log_step "Paso $step_num/14: $step_name"
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

    # Parsear argumentos
    if [ "${1:-}" = "--resume" ] || [ "${1:-}" = "-r" ]; then
        RESUME_MODE=true
        log_step "Modo RESUME activado"
        find_checkpoints
    fi

    if [ "$RESUME_MODE" = false ]; then
        # Ejecución normal: pedir datos al usuario
        validate_dependencies
        get_user_input
        validate_mdp_files
        setup_directory_structure
    else
        validate_dependencies
    fi

    run_step 1  "Copiar archivos"              setup_initial_files
    run_step 2  "Restraints del ligando"        generate_ligand_restraints
    run_step 3  "Actualizar topología"          update_topology
    run_step 4  "Ensamblar complejo"            build_complex
    run_step 5  "Solvatar sistema"              solvate_system
    run_step 6  "Neutralizar sistema"           neutralize_system
    run_step 7  "Crear grupos de índice"        create_index_groups
    run_step 8  "Minimización de energía"       run_minimization
    run_step 9  "Equilibración NVT"             run_nvt_equilibration
    run_step 10 "Equilibración NPT"             run_npt_equilibration
    run_step 11 "Producción"                    run_production
    run_step 12 "Análisis post-producción"      run_analysis
    run_step 13 "Generar resumen"               generate_summary

    # Marcar como completado
    save_checkpoint 14 "COMPLETADO"

    log_step "¡SIMULACIÓN FINALIZADA CON ÉXITO!"
    echo -e "${GREEN}Todos los resultados en: $RUNDIR${NC}"
    echo -e "${GREEN}Para generar gráficas: python3 plot_analysis.py $RUNDIR${NC}\n"
}

# Ejecutar
main "$@"