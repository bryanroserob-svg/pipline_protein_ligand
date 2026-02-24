#!/usr/bin/env python3
"""
Generador automático de gráficas de análisis MD v4.0.
Lee archivos .xvg de una corrida MD y genera un PDF con todas las figuras.

Novedades v4.0:
  - CLI flags: --format, --output-dir, --no-pca, --no-dccm
  - B-factors desde RMSF
  - Tabla resumen de estadísticas en el PDF
  - Evaluación de convergencia (block averaging)
  - Energía de interacción proteína-ligando
  - Stride automático para DCCM en trayectorias largas
  - Contactos nativos

Uso:
    python3 plot_analysis.py                          # Menú interactivo
    python3 plot_analysis.py MD_RUN/mi_corrida/       # Directo
    python3 plot_analysis.py --help                   # Mostrar opciones
"""

import sys
import os
import re
import argparse
import numpy as np
from pathlib import Path
from datetime import datetime

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
except ImportError:
    print("❌ ERROR: matplotlib o numpy no encontrado.")
    print("   Instalar con: pip install matplotlib numpy")
    sys.exit(1)

# ==========================================
# ESTILOS
# ==========================================
plt.rcParams.update({
    'figure.figsize': (10, 5),
    'font.size': 11,
    'font.family': 'sans-serif',
    'axes.titlesize': 13,
    'axes.titleweight': 'bold',
    'axes.labelsize': 11,
    'axes.grid': True,
    'grid.alpha': 0.3,
    'lines.linewidth': 1.2,
    'figure.dpi': 150,
})

C = {
    'blue': '#2563EB', 'red': '#DC2626', 'green': '#059669',
    'orange': '#D97706', 'purple': '#7C3AED', 'teal': '#0D9488',
}


# ==========================================
# PARSER XVG
# ==========================================
def parse_xvg(filepath):
    """Lee un archivo .xvg. Retorna (data_numpy, titulo, xlabel, ylabel, legends)."""
    data, title, xlabel, ylabel, legends = [], "", "", "", []

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#'):
                continue
            elif line.startswith('@'):
                m = re.search(r'"(.+)"', line)
                val = m.group(1) if m else ""
                low = line.lower()
                if 'title' in low and 'sub' not in low:
                    title = val
                elif 'xaxis' in low and 'label' in low:
                    xlabel = val
                elif 'yaxis' in low and 'label' in low:
                    ylabel = val
                elif line.startswith('@ s') and 'legend' in low:
                    legends.append(val)
            elif line and not line.startswith('&'):
                try:
                    data.append([float(x) for x in line.split()])
                except ValueError:
                    continue

    return (np.array(data) if data else None, title, xlabel, ylabel, legends)


# ==========================================
# FUNCIONES DE PLOT
# ==========================================
def plot_timeseries(ax, data, title, xlabel, ylabel, legends=None, color=None):
    if data is None or len(data) == 0:
        return
    x = data[:, 0]
    if data.shape[1] == 2:
        ax.plot(x, data[:, 1], color=color or C['blue'], alpha=0.85)
    else:
        cols = list(C.values())
        for i in range(1, min(data.shape[1], 7)):
            lab = legends[i-1] if legends and len(legends) >= i else f"Col {i}"
            ax.plot(x, data[:, i], color=cols[(i-1) % len(cols)], alpha=0.85, label=lab)
        ax.legend(fontsize=9, framealpha=0.8)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)


def plot_bar(ax, data, title, xlabel, ylabel, color=None):
    if data is None or len(data) == 0:
        return
    x, y = data[:, 0], data[:, 1]
    ax.bar(x, y, width=0.8, color=color or C['blue'], alpha=0.7, edgecolor='none')
    mean_val = np.mean(y)
    ax.axhline(y=mean_val, color=C['red'], linestyle='--', alpha=0.7,
               label=f'Media: {mean_val:.3f}')
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.legend(fontsize=9)


def plot_matrix_from_xpm(ax, filepath):
    """Intenta parsear un .xpm para heatmap. Retorna True si exitoso."""
    values, color_map = [], {}
    try:
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if '"' in line and 'c #' in line:
                    parts = line.split('"')
                    if len(parts) >= 2:
                        char = parts[1][0]
                        vm = re.search(r'/\*\s*"([\d.eE+-]+)"', line)
                        if vm:
                            color_map[char] = float(vm.group(1))
                elif line.startswith('"') and '/*' not in line:
                    row_str = line.strip('",')
                    row = [color_map.get(c, 0.0) for c in row_str if c in color_map]
                    if row:
                        values.append(row)
    except Exception:
        return False

    if not values:
        return False

    matrix = np.array(values)
    im = ax.imshow(matrix, cmap='YlOrRd', aspect='auto', origin='lower')
    ax.set_title('Matriz RMSD')
    ax.set_xlabel('Tiempo (frames)')
    ax.set_ylabel('Tiempo (frames)')
    plt.colorbar(im, ax=ax, label='RMSD (nm)', shrink=0.8)
    return True


def plot_fel(filepath):
    """Genera Free Energy Landscape (FEL) desde proj_2d.xvg. Retorna (fig, True) o (None, False)."""
    data, _, _, _, _ = parse_xvg(filepath)
    if data is None or data.shape[1] < 2:
        return None, False

    pc1 = data[:, 0]
    pc2 = data[:, 1]

    # Histograma 2D
    nbins = 80
    H, xedges, yedges = np.histogram2d(pc1, pc2, bins=nbins)

    # Inversión de Boltzmann: ΔG = -kT * ln(P)
    # T = 300 K, kB = 8.314e-3 kJ/(mol·K)
    kT = 8.314e-3 * 300  # kJ/mol
    H = np.where(H > 0, H, np.nan)  # evitar log(0)
    G = -kT * np.log(H / np.nanmax(H))  # normalizar al mínimo = 0
    G = G.T  # transponer para que el eje Y sea PC2

    fig, ax = plt.subplots(figsize=(8, 7))
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    im = ax.imshow(G, origin='lower', extent=extent, aspect='auto',
                   cmap='RdYlBu_r', interpolation='gaussian')
    cbar = plt.colorbar(im, ax=ax, shrink=0.85, pad=0.02)
    cbar.set_label('ΔG (kJ/mol)', fontsize=12)
    ax.set_title('Free Energy Landscape (FEL)', fontsize=14, fontweight='bold')
    ax.set_xlabel('PC1 (nm)', fontsize=12)
    ax.set_ylabel('PC2 (nm)', fontsize=12)

    # Contornos
    x_centers = 0.5 * (xedges[:-1] + xedges[1:])
    y_centers = 0.5 * (yedges[:-1] + yedges[1:])
    G_contour = np.where(np.isnan(G), np.nanmax(G), G)
    ax.contour(x_centers, y_centers, G_contour, levels=10,
               colors='black', linewidths=0.4, alpha=0.5)

    # Marcar mínimo de energía
    min_idx = np.unravel_index(np.nanargmin(G), G.shape)
    ax.plot(x_centers[min_idx[1]], y_centers[min_idx[0]],
            '*', color='white', markersize=15, markeredgecolor='black',
            markeredgewidth=1.0, label='Mínimo global')
    ax.legend(fontsize=10, loc='upper right')

    fig.tight_layout()
    return fig, True


def plot_dccm(filepath):
    """Genera DCCM desde ca_positions.xvg. Retorna (fig, True) o (None, False)."""
    data, _, _, _, _ = parse_xvg(filepath)
    if data is None or data.shape[1] < 4:
        return None, False

    coords = data[:, 1:]  # Quitar columna de tiempo
    n_atoms = coords.shape[1] // 3
    n_frames = coords.shape[0]

    if n_atoms < 2:
        return None, False

    # Stride automático para trayectorias largas (v4.0)
    if n_frames > 5000:
        stride = n_frames // 5000
        coords = coords[::stride]
        n_frames = coords.shape[0]

    # Reshape a (n_frames, n_atoms, 3)
    coords_3d = coords.reshape(n_frames, n_atoms, 3)

    # Posiciones medias y fluctuaciones
    mean_pos = np.mean(coords_3d, axis=0)
    delta = coords_3d - mean_pos  # (n_frames, n_atoms, 3)

    # Matriz de covarianza: <delta_i · delta_j>
    cov = np.einsum('tix,tjx->ij', delta, delta) / n_frames

    # Normalizar: C_ij = cov_ij / sqrt(cov_ii * cov_jj)
    diag = np.sqrt(np.diag(cov))
    diag[diag == 0] = 1e-10  # evitar división por cero
    dccm = cov / np.outer(diag, diag)

    fig, ax = plt.subplots(figsize=(8, 7))
    im = ax.imshow(dccm, cmap='RdBu_r', vmin=-1, vmax=1,
                   origin='lower', aspect='equal')
    cbar = plt.colorbar(im, ax=ax, shrink=0.85, pad=0.02)
    cbar.set_label('C(i,j)', fontsize=12)
    ax.set_title('Dynamic Cross-Correlation Matrix (DCCM)', fontsize=14, fontweight='bold')
    ax.set_xlabel('Residuo (Cα)', fontsize=12)
    ax.set_ylabel('Residuo (Cα)', fontsize=12)

    # Línea diagonal
    ax.plot([0, n_atoms-1], [0, n_atoms-1], 'k-', linewidth=0.5, alpha=0.3)

    fig.tight_layout()
    return fig, True


# ==========================================
# NUEVAS FUNCIONES v4.0
# ==========================================
def compute_bfactors(rmsf_data):
    """Convierte RMSF (nm) a B-factors (Å²): B = (8π²/3) × RMSF²."""
    if rmsf_data is None or rmsf_data.shape[1] < 2:
        return None
    residues = rmsf_data[:, 0]
    rmsf_nm = rmsf_data[:, 1]
    # RMSF en nm → Å (x10), luego B = (8π²/3) × RMSF²
    rmsf_angstrom = rmsf_nm * 10.0
    bfactors = (8.0 * np.pi**2 / 3.0) * rmsf_angstrom**2
    return np.column_stack([residues, bfactors])


def plot_bfactors(ax, bfactor_data):
    """Plot B-factors calculados desde RMSF."""
    if bfactor_data is None:
        return
    x, y = bfactor_data[:, 0], bfactor_data[:, 1]
    ax.bar(x, y, width=0.8, color=C['teal'], alpha=0.7, edgecolor='none')
    mean_val = np.mean(y)
    ax.axhline(y=mean_val, color=C['red'], linestyle='--', alpha=0.7,
               label=f'Media: {mean_val:.1f} Å²')
    ax.set_title('B-factors (calculados desde RMSF)')
    ax.set_xlabel('Residuo')
    ax.set_ylabel('B-factor (Å²)')
    ax.legend(fontsize=9)


def assess_convergence(rmsd_data):
    """Evalúa convergencia con block averaging del RMSD. Retorna dict con métricas."""
    if rmsd_data is None or rmsd_data.shape[1] < 2:
        return None

    y = rmsd_data[:, 1]
    n = len(y)
    if n < 20:
        return None

    # Dividir en mitades
    first_half = y[:n//2]
    second_half = y[n//2:]

    result = {
        'total_mean': np.mean(y),
        'total_std': np.std(y),
        'first_half_mean': np.mean(first_half),
        'second_half_mean': np.mean(second_half),
        'drift': abs(np.mean(second_half) - np.mean(first_half)),
        'converged': abs(np.mean(second_half) - np.mean(first_half)) < np.std(y),
    }
    return result


def generate_summary_table(stats_dict):
    """Genera una figura con tabla resumen de estadísticas."""
    if not stats_dict:
        return None

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.axis('off')

    rows = []
    for name, vals in stats_dict.items():
        if vals is not None:
            rows.append([name, f"{vals['mean']:.4f}", f"{vals['std']:.4f}",
                        f"{vals['min']:.4f}", f"{vals['max']:.4f}"])

    if not rows:
        return None

    table = ax.table(
        cellText=rows,
        colLabels=['Métrica', 'Media', 'Std', 'Mín', 'Máx'],
        cellLoc='center',
        loc='center',
    )
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 1.5)

    # Estilizar header
    for j in range(5):
        table[0, j].set_facecolor('#2563EB')
        table[0, j].set_text_props(color='white', fontweight='bold')

    # Alternar colores de fila
    for i in range(1, len(rows) + 1):
        color = '#F0F4FF' if i % 2 == 0 else 'white'
        for j in range(5):
            table[i, j].set_facecolor(color)

    ax.set_title('Resumen Estadístico', fontsize=16, fontweight='bold',
                 pad=20, color='#1E3A5F')
    fig.tight_layout()
    return fig


# ==========================================
# CATÁLOGO DE GRÁFICAS
# ==========================================
PLOTS = [
    # (archivo, titulo, tipo, color, seccion)
    ('energy_em.xvg',              'Minimización de Energía',            'ts',  C['red'],    'Minimización'),
    ('temperature_nvt.xvg',        'Temperatura (NVT)',                  'ts',  C['orange'], 'Equilibración'),
    ('pressure_npt.xvg',           'Presión (NPT)',                      'ts',  C['purple'], 'Equilibración'),
    ('density_npt.xvg',            'Densidad (NPT)',                     'ts',  C['green'],  'Equilibración'),
    ('rmsd_backbone.xvg',          'RMSD Backbone',                      'ts',  C['blue'],   'Estabilidad Estructural'),
    ('rmsd_protein.xvg',           'RMSD Proteína Completa',             'ts',  C['green'],  'Estabilidad Estructural'),
    ('gyrate.xvg',                 'Radio de Giro',                      'ts',  C['purple'], 'Estabilidad Estructural'),
    ('rmsd_ligand.xvg',            'RMSD Ligando',                       'ts',  C['orange'], 'Estabilidad del Ligando'),
    ('mindist_prot_lig.xvg',       'Distancia Mínima Proteína-Ligando',  'ts',  C['red'],    'Estabilidad del Ligando'),
    ('native_contacts_dist.xvg',   'Contactos Nativos (d<0.45nm)',       'ts',  C['teal'],   'Estabilidad del Ligando'),
    ('interaction_energy.xvg',     'Energía de Interacción Prot-Lig',    'ts',  C['purple'], 'Estabilidad del Ligando'),
    ('rmsf_residue.xvg',           'Fluctuación por Residuo (RMSF)',     'bar', C['blue'],   'Flexibilidad'),
    ('sasa_protein.xvg',           'Superficie Accesible (SASA)',        'ts',  C['green'],  'Propiedades Superficiales'),
    ('hbond_protein.xvg',          'H-bonds Intra-proteína',             'ts',  C['blue'],   'Interacciones'),
    ('hbond_protein_ligand.xvg',   'H-bonds Proteína-Ligando',           'ts',  C['red'],    'Interacciones'),
    ('contacts_prot_lig.xvg',      'Contactos Proteína-Ligando',         'ts',  C['green'],  'Interacciones'),
    ('eigenval.xvg',               'Eigenvalores (PCA)',                 'ts',  C['purple'], 'PCA'),
    ('proj_pc1.xvg',               'Proyección PC1',                     'ts',  C['blue'],   'PCA'),
    ('proj_pc2.xvg',               'Proyección PC2',                     'ts',  C['red'],    'PCA'),
]


# ==========================================
# GENERADOR DE REPORTE
# ==========================================
def generate_report(rundir, output_dir=None, include_pca=True, include_dccm=True):
    rundir = Path(rundir)
    dirs = [rundir / '04_analysis', rundir / '02_equilibration', rundir / '01_minimization']
    out_dir = Path(output_dir) if output_dir else rundir
    output_pdf = out_dir / 'REPORT_MD.pdf'
    count = 0

    # Colectar estadísticas para tabla resumen
    stats = {}

    print(f"\n{'='*50}")
    print(f"  GENERANDO REPORTE DE ANÁLISIS MD v4.0")
    print(f"{'='*50}")
    print(f"  Corrida: {rundir.name}")
    print(f"  Salida:  {output_pdf}\n")

    with PdfPages(str(output_pdf)) as pdf:
        # Página de título
        fig = plt.figure(figsize=(10, 7))
        fig.text(0.5, 0.65, 'Reporte de Análisis\nDinámica Molecular',
                 ha='center', va='center', fontsize=28, fontweight='bold', color='#1E3A5F')
        fig.text(0.5, 0.45, rundir.name,
                 ha='center', va='center', fontsize=14, color='#555', style='italic')
        fig.text(0.5, 0.35, f'Generado: {datetime.now().strftime("%Y-%m-%d %H:%M")}',
                 ha='center', va='center', fontsize=11, color='#888')
        fig.text(0.5, 0.28, 'Pipeline v4.0',
                 ha='center', va='center', fontsize=10, color='#AAA')
        fig.patch.set_facecolor('white')
        pdf.savefig(fig)
        plt.close(fig)

        current_section = ""

        # Filtrar plots según opciones
        active_plots = PLOTS[:]
        if not include_pca:
            active_plots = [p for p in active_plots if p[4] != 'PCA']

        for filename, title_def, ptype, color, section in active_plots:
            # Buscar archivo en los directorios
            filepath = None
            for d in dirs:
                c = d / filename
                if c.exists():
                    filepath = c
                    break

            if filepath is None:
                print(f"  ⚠ No encontrado: {filename}")
                continue

            data, title_xvg, xlabel, ylabel, legends = parse_xvg(filepath)
            if data is None:
                print(f"  ⚠ Sin datos: {filename}")
                continue

            title = title_xvg or title_def

            # Colectar estadísticas
            if data.shape[1] >= 2:
                stats[title_def] = {
                    'mean': np.mean(data[:, 1]),
                    'std': np.std(data[:, 1]),
                    'min': np.min(data[:, 1]),
                    'max': np.max(data[:, 1]),
                }

            # Separador de sección
            if section != current_section:
                current_section = section
                fs = plt.figure(figsize=(10, 2))
                fs.text(0.5, 0.5, section, ha='center', va='center',
                        fontsize=22, fontweight='bold', color='#2563EB')
                fs.patch.set_facecolor('#F8FAFC')
                pdf.savefig(fs)
                plt.close(fs)

            fig, ax = plt.subplots(figsize=(10, 5))
            if ptype == 'bar':
                plot_bar(ax, data, title, xlabel, ylabel, color)
            else:
                plot_timeseries(ax, data, title, xlabel, ylabel, legends, color)
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)
            count += 1
            print(f"  ✓ {title_def}")

        # B-factors desde RMSF (v4.0)
        rmsf_file = dirs[0] / 'rmsf_residue.xvg'
        if rmsf_file.exists():
            rmsf_data, _, _, _, _ = parse_xvg(rmsf_file)
            bfactor_data = compute_bfactors(rmsf_data)
            if bfactor_data is not None:
                if current_section != 'Flexibilidad':
                    current_section = 'Flexibilidad'
                    fs = plt.figure(figsize=(10, 2))
                    fs.text(0.5, 0.5, 'Flexibilidad', ha='center', va='center',
                            fontsize=22, fontweight='bold', color='#2563EB')
                    fs.patch.set_facecolor('#F8FAFC')
                    pdf.savefig(fs)
                    plt.close(fs)

                fig, ax = plt.subplots(figsize=(10, 5))
                plot_bfactors(ax, bfactor_data)
                fig.tight_layout()
                pdf.savefig(fig)
                plt.close(fig)
                count += 1
                print(f"  ✓ B-factors (desde RMSF)")

        # RMSD Matrix
        xpm = dirs[0] / 'rmsd_matrix.xpm'  # 04_analysis
        if xpm.exists():
            fig, ax = plt.subplots(figsize=(8, 7))
            if plot_matrix_from_xpm(ax, xpm):
                fig.tight_layout()
                pdf.savefig(fig)
                plt.close(fig)
                count += 1
                print(f"  ✓ Matriz RMSD")
            else:
                plt.close(fig)
                print(f"  ⚠ No se pudo parsear rmsd_matrix.xpm")

        # Free Energy Landscape (FEL)
        if include_pca:
            fel_file = dirs[0] / 'proj_2d.xvg'
            if fel_file.exists():
                if current_section != 'Free Energy Landscape':
                    current_section = 'Free Energy Landscape'
                    fs = plt.figure(figsize=(10, 2))
                    fs.text(0.5, 0.5, 'Free Energy Landscape', ha='center', va='center',
                            fontsize=22, fontweight='bold', color='#2563EB')
                    fs.patch.set_facecolor('#F8FAFC')
                    pdf.savefig(fs)
                    plt.close(fs)

                fig_fel, ok = plot_fel(fel_file)
                if ok and fig_fel is not None:
                    pdf.savefig(fig_fel)
                    plt.close(fig_fel)
                    count += 1
                    print(f"  ✓ Free Energy Landscape (FEL)")
                else:
                    print(f"  ⚠ No se pudo generar FEL desde proj_2d.xvg")
            else:
                print(f"  ⚠ No encontrado: proj_2d.xvg (FEL no disponible)")

        # Dynamic Cross-Correlation Matrix (DCCM)
        if include_dccm:
            dccm_file = dirs[0] / 'ca_positions.xvg'
            if dccm_file.exists():
                if current_section != 'Correlación Dinámica':
                    current_section = 'Correlación Dinámica'
                    fs = plt.figure(figsize=(10, 2))
                    fs.text(0.5, 0.5, 'Correlación Dinámica', ha='center', va='center',
                            fontsize=22, fontweight='bold', color='#2563EB')
                    fs.patch.set_facecolor('#F8FAFC')
                    pdf.savefig(fs)
                    plt.close(fs)

                fig_dccm, ok = plot_dccm(dccm_file)
                if ok and fig_dccm is not None:
                    pdf.savefig(fig_dccm)
                    plt.close(fig_dccm)
                    count += 1
                    print(f"  ✓ DCCM (Correlación Dinámica)")
                else:
                    print(f"  ⚠ No se pudo generar DCCM desde ca_positions.xvg")
            else:
                print(f"  ⚠ No encontrado: ca_positions.xvg (DCCM no disponible)")

        # Convergencia (v4.0)
        rmsd_file = dirs[0] / 'rmsd_backbone.xvg'
        if rmsd_file.exists():
            rmsd_data, _, _, _, _ = parse_xvg(rmsd_file)
            conv = assess_convergence(rmsd_data)
            if conv:
                if current_section != 'Convergencia':
                    current_section = 'Convergencia'
                    fs = plt.figure(figsize=(10, 2))
                    fs.text(0.5, 0.5, 'Evaluación de Convergencia', ha='center', va='center',
                            fontsize=22, fontweight='bold', color='#2563EB')
                    fs.patch.set_facecolor('#F8FAFC')
                    pdf.savefig(fs)
                    plt.close(fs)

                fig, ax = plt.subplots(figsize=(10, 4))
                ax.axis('off')
                status = "✓ CONVERGIDO" if conv['converged'] else "⚠ NO CONVERGIDO"
                status_color = '#059669' if conv['converged'] else '#DC2626'
                conv_text = (
                    f"RMSD Backbone - Evaluación de Convergencia\n\n"
                    f"  RMSD total:        {conv['total_mean']:.4f} ± {conv['total_std']:.4f} nm\n"
                    f"  1ª mitad:          {conv['first_half_mean']:.4f} nm\n"
                    f"  2ª mitad:          {conv['second_half_mean']:.4f} nm\n"
                    f"  Drift (|Δ|):       {conv['drift']:.4f} nm\n\n"
                    f"  Estado: {status}"
                )
                ax.text(0.1, 0.5, conv_text, transform=ax.transAxes,
                        fontsize=12, verticalalignment='center', fontfamily='monospace',
                        bbox=dict(boxstyle='round', facecolor='#F8FAFC', edgecolor=status_color,
                                  linewidth=2))
                fig.tight_layout()
                pdf.savefig(fig)
                plt.close(fig)
                count += 1
                print(f"  ✓ Evaluación de convergencia: {status}")

        # Tabla resumen de estadísticas (v4.0)
        fig_table = generate_summary_table(stats)
        if fig_table is not None:
            pdf.savefig(fig_table)
            plt.close(fig_table)
            count += 1
            print(f"  ✓ Tabla resumen estadístico")

    print(f"\n{'='*50}")
    print(f"  ✓ Reporte: {output_pdf}")
    print(f"  ✓ {count} gráficas generadas")
    print(f"{'='*50}\n")


# ==========================================
# CLI
# ==========================================
def build_parser():
    parser = argparse.ArgumentParser(
        description='Generador de reportes de análisis MD v4.0',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Ejemplos:
  python3 plot_analysis.py MD_RUN/mi_corrida/
  python3 plot_analysis.py MD_RUN/mi_corrida/ --output-dir ./reportes/
  python3 plot_analysis.py MD_RUN/mi_corrida/ --no-pca --no-dccm
        """
    )
    parser.add_argument('rundir', nargs='?', default=None,
                       help='Directorio de la corrida MD')
    parser.add_argument('--output-dir', '-o', default=None,
                       help='Directorio de salida para el reporte (default: mismo que rundir)')
    parser.add_argument('--no-pca', action='store_true',
                       help='Omitir análisis PCA y FEL')
    parser.add_argument('--no-dccm', action='store_true',
                       help='Omitir DCCM')
    return parser


if __name__ == '__main__':
    parser = build_parser()

    # Compatibilidad hacia atrás: si no hay args, menú interactivo
    if len(sys.argv) >= 2 and not sys.argv[1].startswith('-'):
        args = parser.parse_args()
        generate_report(args.rundir, output_dir=args.output_dir,
                       include_pca=not args.no_pca, include_dccm=not args.no_dccm)
    elif len(sys.argv) >= 2:
        args = parser.parse_args()
        if args.rundir:
            generate_report(args.rundir, output_dir=args.output_dir,
                           include_pca=not args.no_pca, include_dccm=not args.no_dccm)
        else:
            parser.print_help()
    else:
        md_run = Path('MD_RUN')
        if md_run.exists():
            runs = sorted([d for d in md_run.iterdir()
                          if d.is_dir() and (d / '04_analysis').exists()])
            if runs:
                print("\nCorridas MD con análisis disponible:\n")
                for i, r in enumerate(runs, 1):
                    print(f"  {i}) {r.name}")
                print(f"\nSeleccione la corrida [1-{len(runs)}]: ", end='')
                try:
                    choice = int(input().strip())
                    if 1 <= choice <= len(runs):
                        generate_report(runs[choice - 1])
                    else:
                        print("Selección inválida")
                except (ValueError, KeyboardInterrupt):
                    print("\nCancelado")
            else:
                print("No se encontraron corridas con análisis completado")
                print(f"Uso: python3 {sys.argv[0]} <directorio_corrida_MD>")
        else:
            print(f"Uso: python3 {sys.argv[0]} <directorio_corrida_MD>")
