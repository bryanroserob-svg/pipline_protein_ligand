#!/usr/bin/env python3
"""
Generador automático de gráficas de análisis MD.
Lee archivos .xvg de una corrida MD y genera un PDF con todas las figuras.

Uso:
    python3 plot_analysis.py                          # Menú interactivo
    python3 plot_analysis.py MD_RUN/mi_corrida/       # Directo
"""

import sys
import os
import re
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
def generate_report(rundir):
    rundir = Path(rundir)
    dirs = [rundir / '04_analysis', rundir / '02_equilibration', rundir / '01_minimization']
    output_pdf = rundir / 'REPORT_MD.pdf'
    count = 0

    print(f"\n{'='*50}")
    print(f"  GENERANDO REPORTE DE ANÁLISIS MD")
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
        fig.patch.set_facecolor('white')
        pdf.savefig(fig)
        plt.close(fig)

        current_section = ""

        for filename, title_def, ptype, color, section in PLOTS:
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
        fel_file = dirs[0] / 'proj_2d.xvg'
        if fel_file.exists():
            # Separador de sección FEL
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

    print(f"\n{'='*50}")
    print(f"  ✓ Reporte: {output_pdf}")
    print(f"  ✓ {count} gráficas generadas")
    print(f"{'='*50}\n")


# ==========================================
# CLI
# ==========================================
if __name__ == '__main__':
    if len(sys.argv) >= 2:
        generate_report(sys.argv[1])
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
