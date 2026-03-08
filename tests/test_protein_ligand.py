# -*- coding: utf-8 -*-
"""Tests for pipline_protein_ligand pipeline — XVG parsing, B-factors, convergence.

Run with:  python tests/test_protein_ligand.py
     or:   python -m pytest tests/test_protein_ligand.py -v
"""
from __future__ import annotations

import os
import sys
import tempfile
import types

import numpy as np

# Mock matplotlib before importing plot_analysis (it sys.exit(1) if missing)
if "matplotlib" not in sys.modules:
    _mpl = types.ModuleType("matplotlib")
    _mpl.use = lambda *a, **kw: None

    _plt = types.ModuleType("matplotlib.pyplot")
    _plt.rcParams = {}
    _plt.subplots = lambda *a, **kw: (None, None)
    _plt.figure = lambda *a, **kw: None
    _plt.close = lambda *a, **kw: None
    _plt.colorbar = lambda *a, **kw: None

    _backend = types.ModuleType("matplotlib.backends")
    _pdf = types.ModuleType("matplotlib.backends.backend_pdf")
    _pdf.PdfPages = type("PdfPages", (), {"__init__": lambda s, *a, **kw: None})  # type: ignore

    sys.modules["matplotlib"] = _mpl
    sys.modules["matplotlib.pyplot"] = _plt
    sys.modules["matplotlib.backends"] = _backend
    sys.modules["matplotlib.backends.backend_pdf"] = _pdf

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))


# ---------------------------------------------------------------------------
# plot_analysis tests
# ---------------------------------------------------------------------------

class TestParseXvgProtLig:

    def test_basic_parse(self):
        from plot_analysis import parse_xvg
        with tempfile.NamedTemporaryFile(mode="w", suffix=".xvg",
                                        delete=False) as f:
            f.write('# GROMACS output\n')
            f.write('@ title "RMSD Backbone"\n')
            f.write('@ xaxis label "Time (ps)"\n')
            f.write('@ yaxis label "RMSD (nm)"\n')
            f.write('0.0 0.15\n')
            f.write('10.0 0.18\n')
            name = f.name
        try:
            data, title, xlabel, ylabel, legends = parse_xvg(name)
            assert data is not None
            assert data.shape == (2, 2)
            assert title == 'RMSD Backbone'
            assert xlabel == 'Time (ps)'
        finally:
            os.unlink(name)

    def test_legends_parsed(self):
        from plot_analysis import parse_xvg
        with tempfile.NamedTemporaryFile(mode="w", suffix=".xvg",
                                        delete=False) as f:
            f.write('@ s0 legend "Protein"\n')
            f.write('@ s1 legend "Ligand"\n')
            f.write('0.0 1.0 2.0\n')
            name = f.name
        try:
            data, _, _, _, legends = parse_xvg(name)
            assert len(legends) == 2
            assert 'Protein' in legends
            assert 'Ligand' in legends
        finally:
            os.unlink(name)

    def test_empty_returns_none(self):
        from plot_analysis import parse_xvg
        with tempfile.NamedTemporaryFile(mode="w", suffix=".xvg",
                                        delete=False) as f:
            f.write('# only comments\n')
            name = f.name
        try:
            data, _, _, _, _ = parse_xvg(name)
            assert data is None
        finally:
            os.unlink(name)


class TestComputeBfactorsProtLig:

    def test_conversion(self):
        from plot_analysis import compute_bfactors
        rmsf = np.array([[1.0, 0.1], [2.0, 0.2]])
        result = compute_bfactors(rmsf)
        assert result is not None
        assert result.shape == (2, 2)
        # B = (8π²/3) × (RMSF_Å)²; 0.1 nm = 1 Å
        expected = (8.0 * np.pi**2 / 3.0) * 1.0**2
        assert abs(result[0, 1] - expected) < 0.01

    def test_none_input(self):
        from plot_analysis import compute_bfactors
        assert compute_bfactors(None) is None


class TestAssessConvergenceProtLig:

    def test_converged(self):
        from plot_analysis import assess_convergence
        np.random.seed(42)
        vals = 0.2 + np.random.normal(0, 0.01, 100)
        data = np.column_stack([np.arange(100), vals])
        result = assess_convergence(data)
        assert result is not None
        assert result['converged']

    def test_none(self):
        from plot_analysis import assess_convergence
        assert assess_convergence(None) is None


# ---------------------------------------------------------------------------
# Minimal runner
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    import traceback
    passed = 0
    failed = 0
    for cls in [TestParseXvgProtLig, TestComputeBfactorsProtLig,
                TestAssessConvergenceProtLig]:
        obj = cls()
        for name in sorted(dir(obj)):
            if name.startswith("test_"):
                try:
                    getattr(obj, name)()
                    print(f"  ✓ {cls.__name__}.{name}")
                    passed += 1
                except Exception as e:
                    print(f"  ✗ {cls.__name__}.{name}: {e}")
                    traceback.print_exc()
                    failed += 1
    print(f"\n{passed} passed, {failed} failed")
