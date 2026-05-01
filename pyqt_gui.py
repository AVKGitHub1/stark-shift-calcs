"""PyQt6 front-end for the Stark shift calculator.

All physics lives in :mod:`calcs`; this module is purely GUI plumbing.
"""

from __future__ import annotations

import sys

import matplotlib

matplotlib.use("QtAgg")

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from PyQt6 import QtCore
from PyQt6.QtWidgets import (
    QApplication,
    QComboBox,
    QFileDialog,
    QGridLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMainWindow,
    QMessageBox,
    QPushButton,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

import calcs


_POWER_UNITS = {"uW": 1e-6, "mW": 1e-3, "W": 1.0}
_WAIST_UNITS = {"um": 1e-6, "mm": 1e-3}
_DIST_UNITS = {"nm": 1e-9, "um": 1e-6, "mm": 1e-3, "m": 1.0}


class StarkShiftWindow(QMainWindow):
    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle("Alkali AC Stark Shift Calculator")
        self.resize(1150, 720)
        self._last_result: dict | None = None
        self._build_ui()

    def _build_ui(self) -> None:
        central = QWidget()
        self.setCentralWidget(central)
        outer = QHBoxLayout(central)
        outer.addWidget(self._build_input_panel(), 0)
        outer.addWidget(self._build_output_panel(), 1)

    def _build_input_panel(self) -> QWidget:
        panel = QWidget()
        layout = QVBoxLayout(panel)
        layout.setAlignment(QtCore.Qt.AlignmentFlag.AlignTop)

        atom_group = QGroupBox("Atom")
        atom_layout = QGridLayout(atom_group)
        self.atom_combo = QComboBox()
        for name in calcs.list_atom_species():
            self.atom_combo.addItem(name)
        idx = self.atom_combo.findText("Rubidium87")
        if idx >= 0:
            self.atom_combo.setCurrentIndex(idx)
        atom_layout.addWidget(QLabel("Species:"), 0, 0)
        atom_layout.addWidget(self.atom_combo, 0, 1)
        layout.addWidget(atom_group)

        state_group = QGroupBox("Atomic state")
        sl = QGridLayout(state_group)
        self.n_edit = QLineEdit("5")
        self.l_edit = QLineEdit("0")
        self.j_edit = QLineEdit("0.5")
        self.f_edit = QLineEdit("2")
        self.mf_edit = QLineEdit("")
        self.mf_edit.setPlaceholderText("blank = all m_F")
        for row, (label, w) in enumerate([
            ("n", self.n_edit),
            ("l", self.l_edit),
            ("J", self.j_edit),
            ("F", self.f_edit),
            ("m_F", self.mf_edit),
        ]):
            sl.addWidget(QLabel(label + ":"), row, 0)
            sl.addWidget(w, row, 1)
        layout.addWidget(state_group)

        light_group = QGroupBox("Light")
        ll = QGridLayout(light_group)
        self.wavelength_edit = QLineEdit("1064")
        ll.addWidget(QLabel("Wavelength:"), 0, 0)
        ll.addWidget(self.wavelength_edit, 0, 1)
        ll.addWidget(QLabel("nm"), 0, 2)

        self.sm_edit = QLineEdit("0")
        self.pi_edit = QLineEdit("0")
        self.sp_edit = QLineEdit("1")
        ll.addWidget(QLabel("sigma- power:"), 1, 0)
        ll.addWidget(self.sm_edit, 1, 1)
        ll.addWidget(QLabel("pi power:"), 2, 0)
        ll.addWidget(self.pi_edit, 2, 1)
        ll.addWidget(QLabel("sigma+ power:"), 3, 0)
        ll.addWidget(self.sp_edit, 3, 1)
        ll.addWidget(QLabel("(any units; auto-normalized)"), 4, 0, 1, 3)

        self.power_edit = QLineEdit("1")
        self.power_unit = QComboBox()
        for u in _POWER_UNITS:
            self.power_unit.addItem(u)
        self.power_unit.setCurrentText("mW")
        ll.addWidget(QLabel("Power:"), 5, 0)
        ll.addWidget(self.power_edit, 5, 1)
        ll.addWidget(self.power_unit, 5, 2)

        self.waist_edit = QLineEdit("1")
        self.waist_unit = QComboBox()
        for u in _WAIST_UNITS:
            self.waist_unit.addItem(u)
        ll.addWidget(QLabel("Waist:"), 6, 0)
        ll.addWidget(self.waist_edit, 6, 1)
        ll.addWidget(self.waist_unit, 6, 2)

        self.dist_edit = QLineEdit("0")
        self.dist_unit = QComboBox()
        for u in _DIST_UNITS:
            self.dist_unit.addItem(u)
        ll.addWidget(QLabel("Dist. from waist:"), 7, 0)
        ll.addWidget(self.dist_edit, 7, 1)
        ll.addWidget(self.dist_unit, 7, 2)
        layout.addWidget(light_group)

        self.compute_btn = QPushButton("Compute")
        self.compute_btn.clicked.connect(self._on_compute)
        layout.addWidget(self.compute_btn)

        self.save_btn = QPushButton("Save plot...")
        self.save_btn.clicked.connect(self._on_save_plot)
        self.save_btn.setEnabled(False)
        layout.addWidget(self.save_btn)

        return panel

    def _build_output_panel(self) -> QWidget:
        panel = QWidget()
        layout = QVBoxLayout(panel)

        self.text_output = QTextEdit()
        self.text_output.setReadOnly(True)
        self.text_output.setMaximumHeight(220)
        font = self.text_output.font()
        font.setFamily("Consolas")
        self.text_output.setFont(font)
        layout.addWidget(self.text_output)

        self.figure = Figure(figsize=(7, 5))
        self.canvas = FigureCanvas(self.figure)
        layout.addWidget(self.canvas, 1)

        return panel

    def _gather(self) -> dict:
        try:
            atom_name = self.atom_combo.currentText()
            n = int(self.n_edit.text())
            l = int(self.l_edit.text())
            j = float(self.j_edit.text())
            f = float(self.f_edit.text())
            mf_text = self.mf_edit.text().strip()
            mf = float(mf_text) if mf_text else None
            wavelength_nm = float(self.wavelength_edit.text())
            sm = float(self.sm_edit.text())
            pi_pol = float(self.pi_edit.text())
            sp = float(self.sp_edit.text())
            power = (float(self.power_edit.text())
                     * _POWER_UNITS[self.power_unit.currentText()])
            waist = (float(self.waist_edit.text())
                     * _WAIST_UNITS[self.waist_unit.currentText()])
            distance = (float(self.dist_edit.text())
                        * _DIST_UNITS[self.dist_unit.currentText()])
        except ValueError as exc:
            raise ValueError(f"Could not parse input: {exc}") from exc
        return dict(
            atom_name=atom_name, n=n, l=l, j=j, f=f, mf=mf,
            wavelength=wavelength_nm * 1e-9,
            polarization_powers=(sm, pi_pol, sp),
            power=power, waist=waist, distance=distance,
        )

    def _on_compute(self) -> None:
        try:
            params = self._gather()
        except ValueError as exc:
            QMessageBox.warning(self, "Input error", str(exc))
            return

        self.compute_btn.setEnabled(False)
        self.compute_btn.setText("Computing...")
        QApplication.processEvents()
        try:
            try:
                result = calcs.compute_stark_shifts(**params)
            except calcs.QuantumNumberError as exc:
                QMessageBox.warning(self, "Invalid quantum numbers", str(exc))
                return
            except Exception as exc:
                QMessageBox.critical(self, "Calculation error", repr(exc))
                return
            self._last_result = result
            self._show_result(result)
        finally:
            self.compute_btn.setEnabled(True)
            self.compute_btn.setText("Compute")

    def _show_result(self, result: dict) -> None:
        self.text_output.setPlainText("\n".join(self._format_summary(result)))

        self.figure.clear()
        ax = self.figure.add_subplot(111)
        if result["mf"] is None:
            calcs.plot_level_diagram(ax, result)
            self.figure.tight_layout()
            self.save_btn.setEnabled(True)
        else:
            ax.text(
                0.5, 0.5,
                "Leave m_F blank to draw\nthe full level diagram.",
                ha="center", va="center", transform=ax.transAxes,
                fontsize=12, color="gray",
            )
            ax.set_axis_off()
            self.save_btn.setEnabled(False)
        self.canvas.draw()

    def _format_summary(self, result: dict) -> list[str]:
        pol = result["polarization"]
        lines = [
            f"Atom        : {result['atom_name']}    (I = {result['I']})",
            f"State       : n={result['n']}, l={result['l']}, "
            f"J={result['j']}, F={result['f']}"
            + (f", m_F={result['mf']}" if result["mf"] is not None
               else "    (all m_F)"),
            f"Wavelength  : {result['wavelength'] * 1e9:.3f} nm",
            f"Power       : {result['power'] * 1e3:.4g} mW    "
            f"Waist: {result['waist'] * 1e6:.4g} um    "
            f"Distance: {result['distance'] * 1e6:.4g} um",
            f"Intensity   : {result['intensity']:.3e} W/m^2"
            f"    |E_0|^2 = {result['E_squared']:.3e} V^2/m^2",
            f"Polarization: sigma-={pol['c_minus_sq']:.4f}, "
            f"pi={pol['c_pi_sq']:.4f}, sigma+={pol['c_plus_sq']:.4f}    "
            f"(A={pol['A']:+.4f}, B={pol['B']:+.4f})",
            f"alpha_S={result['alpha_S']:+.3e}, "
            f"alpha_V={result['alpha_V']:+.3e}, "
            f"alpha_T={result['alpha_T']:+.3e}    [Hz m^2/V^2]",
            "",
        ]
        if result["mf"] is not None:
            shifts = result["shifts"][result["mf"]]
            mag = max(abs(s) for s in shifts.values()) or 1.0
            scale, unit = calcs.choose_freq_unit(mag)
            total = sum(shifts.values())
            lines.append(f"Stark shifts for m_F = {result['mf']}:")
            lines.append(f"  Scalar : {shifts['scalar'] / scale:+12.4f} {unit}")
            lines.append(f"  Vector : {shifts['vector'] / scale:+12.4f} {unit}")
            lines.append(f"  Tensor : {shifts['tensor'] / scale:+12.4f} {unit}")
            lines.append(f"  Total  : {total / scale:+12.4f} {unit}")
        else:
            all_vals = [v for d in result["shifts"].values() for v in d.values()]
            mag = max((abs(v) for v in all_vals), default=1.0) or 1.0
            scale, unit = calcs.choose_freq_unit(mag)
            lines.append(f"Stark shifts (per m_F, in {unit}):")
            lines.append(
                f"  {'m_F':>6}  {'scalar':>12}  {'vector':>12}  "
                f"{'tensor':>12}  {'total':>12}"
            )
            for mf in sorted(result["shifts"].keys()):
                s = result["shifts"][mf]
                tot = s["scalar"] + s["vector"] + s["tensor"]
                lines.append(
                    f"  {mf:>+6.1f}  {s['scalar'] / scale:>+12.4f}  "
                    f"{s['vector'] / scale:>+12.4f}  "
                    f"{s['tensor'] / scale:>+12.4f}  "
                    f"{tot / scale:>+12.4f}"
                )
        return lines

    def _on_save_plot(self) -> None:
        path, _ = QFileDialog.getSaveFileName(
            self, "Save level diagram", "stark_shifts.png",
            "PNG image (*.png);;PDF document (*.pdf);;SVG image (*.svg)",
        )
        if not path:
            return
        try:
            self.figure.savefig(path, dpi=200, bbox_inches="tight")
        except Exception as exc:
            QMessageBox.critical(self, "Save failed", str(exc))
            return
        QMessageBox.information(self, "Saved", f"Plot saved to:\n{path}")


def main() -> None:
    app = QApplication(sys.argv)
    win = StarkShiftWindow()
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
