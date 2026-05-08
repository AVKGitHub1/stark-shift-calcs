"""PyQt6 front-end for dynamic scalar/vector/tensor polarizabilities.

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


class PolarizabilityWindow(QMainWindow):
    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle("Alkali Dynamic Polarizability Calculator")
        self.resize(980, 640)
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
        state_layout = QGridLayout(state_group)
        self.n_edit = QLineEdit("5")
        self.l_edit = QLineEdit("0")
        self.j_edit = QLineEdit("0.5")
        self.f_edit = QLineEdit("2")
        self.mf_edit = QLineEdit("")
        self.mf_edit.setPlaceholderText("blank = all m_F")
        for row, (label, widget) in enumerate([
            ("n", self.n_edit),
            ("l", self.l_edit),
            ("J", self.j_edit),
            ("F", self.f_edit),
            ("m_F", self.mf_edit),
        ]):
            state_layout.addWidget(QLabel(label + ":"), row, 0)
            state_layout.addWidget(widget, row, 1)
        layout.addWidget(state_group)

        light_group = QGroupBox("Light")
        light_layout = QGridLayout(light_group)
        self.wavelength_edit = QLineEdit("1064")
        light_layout.addWidget(QLabel("Wavelength:"), 0, 0)
        light_layout.addWidget(self.wavelength_edit, 0, 1)
        light_layout.addWidget(QLabel("nm"), 0, 2)

        self.basis_window_edit = QLineEdit("25")
        light_layout.addWidget(QLabel("Basis window:"), 1, 0)
        light_layout.addWidget(self.basis_window_edit, 1, 1)
        layout.addWidget(light_group)

        self.compute_btn = QPushButton("Compute polarizabilities")
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

        self.figure = Figure(figsize=(7, 4.5))
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
            basis_window = int(self.basis_window_edit.text())
        except ValueError as exc:
            raise ValueError(f"Could not parse input: {exc}") from exc
        if basis_window <= 0:
            raise ValueError("Basis window must be a positive integer.")
        return dict(
            atom_name=atom_name,
            n=n,
            l=l,
            j=j,
            f=f,
            mf=mf,
            wavelength=wavelength_nm * 1e-9,
            basis_window=basis_window,
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
                result = calcs.compute_polarizabilities(**params)
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
            self.compute_btn.setText("Compute polarizabilities")

    def _show_result(self, result: dict) -> None:
        self.text_output.setPlainText("\n".join(self._format_summary(result)))

        self.figure.clear()
        ax = self.figure.add_subplot(111)
        if result["projected_polarizabilities"] is None:
            self._plot_component_bars(
                ax,
                [result["alpha_S"], result["alpha_V"], result["alpha_T"]],
                "J-basis polarizability (Hz m^2/V^2)",
            )
        elif result["mf"] is None:
            self._plot_all_mf_components(ax, result)
        else:
            projected = result["projected_polarizabilities"][result["mf"]]
            self._plot_component_bars(
                ax,
                [projected["scalar"], projected["vector"], projected["tensor"]],
                "Projected polarizability term (Hz m^2/V^2)",
            )
        ax.set_title(self._plot_title(result), fontsize=10)
        self.figure.tight_layout()
        self.save_btn.setEnabled(True)
        self.canvas.draw()

    def _plot_component_bars(self, ax, values: list[float],
                             ylabel: str) -> None:
        labels = ["Scalar", "Vector", "Tensor"]
        colors = ["tab:blue", "tab:green", "tab:red"]
        bars = ax.bar(labels, values, color=colors)
        ax.axhline(0, color="black", lw=0.8)
        ax.set_ylabel(ylabel)
        for bar, value in zip(bars, values):
            y = bar.get_height()
            va = "bottom" if y >= 0 else "top"
            pad = 3 if y >= 0 else -3
            ax.annotate(
                f"{value:+.3e}",
                xy=(bar.get_x() + bar.get_width() / 2, y),
                xytext=(0, pad),
                textcoords="offset points",
                ha="center",
                va=va,
                fontsize=9,
            )

    def _plot_all_mf_components(self, ax, result: dict) -> None:
        projected = result["projected_polarizabilities"]
        mf_values = sorted(projected.keys())
        x_positions = list(range(len(mf_values)))
        width = 0.25
        series = [
            ("Scalar", "scalar", "tab:blue", -width),
            ("Vector", "vector", "tab:green", 0.0),
            ("Tensor", "tensor", "tab:red", width),
        ]
        for label, key, color, offset in series:
            values = [projected[m][key] for m in mf_values]
            xs = [x + offset for x in x_positions]
            ax.bar(xs, values, width=width, label=label, color=color)
        ax.axhline(0, color="black", lw=0.8)
        ax.set_xticks(x_positions)
        ax.set_xticklabels([f"{m:g}" for m in mf_values])
        ax.set_xlabel("m_F")
        ax.set_ylabel("Projected polarizability term (Hz m^2/V^2)")
        ax.legend(fontsize=8)

    def _plot_title(self, result: dict) -> str:
        title = (
            f"{result['atom_name']}  n={result['n']}, l={result['l']}, "
            f"J={result['j']}"
        )
        if result["f"] is not None:
            title += f", F={result['f']}"
            title += (
                f", m_F={result['mf']}" if result["mf"] is not None
                else ", all m_F"
            )
        return f"{title}  at {result['wavelength'] * 1e9:.2f} nm"

    def _format_summary(self, result: dict) -> list[str]:
        return [
            *self._state_summary(result),
            "",
            "J-basis dynamic polarizabilities [Hz m^2/V^2]:",
            f"  Scalar : {result['alpha_S']:+.8e}",
            f"  Vector : {result['alpha_V']:+.8e}",
            f"  Tensor : {result['alpha_T']:+.8e}",
            *self._projected_summary(result),
        ]

    def _state_summary(self, result: dict) -> list[str]:
        state = (
            f"State       : n={result['n']}, l={result['l']}, J={result['j']}"
        )
        if result["f"] is not None:
            state += f", F={result['f']}"
            state += (
                f", m_F={result['mf']}" if result["mf"] is not None
                else "    (all m_F)"
            )
        return [
            f"Atom        : {result['atom_name']}    (I = {result['I']})",
            state,
            f"Wavelength  : {result['wavelength'] * 1e9:.3f} nm",
        ]

    def _projected_summary(self, result: dict) -> list[str]:
        projected = result["projected_polarizabilities"]
        if projected is None:
            return []
        lines = [
            "",
            "Projected hyperfine components [Hz m^2/V^2]:",
            "  These multiply 1, A, and B in the Stark-shift expression.",
        ]
        if result["mf"] is not None:
            p = projected[result["mf"]]
            factors = p["factors"]
            lines.extend([
                f"  Scalar : {p['scalar']:+.8e}",
                f"  Vector : {p['vector']:+.8e}",
                f"  Tensor : {p['tensor']:+.8e}",
                "",
                "Projection factors:",
                f"  Scalar : {factors['scalar']:+.6f}",
                f"  Vector : {factors['vector']:+.6f}",
                f"  Tensor : {factors['tensor']:+.6f}",
            ])
            return lines

        lines.append(
            f"  {'m_F':>6}  {'scalar':>14}  {'vector':>14}  {'tensor':>14}  "
            f"{'v fac':>9}  {'t fac':>9}"
        )
        for mf in sorted(projected.keys()):
            p = projected[mf]
            factors = p["factors"]
            lines.append(
                f"  {mf:>+6.1f}  {p['scalar']:>+14.6e}  "
                f"{p['vector']:>+14.6e}  {p['tensor']:>+14.6e}"
                f"  {factors['vector']:>+9.4f}  {factors['tensor']:>+9.4f}"
            )
        return lines

    def _on_save_plot(self) -> None:
        path, _ = QFileDialog.getSaveFileName(
            self, "Save polarizability plot", "polarizabilities.png",
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
    win = PolarizabilityWindow()
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
