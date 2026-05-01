"""GUI-agnostic AC Stark shift calculations for alkali atoms.

Polarizabilities come from ARC's :class:`arc.DynamicPolarizability` in J basis
(units of ``Hz·m^2/V^2``). The energy shift convention used throughout is

    Delta nu (J, m_J) = -(1/4) |E_0|^2
        * [ alpha_S
            + A * alpha_V * m_J / J
            + B * alpha_T * (3 m_J^2 - J(J+1)) / (J (2J - 1)) ]

with peak field amplitude |E_0|^2 = 2 I / (c epsilon_0), polarization
coefficients A = |c_+|^2 - |c_-|^2 (degree of circular polarization) and
B = (3|c_pi|^2 - 1)/2, where (|c_-|^2, |c_pi|^2, |c_+|^2) are the normalized
power fractions in the spherical basis with the quantization axis chosen
along the magnetic field.

Hyperfine F states are obtained by projecting onto J,m_J:
    Delta nu (F, m_F) = sum_{m_J, m_I = m_F - m_J}
                            |<J m_J; I m_I | F m_F>|^2 * Delta nu (J, m_J).
This is exact in the limit where the laser detuning is large compared with
hyperfine splittings of the connecting upper states (i.e., the polarizabilities
of all F' upper levels of a given J' are approximately equal).
"""

from __future__ import annotations

from typing import Optional

import numpy as np
from scipy.constants import c as C_LIGHT, epsilon_0 as EPS_0
from sympy.physics.quantum.cg import CG

import arc


_ALKALI_NAMES = [
    "Lithium6", "Lithium7",
    "Sodium",
    "Potassium39", "Potassium40", "Potassium41",
    "Rubidium85", "Rubidium87",
    "Caesium",
]


class QuantumNumberError(ValueError):
    """Raised when the requested quantum-number combination is unphysical."""


def list_atom_species() -> list[str]:
    """Names of the alkali atom classes available in the installed ARC."""
    return [name for name in _ALKALI_NAMES if hasattr(arc, name)]


def get_atom(name: str):
    if not hasattr(arc, name):
        raise QuantumNumberError(f"Unknown atom species: {name!r}.")
    return getattr(arc, name)()


def _is_half_integer(x: float) -> bool:
    return abs(2 * x - round(2 * x)) < 1e-9


def _half_integer_range(top: float) -> list[float]:
    n_steps = int(round(2 * top))
    return [-top + k for k in range(n_steps + 1)]


def validate_state(atom, n: int, l: int, j: float, f: float,
                   mf: Optional[float]) -> None:
    """Raise :class:`QuantumNumberError` if the state is not physically valid."""
    if n <= 0:
        raise QuantumNumberError(f"n must be a positive integer (got {n}).")
    if l < 0 or l >= n:
        raise QuantumNumberError(f"l must satisfy 0 <= l < n (got l={l}, n={n}).")
    if not _is_half_integer(j):
        raise QuantumNumberError(f"J must be (half-)integer (got {j}).")
    allowed_j = sorted({abs(l - 0.5), l + 0.5})
    if not any(abs(j - a) < 1e-9 for a in allowed_j):
        raise QuantumNumberError(
            f"J must be l +/- 1/2 = {allowed_j} (got J={j})."
        )

    I = atom.I
    f_min, f_max = abs(j - I), j + I
    if not _is_half_integer(f):
        raise QuantumNumberError(f"F must be (half-)integer (got {f}).")
    if f < f_min - 1e-9 or f > f_max + 1e-9:
        raise QuantumNumberError(
            f"F={f} is outside [|J - I|, J + I] = [{f_min}, {f_max}] "
            f"for J={j}, I={I}."
        )
    if abs((f - f_min) - round(f - f_min)) > 1e-9:
        raise QuantumNumberError(
            f"F={f} is not reachable from |J - I|={f_min} by integer steps."
        )

    if mf is not None:
        if not _is_half_integer(mf):
            raise QuantumNumberError(f"mF must be (half-)integer (got {mf}).")
        if abs(mf) > f + 1e-9:
            raise QuantumNumberError(
                f"|mF| must be <= F (got |mF|={abs(mf)}, F={f})."
            )
        if abs((f - mf) - round(f - mf)) > 1e-9:
            raise QuantumNumberError(
                f"mF={mf} must differ from F={f} by an integer."
            )


def normalize_polarization(p_minus: float, p_pi: float, p_plus: float) -> dict:
    """Normalize the (sigma-, pi, sigma+) power fractions and compute A, B.

    Returns a dict with the normalized intensity fractions ``c_*_sq`` and the
    polarization coefficients ``A`` (vector / circular) and ``B`` (tensor).
    """
    if min(p_minus, p_pi, p_plus) < 0:
        raise ValueError("Polarization powers must be non-negative.")
    total = p_minus + p_pi + p_plus
    if total <= 0:
        raise ValueError("Sum of polarization powers must be positive.")
    c_minus_sq = p_minus / total
    c_pi_sq = p_pi / total
    c_plus_sq = p_plus / total
    return {
        "c_minus_sq": c_minus_sq,
        "c_pi_sq": c_pi_sq,
        "c_plus_sq": c_plus_sq,
        "A": c_plus_sq - c_minus_sq,
        "B": (3 * c_pi_sq - 1) / 2,
    }


def gaussian_intensity_on_axis(power: float, waist: float, distance: float,
                                wavelength: float) -> float:
    """Peak (on-axis) intensity at axial distance ``distance`` from the focus."""
    if waist <= 0:
        raise ValueError("Waist must be positive.")
    if wavelength <= 0:
        raise ValueError("Wavelength must be positive.")
    z_R = np.pi * waist ** 2 / wavelength
    w_z = waist * np.sqrt(1.0 + (distance / z_R) ** 2)
    return 2.0 * power / (np.pi * w_z ** 2)


def get_polarizabilities(atom, n: int, l: int, j: float,
                          wavelength: float,
                          basis_window: int = 25) -> tuple[float, float, float]:
    """Return (alpha_S, alpha_V, alpha_T) in Hz·m^2/V^2 for ``|n, l, J>``."""
    dp = arc.DynamicPolarizability(atom, n, l, j)
    n_min = max(l + 1, n - 5, getattr(atom, "groundStateN", 1))
    n_max = n + basis_window
    dp.defineBasis(n_min, n_max)
    result = dp.getPolarizability(wavelength, units="SI")
    return float(result[0]), float(result[1] or 0.0), float(result[2] or 0.0)


def stark_shift_for_state(alpha_S: float, alpha_V: float, alpha_T: float,
                           E_squared: float, A_pol: float, B_pol: float,
                           j: float, I: float, f: float, mf: float) -> dict:
    """Stark shift contributions (in Hz) for the single state ``|F, mF>``."""
    scalar = vector = tensor = 0.0
    for mj in _half_integer_range(j):
        mi = mf - mj
        if abs(mi) > I + 1e-9:
            continue
        try:
            cg = float(CG(j, mj, I, mi, f, mf).doit())
        except Exception:
            continue
        weight = cg * cg
        if weight < 1e-15:
            continue
        scalar_j = -0.25 * alpha_S * E_squared
        vector_j = (-0.25 * A_pol * alpha_V * E_squared * mj / j
                    if j > 1e-9 else 0.0)
        if j > 0.5 + 1e-9:
            tensor_j = (-0.25 * B_pol * alpha_T * E_squared
                        * (3 * mj * mj - j * (j + 1))
                        / (j * (2 * j - 1)))
        else:
            tensor_j = 0.0
        scalar += weight * scalar_j
        vector += weight * vector_j
        tensor += weight * tensor_j
    return {"scalar": scalar, "vector": vector, "tensor": tensor}


def compute_stark_shifts(*, atom_name: str, n: int, l: int, j: float, f: float,
                          mf: Optional[float], wavelength: float,
                          polarization_powers: tuple[float, float, float],
                          power: float, waist: float, distance: float) -> dict:
    """Top-level entry point. All lengths in m, power in W. Shifts in Hz.

    ``polarization_powers`` is ``(P_sigma_minus, P_pi, P_sigma_plus)``;
    they are normalized internally. If ``mf`` is None, all mF in [-F, F]
    are computed.
    """
    atom = get_atom(atom_name)
    validate_state(atom, n, l, j, f, mf)
    pol = normalize_polarization(*polarization_powers)
    intensity = gaussian_intensity_on_axis(power, waist, distance, wavelength)
    E_squared = 2.0 * intensity / (C_LIGHT * EPS_0)

    alpha_S, alpha_V, alpha_T = get_polarizabilities(atom, n, l, j, wavelength)

    mf_list = _half_integer_range(f) if mf is None else [mf]
    shifts = {
        m: stark_shift_for_state(
            alpha_S, alpha_V, alpha_T,
            E_squared, pol["A"], pol["B"],
            j, atom.I, f, m,
        )
        for m in mf_list
    }

    return {
        "atom_name": atom_name,
        "n": n, "l": l, "j": j, "f": f, "mf": mf,
        "I": atom.I,
        "wavelength": wavelength,
        "polarization": pol,
        "polarization_powers": polarization_powers,
        "power": power, "waist": waist, "distance": distance,
        "intensity": intensity, "E_squared": E_squared,
        "alpha_S": alpha_S, "alpha_V": alpha_V, "alpha_T": alpha_T,
        "shifts": shifts,
    }


_FREQ_UNITS = (("Hz", 1.0), ("kHz", 1e3), ("MHz", 1e6), ("GHz", 1e9))


def choose_freq_unit(magnitude: float) -> tuple[float, str]:
    """Pick a (scale, unit) pair so values display with 1-3 leading digits."""
    if magnitude <= 0 or not np.isfinite(magnitude):
        return 1.0, "Hz"
    scale, unit = 1.0, "Hz"
    for u, s in _FREQ_UNITS:
        if magnitude >= s:
            scale, unit = s, u
    return scale, unit


def _format_half(x: float) -> str:
    twox = 2 * x
    if abs(twox - round(twox)) > 1e-6:
        return f"{x:g}"
    twox_int = int(round(twox))
    if twox_int % 2 == 0:
        return f"{twox_int // 2}"
    sign = "-" if twox_int < 0 else ""
    return f"{sign}{abs(twox_int)}/2"


def plot_level_diagram(ax, result: dict,
                       freq_unit: Optional[str] = None) -> None:
    """Draw a level diagram with arrows for scalar/vector/tensor shifts."""
    from matplotlib.patches import Patch

    f = result["f"]
    mf_values = _half_integer_range(f)
    all_shifts = [v for m in mf_values for v in result["shifts"][m].values()]
    max_abs = max((abs(s) for s in all_shifts), default=1.0) or 1.0

    if freq_unit is None:
        scale, unit = choose_freq_unit(max_abs)
    else:
        scale = dict(_FREQ_UNITS)[freq_unit]
        unit = freq_unit

    colors = {"scalar": "tab:blue", "vector": "tab:green", "tensor": "tab:red"}
    offset = 0.18
    # vertical span in plot data units used to compute label offsets
    span = max_abs / scale if scale else max_abs or 1.0
    # Vertical offset keeps text clear of arrowheads.
    text_vpad = 0.02 * span
    value_labels = []
    arrow_segments = []
    level_segments = []

    for m in mf_values:
        ax.hlines(0, m - 0.4, m + 0.4, color="black", lw=1.2, zorder=1)
        level_segments.append((m - 0.4, 0, m + 0.4, 0))
        shifts = result["shifts"][m]
        for k, key in enumerate(("scalar", "vector", "tensor")):
            value = shifts[key] / scale
            if abs(value) <= 1e-12:
                continue
            x = m + (k - 1) * offset
            arrow_segments.append((x, 0, x, value))
            ax.annotate(
                "", xy=(x, value), xytext=(x, 0),
                arrowprops=dict(arrowstyle="->", color=colors[key], lw=1.6),
                zorder=2,
            )
            # Place the label directly above/below the arrow tip.
            label_y = value + (text_vpad if value >= 0 else -text_vpad)
            label = ax.text(
                x, label_y, f"{value:+.3g}",
                fontsize=7, color=colors[key], ha="center",
                va="bottom" if value >= 0 else "top",
                zorder=3,
            )
            value_labels.append(label)

    ax.axhline(0, color="gray", lw=0.5, zorder=0)
    ax.set_xticks(mf_values)
    ax.set_xticklabels([_format_half(m) for m in mf_values])
    ax.set_xlabel(f"$m_F$ (F = {_format_half(f)})")
    ax.set_ylabel(f"AC Stark shift ({unit})")
    title = (
        f"{result['atom_name']}    "
        f"n={result['n']}, l={result['l']}, "
        f"J={_format_half(result['j'])}, F={_format_half(result['f'])}    "
        f"$\\lambda$ = {result['wavelength'] * 1e9:.2f} nm"
    )
    ax.set_title(title, fontsize=10)

    span = max_abs / scale
    ax.set_xlim(min(mf_values) - 0.6, max(mf_values) + 0.6)
    if span > 0:
        ax.set_ylim(-1.5 * span, 1.5 * span)

    handles = [Patch(color=colors[k], label=f"{k.capitalize()} shift")
               for k in colors]
    _place_legend_smartly(
        ax, handles, value_labels, arrow_segments, level_segments,
        fontsize=8, framealpha=0.9,
    )


_INSIDE_LEGEND_LOCS = (
    "upper right", "upper left", "lower right", "lower left",
    "center right", "center left", "upper center", "lower center", "center",
)


def _place_legend_smartly(ax, handles, value_labels, arrow_segments,
                          level_segments, **legend_kwargs) -> None:
    """Place the legend where it least overlaps rendered plot content."""
    renderer = _get_renderer(ax)
    protected_bboxes = _protected_plot_bboxes(
        ax, renderer, value_labels, arrow_segments, level_segments,
    )

    for loc in _INSIDE_LEGEND_LOCS:
        legend = ax.legend(handles=handles, loc=loc, **legend_kwargs)
        renderer = _get_renderer(ax)
        bbox = legend.get_window_extent(renderer)
        score = sum(_bbox_overlap_area(bbox, other) for other in protected_bboxes)
        legend.remove()
        if score <= 1.0:
            ax.legend(handles=handles, loc=loc, **legend_kwargs)
            return

    ax.legend(
        handles=handles, loc="upper left", bbox_to_anchor=(1.02, 1.0),
        borderaxespad=0.0, **legend_kwargs,
    )


def _get_renderer(ax):
    figure = ax.figure
    if figure.canvas is None:
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        FigureCanvasAgg(figure)
    figure.canvas.draw()
    return figure.canvas.get_renderer()


def _protected_plot_bboxes(ax, renderer, value_labels, arrow_segments,
                           level_segments) -> list:
    from matplotlib.transforms import Bbox

    bboxes = [
        label.get_window_extent(renderer).expanded(1.15, 1.35)
        for label in value_labels
        if label.get_visible()
    ]
    for x0, y0, x1, y1 in arrow_segments:
        bboxes.append(_data_segment_bbox(ax, x0, y0, x1, y1, pad=8))
    for x0, y0, x1, y1 in level_segments:
        bboxes.append(_data_segment_bbox(ax, x0, y0, x1, y1, pad=5))
    return [bbox for bbox in bboxes if isinstance(bbox, Bbox)]


def _data_segment_bbox(ax, x0, y0, x1, y1, pad: float):
    from matplotlib.transforms import Bbox

    p0, p1 = ax.transData.transform([(x0, y0), (x1, y1)])
    bbox = Bbox.from_extents(
        min(p0[0], p1[0]), min(p0[1], p1[1]),
        max(p0[0], p1[0]), max(p0[1], p1[1]),
    )
    return bbox.padded(pad)


def _bbox_overlap_area(a, b) -> float:
    x0 = max(a.x0, b.x0)
    y0 = max(a.y0, b.y0)
    x1 = min(a.x1, b.x1)
    y1 = min(a.y1, b.y1)
    if x1 <= x0 or y1 <= y0:
        return 0.0
    return (x1 - x0) * (y1 - y0)
