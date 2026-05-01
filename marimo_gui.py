"""Marimo front-end for the Stark shift calculator.

Run with::

    marimo edit marimo_gui.py
or::
    marimo run marimo_gui.py

All physics lives in :mod:`calcs`.
"""

import marimo

__generated_with = "0.10.0"
app = marimo.App(width="medium")


@app.cell
def _imports():
    import io

    import matplotlib.pyplot as plt
    import marimo as mo

    import calcs
    return calcs, io, mo, plt


@app.cell
def _title(mo):
    mo.md("# Alkali AC Stark Shift Calculator")
    return


@app.cell
def _atom_input(calcs, mo):
    atom_input = mo.ui.dropdown(
        options=calcs.list_atom_species(),
        value="Rubidium87",
        label="Atom species",
    )
    return (atom_input,)


@app.cell
def _state_inputs(mo):
    n_input = mo.ui.number(value=5, start=1, step=1, label="n")
    l_input = mo.ui.number(value=0, start=0, step=1, label="l")
    j_input = mo.ui.number(value=0.5, step=0.5, label="J")
    f_input = mo.ui.number(value=2, step=0.5, label="F")
    mf_input = mo.ui.text(value="", label="m_F (blank = all)")
    return f_input, j_input, l_input, mf_input, n_input


@app.cell
def _light_inputs(mo):
    wavelength_input = mo.ui.number(
        value=1064.0, step=1.0, label="Wavelength (nm)"
    )
    sm_input = mo.ui.number(value=0.0, step=0.1, label="sigma- power")
    pi_input = mo.ui.number(value=0.0, step=0.1, label="pi power")
    sp_input = mo.ui.number(value=1.0, step=0.1, label="sigma+ power")
    power_input = mo.ui.number(value=1.0, step=0.1, label="Power")
    power_unit = mo.ui.dropdown(
        options=["uW", "mW", "W"], value="mW", label="Power unit"
    )
    waist_input = mo.ui.number(value=1.0, step=0.1, label="Waist")
    waist_unit = mo.ui.dropdown(
        options=["um", "mm"], value="um", label="Waist unit"
    )
    distance_input = mo.ui.number(
        value=0.0, step=0.1, label="Distance from waist"
    )
    distance_unit = mo.ui.dropdown(
        options=["nm", "um", "mm", "m"], value="um", label="Distance unit"
    )
    return (
        distance_input,
        distance_unit,
        pi_input,
        power_input,
        power_unit,
        sm_input,
        sp_input,
        waist_input,
        waist_unit,
        wavelength_input,
    )


@app.cell
def _layout(
    atom_input,
    distance_input,
    distance_unit,
    f_input,
    j_input,
    l_input,
    mf_input,
    mo,
    n_input,
    pi_input,
    power_input,
    power_unit,
    sm_input,
    sp_input,
    waist_input,
    waist_unit,
    wavelength_input,
):
    inputs_panel = mo.vstack([
        mo.md("### Atom and state"),
        atom_input,
        mo.hstack([n_input, l_input, j_input, f_input, mf_input], justify="start"),
        mo.md("### Light"),
        wavelength_input,
        mo.md("**Polarization (relative powers; auto-normalized):**"),
        mo.hstack([sm_input, pi_input, sp_input], justify="start"),
        mo.hstack([power_input, power_unit], justify="start"),
        mo.hstack([waist_input, waist_unit], justify="start"),
        mo.hstack([distance_input, distance_unit], justify="start"),
    ])
    inputs_panel
    return


@app.cell
def _compute_button(mo):
    compute_btn = mo.ui.run_button(label="Compute Stark shifts", kind="success")
    compute_btn
    return (compute_btn,)


@app.cell
def _calculate(
    atom_input,
    calcs,
    compute_btn,
    distance_input,
    distance_unit,
    f_input,
    j_input,
    l_input,
    mf_input,
    n_input,
    pi_input,
    power_input,
    power_unit,
    sm_input,
    sp_input,
    waist_input,
    waist_unit,
    wavelength_input,
):
    _power_factors = {"uW": 1e-6, "mW": 1e-3, "W": 1.0}
    _waist_factors = {"um": 1e-6, "mm": 1e-3}
    _distance_factors = {"nm": 1e-9, "um": 1e-6, "mm": 1e-3, "m": 1.0}

    result = None
    error_msg = None
    if compute_btn.value:
        try:
            _mf_text = (mf_input.value or "").strip()
            _mf_value = float(_mf_text) if _mf_text else None
            result = calcs.compute_stark_shifts(
                atom_name=atom_input.value,
                n=int(n_input.value),
                l=int(l_input.value),
                j=float(j_input.value),
                f=float(f_input.value),
                mf=_mf_value,
                wavelength=float(wavelength_input.value) * 1e-9,
                polarization_powers=(
                    float(sm_input.value),
                    float(pi_input.value),
                    float(sp_input.value),
                ),
                power=float(power_input.value)
                * _power_factors[power_unit.value],
                waist=float(waist_input.value)
                * _waist_factors[waist_unit.value],
                distance=float(distance_input.value)
                * _distance_factors[distance_unit.value],
            )
        except calcs.QuantumNumberError as e:
            error_msg = f"Invalid quantum numbers: {e}"
        except Exception as e:
            error_msg = f"{type(e).__name__}: {e}"
    return error_msg, result


@app.cell
def _summary(calcs, error_msg, mo, result):
    if error_msg:
        summary_view = mo.callout(mo.md(f"**Error:** {error_msg}"), kind="danger")
    elif result is None:
        summary_view = mo.md("*Set inputs above and click* **Compute Stark shifts**.")
    else:
        _pol = result["polarization"]
        _header = (
            f"**{result['atom_name']}** (I = {result['I']}) — "
            f"n={result['n']}, l={result['l']}, "
            f"J={result['j']}, F={result['f']}"
            + (f", m_F={result['mf']}" if result['mf'] is not None else "")
            + f"   |   λ = {result['wavelength']*1e9:.2f} nm"
        )
        _params = (
            f"Power = {result['power']*1e3:.4g} mW, "
            f"waist = {result['waist']*1e6:.4g} um, "
            f"distance = {result['distance']*1e6:.4g} um  →  "
            f"I = {result['intensity']:.3e} W/m²"
        )
        _pol_line = (
            f"σ⁻ = {_pol['c_minus_sq']:.3f}, "
            f"π = {_pol['c_pi_sq']:.3f}, "
            f"σ⁺ = {_pol['c_plus_sq']:.3f} "
            f"(A = {_pol['A']:+.3f}, B = {_pol['B']:+.3f})"
        )
        _alpha = (
            f"α_S = {result['alpha_S']:+.3e}, "
            f"α_V = {result['alpha_V']:+.3e}, "
            f"α_T = {result['alpha_T']:+.3e}  [Hz·m²/V²]"
        )
        if result["mf"] is not None:
            _shifts = result["shifts"][result["mf"]]
            _mag = max(abs(v) for v in _shifts.values()) or 1.0
            _scale, _unit = calcs.choose_freq_unit(_mag)
            _total = sum(_shifts.values())
            _table = (
                f"| | shift ({_unit}) |\n"
                f"|---|---:|\n"
                f"| Scalar | {_shifts['scalar']/_scale:+.4f} |\n"
                f"| Vector | {_shifts['vector']/_scale:+.4f} |\n"
                f"| Tensor | {_shifts['tensor']/_scale:+.4f} |\n"
                f"| **Total** | **{_total/_scale:+.4f}** |\n"
            )
            summary_view = mo.md(
                f"{_header}\n\n{_params}\n\n{_pol_line}\n\n{_alpha}\n\n{_table}"
            )
        else:
            _all = [v for d in result["shifts"].values() for v in d.values()]
            _mag = max((abs(v) for v in _all), default=1.0) or 1.0
            _scale, _unit = calcs.choose_freq_unit(_mag)
            _rows = [f"| m_F | scalar ({_unit}) | vector ({_unit}) | tensor ({_unit}) | total ({_unit}) |"]
            _rows.append("|---:|---:|---:|---:|---:|")
            for _m in sorted(result["shifts"].keys()):
                _s = result["shifts"][_m]
                _t = _s["scalar"] + _s["vector"] + _s["tensor"]
                _rows.append(
                    f"| {_m:+.1f} "
                    f"| {_s['scalar']/_scale:+.4f} "
                    f"| {_s['vector']/_scale:+.4f} "
                    f"| {_s['tensor']/_scale:+.4f} "
                    f"| {_t/_scale:+.4f} |"
                )
            summary_view = mo.md(
                f"{_header}\n\n{_params}\n\n{_pol_line}\n\n{_alpha}\n\n"
                + "\n".join(_rows)
            )
    summary_view
    return


@app.cell
def _plot(calcs, io, mo, plt, result):
    if result is None or result.get("mf") is not None:
        plot_view = mo.md("")
    else:
        fig, ax = plt.subplots(figsize=(9, 5.5), dpi=180)
        calcs.plot_level_diagram(ax, result)
        fig.tight_layout()

        _svg_buf = io.StringIO()
        fig.savefig(_svg_buf, format="svg", bbox_inches="tight")
        _svg = _svg_buf.getvalue()
        _plot_html = mo.Html(
            "<style>"
            ".stark-shift-plot svg {"
            "width: 100%; max-width: 980px; height: auto; display: block;"
            "}"
            "</style>"
            "<div class='stark-shift-plot'>"
            f"{_svg}"
            "</div>"
        )

        _png_buf = io.BytesIO()
        fig.savefig(_png_buf, format="png", dpi=300, bbox_inches="tight")
        _png_buf.seek(0)
        plt.close(fig)

        download_btn = mo.download(
            data=_png_buf.getvalue(),
            filename="stark_shifts.png",
            label="Download plot as PNG",
            mimetype="image/png",
        )
        plot_view = mo.vstack([_plot_html, download_btn])
    plot_view
    return


if __name__ == "__main__":
    app.run()
