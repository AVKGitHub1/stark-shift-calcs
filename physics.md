# Physics of the AC Stark shift calculator

This note documents the formulas and conventions used by [calcs.py](calcs.py).
Quantities and units are SI throughout; energy shifts are reported as
frequencies (`Hz`).

## 1. Setup

An alkali atom in state $|n, l, J, F, m_F\rangle$ sits in a monochromatic
optical field of angular frequency $\omega = 2\pi c / \lambda$, peak electric
field amplitude $E_0$, and polarization decomposed in the spherical basis with
the quantization axis $\hat{z}$ along the magnetic field:

$$
\vec{\varepsilon} \;=\; c_{-1}\,\hat{e}_{-1} + c_{0}\,\hat{e}_{0} + c_{+1}\,\hat{e}_{+1},
\qquad |c_{-1}|^2 + |c_{0}|^2 + |c_{+1}|^2 = 1,
$$

where $\hat{e}_{0} = \hat{z}$, $\hat{e}_{\pm 1} = \mp(\hat{x} \pm i\hat{y})/\sqrt{2}$.
$|c_{+1}|^2$, $|c_0|^2$, $|c_{-1}|^2$ are the fractional powers in
$\sigma^{+}$, $\pi$, $\sigma^{-}$, respectively.

The light is treated as an incoherent mixture of the three components, which
is appropriate for "what fraction of the power is in each polarization"
inputs. The AC Stark shift is linear in intensity, so the sum is exact in
that case.

## 2. Field and intensity

For a Gaussian beam of total power $P$ and waist $w_0$, the on-axis intensity
at axial distance $z$ from the focus is

$$
I(z) = \frac{2P}{\pi\, w(z)^2},
\qquad w(z) = w_0\sqrt{1 + (z/z_R)^2},
\qquad z_R = \pi w_0^2 / \lambda.
$$

The peak field amplitude is then

$$
|E_0|^2 \;=\; \frac{2 I}{c\,\varepsilon_0}.
$$

## 3. Stark shift in the $J$ basis

The standard decomposition of the AC Stark energy shift in terms of
scalar, vector, and tensor dynamic polarizabilities (Manakov,
Ovsiannikov, Rapoport 1986; Le Kien et al. 2013) is

$$
\Delta E(J, m_J)
= -\tfrac{1}{4}\,|E_0|^2 \!\left[
\alpha_S(\omega)
+ A\,\alpha_V(\omega)\,\frac{m_J}{J}
+ B\,\alpha_T(\omega)\,\frac{3 m_J^2 - J(J+1)}{J(2J - 1)}
\right].
$$

The polarization-dependent factors are

$$
A \;=\; |c_{+1}|^2 - |c_{-1}|^2 \quad \text{(degree of circular polarization)},
$$

$$
B \;=\; \frac{3|c_0|^2 - 1}{2} \quad \text{(linear polarization tensor projection)}.
$$

Limits:

- $\sigma^{+}$ pure: $A=+1$, $B=-\tfrac{1}{2}$.
- $\sigma^{-}$ pure: $A=-1$, $B=-\tfrac{1}{2}$.
- $\pi$ pure: $A=0$, $B=+1$.
- The tensor term vanishes for $J\le 1/2$; the vector term vanishes for $J=0$.

The dynamic polarizabilities $\alpha_S$, $\alpha_V$, $\alpha_T$ come directly
from `arc.DynamicPolarizability(atom, n, l, J).getPolarizability(\lambda, units='SI')`.
ARC's `'SI'` mode returns them in `Hz·m²/V²` (i.e., the polarizability divided
by Planck's constant), so the formula above gives a frequency directly when
combined with $|E_0|^2$ in `V²/m²`. No further factor of $h$ is needed.

The basis of intermediate states summed in ARC is set by
`defineBasis(nMin, nMax)` with `nMin = max(l+1, n-5, n_ground)` and
`nMax = n + 25` by default; widen this for highly excited Rydberg states or to
test convergence.

## 4. Hyperfine $F$ states by $J$-basis projection

Rather than constructing $\alpha_S^{F}$, $\alpha_V^{F}$, $\alpha_T^{F}$
separately (which requires the appropriate Wigner-6j Lande factors), this code
projects $|F, m_F\rangle$ onto the $|J, m_J; I, m_I\rangle$ basis and sums the
$J$-basis Stark shifts weighted by the squared Clebsch-Gordan coefficients:

$$
\Delta E(F, m_F)
\;=\; \sum_{m_J,\, m_I = m_F - m_J}
\bigl|\langle J\, m_J;\, I\, m_I\,\big|\, F\, m_F\rangle\bigr|^{2}\;
\Delta E(J, m_J).
$$

This is **exact** in the limit where the laser detuning is large compared
with the hyperfine splittings of the connecting $J'$ manifolds — i.e., the
polarizabilities for different $F'$ levels of a given $J'$ are approximately
equal. This is the standard approximation for far-off-resonant dipole traps.
Near resonance with a hyperfine-resolved line, it does **not** hold and
hyperfine-resolved $\alpha^{F}$ values must be used instead.

The same decomposition is applied separately to the scalar, vector, and
tensor pieces, so the three contributions are reported individually.

Convenient identities for the projection:

$$
\sum_{m_J} \bigl|\langle J m_J;\, I m_I| F m_F\rangle\bigr|^{2} = 1
\quad \Rightarrow \quad \Delta E_S(F, m_F) = -\tfrac{1}{4}\,|E_0|^2\,\alpha_S,
$$

$$
\langle m_J\rangle_{F,m_F} = \sum_{m_J} \bigl|\cdots\bigr|^{2}\, m_J
\;\;\Longrightarrow\;\;
\Delta E_V(F, m_F) = -\tfrac{1}{4}\,|E_0|^2\,A\,\alpha_V\,\frac{\langle m_J\rangle_{F,m_F}}{J},
$$

and analogously for the tensor piece.

For stretched states $|F = J + I, m_F = \pm F\rangle$ the projection collapses
onto a single $|J, m_J\!=\!\pm J\rangle$ component, so the $J$-basis formula is
recovered directly.

## 5. Validation of inputs

[validate_state](calcs.py) raises `QuantumNumberError` for any of:

- $n \le 0$ or $l \ge n$,
- $J \notin \{|l - 1/2|,\, l + 1/2\}$,
- $F < |J - I|$ or $F > J + I$, or $F$ not differing from $|J - I|$ by an integer,
- $|m_F| > F$, or $m_F$ not differing from $F$ by an integer.

If $m_F$ is left blank in the GUI, all $m_F \in \{-F, -F+1, \dots, +F\}$ are
computed and rendered as a level diagram.

## 6. Output

`compute_stark_shifts` returns a dictionary containing the inputs, the
intermediate intensities/fields, the $J$-basis polarizabilities, and a
`shifts` mapping `{m_F: {"scalar": Hz, "vector": Hz, "tensor": Hz}}`. The
total shift is the sum of the three components.

When $m_F$ is `None`, [plot_level_diagram](calcs.py) renders one horizontal
"level" line per $m_F$ at zero detuning, with three coloured arrows per state
showing the magnitude and sign of the scalar (blue), vector (green), and
tensor (red) contributions; the tip of each arrow is annotated with its
numeric value in an auto-selected unit (Hz / kHz / MHz / GHz).

## 7. References and conventions

The conventions here follow:

- N. L. Manakov, V. D. Ovsiannikov, L. P. Rapoport,
  *Phys. Rep.* **141**, 320 (1986).
- F. Le Kien, P. Schneeweiss, A. Rauschenbeutel,
  "Dynamical polarizability of atoms in arbitrary light fields",
  *Eur. Phys. J. D* **67**, 92 (2013).
- J. Mitroy, M. S. Safronova, C. W. Clark,
  "Theory and applications of atomic and ionic polarizabilities",
  *J. Phys. B* **43**, 202001 (2010).
- N. Sibalic et al., "ARC: An open-source library for calculating properties
  of alkali Rydberg atoms", *Comp. Phys. Comm.* **220**, 319 (2017),
  [https://arc-alkali-rydberg-calculator.readthedocs.io](https://arc-alkali-rydberg-calculator.readthedocs.io).

Different references use slightly different normalizations of $\alpha_V$ and
$\alpha_T$ (e.g., some define the vector term with $m_J/(2J)$ instead of
$m_J/J$, absorbing the factor of 2 into $\alpha_V$). The convention used
here matches Le Kien et al. 2013 Eq. (41) and ARC's internal definitions.
