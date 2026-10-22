# Physical and numerical notes

The solver is a modified compPS/compPSc implementation, not the separable Thomson-recursion-plus-Jones approximation. It evolves coupled axisymmetric Stokes I and Q with the thermal Compton scattering redistribution matrix. The original compPS thermal-electron physics and energy-dependent cross-section remain part of the native calculation.

## Polarized source iteration

For each scattering order the source combines four Mueller blocks, schematically

\[
\binom{S_I}{S_Q}(E,\mu,t)=\int dE' d\mu'\,
\begin{pmatrix}R_{II}&R_{IQ}\\R_{QI}&R_{QQ}\end{pmatrix}
\binom{I}{Q}(E',\mu',t).
\]

Both propagation directions and polarization feedback are included. Q is not obtained by multiplying an unpolarized spectrum by a fixed polarization fraction. The native sky-independent reference convention is reversed at the Python boundary to give Q positive parallel to the slab.

## Stable transport

For piecewise-linear source S along a segment, formal transport integrates the exponential attenuation analytically. For constant S and opacity κ, a segment of depth d contributes S(1−exp(−κd))/κ. Positive short-characteristic weights replace a trapezoidal transport approximation that could amplify grazing-angle modes. Unscattered Lambert photons are handled as an incoming boundary intensity, not a delta source in the first depth cell.

## Incident-energy integration

The selected snapshot integrates the redistribution against positive linear hat functions in y = ln E'. If A_j is the legacy energy quadrature weight, it uses

\[
\overline F_{ij}=A_j^{-1}\int F(E_i,e^y)H_j(y)dy,
\]

with eight Gauss points per adjacent interval. All four Mueller entries are integrated; the E/E' factor belongs inside this integration. Ordered energy pairs are evaluated directly, with the native angular reciprocity handling retained. The native cross-section normalization is retained; no output amplitude is fitted to MCRT or to the recursion model.

This is incident-energy product integration with output-energy collocation, not a fully conservative double-bin Galerkin discretization. Small photon-closure error therefore cannot by itself prove full spectral accuracy. Negative-basis cubic experiments are not selected as the production default.

## Validation requirements

1. Verify the scattering-order source increment reaches the stopping tolerance before the cap.
2. Refine depth, internal frequency, and internal angular grids separately.
3. Refine input seed sampling and output angle-window integration.
4. Compare I and signed Q, not only Q/I. Check integrated photon counts and weak tails independently.
5. Use an independent benchmark where possible. Source iteration convergence does not establish physical correctness or mesh convergence, especially in cold, optically thick cases.

References: Poutanen & Svensson (1996), ApJ 470, 249; Nagirner & Poutanen (1993), A&A 275, 325. Retain the upstream convolution implementation attribution to Michal Szanecki.
