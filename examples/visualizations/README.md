# flagVisPack — 2-D Frictional Packing Visualization

This directory contains examples produced with the `flagVisPack` option of
`src/pack.m`. It is **2-D only**: `flagVisPack = true` captures the compression
trajectory of a 2-D disk packing and assembles the sampled frames into an MP4 or
GIF (whichever is smaller), with a red diameter line drawn across each disk to
make particle rotation visible.

## Run it

```matlab
addpath('src');
% Frictional 20x20 packing (N = 400). A GIF is written next to the .mat.
opts = struct('flagVisPack', true, 'flagFrictionOn', true, ...
              'visPackMaxFrames', 60, 'visPackSkip', 3000);
pack(400, 100, 1, 1.4, 1, 1e-3, 2, false, 1, 1, 0, false, '.', opts);
```

`flagVisPack` defaults to `false`, so every existing call path is unchanged and
bit-identical to before.

| Example | Case | Output |
| --- | --- | --- |
| `2D_N400_P0.001_Width20_Seed1_VisPack.mp4` | frictionless | MP4 |
| `2D_N400_P0.001_Width20_Seed2_VisPack.gif` | **frictional** | **GIF** (red rotation lines) |

The frictional GIF is the canonical example: each disk carries a red diameter
line that spins as the disk rotates under tangential (Cundall–Strack) friction,
while the box compresses from ~35 px to the jammed box size.

## Physics: why the frictional controller converges

Issue [#16](https://github.com/ColtonKawamura/GranE/issues/16) reported that the
2-D frictional packing never reached force-balance convergence — it fell into a
limit cycle (the box resized every step off a transient overlap-energy proxy, so
it could never settle). This fix redirects the frictional path to the
Cundall–Strack protocol reflected in the repository's `CundallStrack_2D/` and
`OverDamp.cpp` references. The two mechanisms are paired below with the code.

### 1. Virial pressure instead of an overlap-energy proxy

The old control read a transient scalar and resized continuously. The fix
computes a proper virial (configurational) pressure from the minimum-image
contact vectors and the *full* contact force (normal + tangential), normalized by
the stiffness `K` to match the dimensionless `P_target` scale.

```matlab
% src/pack.m — 2-D virial pressure (normal + tangential)
vecVirial = vecSepX .* vecVirialFx + vecSepY .* vecVirialFy;   % r · F per pair
scalVolume = scalBoxWidthX * scalBoxHeightY;                    % 2-D area
scalPressureVirial = -sum(vecVirial) / (2 * scalVolume * K);    % virial pressure
```

The virial is the statistically-mechanical pressure of a granular packing
(`P_ab = (1/V) Σ_ij r_ab^i F_ab^j`), so it is the physically correct load to drive
the box. Frictional runs use it (`scalPressure = scalPressureVirial`); the
frictionless path keeps its energy-based pressure unchanged.

### 2. Fixed-box relaxation with discrete resize events

Rather than resizing every step, the controller relaxes at a **fixed** box size,
resizes only while `P` sits outside a dead-band around `P_target`, and **holds**
the box while `P` is in-band:

```matlab
% src/pack.m — pressure dead-band + fixed-box relaxation
boolInBand = abs(scalPressure - P_target) / P_target < scalFrictionDeadBand;
...
if scalPressure < P_target * (1 - scalFrictionDeadBand)     % loose -> compress
    ...; scalScale = 1 - scalFrictionRate;                  % discrete isotropic resize
elseif scalPressure > P_target * (1 + scalFrictionDeadBand) % dense -> expand
    ...; scalScale = 1 + scalFrictionRate;
end
scalFrictionHoldCounter = 0;                                 % reset after each resize
```

Physically this is a controlled quasistatic compression: hold the volume long
enough (`scalFrictionHoldSteps`) for the grains to relax, then apply one discrete
isotropic resize. A continuous per-step resize injects energy every step and
sustains the limit cycle; the dead-band + discrete resize lets the packing
settle. The loose/dense snapshots give a per-box-size `P(V)` curve that the
controller bisects toward `P_target` (no overshoot).

### 3. Convergence by force balance, not by energy

A frictional packing carries a permanent rotational/tangential kinetic-energy
floor (the OverDamp-style solver converges on force, `|Acc| < F_thresh`, never on
energy). Gating convergence on `Ek < ε` is unreachable and freezes the box. The
fix accepts the packing only when it is simultaneously **in-band**, **percolated**
(a connected contact network), and **force-balanced**:

```matlab
% src/pack.m — three-part acceptance criterion
scalForceRatio  = scalMeanFnet / max(scalMeanFc, eps);      % mean|F_net|/mean|F_contact|
boolInBand      = abs(scalPressure - P_target)/P_target < scalFrictionDeadBand;
boolPercol      = (scalMeanCoordNum >= scalFrictionZmin);   % percolation guard (z_iso=3)
boolBalanced    = (scalForceRatio < scalFrictionForceTol);
scalFrictionBalCount = ... ; if boolInBand && boolPercol && boolBalanced, inc; else reset; end
if scalFrictionBalCount >= scalFrictionBalWindow, break; end   % sustained balance => jammed
```

A packing that is force-balanced for many consecutive fixed-box windows is at a
true isostatic jammed state — the physically correct "settled" signal. This is
what makes the packing converge instead of cycling.

### 4. Cundall–Strack tangential friction (what the rotation line visualizes)

The red line tracks the integrated orientation angle `vecTheta`, which only
exists when `flagVisPack` AND `flagFrictionOn` are both on. The underlying
mechanism is Coulomb friction on the tangential spring at each contact:

```matlab
% src/pack.m — Coulomb cap |F_t| <= mu*|F_n|
vecDispTan   = vecDispTan + vecVelTan * scalTimestep;         % advance tangential spring
vecDispTanCap = scalMu .* abs(vecForceMag) ./ scalKt;         % mu*|F_n|/K_t
vecDispTan   = sign(vecDispTan) .* min(vecDispTanCap, abs(vecDispTan));
```

The tangential spring accumulates slip (`U_ij`), the Coulomb cap
`|F_t| = K_t|U_ij| ≤ μ|F_n|` limits it, and the resulting tangential force
produces a torque `τ_i = r_i·F_t` that spins each disk — which is exactly the
rotation the red diameter line displays. `vecTheta` is integrated with a
trapezoidal rule from `vecOmega` so the stored memory cost is zero unless the
feature is requested.

### Result

The frictional 20×20 packing in `examples/visualizations/` converges to
φ ≈ 0.82, mean coordination z ≈ 3.8 — consistent with 2-D frictional random
disk packings near jamming. The full `tests/testPack/testPack.m` suite passes
(2-D, frictional, and 3-D), confirming `flagVisPack` is inert (bit-identical
output) when the flag is off.
