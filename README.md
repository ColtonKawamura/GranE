# GranE

**Gran**ular **E**igenvalue 

## Purpose

A new project focusing on the eigenvalue problem of granular packings. 
This project is a continuation of the work in [**Gran**ular **M**echanics **A**coustics](https://github.com/ColtonKawamura/GranMA)


## Roadmap

- [x] 2D packing generation
- [x] 3D packing generation   
- [ ] Hertzian contacts 


## Validating eigenvectors with `simEigenmode.m`

`simEigenmode.m` is the bridge between the analytical eigenvalue problem and the
real packing. The eigenvalue solver predicts a set of **modes** -- eigenpairs of a
springs-and-masses model of the packing (a frequency from the imaginary part of an
eigenvalue, and an attenuation/decay rate from the real part), plus a **mode shape**:
the displacement every particle would carry in that mode. `simEigenmode` does not
trust that prediction on its own; it *recovers* the mode by letting the physical
packing actually move and checks that the observed motion agrees with what was
predicted.

The idea, in words rather than code:

1. **Load and pick a mode.** A saved packing and its predicted eigenpairs are read in
    (or, in the "polynomial/output" path, the eigenpairs are reused directly from a
   saved results structure). A single mode is then chosen, either by a target
   oscillation frequency (the imaginary part of the eigenvalue, $\omega$) or by a
   target attenuation (the real part, $\beta$).

2. **Turn the eigenvector into a starting configuration.** The selected eigenvector
   is reshaped into a per-particle displacement -- this *is* the predicted mode shape.
    Each particle is placed at its equilibrium position offset by a small multiple of
    that shape, and the whole system is held at rest. In other words, the packing is
    *seeded* as the mode and then released.

3. **Let the real packing respond.** A force-based simulation of the actual
   spring/damper contacts (Hertzian spring force plus contact damping), integrated
   with a symplectic Verlet scheme under two-way periodic boundaries, evolves the
   seeded configuration for many timesteps. The simulation runs on the GPU when one
   is available. This is where the check lives: the analytical mode shape is fed in,
   and the *dynamical* response of the same contacts is observed.

4. **Track a probe particle.** One representative particle is watched over time; its
    displacement in $x$ and $y$ is recorded throughout the run, giving a clean signal
    of how the mode oscillates and dies out.

5. **Check the frequency.** The probe's displacement time series is transformed to the
    frequency domain and the dominant oscillation frequency read off. This measured
    frequency is compared against the frequency the eigenvalue solver *predicted* for
    that mode.

6. **Check the attenuation.** The envelope of the displacement magnitude $|\mathbf{u}(t)|$
    is inspected on a log scale; a straight-line fit of $\ln|\mathbf{u}|$ versus time
    yields a measured decay rate, compared against the attenuation predicted by the real
    part of the eigenvalue.

7. **Quantify and report.** The comparison is summarized by an **eigen-frequency
    ratio** (predicted physical frequency over the frequency measured from the
    simulation) -- a value near one means the predicted mode shape really does describe
    the packing's motion. A panel figure plots the probe displacement versus time, the
    decay envelope with its fit, and the frequency-domain spectrum with a reference
    line at the predicted frequency, and the figure is saved for inspection.

Put simply: the solver guesses a mode; `simEigenmode` builds that guess into the
packing, plays it back through a real contact simulation, and confirms that the
packing oscillates at -- and decays at -- the rate the solver claimed. Because the
frequency and decay are re-measured independently from the simulated motion, the
routine is a self-consistency / validation test of the eigenvectors themselves, not a
fresh eigen calculation.
