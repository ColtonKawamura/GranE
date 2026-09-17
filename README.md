# GranE

**Gran**ular **E**igenvalue

## Purpose

A new project focusing on the eigenvalue problem of granular packings.
This project is a continuation of the work in [**Gran**ular **M**echanics **A**coustics](https://github.com/ColtonKawamura/GranMA)

## Documentation

The [project wiki](https://github.com/ColtonKawamura/GranE/wiki) is the main reference: one page per function in `src/` with default syntax, every option, and runnable examples.

## Getting Started

1. Clone the repo and open the project root in MATLAB:

   ```matlab
   addpath('src')
   ```

2. **Generate a packing** with [`pack`](https://github.com/ColtonKawamura/GranE/wiki/pack) (2D or 3D, Hooke or Hertz, optional friction):

   ```matlab
   pack
   ```

3. **Compute damped eigenmodes** with [`processEigenModesDampedPara`](https://github.com/ColtonKawamura/GranE/wiki/processEigenModesDampedPara) (parallel sweep over packing × damping), or build the K, Γ, M matrices yourself with [`matSpringDampMass`](https://github.com/ColtonKawamura/GranE/wiki/matSpringDampMass).

4. **Plot and inspect the results**:
   - Attenuation vs. frequency with [`plotModes`](https://github.com/ColtonKawamura/GranE/wiki/plotModes)
   - Wave speed with [`eigenWavespeed`](https://github.com/ColtonKawamura/GranE/wiki/eigenWavespeed)
   - Saved variables described in the [output files page](https://github.com/ColtonKawamura/GranE/wiki/Output-file)

5. For batch runs on a cluster, see the runner scripts in `hpc/` (see the wiki workflow page for details).

## Roadmap

- [x] 2D packing generation
- [x] 3D packing generation
- [ ] Hertzian contacts
