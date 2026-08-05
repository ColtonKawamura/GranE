Only modify files that already exist in this directory.

## Project Structure

```
GranE/
├── src/              # MATLAB source code
│   ├── pack.m             # 2D/3D packing generation
│   ├── packRepeatTile.m   # Repeat tile packing routine
│   └── cleanRats.m        # Data cleaning utility
├── data/             # Simulation data and results
│   ├── junkyard/          # Scratch/intermediate .mat files
│   └── packings/
│       ├── 2d/            # 2D packing data
│       │   └── hertz/     # Hertzian contact 2D data
│       └── 3d/            # 3D packing data
├── hpc/              # High-performance computing scripts
│   ├── CreatePack.py      # Python script for creating packings
│   ├── RunChainGPU.sh     # GPU job chain runner
│   ├── RunChainJob.sh     # Job chain runner
│   ├── run_one.sh         # Run single simulation
│   ├── commandsPack.txt   # Pack command templates
│   └── localJobs.txt      # Local job configurations
└── README.md         # Project overview and roadmap
```

### Key Directories

- **`src/`** - MATLAB source files for packing generation and data processing
- **`data/`** - Simulation output (`.mat` files). Organized by dimension (`2d`, `3d`) and contact type
- **`hpc/`** - HPC runner scripts and Python utilities for batch job submission

## Workflow

Before starting any task:
1. Create a new git branch: `git checkout -b <short-task-name>`
2. Make all edits on that branch
3. Stage and commit when done: `git add -A && git commit -m "<description>"`

Never work directly on main or master.
