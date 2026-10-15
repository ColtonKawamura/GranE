# AGENTS.md

## Workflow

**Every single edit to any file — no exceptions — must be done on a new branch.**

Steps:
1. FIRST: `git checkout -b <short-task-name>`
2. Make all edits on that branch
3. When done: `git add -A && git commit -m "<description>"`
4. Confirm to the user which branch the changes are on

If you are not on a feature branch, stop and create one before touching any file.
Never work directly on main or master.
Only modify files that already exist in this directory.

## Verification

After any edit to a file in `src/`, run the corresponding test:

```
matlab -batch "run('test/TEST<filename>.m')"
```

Where `<filename>` is the exact name of the file you edited, including the `.m` extension.

**Examples:**
- Edited `src/meetingAug2026.m` → run `matlab -batch "run('test/TESTmeetingAug2026.m')"`
- Edited `src/computeWaveVector.m` → run `matlab -batch "run('test/TESTcomputeWaveVector.m')"`

All assertions must pass before finishing. If the test fails, fix the source file and re-run. Do not modify the test file.

If no corresponding test file exists in `test/`, stop and notify the user before proceeding.

## Key Directories

- **`src/`** - MATLAB source files for packing generation and data processing
- **`data/`** - Simulation output (`.mat` files). Organized by dimension (`2d`, `3d`) and contact type
- **`hpc/`** - HPC runner scripts and Python utilities for batch job submission
- **`test/`** - Agent verification tests. Agent may read but **must not edit** these files.

 # Running tests

 Always run matlab from the project root:
`matlab -batch "cd('/Users/coltonkawamura/repos/GranE'); run('test/TEST<filename>.m')"`
