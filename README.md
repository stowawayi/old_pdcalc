# PDCALC — Physical Damage Calculator

A Fortran 77 implementation of the Probability of Damage (POD) calculation system for nuclear weapon effects analysis, originally developed in 1976 and updated through March 1984.

## Overview

PDCALC computes weapon effects metrics given target characteristics, weapon yield, burst geometry, and weapon accuracy. It supports multiple output modes:

- **POD** — Probability of Damage at a given ground range and height of burst
- **Weapon Radius (WR)** — Effective radius needed to achieve a specified damage level
- **Maximum Distance** — Furthest range at which a specified POD can be achieved
- **Fatality/Casualty estimates** — For personnel targets

### Supported Target Types

| Type | Description |
|------|-------------|
| P    | Overpressure effects |
| Q    | Dynamic pressure effects |
| X    | Personnel/casualty |
| R, S, T | Crater effects |
| A–F, Y, Z | Equivalent Target Area (ETA) types |
| U, L, M, N, O | Specialized target classifications |

## Project History & Modifications

The original 1976 PDCALC code has been modernized:

- Replaced `EQUIVALENCE` statements with combined data arrays (modern Fortran allows >19 continuation lines)
- Improved readability (whitespace, formatting)
- Entire code uses double precision (`real*8`) throughout, including all data statements
- Added Speicher & Brode routines for Peak Dynamic Pressure and Peak Overpressure
- Added routine to find the optimum HOB/ground range given weapon and target characteristics

> **Note:** You must zero out `d`, `wr`, and `pod` before calling `pdcalc` with option 2, otherwise it returns the last "good" value when the computed result is zero.

## Directory Structure

```
old_pdcalc/
├── main.for            # Entry point; dispatches to one of 8 driver programs
├── pdcalc.f            # Core damage calculation engine
├── wrcalc.f            # Weapon radius calculator (P/Q type targets)
├── lncalc.f            # Probability calculations (Gaussian quadrature)
├── wrpers.f            # Personnel casualty calculations
├── etcalc.f            # Equivalent Target Area calculations
├── wrclcy.f            # Cycle-based weapon radius
├── wrcrtr.f            # Crater effects calculations
├── dypres.f            # Dynamic pressure calculations
├── overp.f             # Overpressure calculations
├── ogh.f               # Optional graphics handler (legacy)
├── acon.f              # Mathematical constant initialization
├── errmsg.f            # Error message handling
├── intgf.f             # Integration functions
├── ranf.f              # Random number functions
├── *.h                 # Header files (real8.h, const.h, files.h, basicd.h, cdkpd.h, cdkwr.h)
├── drv[1-8].x          # Driver source files
├── drv[1-8].nml        # Namelist configuration files for each driver
├── gd.dat              # Sample ground range/HOB/yield data (100 kt)
├── ussr.dat            # USSR target/missile data
├── good                # Reformatted RISOP target database
├── types.tgt           # Target type reference data
├── gd.py               # Python script to generate gd.dat test data
├── adjvn.py            # Python utility to adjust vulnerability numbers
├── runv                # Test runner script
├── CMakeLists.txt      # CMake build configuration
├── Makefile            # Traditional make configuration
└── docs/
    ├── PDCALC_April1976.pdf
    ├── PDCALC_March1984.pdf
    ├── MathematicalBackgroundAndProgrammingAidForThePhysicalVulnerabilityHandbook.pdf
    └── b010375_MathematicalBackgroundAndProgrammingAidsForVN.pdf
```

## Building

### Requirements

- A Fortran 77-compatible compiler: `gfortran` (recommended), Intel Fortran (`ifort`), or Absoft
- GNU Make or CMake ≥ 3.10

### Using Make

```bash
make        # Build the main executable
make lib    # Build libpd.a (static library for linking with other codes)
make clean  # Remove build artifacts
```

### Using CMake

```bash
mkdir build && cd build
cmake ..
make
```

CMake automatically selects compiler flags based on the detected compiler (GNU, Intel, or Absoft).

## Running

```bash
./main <driver_number>                  # Run with default namelist (drvX.nml)
./main <driver_number> -i <namelist>    # Run with alternate namelist file
```

There are 8 driver modes, each controlled by a corresponding namelist file (`drv1.nml` through `drv8.nml`).

---

## Driver Reference

### Driver 1 — Ground Range / HOB / Yield File

Reads a file of ground range, height of burst, and yield values and computes damage outputs.

```bash
./main 1
```

**Namelist parameters (`drv1.nml`):**

| Parameter | Description |
|-----------|-------------|
| `cep`     | Circular Error Probable for weapon (km) |
| `r95`     | Radius within which 95% of damage lies (km) |
| `gname`   | Input filename (ground range, altitude, yield) |
| `ivns`    | Starting VN number |
| `ivne`    | Ending VN number |
| `ivnd`    | Step size for VN |
| `jti`     | Target type character (lowercase, e.g. `'p'`) |
| `kfi`     | K-factor |
| `iflg`    | Output mode (see table below) |
| `mode`    | Input format (1 = grn/hob/yld; 2 = grn/hob/yld/iv/jti/kfi) |
| `gamma`   | Reentry angle (degrees) |
| `az`      | Azimuth from DGZ to target (degrees) |

**`iflg` output modes:**

| `iflg` | Output |
|--------|--------|
| 1      | Product POD up to 0.990; `d` must be input |
| 2      | Product POD up to 0.999; `d` must be input |
| 3      | Same as 4 |
| 4      | Weapon radius |
| 5      | Same as 6 |
| 6      | Maximum distance at which a given POD can be achieved; `pod` must be input |
| 7      | Fatality POD and casualty POD (returned in `pod` & `wr`); `d` must be input |
| 8      | Damage sigma input via `pod`; `pod` is output; `d` is input |
| 9      | Damage sigma and weapon radius are input; `pod` is output; `d` is input |
| 10     | Weapon radius input; `pod` is output; `d` is input |

---

### Driver 2 — VN and K-Factor Range Sweep

Exercises the code across ranges of vulnerability number (VN) and k-factor.

```bash
./main 2
```

**Namelist parameters (`drv2.nml`):**

| Parameter | Description |
|-----------|-------------|
| `hobi`    | Height of burst for evaluation (ft) |
| `cep`     | Circular Error Probable (km) |
| `r95`     | 95% damage radius (km) |
| `gname`   | Input data filename |
| `ivnb`    | Starting VN number |
| `ivne`    | Ending VN number |
| `ivns`    | VN step size |
| `jti`     | Target type |
| `kfb`     | K-factor begin value (typically 0; end = min(ivn−1, 9)) |
| `kfs`     | K-factor step (typically 1) |
| `az`      | Azimuth from DGZ to target (degrees) |

---

### Driver 3 — Reformatted RISOP File

Reads a reformatted OPEN-RISOP target file and computes POD for all targets (default: 300 kt yield, 3 nmi CEP).

```bash
./main 3
```

**Namelist parameters (`drv3.nml`):**

| Parameter | Description |
|-----------|-------------|
| `yld`     | Yield (kt) |
| `hobi`    | Height of burst (ft) |
| `cep`     | Circular Error Probable (km) |
| `fname`   | Name of reformatted RISOP file |

**Input file format** (free format, no character fields):

```
line, catcode, xlat, xlon, hgt, r95, vntk, ivn, jti, kfi
```

| Field    | Description |
|----------|-------------|
| `line`   | Line number in original file |
| `xlat`   | Latitude (degrees) |
| `xlon`   | Longitude (degrees) |
| `hgt`    | Target height (ft) |
| `r95`    | 95% damage radius (km) |
| `vntk`   | Full VNTK string (e.g. `10P0`) |
| `ivn`    | VN number portion only |
| `jti`    | Target type (lowercase, e.g. `p`) |
| `kfi`    | K-factor |

---

### Driver 4 — ICBM Parameters with HOB Sweep

Reads a list of available VN numbers with associated CEPs and an ICBM parameter file, then loops over heights of burst to determine POD.

```bash
./main 4
```

**Namelist parameters (`drv4.nml`):**

| Parameter | Description |
|-----------|-------------|
| `iname`   | Filename containing ICBM information |
| `tname`   | Filename containing target types, r95, and VNTK data |
| `hob0`    | Starting HOB (ft) |
| `hob1`    | Ending HOB (ft) |
| `dhob`    | HOB step size (ft) |

**ICBM file format** (`a22,1x,a20,1x,i5,1x,i5,1x,f12.5,1x,f12.5,1x,f12.5,1x,i3`):

| Field   | Description |
|---------|-------------|
| `tnato` | NATO name for missile |
| `tussr` | USSR name for missile |
| `nsoft` | Number of missiles on soft launchers (future use) |
| `nhard` | Number of missiles on hard launchers (future use) |
| `yld`   | Yield (kt) |
| `rmax`  | Maximum range (km) (future use) |
| `cep`   | CEP of missile (ft) |
| `numwh` | Number of warheads per missile (future use) |

**Target file format** (`i5,1x,a20,1x,f17.5,1x,i2,1x,a,1x,a`):

| Field   | Description |
|---------|-------------|
| `num`   | Line number from original file |
| `akey`  | VNTK-R95 key string |
| `r95`   | 95% damage radius |
| `ivn`   | VN number portion |
| `jtinp` | Target type (uppercase; code converts to lowercase internally) |
| `kfi`   | K-factor |

---

### Drivers 5–8

Additional specialized drivers. See `drv5.nml` through `drv8.nml` for their respective namelist parameters.

---

## Key Concepts

**CEP (Circular Error Probable):** The radius within which 50% of weapons land. A smaller CEP means greater accuracy.

**Vulnerability Number (VN):** A standardized measure of a target's resistance to weapon effects. Combined with target type and k-factor to form the full VNTK descriptor (e.g., `10P0`).

**K-Factor:** An adjustment factor applied to the vulnerability number to account for specific target hardening or exposure conditions.

**Height of Burst (HOB):** The altitude at detonation. Optimum HOB maximizes damage area for a given yield and target type.

**Ground Range (GR):** Horizontal distance from the detonation point directly above ground zero (DGZ) to the target.

---

## Covariance Delivery Error Mode

In addition to the original scalar CEP input, a 3×3 delivery error covariance mode is available. It models anisotropic weapon accuracy and correlated HOB error without changing any existing input files.

### Overview

The covariance mode is activated by supplying an optional `/covlst/` namelist **after** the standard `/plst/` namelist in the input file. If `/covlst/` is absent, or if `sig_cr` and `sig_dr` are both zero, the driver falls back to the original CEP path (`pdcalc` → `lncalc`) with no change in behavior.

When `sig_cr > 0` or `sig_dr > 0`, the driver calls `pdcov` instead of `pdcalc`. `pdcov` dispatches to `lncov`, which performs a 125-point (5×5×5) Gauss-Hermite quadrature over the full 3×3 covariance.

### Coordinate frame

Errors are specified in the **weapon delivery frame**:

| Axis | Description |
|------|-------------|
| Cross-range (`sig_cr`) | Perpendicular to the flight path in the horizontal plane |
| Down-range (`sig_dr`) | Along the flight path in the horizontal plane |
| HOB (`sig_hob`) | Vertical (altitude) |

For ICBMs/SLBMs this frame is aligned with atmospheric reentry at approximately 300,000 ft.

### `/covlst/` namelist parameters

| Parameter | Units | Description |
|-----------|-------|-------------|
| `sig_cr` | km | Cross-range delivery 1-sigma (converted to feet internally) |
| `sig_dr` | km | Down-range delivery 1-sigma (converted to feet internally) |
| `sig_hob` | feet | HOB delivery 1-sigma |
| `rho_cr_dr` | — | Cross-range / down-range correlation (−1 to +1) |
| `rho_cr_hob` | — | Cross-range / HOB correlation (−1 to +1) |
| `rho_dr_hob` | — | Down-range / HOB correlation (−1 to +1) |

**Driver 4 exception:** `sig_cr` is read from the missile data file (one column per missile). The `/covlst/` namelist uses `sig_dr_ratio` instead of `sig_cr`:

| Parameter | Units | Description |
|-----------|-------|-------------|
| `sig_dr_ratio` | — | `sig_dr = sig_cr × sig_dr_ratio`; set > 0 to activate covariance path |
| `sig_hob` | feet | HOB delivery 1-sigma |
| `rho_cr_dr` | — | Cross-range / down-range correlation |
| `rho_cr_hob` | — | Cross-range / HOB correlation |
| `rho_dr_hob` | — | Down-range / HOB correlation |

### Example input file

```fortran
 &plst
  ivn=10, jti='p', kfi='0',
  yld=550.0, r95=0.1, cep=0.0,
  hob0=0.0, hob1=5000.0, dhob=100.0,
  d0=0.0, d1=1.5, nd=64, az=0.0
 /
 &covlst
  sig_cr=0.150, sig_dr=0.100, sig_hob=150.0,
  rho_cr_dr=0.0, rho_cr_hob=0.0, rho_dr_hob=0.0
 /
```

`sig_cr` and `sig_dr` are in km; the driver converts them to feet. `sig_hob` is in feet. All correlation values default to 0 if omitted. The `/plst/` `cep` value is ignored when the covariance path is active.

### Algorithm summary

1. HOB sensitivity `α = d(WR)/d(HOB)` computed by central difference of two `wrcalc` calls
2. Conditional horizontal covariance via Schur complement: `Σ_h|z = Σ_h − b bᵀ / σ²_hob`
3. Eigendecompose 2×2 conditional covariance into independent principal axes
4. Triple 5-point Gauss-Hermite quadrature: outer = HOB (or single pass if `sig_hob=0`), inner = 2D horizontal eigenbasis
5. Target size `r95` added isotropically: `0.231 × r95²` on both horizontal diagonal terms

---

## Documentation

Reference documentation is included in the repository:

- `PDCALC_April1976.pdf` — Original PDCALC system documentation
- `PDCALC_March1984.pdf` — Updated documentation with enhancements
- `MathematicalBackgroundAndProgrammingAidForThePhysicalVulnerabilityHandbook.pdf`
- `b010375_MathematicalBackgroundAndProgrammingAidsForVN.pdf`
