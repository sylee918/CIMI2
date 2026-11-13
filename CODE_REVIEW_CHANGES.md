# CIMI Code Review and Optimization Summary

- **Date**: 2026-11-13
- **Reviewer**: Sang-Yun Lee // Antigravity Gemini 4 Argon
- **Repository**: `c:\Users\slee122\Cursor\CIMI`
- **Build Verification**: `ifx -r8 -132 -save -O3 -xHost -qopenmp -mcmodel=medium -Wl,--no-relax` (`make EXE` exited with code 0, no errors)

---

## Modified Files Summary (`git diff --stat`)

| File | Diff Size | Key Areas |
| :--- | :--- | :--- |
| `src/cimi.f90` | 300 lines (+/-) | Physics formula typo, `NaN` prevention, flag logic, OpenMP parallelization, search/loop optimization |
| `src/ModMate.f90` | 138 lines (+/-) | 4D array memory layout reorder, time/longitude boundary handling, `parameter` conversion |
| `src/ModIonDiff.f90` | 84 lines (+/-) | Remove `contiguous`, uninitialized variable fix, binary search hoisting |
| `src/ModTwinsHden.f90` | 71 lines (+/-) | Thread safety (`recursive`, local arrays), `NumHour`/`acos` bounds guards |
| `src/Lstar2.f90` | 47 lines (+/-) | Eliminate inner-loop `allocate`/`log10`, index bounds guards |
| `src/ModCurvScatt.f90` | 3 lines (+/-) | Prevent `NaN` in real-exponent power of curvature derivatives |
| `src/w05sc_fok.f90` | 3 lines (+/-) | Initialize `skip` variable under `-save` |
| `src/rcm6.f90` | 1 line (+) | Initialize `lsflag` array under `-save` |

---

## 1. Critical Bug Fixes and Logic Errors

### 1.1 Physics Formula Typo & `NaN` Prevention
* **Dipole Magnetic Field Formula Typo (`src/cimi.f90` — `fieldpara`, L1899)**
  * Fixed the ionospheric reference-altitude (`rc`) dipole field magnitude formula from `1.+3.*sinlat` to `1.+3.*sinlat*sinlat` ($\sqrt{1 + 3\sin^2\lambda}$).
* **Low-Latitude Precipitation Flux `NaN` Prevention (`src/cimi.f90` — `p_preci`, L4188–4205)**
  * At the lowest latitude grid point (`i=1`, $r \approx 1.097\,R_E$), which lies below the 800 km altitude shell ($r_{800} = 1.1256\,R_E$), `bo(i,j)/B800 > 1` caused `sqrt(negative)` to write `NaN` into `SA800`, `prePi`, and `prePe`. Clamped `max(0., 1.-bo(i,j)/Brc)` and `max(0., 1.-bo(i,j)/B800)` and guarded with `if (B800.gt.0..and.SArc.gt.0.)`.
* **Curvature Scattering Power-Law `NaN` Guard (`src/ModCurvScatt.f90` — L297)**
  * Protected `zeta1_flc(i,j)**a1 * zeta2_flc(i,j)**a2` with `max(abs(...), 1.e-6)` to prevent `NaN` when curvature derivatives become negative or zero in stretched magnetotail configurations.
* **`acos` Domain Clamping (`src/ModTwinsHden.f90`, `src/cimi.f90` — `hodges_geocorona`)**
  * Clamped `cos_phi` and `cos_theta` to `[-1., 1.]` before `acos` and added `radxy > 0.`, `rad > 0.` guards to prevent `NaN` from floating-point rounding or origin evaluation.

### 1.2 Input Flag and Branching Logic Fixes
* **Wave Diffusion Flag Logic (`src/cimi.f90` — L524–530)**
  * Fixed `if ((ichor.ne.0.or.ihiss.ne.0).and.n_e.gt.0) iwave=1` and moved `call readDiffCoef(ichor,ihiss)` inside `if (iwave.eq.1)` so ion-only (`n_e=0`) or EMIC-only runs do not access `ekev(0,:)` out of bounds or trigger electron VLF diffusion.

### 1.3 Uninitialized Variables & Stale State Under `-save`
* **`lsflag` Reset (`src/rcm6.f90` — `aphi`, L606)**: Added `lsflag(j) = .false.` at the start of `do j=2,jdim-1` so boundary flags do not remain latched to `.true.` across timesteps.
* **`skip` Reset (`src/w05sc_fok.f90` — `epotval`, `mpfac`, `scplm`)**: Initialized `skip = 0` before `do j=1,csize`.
* **`EmicOpowerInPs` / `EmicOpowerOutPs` Initialization (`src/ModIonDiff.f90` — `diffuse_UV_ions`)**: Initialized both to `0.` before the early-exit `goto 9999` path.
* **`find_ceSigma` Initialization & Electron Handling (`src/cimi.f90` — L2521)**: Initialized `d = -99.` and explicitly set `ceSigma(n,k,1) = 0.` for electrons (`js(n) >= 5`).
* **`lossconeN` / `lossconeS` & `mapgrid` (`src/cimi.f90`)**: Initialized `dkc=1., dko=0.` when `mc == 0` and initialized `ibp1=1` in `mapgrid`.

### 1.4 Array Bounds (Off-by-One) and Interpolation Boundary Handling
* **`ModMate` Time/Longitude/Search Bounds (`src/ModMate.f90`)**
  * `mate_geocorona`: Added boundary clamping for `nHourMate <= 1`, `t <= TimeMate(1)`, and `t >= TimeMate(nHourMate)`.
  * `lintp3`: Added periodic longitude wrap-around interpolation (`z(nz) -> z(1)+360°`) when `z(nz) < 360.0` and `z1 > z(nz)`.
  * `locate_fast`: Returned `j = max(1, n - 1)` when `x >= xx(n)`.
  * `lintp4`: Fixed swapped `mint`/`maxt` and `minx`/`maxx` bounds variables.
* **`Lstar2.f90` Bounds Guards**: Added `if (ii.gt.1)` (L100), `if (ii1.lt.iba(j)...)` (L135), and `if (i.gt.1)` in `locateB` (L186).
* **`EquatorialXing` & `calc_flc_para` (`src/cimi.f90`)**: Added `minloc(bs(1:npf),dim=1)` fallback when no strict interior minimum exists, and guarded `calc_flc_para` with `if (imn.ge.4 .and. imn+3.le.n)`.
* **`NumHour` Clamping (`src/cimi.f90` L1740, `src/ModTwinsHden.f90`)**: Clamped `NumHour = min(max(int(t/3600.), 0), 71)` to prevent out-of-bounds table lookups.

---

## 2. Performance, Memory, and OpenMP Optimization

* **Column-Major 4D Array Layout in `ModMate` (`src/ModMate.f90`)**
  * Reordered `logHdenMate` from `(nHourMate, nRMate, nLatMate, nLonMate)` to `(nRMate, nLatMate, nLonMate, nHourMate)` so `logHdenMate(:,:,:,iTime)` passed to `lintp3` is contiguous in memory, eliminating 3D temporary array copies on every call.
* **Eliminate Non-Contiguous Slice Copies in `diffuse_UV_ions` (`src/ModIonDiff.f90`)**
  * Removed `contiguous` from dummy arguments (`Gjac0`, `y0`, `tya0`, `f2`, etc.) so passing non-contiguous slices from `cimi.f90` inside OpenMP loops no longer forces temporary heap allocation and copy-in/copy-out per `(i,j)` iteration.
* **Eliminate Inner-Loop Heap Allocation in `calc_Lstar2` (`src/Lstar2.f90`)**
  * Removed millions of inner-loop `allocate(logbm_i)` / `log10(bmm)` / `deallocate(logbm_i)` calls by precomputing `logbmm(1:iba(j),j)` once per `m` at the top of the `m` loop.
* **OpenMP Parallelization of `sume` & `sum3band` + $O(1)$ Energy-Band Indexing (`src/cimi.f90`)**
  * Added `!$omp parallel do` to serial bottlenecks `sume` and `sum3band` (with `reduction(+:rbsum_loc,rcsum_loc)` in `sum3band`).
  * Replaced the 24-iteration `kkloop` linear search in `sume`, `sum3band`, `Bdrift`, and `Edrift` with direct $O(1)$ indexing `kk = (k + 1) / 2` (leveraging `iw = 2*je`).
* **Hoist Loop-Invariant Pitch-Angle Binary Search (`src/cimi.f90` — `diffuse_VLF`, `src/ModIonDiff.f90` — `diffuse_UV_ions`)**
  * Moved `locate1(wPA, mpa, ao, jpa)` out of the `k=0,iw+1` loop (to right after the `m` loop), reducing search calls by 150x.
* **Remove Redundant $O(N)$ Scans in Interpolation/Search Routines (`src/cimi.f90`, `src/ModIonDiff.f90`, `src/ModMate.f90`)**
  * Removed the per-call `monoCheck` loop in `locate1` and replaced full-array `minval`/`maxval` scans in `lintp`, `lintp2`, and `lintp2b` with $O(1)$ endpoint checks (`min(x(1),x(n))`, `max(x(1),x(n))`).
* **Column-Major Loop Ordering (`src/cimi.f90` — `VdriftB`, `VdriftE`, `p_result`)**
  * Reordered nested loops to match Fortran column-major memory layout (`j` outer / `i` inner, `m` outer / `k` inner) and added `i > 1` guards.

---

## 3. Readability, Thread Safety, and Maintainability

* **Convert Implicitly `SAVE`d Constants to `parameter` (`src/ModMate.f90`, `src/ModIonDiff.f90`)**
  * Changed declaration-initialized constants (`cRadToDeg`, `cDegToRad`, `cDayToSec`, `cHourToSec`, `cSecToHour`, `iDaa`, `iDpp`) to `real, parameter`.
* **Thread-Safe `twins_geocorona` (`src/ModTwinsHden.f90`)**
  * Added `recursive` and `implicit none`, and replaced shared module buffer writes (`param2beUsed`, `Tmatrix2beUsed`) with thread-local arrays (`p`, `tm`).
