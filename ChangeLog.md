<!--
SPDX-FileCopyrightText: 2022-2026 Helmholtz-Zentrum hereon GmbH
SPDX-License-Identifier: CC0-1.0
-->

# Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- `omexdia_c`: OMExDia extended with methane and sulfate cycling
- `omexdia_c_bottom`: trait-based, bottom-only variant of the methane/sulfate OMExDia model
- Zenodo metadata (`.zenodo.json`) with full creator/contributor list, ORCID identifiers, and funding acknowledgements
- Pre-commit hooks enforcing REUSE license compliance and validating `.zenodo.json` against the Zenodo metadata schema

### Changed
- Updated Readme with nope
- Coupled `omexdia_c` via fluxes rather than states
- Used salinity as a proxy for sulfate in `omexdia_c`
- Added porosity to `omexdia_c`
- Replaced Monod kinetics with exponential functions in `omexdia_c`

### Fixed
- Added missing SPDX license headers to `omexdia_c.F90` and `omexdia_c_bottom.F90`

## [0.3.0] - 2026-08-17

### Fixed
- Completed licensing fix
- Started SPDX cleanup

### Changed
- Fixed `dneitrification_scale_factor` name

### Added
- Added denit/nitri scale factors to parameters

## [0.2.0] - 2025-07-11

### Changed
- Adjusted model name
- Finalized comments
- Corrected k calculation
- Cleaned up folders
- Modified Tang equation
- Improved comments
- Adjusted filenames
- Removed unnecessary parts
- Removed unused code
- Changed name
- Adjusted N2O yield calculation

### Added
- Sediment oxygen uptake diagnostic
- Option to use different emission factors
- Omexdia_n2o_nope and nope modules

### Fixed
- Debugged sediment oxygen uptake diagnostic

## [0.1.0] - 2024-12-12

### Added
- Initial project structure
- First working version
- Monthly detritus fluxes
- Total nitrogen calculations
- Diagnostics
- Omexdia bottom model (initial version)
- GOTM, 0D, and Python instructions
- FABM configuration file
- Omexdia_p modernization
- Extended light model with horizontally varying parameters
- Hereon structure
- PyFABM test case

### Changed
- Removed GOTM and oxypom submodels
- Fixed naming conflicts
- Adjusted directory structure
- Renamed HZG to Hereon

### Removed
- Removed all unnecessary code
- Deleted old omexdia_n2o files
- Moved old fabm.nml testcases

[Unreleased]: https://github.com/username/project/compare/v0.3.0...HEAD
[0.3.0]: https://github.com/username/project/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/username/project/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/username/project/commits/v0.1.0