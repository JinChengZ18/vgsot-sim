# vgsot-sim documentation

Documentation for the `vgsot-sim` VGSOT/SOT-MTJ switching simulator. Start with
the root [`README.md`](../README.md) for installation and quick-start examples;
the pages below go deeper.

## User guide

- [Simulation cases](cases.md) — the four CLI cases and the physics each one demonstrates
- [Default parameters](parameters.md) — per-case parameter tables and what they control
- [Parameter validation](parameter_validation.md) — how the calibrated defaults are checked against the device data

## Reference

- [API reference](api.md) — configs, result dataclasses, high-level cases, and low-level kernels
- [Project structure](structure.md) — package layout, information flow, and public API surface
- [Technical details](technical_details.md) — physics derivations, calibration notes, and figure provenance
- [Reproducing the Chapter 2 figures and tables](reproducing_figures.md) — every figure and table mapped to its script, command, input data, and runtime

## Development

Maintainer-facing notes, kept separate from the user documentation above:

- [Implementation status](maintenance/IMPLEMENTATION_STATUS.md) — physics-feature coverage map (implemented / pending / known caveats)
- [Version notes](maintenance/version_notes.md) — release-by-release physics changes
- [Claims backing matrix](../CLAIMS.md) — every documented claim mapped to its test or committed artifact
