# scdown 0.4.0

- added `scope_table()` to make the intended boundary of `scdown` explicit
- added `export_for_handoff()` to prepare objects for `muscat`, `dreamlet`, `miloR`, and CellChat workflows
- expanded `build_scdown_report()` with scope guidance plus marker and communication test sections
- strengthened README, pkgdown, and the recommended pipeline vignette around exploratory-vs-specialist handoff

# scdown 0.3.0

- clarified `scdown` as an annotation-after downstream hub rather than a full upstream workflow
- added a recommended pipeline helper and a handoff vignette for PBMC, tumor immune, and patient comparison use cases
- updated README and pkgdown site to explain where `scdown` fits and which specialist tools to use before or after it

# scdown 0.2.0

- switched to a sparse-first internal matrix representation
- added custom signatures, marker panels, lineage mapping, and ligand-receptor pairs
- split exploratory functions from sample-aware testing functions
- added CI workflows and build hygiene files

# scdown 0.1.0

- Initial release.
