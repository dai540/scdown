# Changelog

## scdown 0.5.1

- removed non-package GitHub issue, PR, and release automation files to
  keep the repository focused on the R package itself
- simplified `.Rbuildignore` and `.gitignore`
- rebuilt pkgdown after pruning repository-only files

## scdown 0.5.0

- simplified the package to annotation-after, within-dataset downstream
  analysis
- removed group-aware tests, handoff helpers, composition summaries, and
  pipeline guidance
- focused the API on maps, markers, marker-panel annotation checks,
  average expression, cluster similarity, signature scoring, exploratory
  communication, and compact reporting
- cleaned the documentation, vignette, tests, and generated site around
  the smaller scope

## scdown 0.4.2

- fixed the GitHub release workflow by installing vignette-building
  dependencies before `R CMD build`
- added manual dispatch support to the release workflow so a release can
  be rerun without creating another tag first

## scdown 0.4.1

- added a release workflow that can publish GitHub releases from pushed
  version tags
- expanded CI from a single Linux check to a multi-OS matrix with an
  extra devel check on Linux
- split and strengthened tests around sparse inputs, failure modes,
  report output, and handoff schemas
- added repository-facing metadata files such as `CITATION.cff` and
  issue templates to improve public package trust signals
- clarified in the HTML report that built-in tests are
  screening-oriented and should hand off to specialist methods for
  primary claims

## scdown 0.4.0

- added `scope_table()` to make the intended boundary of `scdown`
  explicit
- added `export_for_handoff()` to prepare objects for `muscat`,
  `dreamlet`, `miloR`, and CellChat workflows
- expanded
  [`build_scdown_report()`](https://dai540.github.io/scdown/reference/build_scdown_report.md)
  with scope guidance plus marker and communication test sections
- strengthened README, pkgdown, and the recommended pipeline vignette
  around exploratory-vs-specialist handoff

## scdown 0.3.0

- clarified `scdown` as an annotation-after downstream hub rather than a
  full upstream workflow
- added a recommended pipeline helper and a handoff vignette for PBMC,
  tumor immune, and patient comparison use cases
- updated README and pkgdown site to explain where `scdown` fits and
  which specialist tools to use before or after it

## scdown 0.2.0

- switched to a sparse-first internal matrix representation
- added custom signatures, marker panels, lineage mapping, and
  ligand-receptor pairs
- split exploratory functions from sample-aware testing functions
- added CI workflows and build hygiene files

## scdown 0.1.0

- Initial release.
