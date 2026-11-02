# Shared numerical fixtures

JSON fixtures shared with the Python package UnitRoot-Tools. Both
implementations must reproduce every fixture within its stated tolerance
(see the Validation standard in `spec/specification.qmd`).

Each file contains:

- `id`, `description`
- `source`: provenance of the data and of the expected values, with a hash of the data source
- `data`: the series
- `settings`: test or strategy options
- `expected`: statistics, critical values, decisions, history
- `tolerance`: `rtol` and `atol`, with a justification when looser than the default
