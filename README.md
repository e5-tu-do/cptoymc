# CPToyGen

A generator of pseudo-datasets of decays with time-dependent CP violation.

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.20356.svg)](http://dx.doi.org/10.5281/zenodo.20356)

## Remarks for developers
- The `master` branch should always correspond to a tagged version.
- All features should be merged into the `develop` branch.
- Features should be developed in a `feature-xyz` branch, where `xyz` should be replaced with some meaningful name. Feature branches should have a limited lifetime. After finishing the development of the feature, merge its branch to `develop` and delete it (both local and remote).
- Adaptations for specific analyses should be developed in a `dev-analysisname` branch.
- Release branches have the form `release-x.y` and are based on the develop branch. After preparation, the release is merged into the master, the master is tagged, and the release branch is deleted (local and remote).
