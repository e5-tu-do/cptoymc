# CPToyGen

A ToyMC generator for CP studies with P2VP decays (like Bd2JpsiKS and Bs2JpsiKS).

## Remarks for developers
- The `master` branch should always correspond to a tagged version.
- All features should be merged into the `develop` branch.
- Features should be developped in a `feature-xyz` branch, where `xyz` should be replaced with some meaningful name. Feature branches should have a limited lifetime. After finishing the development of the feature, merge its branch to `develop` and delete it (both local and remote).
- Adaptations for specific analyses should be developped in a `dev-analysisname` branch.
- Release branches have the form `release-x.y` and are based on the develop branch. After preparation, the release is merged into the master, the master is tagged, and the release branch is deleted (local and remote).
