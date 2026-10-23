## RPesto CRAN submission

This submission is for fixing the R check, specifically avoiding the non-API calls

I have done an R CMD check on the following platforms
- Ubuntu 24.04.5 LTS-R (R-devel)
- Mac M1 (R-devel)
- Mac M1 (R-4.6.1)
- Windows Server 2025 x64 (R-devel)

RPesto includes rust code, and so it depends on cargo and rustc. I have bundled the rust dependencies in an archive src/rust/vendor.tar.xz such that the installation script does not download any external crates.

I have also used the CRAN win-builder service, with the result 

Installation time in seconds: 126
Check time in seconds: 118
Status: OK
R Under development (unstable) (2026-10-22 r90748 ucrt)