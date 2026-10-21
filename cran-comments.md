## RPesto CRAN submission

This submission is for fixing the R check, specifically avoiding the non-API calls

I have done an R CMD check on the following platforms
- Ubuntu 24.04.5 LTS-R (R-devel)
- Mac M1 (R-devel)
- Mac M1 (R-4.6.1)
- Windows Server 2025 x64 (R-devel)

On Windows, there was one warning:

    Warning: Warning: file 'RPesto/configure' did not have execute permissions: corrected

I am not sure how to correct this, but on the other platforms, there were no issues 

RPesto includes rust code, and so it depends on cargo and rustc. I have bundled the rust dependencies in an archive src/rust/vendor.tar.xz such that the installation script does not download any external crates.

I have also used the CRAN win-builder service, with the result 
