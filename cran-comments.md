## RPesto initial CRAN submission

This is the initial submission of the RPesto package.

I have done an R CMD check on the following platforms
- Ubuntu 24.04.3 LTS-R (R-devel)
- openSUSE Leap 15.6 (R 4.5.1)
- macOS (R 4.5.1)
with the following result

0 errors | 0 warnings | 0 note

RPesto includes rust code, and so it depends on cargo and rustc. I have bundled the rust dependencies in an archive src/rust/vendor.tar.xz such that the installation script does not download any external crates.

The CRAN check has one note, 

0 errors | 0 warnings | 1 note

about the spelling of words in the description. The following words are not mis-spelled:
  Hoehna (7:197)
  Kopperud (7:186)
