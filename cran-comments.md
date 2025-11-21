## RPesto second CRAN submission

This is the second submission of the RPesto package.

I have fixed the .Rd files as per the previous instructions, thanks.

In the previous message, I was told that writing to files in the package directory is not allowed. I understand this makes sense for user functions. However, the "Writing R Extensions" document specifically permits for a `configure` script to write/create a `src/Makevars` file (see https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Configure-example). Thus, is the flagging of my `tools/config.R` a false positive? See also the discussion in this github issue: https://github.com/extendr/rextendr/issues/453

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
