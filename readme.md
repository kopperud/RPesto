## RPesto

`RPesto` is an alternative implementation of the `Pesto` software in R. It can estimate branch-specific diversification rates for large trees.

## Installation

The backend of `RPesto` is implemented in `rust`, therefore you need to have rust and cargo installed. To install this, I recommend to follow the instructions on [https://rustup.rs/](https://rustup.rs/).

Next, we will use the R-package `remotes` to install this package

```R
install.packages("remotes")

library(remotes)
install_github("kopperud/RPesto")
```

Now you can load `RPesto` and fit the model to your tree

```R
library(RPesto)
library(ggtree)

data("primates")

sampling_fraction <- 0.635

phy <- fit_bds(primates, sampling_fraction)

ggtree(phy, aes(color = netdiv))
```
