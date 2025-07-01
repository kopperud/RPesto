## RPesto

`RPesto` is an alternative implementation of the `Pesto` software in R. It can estimate branch-specific diversification rates for large trees.

## Installation

The backend of `RPesto` is implemented in `rust`, therefore you need to have rust and cargo installed. To install this, I recommend to follow the instructions on [https://rustup.rs/](https://rustup.rs/).

Next, we will use the R-package `remotes` to install necessary R package dependencies

```R
install.packages("remotes")

library(remotes)

install_github("YuLab-SMU/treeio")
install_github("YuLab-SMU/tidytree")
```

Next we install `RPesto`

```R
install_github("kopperud/RPesto")
```

Now you can load `RPesto` and fit the model to your tree, where we have to assume some taxon sampling probability.

```R
library(RPesto)

data("primates")

sampling_fraction <- 0.635

tree <- fit_bds(primates, sampling_fraction)
```

We can also use `ggtree` to plot the results
```R
library(ggtree)

ggtree(tree, aes(color = netdiv))
```
