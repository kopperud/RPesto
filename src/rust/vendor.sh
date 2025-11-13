#!/bin/sh -e

cargo vendor

# c.f. https://reproducible-builds.org/docs/archives/
tar \
  --owner=0 \
  --group=0 \
  --numeric-owner \
  --xz \
  --create \
  --file=vendor.tar.xz \
  vendor

echo
echo
echo "#############################################"
echo "#                                           #"
echo "#  UPDATING src/rust/vendor-config.toml !!!  #"
echo "#                                           #"
echo "#############################################"

echo  "[source.crates-io]" > vendor-config.toml
echo  "replace-with = \"vendored-sources\"" >> vendor-config.toml
echo  "" >> vendor-config.toml
echo  "[source.vendored-sources]" >> vendor-config.toml
echo  "directory = \"vendor\"" >> vendor-config.toml
