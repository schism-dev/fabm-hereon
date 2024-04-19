# FABM-hereon

This is a collection of [FABM](https://fabm.net) ports of the GOTM light model with horizontal parameter variation and the OmexDia model.

## How to build

This code must be compiled together with FABM. To do this, provide the following additional arguments to cmake [when you build FABM](https://github.com/fabm-model/fabm/wiki/Building-and-installing): `-DFABM_INSTITUTES=hereon -DFABM_HEREON_BASE=</path/to/fabm-hereon>`

Here, `</path/to/fabm-hereon>` is the directory with the FABM-hereon code (the same directory that contains this ReadMe file). Note that `-DFABM_INSTITUTES=hereon` will make FABM compile OmexDia and light as the *only* available biogeochemical models. If you additionally want to have access to other biogeochemical models included with FABM, you can set `FABM_INSTITUTES` to a semi-colon separated list, e.g., `-DFABM_INSTITUTES="hereon;iow"` (to prevent the shell from interpreting the semi-colons, you typically have to enclose this list with quotes).

For instance, to use hereon models with the latest stable release of the Semi-implicit Cross-scale Hydroscience Integrated System Model (SCHISM), do the following:

```
git clone --recurse-submodules -b v5.11.1 https://github.com/schism-dev/schism.git schism
git clone https://github.com/fabm-model/fabm.git
git clone https://github.com/schism-dev/fabm-hereon.git
mkdir build
cd build
cmake ../schism/src -DBLD_STANDALONE=ON -DUSE_FABM=ON -DFABM_BASE=../fabm -DFABM_INSTITUTES=hereon -DFABM_HEREON_BASE=../fabm-hereon
make pschism
```

This will create the `pschism` executable with support for FABM light and omexdia.
