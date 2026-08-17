<!--
SPDX-FileCopyrightText: 2022-2026 Helmholtz-Zentrum hereon GmbH
SPDX-License-Identifier: CC0-1.0
-->


# FABM-hereon

This is a collection of [FABM](https://fabm.net) models developed at Helmholtz-Zentrum hereon GmbH by the ecosystem modeling group.  It currently includes
1. A `light` model, based on the GOTM `light` implementation but with more exact layer averaging
2. The **OMExDia** model in several variations. OMExDia is a sediment biogeochemical model
  *  OMExDia with phosphorus (`omexdia_p`)
  *  OMExDia with methane (`omexdia_c`)
1. The **NOPE** model for Nitrous Oxide Production and Emission
  
## How to build

This code must be compiled together with FABM. To do this, provide the following additional arguments to cmake [when you build FABM](https://github.com/fabm-model/fabm/wiki/Building-and-installing): `-DFABM_INSTITUTES=hereon -DFABM_HEREON_BASE=$FABM_HEREON_BASE`

Here, `$FABM_HEREON_BASE` is an environment variable pointing to the directory with the FABM-hereon code (the `src` subdirectory as seen from the directory containing this `ReadMe` file). Note that `-DFABM_INSTITUTES=hereon` will make FABM compile NOPE, OMExDia and light as the *only* available biogeochemical models. If you additionally want to have access to other biogeochemical models included with FABM, you can set `FABM_INSTITUTES` to a semi-colon separated list, e.g., `-DFABM_INSTITUTES="hereon;iow"` (to prevent the shell from interpreting the semi-colons, you typically have to enclose this list with quotes).

### Building with GOTM
To use hereon models with the General Ocean Turbulence Model (GOTM), do the following, after defining GOTM_BASE and FABM_HEREON_BASE

```
git clone --recursive https://github.com/gotm-model/code.git $GOTM_BASE
git clone https://github.com/schism-dev/fabm-hereon.git $FABM_HEREON_BASE

cd $GOTM_BASE && git submodule update --init --recursive
mkdir -p $GOTM_BASE/build
cd $GOTM_BASE/build 
cmake -B $GOTM_BASE/build -S $GOTM_BASE  -DGOTM_USE_FABM=ON -DFABM_INSTITUTES=hereon -DFABM_HEREON_BASE=$FABM_HEREON_BASE
make
```

This will create the `gotm` executable with support for FABM light and omexdia.

### Building with SCHISM
To use hereon models with the latest stable release of the Semi-implicit Cross-scale Hydroscience Integrated System Model (SCHISM), do the following:

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

## Building the 0D driver

To use hereon models in a 0d setup, do the following, after defining FABM_BASE, GOTM_BASE and FABM_HEREON_BASE.

```
mkdir $FABM_BASE/build-0d
cmake -B $FABM_BASE/build-0d -S $FABM_BASE/src/drivers/0d -DGOTM_BASE=$GOTM_BASE -DFABM_INSTITUTES=hereon -DFABM_HEREON_BASE=$FABM_HEREON_BASE
make
```


### Building pyfabm

To use hereon models with pyfabm, do the following, after defining FABM_BASE and FABM_HEREON_BASE. See also the detailed instructions at https://github.com/fabm-model/fabm/wiki/python.

```
cat <<EOT > $FABM_BASE/setup.cfg
[build_ext]
cmake_opts=-DFABM_EXTRA_INSTITUTES=hereon -DFABM_HEREON_BASE="$FABM_HEREON_BASE"
force=1
debug=1
EOT

python -m pip install $FABM_BASE/setup.cfg
```

## How to couple

FABM provides several built-in models that you can use to provide source and sink terms, as well as boundary conditions. The following models are available

* `interior_constant`, `horizontal_constant`
* `surface_flux`, `constant_surface_flux`, `external_surface_flux`
* `external_bottom_flux`
* `interior_source`
* `bottom_source`
* `interior_relaxation`
* `column_projection`
* `weighted_sum`
* `horizontal_weighted_sum`
* `bottom_field`
 
 ### Interior sources

 You can convert the flux from one model `source` to another model `target`  by specifying an `interior_source` model:

 ```
instances:
  interior_source:
    model: interior_source
    coupling:
      source: source/flux
      target: target/c
```

### Horizontal constant

You can specify a new horizontally constant property, that can be used in the `coupling:` section of onother model as `horizontal_constant/data`:

```
instances:
  horizontal_constant:
    model: horizontal_constant
    parameters:
      value: 0.0001                   # value

```
