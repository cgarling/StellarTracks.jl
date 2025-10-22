# [Using StellarTracks.jl from Python](@id python)

If you are mainly a Python user and wish to use StellarTracks.jl and/or [BolometricCorrections.jl](@extref BolometricCorrections overview) *from Python*, this guide will describe how to do so. 

Julia can be accessed from Python through the [`juliacall` package](https://juliapy.github.io/PythonCall.jl/stable/juliacall/). You should also install the [`juliapkg` package](https://github.com/JuliaPy/PyJuliaPkg) to manage your Julia dependencies. In most configurations these can be easily installed with `pip install juliacall` and `pip install juliapkg`. 

!!! info
    The use of a Python virtual environment is recommended as the total installed data volume (Julia, packages, stellar track and bolometric correction files) may be large, and virtual environments can be removed easily when you are finished with them. We recommend the IPython shell as an extension will be loaded alongside `juliacall` to better display Julia objects.

If you don't already have Julia installed, `juliacall` can install the Julia runtime when its loaded for the first time. You should do this before moving on and installing packages; from Python, run

```python
from juliacall import Main as jl
```

and a Julia runtime should be installed if you do not already have Julia installed.

To install StellarTracks.jl and BolometricCorrections.jl into your `juliapkg`, environment, execute the following from Python:

```python
>>> import juliapkg
>>> juliapkg.add("StellarTracks")
>>> juliapkg.add("BolometricCorrections")
```

Now you should be able to load the packages from Python with

```python
>>> from juliacall import Main as jl
>>> jl.seval("using StellarTracks, BolometricCorrections")
```

after which the contents of the package can be accessed as

```python
>>> jl.StellarTracks.PARSEC.eep_idxs
Julia: (PMS_BEG = 1, MS_BEG = 201, MS_TMIN = 401, MS_TO = 601, RG_TIP = 1101, HE_BEG = 1131, END_CHEB = 1631, TPAGB_BEG = 1731)
>>> jl.StellarTracks.PARSEC.eep_idxs.RG_TIP
1101
>>> jl.BolometricCorrections.YBC.ybc_url
'https://gitlab.com/cycyustc/ybc_tables.git'
```

!!! warning "Multi-Threading"
    Some methods in these packages make use of multi-threading in Julia. When executing `from juliacall import Main as jl` from Python, you may see a warning like this one printed in the console: 

    Julia was started with multiple threads but multithreading support is experimental in JuliaCall. It is recommended to restart Python with the environment variable PYTHON\_JULIACALL\_HANDLE\_SIGNALS=yes set, otherwise you may experience segfaults or other crashes. Note however that this interferes with Python's own signal handling, so for example Ctrl-C will not raise KeyboardInterrupt. See https://juliapy.github.io/PythonCall.jl/stable/faq/#Is-PythonCall/JuliaCall-thread-safe? for further information. You can suppress this warning by setting PYTHON\_JULIACALL\_HANDLE\_SIGNALS=no.

    In this case, starting Python with the environment variable `PYTHON_JULIACALL_HANDLE_SIGNALS=yes` is recommended. From a terminal this can be accomplished with `PYTHON_JULIACALL_HANDLE_SIGNALS=yes python`. Without this setting the multi-threaded functions may cause segmentation faults or other errors.

Some self-enclosed examples are given below.

## Example Isochrone Generation

```python
>>> from juliacall import Main as jl
>>> import math
>>> jl.seval("using StellarTracks, BolometricCorrections")
>>> tracklib = jl.PARSECLibrary()
# We will choose to use Gaia eDR3 bolometric corrections
>>> bcg = jl.YBCGrid("gaiaEDR3")
>>> age = 1e9 # Age [yr], 
>>> MH = -1.7 # [M/H] metallicity
>>> Av = 0.02 # Reddening (A_v = 0.02 mag)
>>> iso = jl.isochrone(tracklib, bcg, math.log10(age), MH, Av)
Table with 9 columns and 1709 rows:
      eep  m_ini     logTe    Mbol     logg     C_O  G        G_BP     G_RP
    ┌─────────────────────────────────────────────────────────────────────────
 1  │ 200  0.10134   3.53186  11.8705  5.35728  0.0  12.5144  13.6954  11.4482
 2  │ 201  0.106499  3.54134  11.6584  5.33164  0.0  12.2518  13.3588  11.2163
 3  │ 202  0.107108  3.5424   11.6346  5.32877  0.0  12.2225  13.3218  11.1903
 4  │ 203  0.107746  3.54349  11.6099  5.3258   0.0  12.1921  13.2836  11.1634
 5  │ 204  0.108413  3.54461  11.5844  5.32274  0.0  12.1608  13.2444  11.1356
 6  │ 205  0.10911   3.54575  11.5582  5.31959  0.0  12.1286  13.2041  11.1069
 7  │ 206  0.109835  3.54693  11.5312  5.31636  0.0  12.0957  13.1631  11.0775
 8  │ 207  0.11059   3.54812  11.5036  5.31306  0.0  12.062   13.1213  11.0475
 9  │ 208  0.111373  3.54933  11.4755  5.30968  0.0  12.0276  13.0788  11.0168
 10 │ 209  0.112186  3.55055  11.4468  5.30625  0.0  11.9927  13.0358  10.9856
 ⋮  │  ⋮      ⋮         ⋮        ⋮        ⋮      ⋮      ⋮        ⋮        ⋮
```

The returned isochrone `iso` is a Julia `TypedTables.Table` object. The columns are Julia vectors and can be accessed as `iso.eep`, for example. The columns can be interpreted as NumPy arrays:

```python
>>> import numpy as np
>>> np.array(iso.m_ini)
array([0.10133962, 0.10649858, 0.10710792, ..., 1.85926778, 1.85927044,
       1.8592731 ], shape=(1709,))
```

# Example Track Generation

Now we look at interpolating a stellar track at a given initial stellar mass. 

```python
>>> from juliacall import Main as jl
>>> import math
>>> jl.seval("using StellarTracks, BolometricCorrections")
>>> tracklib = jl.PARSECLibrary()
# We will choose to use Gaia eDR3 bolometric corrections
>>> bcg = jl.YBCGrid("gaiaEDR3")
>>> age = 1e9 # Age [yr], 
>>> MH = -1.7 # [M/H] metallicity
>>> Av = 0.02 # Reddening (A_v = 0.02 mag)
>>> m_ini = 0.81 # Initial stellar mass ([solar masses])
>>> track = tracklib(MH, m_ini) # Interpolate a stellar track with correct properties
InterpolatedTrack with M_ini=0.81, MH=-1.6999999999999997, Z=0.0003100162055070581, Y=0.24905182884580257, X=0.7506381549486904.
```

`track` is now a stellar track that has been interpolated from the initial PARSEC grid to our desired metallicity and initial stellar mass. Calling `track(log10(age [yr]))` will now return interpolated quantities at our desired age:

```python
>>> values = track(math.log10(age)) # Interpolate the track at the desired age
(logTe = 3.7973512914337046, Mbol = 5.153244139448194, logg = 4.635189907463032, C_O = 0.0)
```

Now, to generate mock photometry in our desired Gaia eDR3 filters, we do

```python
# Interpolate the bolometric correction grid to our desired metallicity, reddening
>>> bcg_interp = bcg(MH, Av)
>>> phot = bcg_interp(values) # Call bolometric correction grid with track values to get mock photometry
3-element StaticArraysCore.SVector{3, Float32} with indices SOneTo(3):
 -0.027466586
 -0.27473205
  0.38811892
```

The Julia `SVector` type of `phot` is most akin to a tuple, but can also be converted to a NumPy array

```python
>>> tuple(phot)
(-0.027466585859656334, -0.2747320532798767, 0.3881189227104187)
>>> np.array(phot)
array([-0.02746659, -0.27473205,  0.38811892], dtype=float32)
```

These are now magnitudes corresponding to the filter names contained in `bcg_interp`:

```python
>>> [jl.String(s) for s in jl.filternames(bcg_interp)]
['G', 'G_BP', 'G_RP']
```
