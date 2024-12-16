Ensemble simulations and Stochastic processes
=========================================================

To use the ensemble and stochastic modules in CROCO, define the following CPP keys in **cppdefs.h**:

================= ===============================================================
STOGEN            Activate the option to generate stochastic processes
ENSEMBLE          Activate the option to generate parallel ensemble simulations
================= ===============================================================

STOGEN presentation
===================

This tool generate maps of stochastic processes to pertubate various model parameters or variables:

- with various space and time correlation structures,
- with various marginal probability distributions (at a given time and location).

Space correlation
-----------------

Spatial correlation can be obtained by any of the following algorithms:

- the recursive application of a diffusion operator on a white noise (more efficient for short-range correlation scales), in module ``stodiff``,
- the convolution of a white noise with a filtering kernel, with a Quasi Monte Carlo approximation of the convolution integral (more efficient for long-range correlation scales), in module ``stokernel``.

Other algorithms can easily be introduced in the code as a new possible option (in routine ``sto_new_2d`` in module ``stopar``).

Time correlation
----------------

Time correlation is obtained by autoregressive processes (in module ``stopar``):

- a new map of spatially correlated noise is used to feed the autoregressive equation at each time step,
- different time correlation spectra can be produced by using autoregressive processes of various orders (the higher the order, the smoother the process),
- for long correlation time scales (as compared to time step) the process can be only updated every n time steps, and then interpolated in time in between. This can save a substantial computation time if the generation of the spatially correlated random noise is expensive.

In the current version, only a linear time interpolation is available, but the code allows for storing an arbitrary number of time slices in memory to allow for higher order interpolation schemes.

Marginal distribution
---------------------

By default, the marginal distribution of the stochastic fields is a normal distribution (with zero mean and unit standard deviation). Other distributions can be obtained by applying a transformation to the reference normal processes (in module ``stomarginal``), including:

- a bounded normal distribution: all values outside the bound are reset to the bound,
- a wrapped normal distribution, which is useful for cyclic variables (like angles),
- a lognormal distribution, which is useful where positive random numbers are requested,
- a bounded distribution in a given interval, which is here obtained by a transformation by the arctangent function.

New options can easily be included.

Random number generator
-----------------------

The code includes three random number generators (shr3, kiss32, and kiss64) and two algorithms to transform the integer pseudo-random sequence into Gaussian numbers (the polar method and the ziggurat method), all included in the modules ``storng_kiss`` and ``storng_ziggurat``. The default option is kiss32 with ziggurat, but the user is free to use any combination of them. Other random number generators can also easily be included. A tool to compare their computational performance is provided in module `storng_check`.

The result I obtain from this comparison is:

.. code-block:: text

   We test 3 random number generators:
     -kiss64 (period ~ 2^250, with four 64-bit seeds),
     -kiss32 (period ~ 2^123, with four 32-bit seeds),
     -shr3   (period ~ 2^64,  with one  32-bit seed)
   together with 2 methods to generate normal numbers:
     -the polar method (Marsaglia, 1964),
     -the ziggurat method (Marsaglia, 2000).

   Typical results obtained (in seconds), to generate 10^8 numbers:
   polar    (Marsaglia, 1964) + kiss64 (Marsaglia, 2009)   1.89253000000000
   polar    (Marsaglia, 1964) + kiss32 (Marsaglia, 1992)   1.50936100000000
   polar    (Marsaglia, 1964) + shr3 (simple & fast rng)   1.14710500000000
   ziggurat (Marsaglia, 2000) + kiss64 (Marsaglia, 2009)   1.26209300000000
   ziggurat (Marsaglia, 2000) + kiss32 (Marsaglia, 1992)  0.693772999999999
   ziggurat (Marsaglia, 2000) + shr3 (simple & fast rng)  0.523565000000000
   (the last 4 are obtained with real compuations in single precision)

   In view of these results, the default choice used in the code is:
   ziggurat (Marsaglia, 2000) + kiss32 (Marsaglia, 1992) which combines a long enough period with good efficiency.
   Options are included in the namelist to modify the default.



Descritpion of the STOGEN code
==============================

There are two kinds of modules:

- the external code, which mimics the time iteration of a geohysical model (like NEMO), and illustrates how to make requests to the stochastic modules;
- the stochastic modules, which receive requests from the external code and produce the stochastic processes with the requested properties.

External code
-------------

- **stogen** : Main program with empty model illustrating the use of the stochastic modules, including:

    - initialization of model parameters (grid, number of time steps, restart options),
    - initialization of stochastic code (call to ``sto_mod_init``),
    - a loop on time steps to update the stochastic fields (call to `sto_mod`) and store them in files (call to ``sto_write``).

- **stomod** : Main stochastic module (model dependent), embedding all dynamical stochastic parameterizations:

    - initialization phase (routine ``sto_mod_init``):

      - initialization of every dynamical stochastic parameterizations (here only ``sto_template_init``),
      - initialization of the structure of the stochastic arrays (call to ``sto_array_init``),
      - initialization of the time iteration of the stochastic arrays (call to ``sto_par_init``);

    - time update (routine ``sto_mod``):

      - update stochastic fields (call to ``sto_par``),
      - apply dynamical stochastic parameterization (call to ``sto_template``).

    The routines may need to be organized differently depending on where the stochastic parameterization code must be used in the geohysical model.
    
- **stotemplate** : Template for including a new dynamical stochastic parameterization in the geohysical model. This illustrates how to make requests for stochastic fields with user-defined fetaures and how to use the resulting stochastic fields in the model.

    - initialization phase (routine ``sto_template_init``):

      - request index for a new stochastic field (call to ``sto_array_request_new``),
      - define the features of the stochastic field with the corresponding index (by filling parameters like ``stofields(index)%type_xy`` specifying the requested type of xy correlation strcuture),

    - time update (routine ``sto_template``):

      - make use of the stochastic field in the model (``stofields(index)%sto2d``, ``stofields(index)%sto3d``, or ``stofields(index)%sto0d``, depending on the requested dimension of the stochastic field, stored in ``stofields(index)%dim``).

- **stowrite** : Write the resulting stochastic fields in a NetCDF file.

- **stoexternal** : This module is used by the stochastic modules below to get all information they need from the geohysical model: type of variables, description of the model grid,ensemble parameters, lateral boundary conditions (or connection between subdomains). This is the only place where model data go to the stochastic modules, so that this can be easily identified and possibly upgraded. This is model dependent.

Stochastic modules
------------------

- **stoarray** : This is the data module, where all stochastic fields are defined and stored:

    - it receives the requests from the users (with routines ``sto_array_request_size`` and ``sto_array_request_new``),
    - it allocates the required arrays according to requests, and check the consistency of the requested features (with routine ``sto_array_init``).

- **stopar** : This is the time evolution module, where all stochastic fields are updated at each time step:

    - initialization phase (routine ``sto_par_init``):

      - seed random number generator (according to subdomain and member indices),,
      - initialize methods to generate spatially correlated random fields (calls to ``sto_diff_init``, ``sto_kernel_init``,...),
      - initialize transformations to requested marginal distributions (call to ``sto_marginal_init``),
      - initialize parameters of autoregressive processes,
      - initialize random fields (from restart or from the requested method to generate spatially correlated random fields: ``sto_diff``, ``sto_kernel``,...);

    - time update (routine ``sto_par``):

      - forward the autoregressive process in time  (or interpolate between a past and future state of the autoregressive process),
      - perform the transformation to the requested marginal distribution.

- **stowhite** : Generate a map of Gaussian white noise, with zero mean and unit standard deviation.

- **stodiff** : Generate a map of spatiallye correlated noise, with zero mean and unit standard deviation, using the recursive application of a diffusion operator on a white noise.

- **stokernel** : Generate a map of spatiallye correlated noise, with zero mean and unit standard deviation, using the convolution of a white noise with a filtering kernel.

    The convolution integral is computed using a Quasi Monte Carlo approximation, by summing over a limited number of kernel locations.

    The Quasi Monte Carlo sequence of kernel locations is obtained from a 2D random Sobol sequence (with module ``stosobolseq``).

    Options for the filtering kernel include: Gaussian kernel, Laplacian kernel, box kernel, triangle kernel, Mexican hat wavelet (Ricker wavelet), Morlet wavelet a (with specific choice of frequency, adjust if needed).

    Options for computing distances include: grid coordinates, Cartesian coordinates, spherical coordinates (more expensive).

- **stosobolseq** : Module to generate mutlidimensional Sobol sequences (obtained from https://github.com/DaanVanVugt).

- **stomarginal** : Transform the Gaussian process to the requested marginal distribution.

- **storng_kiss** : Random number generator. This includes the kiss32 and kiss64 random number generators and the polar method to transform the integer sequence into Gaussian numbers.

- **storng_ziggurat** : Random number generator. This includes the shr3 random number generator and the ziggurat method to transform the integer sequence into Gaussian numbers.

- **storng_check** : Check relative performance of random number generators.

- **storst** : Read and write from restart file (not yet implemented).



