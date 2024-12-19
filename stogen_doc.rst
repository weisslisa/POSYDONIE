=========================================================
Ensemble simulations and Stochastic processes
=========================================================

The ``croco/STOGEN`` source directory contains modules to introduce stochastic parameterizations in CROCO (based on the original tool available at https://github.com/brankart/stogen.git).

It includes the modules to generate stochastic processes with the required statistical properties, the interface with CROCO, and the stochastic parameterizations (described below in the "STOGEN code descritpion" section).

To use the ensemble and stochastic modules in CROCO, define the following CPP keys in **cppdefs.h**:

================= ===============================================================
STOGEN            Activate the option to generate stochastic processes
ENSEMBLE          Activate the option to generate parallel ensemble simulations
================= ===============================================================

Deactivating this keys should remove all ensemble and stochastic features.

Verify that the namelist files ``namelist_sto_cfg`` and ``namelist_sto_ref`` are copied in the running directory (from the ``croco/STOGEN`` source directory).

Edit both **namelist_sto** to modify the parameters of the stochastic fields (parameters are described hereafter).



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

The code includes three random number generators (shr3, kiss32, and kiss64) and two algorithms to transform the integer pseudo-random sequence into Gaussian numbers (the polar method and the ziggurat method), all included in the modules ``storng_kiss`` and ``storng_ziggurat``. The default option is kiss32 with ziggurat, but the user is free to use any combination of them. Other random number generators can also easily be included. 



STOGEN code description
=======================

- the stochastic modules
- the interface modules between STOGEN and CROCO
- the stochastic parametereization modules in CROCO 


Stochastic modules
------------------

These modules generate maps of stochastic processes with the requested statistical properties (time and space correlation, marginal distribution). 
They are fully independant of CROCO and can follow updates of the stogen package without interference.

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


Interface with CROCO
--------------------

- **stomod** : Main stochastic module (model dependent), embedding all dynamical stochastic parameterizations:

    - initialization phase (routine ``sto_mod_init``):

      - initialization of every dynamical stochastic parameterizations (here only ``sto_template_init``),
      - initialization of the structure of the stochastic arrays (call to ``sto_array_init``),
      - initialization of the time iteration of the stochastic arrays (call to ``sto_par_init``);

    - time update (routine ``sto_mod``):

      - update stochastic fields (call to ``sto_par``),
      - apply dynamical stochastic parameterization (call to ``sto_template``).

    The routines may need to be organized differently depending on where the stochastic parameterization code must be used in CROCO.

- **stoexternal** : This module is used by the stochastic code to collect all information required from CROCO: type of variables, description of the model grid,ensemble parameters, lateral boundary conditions (or connection between subdomains). This is the only place where model data go to the stochastic modules, so that this can be easily identified and possibly upgraded. This is model dependent.

- **storst** : Read and write from restart file.

- **Intervention in the CROCO code** : 

   - definition of the CPP keys STOGEN and ENSEMBLE in ``cppdefs.h``
   - initialization of the stochastic modules in ``main.F``
   - update the stochastic processes in ``step.F``
   - use of the stochastic arrays if and where needed 


Stochastic parametereizations in CROCO 
--------------------------------------

There must be an additional module for each type of paramteterization introduced in CROCO.

At this stage, we have introduced only one option:

- **stobulk**: Stochastic parameterization of the bulk formulation for the air-sea fluxes. At this stage, we have only introduced simple perturbations of the drag coefficient (multiplicative noise). 

**Time correlation**

.. code-block:: text

   stofields(jsto)%type_t='arn'  ! Use autorgeressive processes
   stofields(jsto)%nar_order=2   ! of order 2
   stofields(jsto)%corr_t=10.0   ! with a correlation time scale of 5 time steps
   stofields(jsto)%nar_update=5  ! with update every 5 time steps (and interpolation inbetween)

Other possibilities for ``type_t`` are ``white`` (default) and `constant`.
The cost is proportional to ``nar_order`` (default=1). The cost decreases with ``nar_update`` (default=1), especially if the scheme to generate the spatially correlated noise is expensive.
``nar_update`` should not be larger than ``corr_t``.

**Space correlation**

With diffusive operator

.. code-block:: text

   stofields(jsto)%type_xy='diffusive'  ! Use diffusive operator to obtain space correlation
   stofields(jsto)%diff_passes=50       ! number of passes of the diffusion operator
   stofields(jsto)%diff_type=1          ! type of diffusion operator (1=laplacian with mask)

The cost is proportional to ``diff_passes``, and may become quite large.
``diff_type=0`` : Laplcian diffusion without accounting for the land mask (default).
``diff_type=1`` : the land mask boundary conditions are applied to the diffusion operator.

With kernel convolution

.. code-block:: text

   stofields(jsto)%type_xy='kernel'     ! Use kernel convolution to obtain space correlation
   stofields(jsto)%corr_xy=10.0         ! correlation length scale (in degrees, in spherical coordinates)
   stofields(jsto)%ker_type=0           ! type of kernel (0=gaussian, default)
   stofields(jsto)%ker_coord=2          ! type of grid coordinates (2=spherical coordinates)

The type of kernel (`ker_type`) can be: Gaussian (0), Laplacian (1), Box kernel (2), Triangle kernel (3), Mexican hat wavelet (4), or Morlet wavelet (5).
The type of grid coordinates can be: grid (0, default), Cartesian (1), spherical (2). Options 0 and 1 are much less expensive, but, if used on the sphere, the behaviour at the poles becomes grossly inconsistent and the periodicity conditions are lost.

**Marginal distribution**

.. code-block:: text

   stofields(jsto)%type_variate='lognormal'   ! Use a lognormal marginal distribution
   stofields(jsto)%ave=1.0                    ! with mean equal to 1.0
   stofields(jsto)%std=0.5                    ! with standard deviation equal to 0.5

Or

.. code-block:: text

   stofields(jsto%type_variate='bounded_atan' ! Use a bounded distribution with arctangent transformation
   stofields(jsto)%min=0.                     ! with the lower bound equal to 0.
   stofields(jsto)%max=1.                     ! with the upper lower bound equal to 1.
   stofields(jsto)%ave=0.7                    ! with a transformed value of the normal mean equal to 0.7
   stofields(jsto)%std=0.3                    ! with a measure of the spread around the mean equal to 0.3

Currently available options are:

- a bounded normal distribution: all values outside the bound are reset to the bound, with ``type_variate='normal'`` and two different values for ``min`` and ``max``;
- a wrapped normal distribution, which is useful for cyclic variables (like angles), with ``type_variate='wrapped_normal'`` and two different values for ``min`` and ``max``;
- a lognormal distribution, which is useful where positive random numbers are requested with ``type_variate='lognormal'`` and explicitly specified mean and standard deviation;
- a bounded distribution in a given interval, which is here obtained by a transformation by the arctangent function, with ``type_variate='bounded_atan'``.



ENSEMBLE code description
=======================

- **ensmpi** : Run ensemble simulations, i.e. several members in parallel using MPI

This module requires to create as many croco.in files as there are members in the desired ensemble, named ``001croco.in``, ``002croco.in`` ... ``0nncroco.in``.

In the ``*croco.in`` file, output names should also be updated such as for example ``001config-name_rst.nc``, ``001config-name_his.nc``, ``001config-name_avg.nc``.



Parameters of the stochastic fields
===================================

Edit both ``namelist_sto`` to modify the parameters of the stochastic fields and ensemble generation:

**general parameters for ensemble simulations**

.. code-block:: text

   nn_ens_size = 1         ! ensemble size (number of members to be run in parallel)
   nn_ens_start = 1        ! index of the first ensemble member
   ln_ens_rst_in = .FALSE. ! use ensemble (T) or single (F) input restart file

**general parameters for stochastic parameterization**

.. code-block:: text

   ln_rststo   = .false.                ! start from mean parameter (F) or from restart file (T)
   ln_rstseed  = .false.                ! read seed of RNG from restart file
   cn_storst_in  = "restart_sto_in.nc"  ! input restart file for stochastic parameters
   cn_storst_out = "restart_sto_out.nc" ! output restart file for stochastic parameters

**parameters for stochastic bulk formulation** (alternatively parameters can also be modify in the **stobulk** module)

.. code-block:: text

   std  = 0.1    ! standard deviation of the multiplicative noise
   tcor = 10.0   ! time correlation (in days)
   npasses = 100 ! number of passes of the horizontal Laplacian filter
   arorder = 2   ! order of the autorgeressive processes
   nupdate = 6   ! udpate frequency of the autorgeressive processes (in time steps)
