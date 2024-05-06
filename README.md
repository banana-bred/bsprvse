# BSPRVSE

Use a basis of **B**-**sp**lines to solve the (**r**o)**v**ibrational **S**chrödinger **e**quation,
$$H_{j}(R) \psi_\nu(R) \equiv \left[ -\frac{1}{2  μ}  \frac{d^2}{dR^2} + V(R) + \frac{j(j+1)}{2  μ  R^2} \right] \psi_\nu(R) = E_\nu \psi_\nu(R),$$
where $R$ is the internuclear distance, $\mu$ is the reduced mass of the molecule, $j$ is the rotational quantum number of the molecule, and $\nu$ is the vibrational quantum number of the molecule.
The potential $V(R)$ can be purely real, or can have a purely imaginary potential added to the tail end, resulting in a complex absorbing potential (CAP).
Using a purely real potential will result in calculating real-valued bound-state wavefunctions with real energies $E_\nu$, although the output variables are actually complex.
Using a CAP, the wavefunctions and energies become complex (the Hamiltonian is no longer Hermitian, so the eigenenergies no longer have to be real), although the complex part of $E_\nu$ for bound states should be negligible.
The available imaginary potentials are as follows :

<center>

| type         | form                                        |
| -------------|---------------------------------------------|
| exponential  | $V_\text{CAP}(R) = V(R) - i A N e^{-2 / x}$ |
| linear       | $V_\text{CAP}(R) = V(R) - i A x$            |
| quadratic    | $V_\text{CAP}(R) = V(R) - i3Ax^2/2$         |
| cubic        | $V_\text{CAP}(R) = V(R) - i 2A x^3$         |
| quartic      | $V_\text{CAP}(R) = V(R) - i 5A x^4/2$       |

</center>

where $x = (R - R_0) / L$, $L$ is the length of the imaginary potential, and $R_0$ is the starting value of the imaginary potential.
The imaginary potential starts at $R_0$ and increases monotonically until the largest available value of $R$ given by the internuclear potential.
See [[1]](#1) for more details on absorbing potentials as they're implemented by this code.

## Dependencies
- A working installation of `LAPACK`
- [*optional*] The [Fortran Package Manager (fpm)](https://github.com/fortran-lang/fpm)

## Building with [fpm](https://github.com/fortran-lang/fpm)
In the package directory, just run

    $ fpm build --profile release

The archive file `libbsprvse.a` and `.mod` files will be placed in the generated `build` subdirectory.
These files are necessary to compile another program that uses this library.

## Building without [fpm](https://github.com/fortran-lang/fpm)
In the package directory, just run the compile script

    $ ./compile

The archive file `libbsprvse.a` will be placed in the `lib` subdirectory and is necessary to link against when compiling another
program that uses this library. The `bsprvse` executable will be placed in the `bin` subdirectory.

## Testing
Test suite not yet implemented...

## Using `bsprvse`

### The `bsprvse` executable

Execution of the `program` by calling the `bsprvse` executable is controlled by `namelist` variables read via `stdin`, which can be supplied via `stdin` redirection, i.e. :

    $ ./bin/bsprvse < input/TEMPLATE.namelist

and with [fpm](https://github.com/fortran-lang/fpm),

    $ fpm run < input/TEMPLATE.namelist

The input parameters, defined in the `$input_parameters` namelist, are the following :

```
$input_parameters

j = 0
  !! The rotational quantum number $j$

nwf = 65
  !! The number of wavefunctions to calculate, from ν = 0 to ν = nwf - 1
nR_wf = 1000
  !! The number of $R$-values for the wavefunctions
reduced_mass = 7.3545994295895865
  !! The reduced mass of the molecule, in atomic mass units, i.e., the mass of
  !! carbon-12 is 12 atomic mass units.

np = 400
  !! The number of B-spline intervals
order = 5
  !! The order of the B-splines
legpoints = 10
  !! The number of Gauss-Legendre quadrature points used to calcuate integrals

R_max = 10.0
  !! The largest value to consider for the internuclear potential

CAP_exists = .true.
  !! Use a comlex absorbing potential (CAP) ?
  !! yes -> CAP_exists = .true.
  !! no  -> CAP_exists = .false.
CAP_type = 0
  !! The CAP type:
  !!   0: exponential
  !!   1: linear
  !!   2: quadratic
  !!   3: cubic
  !!   4: quartic
CAP_length = 4.0
  !! The length of the CAP, in atomic units. The CAP starts at R_max - CAP_length
  !! and continues to R_max, where it takes its maximal value. R_max is the
  !! largest value at which the internuclear potential was calcualted
CAP_strength = 0.01
  !! The CAP_strength $A$, in atomic units.

potential_file = "data/CFp_pot_30.dat"
  !! The location for the input internuclear potential in the format
  !!   R1 V1
  !!   R2 V2
  !!   R3 V3
  !!    .  .
  !!    .  .
  !!    .  .
  !!
  !! where R and V are both in atomic units

output_directory = "output"
  !! The directory in which to write the wavefunctions and energies.
  !! This directory must already exist


/
```

### `bsprvse` in your [fpm](https://github.com/fortran-lang/fpm) project

To use this project within your [fpm](https://github.com/fortran-lang/fpm) project, add the following to your `fpm.toml` file:

    [dependencies]
    bsprvse = { git = "https://github.com/banana-bred/bsprvse" }

<!-- or -->

<!--     [dependencies] -->
<!--     bsprvse = {"namespace" = "..."} -->

The module `bsprvse` contains the public interface `solve_RVSE()` :

which may be invoked as
```
    solve_RVSE(R_vals, V_vals, j, reduced_mass, nwf, nR_wf, R_wf, wf, wf_nrg, &
               np, legpoints, order                                           &
    )
```
or

```
    solve_RVSE(R_vals, V_vals, j, reduced_mass, nwf, nR_wf, R_wf, wf, wf_nrg, &
               np, legpoints, order,                                          &
               CAP_exists, CAP_length, CAP_type, CAP_strength                 &
    )
```

The input variables are as follows :

<center>

| variable         | type                          | explanation |
| ---------------- | ----------------------------- | ----------------------------------------------- |
| R_vals(:)        | real(dp), intent(in)          | Array containing the values of the internuclear distance in atomic units|
| V_vals(:)        | real(dp), intent(in)          | Array containing the values of the intermolecular potential in atomic units|
| j                | integer, intent(in)           |  The value $j$ in the Schrödinger equation|
| reduced_mass     | real(dp), intent(in)          | The system's reduced mass in atomic units (electron mass = 1; not atomic mass units)|
| nwf              | integer, intent(in)           | The number of wavefunctions to calculate|
| nR_wf            | integer, intent(in)           |  The number of $R$-grid points on which to evaluate the wavefunctions|
| R_wf(:)          | real(dp), intent(in)          | The $R$-grid on which wavefunctions will be calculated
| wf(:,:)          | complex(dp), intent(out)      |  Array containing the values of the wavefunctions indexed as `(iR, iv)`, where `iR` runs over the|
|                  |                               |  internuclear distances and iv runs over the vibrational quantum number ν|
| wf_nrg(:)        | complex(dp), intent(out)      | The energies of the wavefunctions (in atomic units)|
| np               | integer, intent(in)           | The number of B-spline intervals|
| legpoints        | integer, intent(in)           | The number of Gauss-Legendre quadrature points used to calculate integrals|
| order            | integer, intent(in)           | The order of the B-splines|
| CAP_exists       | logical, intent(in)           | Add an imaginary potential to the real internuclear potential ?|
|                  |                               |   yes -> `.true.` |
|                  |                               |   no  -> `.false.` |
| CAP_length       | real(dp), intent(in)          | The length of the imaginary potential in atomic units. Has no effect if `CAP_exists` is `false`|
| CAP_strength     | real(dp), intent(in)          | The strength ($A$) of the CAP in atomic units|
| CAP_type         | integer, intent(in)           | The type of CAP : |
|                  |                               |   0 : exponential |
|                  |                               |   1 : linear |
|                  |                               |   2 : quadratic |
|                  |                               |   3 : cubic |
|                  |                               |   4 : quartic |

</center>

where `dp` represents double precision, as defined by `real64` in the `iso_fortran_env` intrinsic module.
The former call (without the CAP variables) is essentially the same as calling the latter call where `CAP_exists` is `.false.`.


## Reference(s)

<a id="1">[1]</a>
Á. Vibók and G. G. Balint-Kurti
*Parameterization of Complex Absorbing Potentials for Time-Dependent Quantum Dynamics*,
J. Phys. Chem. 1992, 96, 8712-8719
URL: [https://pubs.acs.org/doi/pdf/10.1021/j100201a012](https://pubs.acs.org/doi/pdf/10.1021/j100201a012)
