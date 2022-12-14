
<!-- README.md is generated from README.Rmd. Please edit that file -->

# native-wavefunction

<!-- badges: start -->
<!-- badges: end -->

The **native-wavefunction** project is a collection of functions that
calculate and plot the hydrogen wave function using only native R. All
functions were implemented by myself using the references below except
for the associated Laguerre polynomials which use the
[mpoly](https://github.com/dkahle/mpoly) implementation. A Monte Carlo
simulation is then used to choose points in space where an electron
would be likely to exist. These points can then plotted using
[plotly](https://github.com/plotly/plotly.R) to create pointillist
representations of the hydrogen orbitals.

The functions in **native-wavefunction** are by no means fast. This
project was a proof of concept and takes a significant amount of time to
simulate more complex orbitals. An improved wave function implementation
can be found in my
[**wavefunction**](https://github.com/hmlea/wavefunction) package. The
**wavefunction** package uses the [Boost
implementation](https://www.boost.org/doc/libs/1_80_0/libs/math/doc/html/math_toolkit/sf_poly/sph_harm.html)
of the spherical harmonics and is significantly faster.

This project was greatly inspired by Alvin Q. Meng and his own project,
[Evanescence](https://al2me6.github.io/evanescence)
([GitHub](https://github.com/al2me6/evanescence)). I would not have been
able to create this project without the insight into his project that he
enthusiastically shared with me. Please check out his project at the
link above.

## Usage

To use this project, download it and source `monte_carlo.R`; this will
provide access to all the main functions below:

### Simulating Orbitals

    sim(n, l, m, max_psi=0.0026, num_points=2048, r_lim=35,
        pos_col="#d13010", neg_col="#2d709f", bg_col="gray",
        benchmark=TRUE, plot_orb=TRUE)

-   `n`   The principle quantum number, n
-   `l`   The azimuthal quantum number, l
-   `m`   The magnetic quantum number, m
-   `max_psi`   The maximum value of $\Psi$ for the orbital
-   `num_points`   The number of points to simulate for the orbital
-   `pos_col`   The color of points with a positive $\Psi$ value
-   `neg_col`   The color of points with a negative $\Psi$ value
-   `bg_col`   The color of the background of the plot
-   `benchmark`   Benchmark the simulation time and print it to the
    console
-   `plot_orb`   Plot the orbital; if false, a data frame containing
    coordinates and $\Psi$ will be returned

### Plotting Orbitals

    plot_data(coords, pos_col="#d13010", neg_col="#2d709f", bg_col="gray")

-   `coords`   A data frame containing coordinates and $\Psi$ from the
    `sim` function with `plot_orb = FALSE`
-   `pos_col`   The color of points with a positive $\Psi$ value
-   `neg_col`   The color of points with a negative $\Psi$ value
-   `bg_col`   The color of the background of the plot

### Calculating $\Psi$

    psi(n, l, m, r, theta, phi)

-   `n`   The principle quantum number, n
-   `l`   The azimuthal quantum number, l
-   `m`   The magnetic quantum number, m
-   `r`   The radius in spherical coordinates
-   `theta`   The angle $\Theta$ in spherical coordinates
-   `phi`   The angle $\Phi$ in spherical coordinates

### Calculating the Radial Component

    radial_comp(n, l, r)

-   `n`   The principle quantum number, n
-   `l`   The azimuthal quantum number, l
-   `r`   The radius in spherical coordinates

### Calculating the Real Spherical Harmonics

    real_sph_harm(l, m, theta, phi)

-   `l`   The azimuthal quantum number, l
-   `m`   The magnetic quantum number, m
-   `r`   The radius in spherical coordinates
-   `phi`   The angle $\Phi$ in spherical coordinates

### Testing

    source("tests/tests.R")

-   Runs a test file that compares the values of the
    **native-wavefunction** project to the **wavefunction** package

Tests can also be manually created using:

    source("tests/core.R")
    t = gen_tests(num_tests, seed, max_n, max_r)

-   `num_tests`   The number of tests to create; the default is 1000
-   `seed`   The seed; the default is the current UNIX timestamp
-   `max_n`   The maximum principle quantum number, n, for the tests
-   `max_r`   The maximum radius, r, of the spherical coordinates for
    the tests

Tests can be summarized using:

    tol = summarize_tests(tests, tolerance, omit_index)

-   `tests`   A test to summarize
-   `tolerance`   The maximum relative error for values to be considered
    equal; the default is 10<sup>-7</sup>
-   `omit_index`   Omit the index of values that aren’t equal; the
    default is `TRUE`

`tol` then stores the minimum tolerance for all values to be considered
equal. This can be confirmed by using `tol` with `summarize_tests()`.

## Examples

Simulate and plot the the 2p<sub>x</sub> orbital using electron density:

    sim(2, 1, -1)

Simulate the electron density of the 2s orbital; then plot the orbital
on a white background and color positive values of $\Psi$ green and
negative values of $\Psi$ orange:

    coordinates = sim(2, 0, 0, plot_orb=FALSE)
    plot_data(coordinates, pos_col="green", neg_col="orange", bg_col="white")

Simulate and store a plot of the 1s orbital with 4000 points on a cream
colored background without benchmarking; then display the plot:

    fig = sim(1, 0, 0, num_points=4000, bg_col="#fffdd0", benchmark=FALSE)
    fig

Calculate and store $\Psi$ for an orbital with quantum numbers $n=6$,
$l=4$, and $m=-3$ at the spherical coordinates $r=0.9823$,
$\Theta=1.65$, and $\Phi=5.392$; then print the $\Psi$ and $|\Psi|^2$
values with labels:

    psi_val = psi(6, 4, -3, 0.9823, 1.65, 5.392)
    print(paste("Psi:", psi_val))
    print(paste("|Psi|^2:", abs(psi_val)**2))

Test to see how equal the **native-wavefunction** project is compared to
the **wavefunction** package using 2000 different $\Psi$ values:

    source("tests/core.R")
    test = gen_tests(2000)
    tol = summarize_tests(test)
    summarize_tests(test, tolerance=tol)

## To Do

In the future, I would like to:

-   Turn this project into a full package
-   Clean up and refactor existing code
-   Create my own implementation of the associated Laguerre polynomials

## References

1.  Meng 2022, [Evanescence](https://al2me6.github.io/evanescence)
    ([GitHub](https://github.com/al2me6/evanescence))
2.  Tully *et al.* 2013, [Interactive Web-Based Pointillist
    Visualization of Hydrogenic Orbitals Using
    Jmol](https://doi.org/10.1021/ed300393s)
3.  DeBruyne 2003, [Bra, Ket, Dirac, and all
    that…](https://faculty.washington.edu/seattle/physics441/441xxxindex.html)
4.  Lobo *et al.* 2019, [A smooth path to plot hydrogen atom via Monte
    Carlo method](https://doi.org/10.1590/1806-9126-RBEF-2019-0073)
5.  Blanco *et al.* 1997, [Evaluation of the rotation matrices in the
    basis of real spherical
    harmonics](https://doi.org/10.1016/S0166-1280(97)00185-1)
6.  Kahle 2013,
    [mpoly](https://journal.r-project.org/archive/2013-1/kahle.pdf)
    ([GitHub](https://github.com/dkahle/mpoly))
7.  Sievert 2020, [Interactive Web-Based Data Visualization with R,
    plotly, and shiny](https://plotly-r.com/)
    ([GitHub](https://github.com/plotly/plotly.R))
