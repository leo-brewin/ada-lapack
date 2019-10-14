# Ada Lapack

This is a port of a subset of [Lapack-3.4.1][1] to Ada2012.

The codes provided here are a direct translation (using [f2a.pl][2] from [Oliver
Kellogg][3]) of the Fortran source to Ada. The Ada source is written in the form of a
generic package and covers the Lapack routines for

* Matrix determinants and inverses for general matrices,
* Eigenvalues and eigenvectors of general, real and hermitian symmetric matrices,
* Solutions of systems of equations for general, real and hermitian symmetric coefficient matrices,
* Singular value decomposition for general matrices

Also included is a set of functions and procedures that provide a more familiar Ada
style interface to the Lapack routines.

This GitHub site is copy of the original [Ada Lapack][4] site on Sourceforge.

See also [AdaLAPACK][5], a binding to the Fortran LAPACK.

__Note__ that the codes provided here requires an Ada 2012 compiler. They will not
compile under Ada95 or earlier (see `RATIONALE.txt` for an explanation).

## Compiling and testing

The package and tests can be run using the supplied Makefile.
To compile the package and run the tests simply type

    make

The compilation and testing can be done separately using

    make lib
    make tests

The compiler will report warnings of the form

    ada_lapack.adb:14444:10: info: code between label and backwards goto rewritten as loop
    ada_lapack.adb:14470:10: info: code between label and backwards goto rewritten as loop

These can be safely ignored.

Each test, one for each of the Lapack routines, compares the output for that test (in
`tests/output`) against the expected output (in `tests/expected`). If the codes are
working correctly (as they should) then you should see no output from the tests. To be
doubly sure you can run `make diff` from the `tests` directory. The `.diff` files will be
non-empty only when there was a problem with that particular test.

[Simon Wright](https://github.com/simonjwright) has kindly provided a set of GNAT project files that can
be used to compile the package and any project that uses the package.

To compile the package type

    gprbuild -p -P ada_lapack.gpr

To compile one of the test codes (say tdgeev) type

    cd test
    gprbuild -p -P tests.gpr -o tdgeev tdgeev.adb

It should be a simple matter to adapt these project files to your own needs.

## License

The main Ada Lapack files `source/ada_lapack.adb` and `source/ada_lapack.ads` carry the
standard [Lapack copyright notice][6]. All other files in this collection are
distributed under the [MIT][7] license. See the file LICENSE.txt for the full details.

 [1]: http://www.netlib.org/lapack/#_lapack_version_3_4_1
 [2]: http://www.okellogg.de/for2ada95-1.4.tar.gz
 [3]: http://www.okellogg.de/x.html
 [4]: https://sourceforge.net/projects/ada-lapack/
 [5]: https://sourceforge.net/projects/adalapack/
 [6]: http://www.netlib.org/lapack/LICENSE.txt
 [7]: https://opensource.org/licenses/MIT
