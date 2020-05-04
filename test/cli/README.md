# CLI Test

Here are test files for command line interface tests, i.e. the app is executed with defined input files and parameters.
The test then validates whether the output is correct.

Each test fixture should be inherited from the `cli_test` class: It provides the functionality of executing the app,
finding the input files, capturing the output and creating individual test directories.
The test output files are stored in the directory `test/output`

Attention: The default `make` target does not build tests.
Please invoke the build with `make cli_test` or use `make test` to build and run all kinds of tests.
