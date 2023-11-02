# Source Directory: `/src`

This directory contains the essential source files for the CS279 Problem Set 2.

## Files

- **`makefile`**: The build file used for compiling the source code in this directory.
- **`molecule.cpp`**: Source file containing implementations for molecular structure and computations.
- **`eht_matrices.cpp`**: Source file dedicated to calculations involving extended Hückel theory matrices.

## Building the Source Code

To compile the source code and produce the necessary object files, navigate to the `/src` directory and execute:

```bash
make all
```

This command will compile `molecule.cpp` and `eht_matrices.cpp`, generating corresponding object files in the `/lib` directory.

## Usage

The source files `molecule.cpp` and `eht_matrices.cpp` encapsulate the primary functionalities required for the problem set. Include the relevant headers from the `/inc` directory in your main application or tests to use these functionalities.

## Notes

- Before compiling, ensure all dependencies, such as the Armadillo linear algebra library, are correctly installed.
- Update the `makefile` as necessary if any modifications are made to the source files or if additional source files are added.

## Contributing

Please refer to the main repository's guidelines for modifications, suggestions, or adding new features. Direct any discussions to the project maintainers.
