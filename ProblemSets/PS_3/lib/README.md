# `/lib` Directory

Welcome to the `lib` directory of our project!

## Overview

The `lib` directory is a key component of the project's structure, designated for storing object files and potentially static or shared libraries. These are created during the compilation of source files. Object files, denoted with the `.o` extension, contain machine code but aren't executable on their own. They play a crucial role in the building phase of the project, as they are linked together to create the final executable or libraries.

## Structure

Initially, this directory might appear empty. It gets populated with object files following the project compilation. Common examples include:

- `molecule.o`: Resulting from the compilation of `molecule.cpp`.
- `eht_matrices.o`: Generated by compiling `eht_matrices.cpp`.

## Usage

### Compilation
When using a compiler such as `g++`, the output object files should be directed to this `lib` directory. This keeps the compiled outputs organized and separate from the source files and executables.

### Linking
Object files in this directory are linked in the final step of the compilation to form the executable file, typically located in the `bin` directory.

### Clean Builds
For clean builds or troubleshooting issues, the contents of this directory (object files) can be removed. It's a common practice during full project rebuilds or resolving certain compilation and linking errors.

## Best Practices

- **Organizational Clarity**: Keeping object files in the `lib` directory helps in distinguishing between different types of files in the project (source code, intermediate files, executables) and enhances the manageability of the build process.

- **Version Control**: It is advisable not to track object files in version control systems (e.g., Git). Ensure your `.gitignore` file excludes these `.o` files.

- **Naming Conventions**: Consistently name object files to reflect the source files they are compiled from, aiding in clarity and maintainability.

## Additional Notes

- This directory can also be used to store static or shared libraries (such as `.a` or `.so` files in Unix-like systems), if your project expands to include them.

- Regularly reviewing the contents of this directory can be beneficial, particularly when diagnosing build or linking issues.

For more details or queries related to the build process or project organization, please refer to the main project documentation or reach out to the project maintainers.
