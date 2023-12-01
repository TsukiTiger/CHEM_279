# `/lib` Directory

Welcome to the `/lib` directory for the CS279 Problem Set 4!

## Overview

The `/lib` directory plays an essential role in the project's structure, dedicated to storing the object files and library files generated during the compilation process. Object files, typically with `.o` extensions, contain compiled code but are not independently executable. They are instrumental in the build phase, as they are combined to create the final executable or library files for the project.

## Structure

- Initially, this directory is empty.
- Gets populated with object files after compiling the source files from the `/src` directory.
- Common object files include:
    - `molecule.o`: Compiled from `molecule.cpp`.
    - `AO.o`: Compiled from `AO.cpp`.
    - `CNDO.o`: Compiled from `CNDO.cpp`.

## Usage

- **Compilation**: During compilation, compilers like `g++` output object files into this directory, maintaining an organized structure separating compiled outputs from source files and executables.
- **Linking**: The object files are linked in the final step of the compilation to form the project's executable file, usually found in the `/bin` directory.
- **Clean Builds**: For fresh builds or troubleshooting, contents of this directory can be cleared, a common practice during complete rebuilds of the project or when resolving compilation/linking issues.

## Best Practices

- **Organization**: Keeping object files in the `lib` directory aids in distinguishing between source files, intermediate files, and executables, streamlining the build process.
- **Version Control**: Object files should typically not be tracked in version control systems like Git. Ensure `.gitignore` excludes `.o` files.
- **Naming Conventions**: Consistently name object files to reflect their source files, enhancing clarity and maintainability.

## Additional Notes

- The directory can also host static or shared libraries (`.a` or `.so` files) if the project includes such components.
- Regular review of this directory's contents can be helpful, especially when diagnosing build or linking issues.

For more information or inquiries related to the build process or overall project organization, please consult the main project documentation or contact the project maintainers.