# `/lib` Directory

Welcome to the `lib` directory!

## Overview

The `lib` directory is specifically designated for housing object files and, if applicable, static or shared libraries generated during the compilation process of the project. Object files are intermediate files that are produced when source files are compiled. These files contain machine code but are not executable by themselves. They are linked together to create the final executable or library.

## Structure

Initially, this directory is empty, but it will be populated with object files as you compile the project. For instance:

- `PS_2.o`: This object file is the result of the compilation of `PS_2.cpp`.

## Usage

1. **Compilation**: When compiling a source file using a compiler like `g++`, specify the output directory for the object file. Direct the compiler to place the object file in this `lib` directory.

2. **Linking**: After the compilation of all source files, the object files located in the `lib` directory can be linked together to form the final executable, usually located in the `bin` directory.

3. **Clean Builds**: Occasionally, you might find it useful to remove all object files and start a clean build. In such scenarios, the contents of the `lib` directory can be cleaned or deleted.

## Best Practices

- **Separation of Concerns**: Maintaining object files in the `lib` directory aids in separating source code, intermediate files, and executables. This organization simplifies the build process and makes it more manageable.

- **Version Control**: It’s conventional not to store object files (such as `.o` files) in version control systems like Git. Ensure your `.gitignore` is configured to exclude these files from your repository.

- **Consistency**: Maintain consistent naming conventions for object files, ideally reflecting the name of the source file from which they originate.

## Notes

- If you decide to generate static or shared libraries (like `.a` or `.so` files on Unix-like systems) in the future, this directory can also host them.

- Periodically reviewing this directory is good practice, especially when resolving build or linking issues.

For any additional details or queries regarding the build process or the project's organizational structure, please refer to the main project documentation or consult with the project maintainers.
