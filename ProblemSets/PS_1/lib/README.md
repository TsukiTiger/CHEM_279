# `/lib` Directory

Welcome to the `lib` directory!

## Overview

The `lib` directory is designated for storing object files and, potentially, static or shared libraries generated during the compilation process. Object files are intermediate files produced when source files are compiled. They contain machine code, but aren't executable on their own. Instead, they are linked together to create the final executable or library.

## Structure

Currently, this directory is empty. However, as you compile the project, it will be populated with object files like:

- `PS_1.o`: The object file resulting from the compilation of `PS_1.cpp`.

## Usage

1. **Compilation**: When you compile a source file (e.g., with `g++` or another compiler), you can specify an output directory for the object file. Ensure you direct the compiler to place the object file in this `lib` directory.

2. **Linking**: Once all the source files are compiled, the resulting object files in the `lib` directory can be linked together to produce the final executable. This executable will typically reside in the `bin` directory.

3. **Clean Builds**: Sometimes, it's useful to remove all object files and start a fresh build. In such cases, you might delete or clean the contents of the `lib` directory.

## Best Practices

- **Separation of Concerns**: Keeping object files in the `lib` directory ensures a clean separation between source code, intermediate files, and executables. This makes the build process more organized and easier to manage.

- **Version Control**: Typically, object files (like `.o`) are not stored in version control systems like Git. Ensure you have appropriate `.gitignore` rules set up to exclude these files from your repository.

- **Consistency**: Always ensure consistent naming conventions for object files, typically reflecting the name of the source file they originate from.

## Notes

- If, in the future, you decide to create static or shared libraries (`.a` or `.so` files on Unix-like systems), they can also be stored in this directory.

- It's good practice to periodically check this directory, especially when troubleshooting build or linking issues.

For further information or queries about the build process or organization of the project, refer to the main project documentation or consult with the project maintainers.