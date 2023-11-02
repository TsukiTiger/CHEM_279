# `/inc` Directory

Welcome to the `inc` directory!

## Overview

The `inc` directory is a crucial component of the project, designated for storing header files. Header files are instrumental in a modular code architecture, containing declarations of functions, classes, structures, and other types. These declarations are essential for enabling code reuse and improving maintainability across various parts of the project.

## Structure

This directory currently houses the following header files:

- `eht_matrices.h`: This header file contains declarations related to the Extended Hückel Theory (EHT) matrices. It is key to the computational aspects dealing with molecular orbitals and electron interactions.

- `molecule.h`: In this header file, the `Molecule` class and its associated properties and methods are declared. This class plays a central role in representing and manipulating molecular structures within the project.

## Usage

1. **Including Headers**: To utilize the declarations in these header files within your source files (`*.cpp`), include them at the beginning of your source file using the `#include` directive, like so:
   ```cpp
   #include "eht_matrices.h"
   #include "molecule.h"
   ```

2. **Extending the Project**: As the project expands, new header files might be introduced. Always place new header files in this `inc` directory to maintain an orderly structure.

3. **Consistency and Clarity**: Maintain uniform naming conventions and a consistent coding style across all header files. Clear, consistent commenting and documentation are also essential for readability and maintainability.

## Best Practices

- **Header Guards**: Use header guards in each header file to avoid issues with double inclusion. This is commonly implemented using the `#ifndef`, `#define`, and `#endif` preprocessor directives.

- **Documentation**: Thoroughly document the purpose and usage of the declarations in your header files. Well-documented code aids comprehension and facilitates ease of use by other developers.

- **Minimize Dependencies**: Aim to reduce dependencies within your header files to avert circular dependencies and streamline the compilation process.

## Notes

- Header files should primarily consist of declarations, not definitions (excluding inline functions and templates). They are included in source files that are compiled.

- Regularly review and update these header files to ensure they remain relevant and accurate, aligning with the ongoing developments and requirements of the project.

For any additional information or assistance regarding the header files and their management, please refer to the main project documentation or contact the project maintainers.