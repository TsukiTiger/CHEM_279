# `/Inc` Directory

Welcome to the `inc` directory!

## Overview

The `inc` directory is reserved for storing header files for the project. Header files typically contain declarations of functions, classes, structures, and other types that can be used across multiple source files. They are essential for facilitating modular code design.

## Structure

Currently, this directory contains the following header files:

- `PS_1.h`: This header file pertains to Problem Set 1. It includes necessary declarations that are used in its associated source files.

## Usage

1. **Including Headers**: When you need to use any functionality declared in a header file within a source file, include it using the `#include` directive. For instance:
   ```cpp
   #include "PS_1.h"
   ```

2. **Extending the Project**: As the project grows, you might need to add more header files. Always place them in this `Inc` directory to maintain a clean and organized project structure.

3. **Consistency**: Ensure consistency in naming conventions, coding style, and comments within header files to make the codebase maintainable and readable.

## Best Practices

- **Header Guards**: Ensure that each header file has header guards to prevent double inclusion. This is typically done using `#ifndef`, `#define`, and `#endif` preprocessor directives.

- **Documentation**: Always document your declarations in the header files. This helps other developers understand the purpose and usage of functions, classes, or other types you've declared.

- **Dependencies**: Minimize the dependencies within header files. This reduces the chances of circular dependencies and compilation issues.

## Notes

- Remember that header files are not typically compiled, but they are included in source files that are. Ensure that your header files only contain declarations and not definitions (with some exceptions like `inline` functions or template definitions).

- It's a good practice to periodically review the header files to ensure they are up-to-date and only contain relevant information.

If you have further queries or face any issues related to the header files, please refer to the main project documentation or reach out to the project maintainers.