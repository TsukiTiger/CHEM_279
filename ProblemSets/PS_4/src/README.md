# Source Directory: `/src`

Welcome to the `/src` directory for CS279 Problem Set 4!

## Overview

This directory hosts the core source files that are integral to the computational tasks and algorithms presented in Problem Set 4.

## Directory Contents

- **`CNDO.cpp`**: Contains the implementation of the Complete Neglect of Differential Overlap (CNDO) method, which is central to the electronic structure calculations in this problem set.

- **`molecule.cpp`**: Manages the molecular data structures and provides functionalities for molecular property calculations and manipulations.

- **`makefile`**: A script to automate the compilation of these source files, ensuring smooth integration and building of the project components.

- **`README.md`**: This file, providing a comprehensive guide and information about the source files and their roles in the project.

## Building the Source Code

To compile the source code and generate the corresponding object files, follow these steps:

1. Navigate to the `/src` directory in the terminal.
2. Run the command:

   ```bash
   make all
   ```

   This will compile all the source files and place the object files in the `/lib` directory.

## Usage

The source files, particularly `CNDO.cpp`, and `molecule.cpp`, encapsulate the core functionalities required for the computational aspects of the problem set. Include the relevant header files from the `/inc` directory in your main application or test scripts to utilize these functionalities.

## Best Practices

- Ensure all dependencies, such as external libraries or tools, are properly installed and configured before compiling the source code.
- Keep the source code well-organized and documented to facilitate easy understanding and modifications.
- Regularly update and maintain the `makefile` to reflect any changes in the source files or the addition of new files.

## Contributing

Contributions to enhance or extend the functionalities of these source files are always welcome. For guidelines on contributing, please refer to the main repository's documentation. For any discussions, suggestions, or queries, feel free to contact the project maintainers.