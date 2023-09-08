# `/utils` Directory

## Overview

Welcome to the `utils` directory! This directory contains utility header files that are essential for the various functionalities of the project. These header files define classes, structures, and constants that are used throughout the codebase to maintain modularity and reusability.

## Directory Contents

- `atom.h`: This header file contains the Atom class, used for representing atomic entities in the project. The Atom class is central to several algorithms and computations.

- `vector3d.h`: This header file includes the Vector3D class, a representation of 3D vectors. This class provides essential methods for vector manipulations, like addition, subtraction, and dot product calculations.

- `constants.h`: This header file contains constants that are used throughout the project, such as physical constants or specific values that shouldn't be changed during the execution of the program.

## Usage

### How to Include Utilities

To make use of these utility header files in your source code, include them at the top of your `.cpp` files using the `#include` directive. For example:

```cpp
#include "utils/atom.h"
#include "utils/vector3d.h"
#include "utils/constants.h"
```

### Best Practices

1. **Don't Modify**: These files are meant to be reusable and should not be modified unless absolutely necessary. Any changes here could have far-reaching impacts on the rest of the project.

2. **Documentation**: Always keep the documentation in these files up to date. If you make a change, update the corresponding comments to reflect that change.

3. **Single Responsibility**: Each utility should have a single responsibility. For example, the `atom.h` file should only contain matters related to the Atom class.

4. **Version Control**: Keep track of any changes made to these files, as they can affect multiple parts of the project.

## Contributing

If you feel the need to add a new utility, please ensure that it adheres to the standards and patterns already set in the existing code. New utilities should be general enough to be reusable and should be well-documented.

## Questions or Issues

If you encounter any difficulties or have questions about how to use these utilities, please refer to the main project README or consult with the project maintainers.

Thank you for visiting the `utils` directory. Happy coding!