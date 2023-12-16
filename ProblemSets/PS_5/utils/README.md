# `/utils` Directory

## Overview

The `utils` directory is a vital part of our project, serving as the hub for various utility header files. These headers play a crucial role in the project by defining classes, functions, structures, and constants. Their primary purpose is to ensure modularity, reusability, and a clean organization of the codebase.

## Directory Contents

- `atom.h`: Defines the Atom class, crucial for representing atoms within our simulations or calculations. The Atom class provides the backbone for numerous algorithms and computations in the project.

- `vector3d.h`: Contains the Vector3D class, designed for representing and manipulating 3-dimensional vectors. It includes methods for basic vector operations such as addition, subtraction, and scalar products.

- `constants.h`: Houses essential constants used across the project, encompassing physical constants and immutable values crucial for various calculations and operations.

- `factorial.h`: Provides functionality for computing factorials, often used in mathematical and scientific computations within the project. This may include functions like `factorial(int)` or `doubleFactorial(int)`.

## Usage

### Including Utilities

To incorporate these utility headers in your code, use the `#include` directive in your `.cpp` files:

```cpp
#include "utils/atom.h"
#include "utils/vector3d.h"
#include "utils/constants.h"
#include "utils/factorial.h"
```

### Best Practices

1. **Modification**: Avoid altering these files unless necessary. Changes can lead to widespread implications across the project.

2. **Documentation**: Keep the documentation within these files current. Any updates or modifications should be accurately reflected in their comments.

3. **Focused Responsibility**: Ensure each utility file focuses on a single aspect or functionality (e.g., `factorial.h` should solely deal with factorial calculations).

4. **Version Control**: Monitor changes to these utilities closely, as they can impact various components of the project.

## Contributing

Introducing new utilities should be done with careful consideration, ensuring they align with existing standards and patterns. They should be broadly applicable, well-documented, and maintain the codebase's integrity and readability.

## Support

For any assistance, queries, or discussion about these utilities, consult the main project README or reach out to the project maintainers.

Thank you for exploring the `utils` directory — your interest and contributions help in shaping a robust and efficient project!