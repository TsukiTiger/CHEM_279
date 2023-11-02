# Root makefile for the project

# Phony targets are those that do not correspond to file names.
.PHONY: object test run clean

# Default target: when 'make' is called without arguments, it builds 'object'
default: object

# Build object files by delegating to the src makefile
object:
	@echo "Building object files..."
	@$(MAKE) -C src

# Build and run tests by delegating to the test makefile
test:
	@echo "Building tests executable..."
	@$(MAKE) -C test

# Run test executables
run:
	@echo "Running tests..."
	@echo "H2: "
	@./bin/test ./data/H2.txt
	@echo "C2H2: "
	@./bin/test ./data/C2H2.txt
	@echo "Calculating Energy Difference: "
	@./bin/energy_diff



# Clean the project by delegating to the src and test makefiles
clean:
	@echo "Cleaning up object files in src..."
	@$(MAKE) -C src clean
	@echo "Cleaning up test binaries and objects..."
	@$(MAKE) -C test clean
