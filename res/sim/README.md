# Simulator Configuration Permutations

This directory contains configuration files defining which simulator template 
permutations to build.

## File Format

Each line defines one configuration as comma-separated values:
```
dimension,interpolation_order,refined_particle_count
```

Example:
```
1,2,3
```
means: 1D simulation, 2nd order interpolation, 3 refined particles per cell

## Files

- **all.txt**: Complete set of all supported configurations (default)
- Custom permutation files can be specified via `-DPHARE_PERMUTATIONS=/path/to/file.txt`

## Adding New Configurations

1. Ensure the configuration is scientifically valid
2. Add line to all.txt (or your custom file)
3. Rebuild (this will compile a new Python module for that configuration)

Each configuration results in a separate compiled Python module (e.g., `cpp_1_2_3.so`).

## Performance Considerations

Each configuration results in a separate compiled Python module (~5-20MB each).
Building all 23 configurations takes significant time and disk space.

For development, consider creating a minimal set with just the configurations you need:

```bash
# Create a file with only the configurations you need
echo "1,2,3" > my_configs.txt
echo "2,2,4" >> my_configs.txt

# Build only those configurations
cmake -DPHARE_PERMUTATIONS=/path/to/my_configs.txt ..
```

## How It Works

At build time, CMake reads this file and generates one Python module per line:
- Configuration `1,2,3` becomes module `pybindlibs.cpp_1_2_3`
- Configuration `2,2,4` becomes module `pybindlibs.cpp_2_2_4`
- etc.

When you create a Simulator with specific parameters in Python, PHARE automatically
loads the appropriate module for that configuration.
