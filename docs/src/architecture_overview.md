# Architecture

Agate separates model-family definitions from the generic machinery that constructs and runs them.

## Model workflow

A model build has four broad stages:

1. **Model definition** (`Models/`, `Factories/`, and `Library/`).
   A model family provides its community, parameter definitions, tracer equations, and other scientific defaults. External packages can provide families through the same factory and recipe interfaces.

2. **Configuration and construction** (`Configuration/` and `Construction/`).
   Agate normalizes community and interaction inputs, resolves parameter defaults and overrides, and assembles the concrete biogeochemistry object. Recipes capture the authored scientific definition; manifests record the resolved scientific state.

3. **Equation assembly** (`Equations/` and `Library/`).
   Reusable scientific components and compiled equations are combined into tracer tendencies for the constructed model.

4. **Runtime and inspection** (`Runtime/`, `Diagnostics/`, and `Introspection.jl`).
   Runtime utilities provide GPU-safe tracer access and indexing for Oceananigans kernels, while diagnostics and introspection tools expose the structure of the assembled model.

## Source tree

```text
src/
∟ Models/                 # model-family definitions and user-facing constructors
∟ Factories/              # parameter definitions and family defaults
∟ Configuration/          # community and interaction normalization
∟ Construction/           # recipes, manifests, and model assembly
∟ Equations/              # compiled equation wrappers
∟ Library/                # reusable scientific components
∟ Runtime/                # runtime tracer access and indexing
∟ Diagnostics/            # model checks and diagnostics
∟ Introspection.jl        # model inspection utilities

examples/                 # runnable examples
test/                     # behavior and integration tests
paper/                    # paper-related workflows
```
