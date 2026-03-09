# Architecture
## Overview

At a high level, Agate is organized into four broad layers:

1. **Model constructors** in `Models/`, such as [`NiPiZD`](@ref NiPiZD) and [`DARWIN`](@ref DARWIN), which provide the user-facing constructors and define model-specific variants.

2. **Reusable model building blocks** in `Library/`, which provide shared components for processes such as photosynthesis, predation, remineralization, light limitation, mortality, and allometry.

3. **Construction-time machinery** in `Factories/`, `Configuration/`, `Construction/`, and `Equations/`, which resolves defaults, normalizes configuration inputs, and assembles a concrete biogeochemistry object.

4. **Runtime and inspection utilities** in `Runtime/`, `Diagnostics/`, and `Introspection.jl`, which support tracer access, kernel use, diagnostics, and model summary after construction.

## Model workflow

A typical model build in Agate has seven components.

1. **Model entry point** (`Models/`).
   Each model is defined using a user-facing constructor such as `Agate.Models.<Model>.construct(...)`. The `Models/` layer defines the public entry points for the built-in model families.

2. **Model defaults and community selection** (`Models/` and `Factories/`).
   For each model constructor, Agate selects the factory, parameter definitions, and default community information needed for the build. `Models/` provides model-family-specific definitions, while `Factories/` provides shared defaults and construction metadata.

3. **Configuration normalization** (`Configuration/`).
   Community structure, trophic roles, and interaction inputs are converted into a consistent internal representation. The `Configuration/` layer is responsible for parsing and normalizing these inputs before assembly begins.

4. **Parameter resolution** (`Factories/`, `Configuration/`, and `Construction/`).
   After the configuration has been normalized, Agate resolves explicit inputs, applies defaults, and constructs any derived matrices required by the model. This step draws on shared metadata, normalized configuration, and assembly logic.

5. **Model assembly** (`Construction/`, `Equations/`, and `Library/`).
   Tracer equations and supporting components are assembled into a concrete Agate model object. `Construction/` contains the core assembly pipeline, while `Equations/` and `Library/` provide equation wrappers and reusable scientific components.

6. **Runtime integration** (`Runtime/`).
   Once assembled, the model exposes tracer accessors, indexing utilities, and other helpers needed at runtime. The `Runtime/` layer provides GPU-safe tracer access and related support for Oceananigans kernels.

7. **Inspection and diagnostics** (`Diagnostics/` and `Introspection.jl`).
   After construction, helper tools are available to inspect the generated model and evaluate its structure or behavior. `Diagnostics/` and `Introspection.jl` provide utilities for understanding, checking, and summarizing the assembled model.

## Source tree

- `Models/` defines the user-facing entry points for built-in model families.
- `Factories/` defines shared defaults, parameter metadata, and related construction settings.
- `Configuration/` parses communities, trophic roles, and interaction structure into a normalized internal form.
- `Construction/` assembles the concrete model object from validated inputs, defaults, and derived quantities.
- `Equations/` provides compiled equation wrappers used during model assembly.
- `Library/` contains reusable model components shared across model families.
- `Runtime/` provides tracer accessors, indexing utilities, and other runtime helpers for Oceananigans kernels.
- `Diagnostics/` and `Introspection.jl` support inspection, checking, and summary of the generated model.

At the repository level, the main directories are:

```text
src/
∟ Agate.jl                # top-level module wiring
∟ Models/                 # user-facing constructors and model variants
∟ Factories/              # shared defaults and parameter metadata
∟ Configuration/          # community parsing and interaction handling
∟ Construction/           # shared model assembly pipeline
∟ Equations/              # compiled equation wrappers
∟ Library/                # reusable model components
∟ Runtime/                # tracer indexing and accessor machinery
∟ Diagnostics/            # diagnostics and consistency checks
∟ Introspection.jl        # tools for model summary and inspection

docs/
∟ src/                    # documentation source files
∟ make.jl                 # Documenter build configuration and navigation

examples/                 # example scripts

test/                     # integration and behavior checks

paper/                    # paper-related and GPU-oriented workflows
```
