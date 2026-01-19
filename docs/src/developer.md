# Developer guide {#developer-guide}

This section is for contributors who want to extend Agate: add new models, add new parameters, or change the compiled dynamics.
It assumes you are comfortable reading Julia code and running the test suite.

## Key design ideas

- **Explicit model boundaries**: each model lives under `src/Models/<ModelName>/` and exposes a small, keyword-driven `construct` API.
- **Parameter metadata**: model parameters are declared in a `parameter_directory` so the constructor can validate shapes and provide clearer errors.
- **Role-aware interactions**: consumer-by-prey interactions are represented *canonically* as rectangular matrices, and models consume these directly.
- **GPU and parametric `FT`**: code paths avoid dynamic dispatch and allocate arrays using the chosen floating-point type and architecture.
- **Variants**: paper- or project-specific configurations live as lightweight variant specs (instead of copying entire model modules).

## Where things live

- `src/Constructor/`: the general constructor pipeline and input validation.
- `src/Utils/`: shared infrastructure such as interaction contexts, rectangular interaction storage, and small utilities.
- `src/Models/`: model-specific implementations.
- `src/Library/`: reusable biogeochemical building blocks and parameterizations.

## Start here

- Read [Adding a model](@ref "Adding a model") for the minimum set of files and functions to implement.
- Read [Variants](@ref "Variants") for the recommended way to manage manuscript / experiment configurations without repo bloat.
- Reference docs:
  - [Constructor API](@ref "Constructor API")
  - [Parameters and interaction matrices](@ref "Parameters and interaction matrices")
  - [Callable dynamics API](@ref "Callable dynamics API")
  - [API reference](@ref "API reference")

## Development workflow

- Run the tests in `test/` regularly while refactoring.
- Prefer small, explicit functions with clear type signatures.
- Keep docstrings aligned with the current implementation (avoid describing removed behavior).
