# Architecture

Agate separates scientific model definition from setup-time compilation and runtime execution.

## Model workflow

A process-defined model moves through five stages:

1. **Definition** (`Configuration/`, `Processes/`, and `Factories/parameter_directory.jl`).
   `Population` and `Pool` components describe state identity and structure. Named processes describe scientific transformations, and parameter definitions bind defaults to formulation-owned semantic requirements.

2. **Normalization and realization** (`Processes/` and `Configuration/`).
   Agate canonicalizes process identity and factor order, realizes structured components into concrete tracer classes, discovers process drivers, and resolves process-local participant axes.

3. **Contribution compilation** (`Compilation/`).
   Each named process produces typed contributions to the concrete tracers it affects. Contributions are grouped by target tracer and lowered during setup into static compiled equations.

4. **Construction and replay** (`Construction/`).
   Direct `ModelDefinition` construction resolves defaults and overrides from the model definition. Named model families such as NiPiZD and DARWIN use the same definition-driven core while adding durable recipe/replay identity.

5. **Runtime and inspection** (`Runtime/`, `Diagnostics/`, and `Introspection.jl`).
   Runtime kernels evaluate lean compiled terms with resolved tracer and parameter indices. Diagnostics and introspection expose the realized model without reinterpreting the model definition.

The central data flow is:

```text
components
+ named processes
+ parameter definitions
        |
        v
normalization
        |
        v
process-local topology + parameter bindings
        |
        v
typed process contributions
        |
        v
setup-time grouping / static lowering
        |
        v
one compiled equation per concrete tracer
        |
        v
lean runtime
```

## Scientific boundaries

Components describe structure rather than ecological role. A mixotroph is an ordinary population participating in both growth and grazing. Bacterioplankton are ordinary populations that may consume POM through `Consumption` and be consumed as living prey through `Grazing`. Structured material pools use the same component-layout machinery as structured populations.

Named factors are multiplicative within a process, while independent named processes add through their contributions to a tracer equation. Routing and stoichiometry map process rates into affected material and currency pools.

## Source tree

```text
src/
|-- Configuration/         # components, realization, interactions
|-- Processes/             # process definitions and normalization
|-- Compilation/           # topology, contributions, static lowering
|-- Factories/             # parameter defaults and named-family adapters
|-- Construction/          # direct construction, recipes, manifests, replay
|-- Equations/             # compiled equation wrappers
|-- Library/               # reusable scientific formulations
|-- Runtime/               # runtime tracer access and indexing
|-- Diagnostics/           # model checks and diagnostics
|-- Models/                # bundled named model families
`-- Introspection.jl       # model inspection utilities

examples/                  # runnable examples
test/                      # behavior and integration tests
```
