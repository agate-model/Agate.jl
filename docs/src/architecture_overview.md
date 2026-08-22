# Architecture

Agate separates scientific model definition from setup-time compilation and runtime execution.

## Model workflow

A process-defined model moves through five stages:

1. **Definition** (`Configuration/`, `Processes/`, and `Parameters/`).
   `Population` and `Pool` components describe state identity and structure. A population owns ecological classes and may carry one or more aligned prognostic states. Named processes describe scientific transformations, and parameter definitions bind defaults to formulation-owned semantic requirements.

2. **Normalization and realization** (`Processes/` and `Configuration/`).
   Agate canonicalizes process identity and factor order, realizes ecological classes separately from concrete prognostic tracers, discovers process drivers, and resolves process-local participant axes.

3. **Flux compilation** (`Compilation/`).
   Each named process produces generic target + rate + weight flux specifications. Fluxes are grouped by target tracer and lowered during setup into static compiled equations with no target symbols or process metadata in runtime terms.

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
generic target + rate + weight fluxes
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

Components describe structure rather than ecological role. A `Population` owns ecological classes; prognostic states are aligned inventories carried by those classes. A one-state population preserves the established tracer identity (`P_1`, `P_2`, ...), while a multi-state population realizes class-qualified state tracers such as `P_1_carbon` and `P_1_nitrogen`. Interaction and parameter axes remain ecological-class axes and therefore do not grow with state multiplicity. `population_state(:P, :nitrogen)` identifies one aligned state at setup time; compilation resolves that ecological state reference to a concrete static tracer operand, so no state-name lookup enters runtime kernels. Derived quantities such as N:C can therefore be expressed from prognostic state operands while preserving concrete/isbits compiled equations.

A mixotroph is an ordinary population participating in both growth and living-prey consumption. `Consumption` is the single consumer-resource process for living prey, bacterivory, mixotrophy, and material-pool consumption. Product destinations are expressed through `ProductRouting`, so routing has one representation across process types. Collection-valued participant roles use plural keywords such as `populations=`, `consumers=`, `resources=`, `sources=`, and `destinations=`; each accepts either one `Symbol` or a tuple and is canonicalized to a tuple during authoring. Bacterioplankton may consume POM and be consumed as living prey through the same consumer-resource machinery. Structured material pools use the same component-layout machinery as structured populations.

Named factors are multiplicative within a process, while independent named processes add through their fluxes to a tracer equation. Routing and stoichiometry map process rates into affected material and currency pools.

## Source tree

```text
src/
|-- Configuration/         # components, realization, interactions
|-- Processes/             # process definitions and normalization
|-- Compilation/           # topology, fluxes, static lowering
|-- ModelFamilies/         # registered family identity and canonical definitions
|-- Parameters/            # parameter definitions, defaults, and provisions
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
