# Architecture

Agate separates scientific model definition from setup-time compilation and runtime execution.

## Model workflow

A process-defined model moves through five stages:

1. **Definition** (`Configuration/`, `Processes/`, and `Parameters/`).
   `Population` and `Pool` components describe state identity and structure. A population owns ecological classes and may carry one or more aligned prognostic states. Named processes describe scientific transformations and bind formulation-local parameter slots directly to stable model-parameter names. The keyed parameter block owns each parameter's default and storage policy.

2. **Validation, canonicalization, and realization** (`Processes/` and `Configuration/`).
   Agate validates scientific references and structural compatibility, canonicalizes process identity and factor order, canonicalizes intrinsic or family population realization once before layout construction, realizes ecological classes separately from concrete prognostic tracers, discovers process drivers, and resolves process-local participant axes.

3. **Flux compilation** (`Compilation/`).
   A setup-time compile context carries the canonical definition, realized layout, and parameter plan through lowering. Parameter operands resolve ecological class identity directly against realized storage labels during setup, so only final static indices enter runtime terms. Validated Growth routing is represented by typed canonical routing variants rather than symbolic mode flags. Each named process produces generic target + rate + weight flux specifications, which are grouped by target tracer and lowered into static compiled equations with no target symbols or process metadata in runtime terms.

4. **Construction and replay** (`Construction/`).
   Direct `ModelDefinition` construction resolves defaults and overrides from the model definition. Named model families such as NiPiZD and DARWIN use the same definition-driven core while adding durable recipe/replay identity. Recipes record a registered family, its exact scientific definition version, and the canonical realization inputs required by the family constructor; replay then follows the normal family construction path.

5. **Runtime and inspection** (`Runtime/`, `Diagnostics/`, and `Introspection.jl`).
   Runtime kernels evaluate lean compiled terms with resolved tracer and parameter indices. Diagnostics and introspection expose the realized model without reinterpreting the model definition.

The central data flow is:

```text
components
+ named processes
+ keyed parameters
        |
        v
validation + canonicalization
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

## Custom process extension boundary

Custom process implementations that need new flux topology can extend
`Processes.process_facts` to attach setup-validated facts to an `AbstractProcess` and
`Compilation.process_fluxes` to lower the resulting `NamedProcess` using the shared
`CompileContext`. Custom processes otherwise follow the same model-definition, validation,
canonicalization, parameter-planning, construction, and runtime pipeline as built-in processes. Formulations and
factors use the same method-based slot/input/rate protocol.

## Scientific boundaries

Components describe structure rather than ecological role. A `Population` owns ecological classes; prognostic states are aligned inventories carried by those classes. A one-state population preserves the established tracer identity (`P_1`, `P_2`, ...), while a multi-state population realizes class-qualified state tracers such as `P_1_carbon` and `P_1_nitrogen`. Interaction and parameter axes remain ecological-class axes and therefore do not grow with state multiplicity. `population_state(:P, :nitrogen)` identifies one aligned state at setup time; compilation resolves that ecological state reference to a concrete static tracer operand, so no state-name lookup enters runtime kernels. Derived quantities such as N:C can therefore be expressed from prognostic state operands while preserving concrete/isbits compiled equations.

A mixotroph is an ordinary population participating in both growth and living-prey consumption. `Consumption` is the single consumer-resource process for living prey, bacterivory, mixotrophy, and material-pool consumption. Process products are expressed directly through `products=` or `unassimilated_products=`. `Products` provides conservative named allocation when one process flux has multiple destinations. Collection-valued participant roles use plural keywords such as `populations=`, `consumers=`, `resources=`, and `sources=`; each accepts either one `Symbol` or a tuple and is canonicalized to a tuple during authoring. Remineralization has one `destination=` because its current scientific contract is many sources to one destination. Bacterioplankton may consume POM and be consumed as living prey through the same consumer-resource machinery. Structured material pools use the same component-layout machinery as structured populations.

Named factors are multiplicative within a process, while independent named processes add through their fluxes to a tracer equation. Products and stoichiometry map process rates into affected material and currency pools. A product may target one pool directly or a currency-to-pool mapping under `FixedStoichiometry`; setup canonicalizes both authoring forms to one product-to-currency-to-component target shape before lowering, so compilation uses the same allocation machinery for both.

## Source tree

```text
src/
|-- Configuration/         # components, realization, interactions
|-- Processes/             # process authoring, validation, canonicalization, parameter schema
|-- Compilation/           # runtime IR, process flux lowering, static compilation
|-- ModelFamilies/         # registered family identity and canonical definitions
|-- Parameters/            # keyed parameters, defaults, and storage policy
|-- Construction/          # direct construction, recipes, manifests, replay
|-- Library/               # reusable scientific formulations
|-- Runtime/               # active parameters and box-ODE utilities
|-- Diagnostics/           # model checks and diagnostics
|-- Models/                # bundled named model families
`-- Introspection.jl       # model inspection utilities

examples/                  # runnable examples
test/                      # behavior and integration tests
```
