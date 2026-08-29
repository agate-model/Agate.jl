# Architecture

Agate separates scientific model definition from setup-time compilation and runtime execution.

## Model workflow

A process-defined model moves through five stages:

1. **Definition** (`Configuration/`, `Processes/`, and `Parameters/`).
   `Population` and `Pool` components describe state identity and structure. A population owns ecological classes and may carry one or more aligned prognostic states. Named processes describe scientific transformations and bind formulation-local parameter slots directly to stable model-parameter names. The keyed parameter block owns each parameter's default and storage policy.

2. **Validation, canonicalization, and realization** (`Processes/` and `Configuration/`).
   Agate validates scientific references and structural compatibility, canonicalizes process identity and factor order, canonicalizes intrinsic or family population realization once before layout construction, realizes ecological classes separately from concrete prognostic tracers, discovers process drivers, and resolves process-local participant axes.

3. **Flux compilation** (`Compilation/`).
   A setup-time compile context carries the canonical definition, realized layout, and parameter plan through lowering. Parameter operands resolve ecological class identity directly against realized storage labels during setup, so only final static indices enter runtime terms. Participant states and one-tracer-per-class components are realized through one shared tracer + ecological-position traversal, which is reused by one-axis and two-axis process lowering. Growth canonicalization resolves one reference resource draw plus any additional stoichiometric resource draws, so lowering follows one resource-transfer rule regardless of the authored nutrient-response form. Growth also owns its maximum-rate binding; Smith and Geider light factors consume that resolved process scale rather than declaring a second parameter slot. Each named process produces generic target + rate + weight flux specifications, which are grouped by target tracer and lowered into static compiled equations with no target symbols or process metadata in runtime terms.

4. **Construction and replay** (`Construction/`).
   Direct `ModelDefinition` construction resolves defaults and overrides from the model
   definition. Bundled families such as NiPiZD and external registered families such as
   DARWIN in AgateRegistry.jl use the same definition-driven core while adding durable
   recipe/replay identity. Recipes record a registered family, its exact scientific
   definition version, and the
   canonical realization inputs required by the family constructor; replay then follows the
   normal family construction path.

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
`CompileContext`. Parameterized custom processes use `Compilation.process_parameter_operands`
to receive process-owned static operands without depending on dense binding references. Custom
processes otherwise follow the same model-definition, validation, canonicalization,
parameter-planning, construction, and runtime pipeline as built-in processes. Formulations and
factors use the same method-based slot/input/rate protocol.

## Scientific boundaries

Components describe structure rather than ecological role. A `Population` owns ecological classes and declares its prognostic `states` plus one `reference_state`; `variable_states` is the derived set of all remaining states. Prognostic states are aligned inventories carried by those classes. A one-state population uses class tracer identities such as `P_1` and `P_2`, while a multi-state population realizes class-qualified state tracers such as `P_1_carbon` and `P_1_nitrogen`. Interaction and parameter axes remain ecological-class axes and therefore do not grow with state multiplicity. Setup resolves authored state identities to internal population-state references and compilation lowers those references to concrete static tracer operands, so no state-name lookup enters runtime kernels. `state_element` centrally interprets the optional conserved-element bookkeeping identity for each state: canonical elemental state names map to themselves and non-elemental states such as chlorophyll map to `nothing`. Derived quantities such as N:C can then be expressed from prognostic state operands.

A mixotroph is an ordinary population participating in both growth and living-prey consumption. `Consumption` is the single consumer-resource process for living prey, bacterivory, mixotrophy, and material-pool consumption. Process products are expressed directly through `products=` or `unassimilated_products=`. `Products` provides conservative named allocation when one process flux has multiple destinations. Collection-valued participant roles use plural keywords such as `populations=`, `consumers=`, `resources=`, and `sources=`; each accepts either one `Symbol` or a tuple and is canonicalized to a tuple during authoring. Remineralization maps many `sources=` to one `destination=`. Bacterioplankton may consume POM and be consumed as living prey through the same consumer-resource machinery. Structured material pools use the same component-layout machinery as structured populations.

Named factors are multiplicative within a process, while independent named processes add through their fluxes to a tracer equation. Products and stoichiometry map process rates into affected material and currency pools. A product may target one pool directly or a currency-to-pool mapping under `FixedStoichiometry`; setup canonicalizes both authoring forms to one product-to-currency-to-component target shape before lowering, so compilation uses the same allocation machinery for both.

Parameter slots are resolved once during canonicalization. Each canonical process carries dense references into the single ordered `ParameterBinding` tuple, arranged alongside its process, factor, product, and stoichiometry structure. Lowering consumes those references directly; scientific tree paths remain setup-time validation/error context rather than a second lookup key representation.

## Source tree

```text
src/
|-- Configuration/         # components, diameter specs, realized layout, interactions
|-- Processes/             # authoring, canonicalization, parameter planning, authored validation, parameter validation
|-- Compilation/           # runtime IR, process flux lowering, static compilation
|-- ModelFamilies/         # named family interface and definition/version hooks
|-- Parameters/            # keyed parameters, defaults, and storage policy
|-- Construction/          # parameter realization, construction, recipes, manifests, replay
|-- Library/               # scientific kernels, forcings, and allometric utilities
|-- Runtime/               # active parameters and box-ODE utilities
|-- Diagnostics/           # model checks and diagnostics
|-- Models/                # bundled named model families
`-- Introspection.jl       # model inspection utilities

examples/                  # runnable examples
test/                      # behavior and integration tests
```
