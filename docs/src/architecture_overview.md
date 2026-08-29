# Architecture

Agate separates scientific model definition from setup-time compilation and runtime execution.

## Model workflow

A process-defined model moves through five stages:

1. **Definition** (`Configuration/`, `Processes/`, and `Parameters/`).
   `Plankton` describes prognostic state identity and structural realization, while `Pool` describes one scalar elemental material inventory. A plankton realizes PFTs, each optionally subdivided into SizeClasses, and carries one or more aligned prognostic states per realized PFT entity. Named processes describe scientific transformations and bind formulation-local parameter slots directly to stable model-parameter names. The keyed parameter block owns each parameter's default and storage policy.

2. **Validation, canonicalization, and realization** (`Processes/` and `Configuration/`).
   Agate validates scientific references and structural compatibility, canonicalizes process identity and factor order, canonicalizes intrinsic or family plankton realization once before layout construction, realizes PFT entities and optional SizeClasses separately from concrete prognostic tracers, discovers process drivers, and resolves process-local participant axes.

3. **Flux compilation** (`Compilation/`).
   A setup-time compile context carries the canonical definition, realized layout, and parameter plan through lowering. Parameter operands resolve realized entity identity directly against storage labels during setup, so only final static indices enter runtime terms. Participant states and one-tracer-per-entity components are realized through one shared tracer + entity-position traversal, which is reused by one-axis and two-axis process lowering. Growth authoring declares `reference_resource` and any element-keyed `additional_resources` explicitly, so canonicalization no longer infers material transfer from nutrient-factor topology. Growth may have zero or more multiplicative factors; with none, its rate is simply the maximum-rate scale times biomass. Additional stoichiometric draws are valid only for Elements that the Growth plankton does not already carry as explicit prognostic States. `NutrientLimitation` may therefore mix external `NutrientResponse` and internal `QuotaResponse` subfactors while material transfer remains process-owned. Growth also owns its maximum-rate binding; Smith and Geider light factors consume that resolved process scale rather than declaring a second parameter slot. Each named process produces generic target + rate + weight flux specifications, which are grouped by target tracer and lowered into static compiled equations with no target symbols or process metadata in runtime terms.

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

Components describe structure rather than ecological role. A `Plankton` declares prognostic `states`, one `reference_state`, and realizes one or more PFTs. Each PFT may optionally be subdivided into SizeClasses. `variable_states` is the derived set of all remaining states. Each realized PFT entity carries one tracer per prognostic state: a one-state plankton uses identities such as `P_1` and `P_2`, while a multi-state plankton realizes state-qualified tracers such as `P_1_carbon` and `P_1_nitrogen`. Interaction and plankton parameter axes remain realized-plankton-entity axes and therefore do not grow with state multiplicity. Setup resolves authored state identities to internal plankton-state references and compilation lowers those references to concrete static tracer operands, so no state-name lookup enters runtime kernels. `state_element` centrally interprets the optional conserved-element bookkeeping identity for each state: canonical elemental state names map to themselves and non-elemental states such as chlorophyll map to `nothing`. Derived quantities such as N:C can then be expressed from prognostic state operands.

`plankton_pfts` is the authority for PFT ordering. `pft_size_structures` is keyed metadata and is canonicalized to that hierarchy before layout realization, so its NamedTuple key order cannot independently change tracer, parameter, interaction, or manifest order. Physical layout uses a structural component order: all scalar Pool tracers are realized first, followed by Plankton in authored component order and each Plankton's PFT/SizeClass order.

A PFT is defined by its functional parameter identity, not by size. Size structure is optional metadata within a PFT: an unsized PFT realizes one entity named by the PFT, while a sized PFT realizes one or more `<pft>_<index>` SizeClasses. `n=0` is only a high-level named-model constructor shorthand for an unsized PFT and is normalized away before core realization; `n=1` denotes one actual SizeClass with a real diameter.

A mixotroph is an ordinary plankton participating in both growth and living-prey consumption. `Consumption` is the single consumer-resource process for living prey, bacterivory, mixotrophy, and material-pool consumption. Process products are expressed directly through `products=` or `unassimilated_products=`. `Products` provides conservative named allocation when one process flux has multiple destinations. Collection-valued participant roles use plural keywords such as `plankton=`, `consumers=`, `resources=`, and `sources=`; each accepts either one `Symbol` or a tuple and is canonicalized to a tuple during authoring. Remineralization maps many `sources=` to one `destination=`. Bacterioplankton may consume POM and be consumed as living prey through the same consumer-resource machinery. Pools are scalar material inventories; size/PFT realization is owned by `Plankton`.

Named factors are multiplicative within a process, while independent named processes add through their fluxes to a tracer equation. Growth bookkeeping is validated per Element: an Element represented by an explicit prognostic plankton State is updated by an explicit state-changing process such as `NutrientUptake` and cannot also be supplied implicitly through Growth stoichiometry, while implicit Elements may be coupled through `FixedStoichiometry`. Products and stoichiometry map process rates into affected material and element pools. A product may target one pool directly, derive several elemental products from a one-element source through `FixedStoichiometry`, or route the actual elemental inventories of a multi-state plankton through an element-to-pool mapping. Multi-state mortality and living-prey consumption use the reference state to define the shared specific loss intensity, apply that intensity to every prognostic state, and route only states with an Element into elemental products. Non-elemental states such as chlorophyll are removed proportionally without creating a second elemental inventory.

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
