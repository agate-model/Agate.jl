# Architecture

Agate separates scientific model definition from setup-time compilation and runtime execution.

## Model workflow

A process-defined model moves through five stages:

1. **Definition** (`Configuration/`, `Processes/`, and `Parameters/`).
   `Plankton` describes prognostic state identity and structural realization, while `Pool` describes one scalar elemental material inventory. A plankton realizes PFTs, each with one or more SizeClasses, and carries one or more aligned prognostic states per realized SizeClass. Named processes describe scientific transformations and bind formulation-local parameter slots directly to stable model-parameter names. The keyed parameter block owns each parameter's default and storage policy.

2. **Validation, canonicalization, and realization** (`Processes/` and `Configuration/`).
   Agate validates scientific references and structural compatibility, canonicalizes process identity and factor order, canonicalizes intrinsic or family plankton realization once before layout construction, realizes each PFT into one or more SizeClasses separately from concrete prognostic tracers, discovers process drivers, and resolves process-local participant axes.

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
                       ModelDefinition
        authored components / processes / parameters
                              |
                              v
                 CanonicalModelDefinition
            validated scientific semantics
                              |
                              v
                       ModelLayout
            realized entities / tracers / topology
                              |
                              v
                      ParameterPlan
          axes / labels / parameter storage plan
                       /             \
                      /               \
                     v                 v
       resolved parameter values    CompileContext
                     |                 |
                     v                 v
       runtime parameter values     FluxSpec tuple
                     |                 |
                     |                 v
                     |          static tendency IR
                     |                 |
                      \               /
                       \             /
                        v           v
                           AgateBGC
```

Canonicalization establishes scientific meaning. Model realization establishes named topology.
Parameter planning establishes labelled storage. Runtime indices and array positions implement
those named structures; they do not define scientific identity.

## Custom process extension boundary

Custom process implementations that need new flux topology can extend
`Processes.process_facts` to attach setup-validated semantic facts to an `AbstractProcess` and
`Compilation.process_fluxes` to lower the resulting `CanonicalProcess` using the shared
`CompileContext`. Parameterized custom processes use `Compilation.process_parameter_operands`
to receive process-owned static operands without depending on dense binding references. Custom
processes otherwise follow the same model-definition, validation, canonicalization,
parameter-planning, construction, and runtime pipeline as built-in processes. Formulations and
factors use the same method-based slot/input/rate protocol.

## Scientific boundaries

Components describe structure rather than ecological role. A `Plankton` declares prognostic `states`, one `reference_state`, and realizes one or more PFTs. Every PFT realizes one or more SizeClasses. `variable_states` is the derived set of all remaining states. Each realized SizeClass carries one tracer per prognostic state: a one-state plankton uses identities such as `P_1` and `P_2`, while a multi-state plankton realizes state-qualified tracers such as `P_1_carbon` and `P_1_nitrogen`. Interaction and plankton parameter axes are realized per parameter from its actual process applicability and therefore do not grow with state multiplicity. Setup resolves authored state identities to internal plankton-state references and compilation lowers those references to concrete static tracer operands, so no state-name lookup enters runtime kernels. `state_element` centrally interprets the optional conserved-element bookkeeping identity for each state: canonical elemental state names map to themselves and non-elemental states such as chlorophyll map to `nothing`. Derived quantities such as N:C can then be expressed from prognostic state operands.

`plankton_pfts` maps each logical Plankton directly to its named PFT size specifications. PFT mapping order is canonicalized by identity before layout realization, so authored key order does not change tracer, parameter, interaction, or manifest order. Explicit diameter vectors are likewise canonicalized by ascending diameter before SizeClass identities are assigned. Physical layout uses a structural component order: all scalar Pool tracers are realized first, followed by Plankton in component order and each Plankton's canonical PFT/SizeClass order.

A PFT is defined by its functional parameter identity, not by size. Every PFT realizes at least one SizeClass. Without an explicit size structure it realizes one implicit singleton SizeClass named by the PFT, with no diameter metadata; an explicit size structure realizes one or more `<pft>_<index>` SizeClasses with physical diameters. `n=0` is only a high-level named-model constructor shorthand for that implicit singleton and is normalized to `nothing` before core realization; `n=1` denotes one explicit SizeClass with a real diameter.

A mixotroph is an ordinary plankton participating in both growth and living-prey consumption. `Consumption` is the single consumer-resource process for living prey, bacterivory, mixotrophy, and material-pool consumption. Process products are expressed directly through `products=` or `unassimilated_products=`. `Products` provides conservative named allocation when one process flux has multiple destinations. Collection-valued participant roles use plural keywords such as `plankton=`, `consumers=`, `resources=`, and `sources=`; each accepts either one `Symbol` or a tuple and is canonicalized to a tuple during authoring. Remineralization maps many `sources=` to one `destination=`. Bacterioplankton may consume POM and be consumed as living prey through the same consumer-resource machinery. Pools are scalar material inventories; size/PFT realization is owned by `Plankton`.

Named factors are multiplicative within a process, while independent named processes add through their fluxes to a tracer equation. Growth bookkeeping is validated per Element: an Element represented by an explicit prognostic plankton State is updated by an explicit state-changing process such as `NutrientUptake` and cannot also be supplied implicitly through Growth stoichiometry, while implicit Elements may be coupled through `FixedStoichiometry`. `NutrientLimitation.responses` is keyed by Element identity, independently of whether each response reads an external Pool or an internal quota State. Built-in Growth and `NutrientUptake` do not synthesize arbitrary non-elemental prognostic states; model families that require such synthesis can currently provide it through the custom-process extension boundary. Products and stoichiometry map process rates into affected material and element pools. A product may target one pool directly, derive several elemental products from a one-element source through `FixedStoichiometry`, or route the actual elemental inventories of a multi-state plankton through an element-to-pool mapping. Multi-state mortality and living-prey consumption use the reference state to define the shared specific loss intensity, apply that intensity to every prognostic state, and route only states with an Element into elemental products. Non-elemental states such as chlorophyll are removed proportionally without creating a second elemental inventory. When one mortality or living-prey consumption process routes multi-element products from multiple source plankton, those sources currently must expose the same prognostic Element set; heterogeneous source Element sets are a deferred extension rather than part of the v0.12 contract.

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
