# Biological interpretation of Wind and Current metrics
# Extracted from eReefs (monthly, GBR1)

---

## 📊 Metric Categories

### 1. Basic physical metrics
### 2. Directional metrics  
### 3. GBR-specific metrics
### 4. Seasonal metrics
### 5. Biological proxy metrics
### 6. Larval transport metrics

---

## 1️⃣ BASIC PHYSICAL METRICS

### `mean_wind_speed` / `mean_current_speed` (m/s)
**What it measures:** Average velocity magnitude  
- **High values**: Greater physical stress on corals, higher wave energy
- **Low values**: Calmer conditions, less mechanical stress
- **Relevance**: Selection for morphological traits (branching vs massive), tissue strength, attachment, pop size

### `sd_wind_speed` / `sd_current_speed` (m/s)
**What it measures:** Temporal variability in speed  
- **High SD**: Unpredictable environment → selection for plasticity
- **Low SD**: Stable environment → selection for specialization
- **Relevance**: Environmental predictability, disturbance regime

### `cv_wind_speed` / `cv_current_speed` (coefficient of variation)
**What it measures:** Relative variability (SD/mean)  
- **High CV**: Proportionally variable conditions → bet-hedging strategies?
- **Low CV**: Consistent conditions → specialists?
- **Relevance**: Environmental stochasticity

### `wind_speed_range` / `current_speed_range` (m/s)
**What it measures:** Difference between calm and extreme conditions  
- **Large range**: Sites experience both extremes → need for wide environmental tolerance
- **Small range**: Narrow environmental envelope → specialisation possible

---

## 2️⃣ DIRECTIONAL METRICS

### `mean_wind_direction` / `mean_current_direction` (degrees)
**What it measures:** Predominant flow direction (0° = North, clockwise)  
- Determines fetch (distance over which wind/waves build)
- Affects which direction larvae are transported
- Determines shelter/exposure patterns
- ***Relevance**: Direction of larval transport, isolation-by-distance patterns

**Example:**
- 135° (SE) = trade wind alignment on GBR
- 315° (NW) = winter storm direction

### `directional_constancy` (0-1)
**Formula:** vector_magnitude / mean_speed 
**What it measures:** How consistently flow comes from same direction  
- **1.0**: Wind/current always from same direction → predictable transport
- **0.5**: Moderate variability → mixed transport patterns
- **0.0**: Completely variable direction → local retention?

**Critical for:**
- **Larval dispersal directionality**: High constancy = unidirectional gene flow
- **Connectivity modeling**: Predictable vs chaotic dispersal
- **IBD vs IBE**: Constant flow = distance matters, variable = environment matters

### `directional_variability` (0-1)
**What it measures:** Inverse of constancy (1 - constancy)  
**Biological meaning:**
- **High**: Variable flow → larvae mixed locally → panmixia or local structure
- **Low**: Constant flow → larvae swept away → long-distance dispersal

### `predominant_wind_sector` / `predominant_current_sector`
**What it measures:** Cardinal/intercardinal direction (N, NE, E, SE, S, SW, W, NW)  
- Categorical version of mean direction

---

## 3️⃣ GBR-SPECIFIC METRICS

# Would use & interpret the following with caution

### `trade_wind_component` (m/s)
**Formula:** u·cos(135°) + v·sin(135°) 
**What it measures:** Wind/current component aligned with SE trade winds (135°)  
**Biological meaning:**
- **Positive**: Aligned with SE trades (typical for GBR)
- **Negative**: Against trades (unusual, sheltered)
- **Zero**: Perpendicular to trades

### `alongshore_wind_component` / `alongshore_current_component` (m/s)
**What it measures:** North-South component (v)  
- **Positive**: Northward flow → connectivity to northern populations
- **Negative**: Southward flow → connectivity to southern populations
- **Critical for:** N-S gene flow along GBR

### `crossshore_wind_component` / `crossshore_current_component` (m/s)
**What it measures:** East-West component (u)  
- **Positive**: Offshore (eastward) → larvae exported to open ocean
- **Negative**: Onshore (westward) → larvae retained on reef

### Exposure to Cardinal Directions
- `easterly_exposure`, `westerly_exposure`
- `northerly_exposure`, `southerly_exposure`

**What they measure:** Mean flow FROM each direction  
- **Easterly**: Exposure to offshore winds/waves
- **Westerly**: Sheltered from open ocean, but exposed to land effects
- **Northerly**: Potential connectivity from northern reefs
- **Southerly**: Potential connectivity from southern reefs

---

## 4️⃣ SEASONAL METRICS

### Seasonal Speeds
- `summer_wind_speed`, `winter_wind_speed`, etc.
- **Summer (Dec-Feb)**: Spawning season
- **Winter (Jun-Aug)**: SE trade wind peak, cooler temps
- **Spring (Sep-Nov)**: Pre-spawning season
- **Autumn (Mar-May)**: Post-spawning

### `seasonal_wind_range` / `seasonal_current_range` (m/s)
**What it measures:** Difference between strongest and weakest season  
- **High range (>2 m/s)**: Strong seasonal contrast
- **Low range (<1 m/s)**: Weak seasonality

### `dominant_wind_season` / `dominant_current_season`
**What it measures:** Which season has strongest conditions  

### `summer_winter_diff` (m/s)
**What it measures:** Absolute difference between summer and winter  
- Simplified seasonal variability metric

---

## 5️⃣ BIOLOGICAL PROXY METRICS

### `mean_wave_energy_proxy` (m³/s³)
**What it measures:** Wind speed cubed (wave energy ∝ wind³)  
- Approximates wave energy without wave data
- **High**: Strong wave action → mechanical stress
- **Low**: Calm waters → less mechanical damage

### `extreme_wind_frequency` / `extreme_current_frequency` (%)
**What it measures:** Percentage of months with >95th percentile conditions  
- **High**: Frequent disturbances
- **Low**: Rare disturbances

### `calm_wind_frequency` / `calm_current_frequency` (%)
**What it measures:** Percentage of months with very low flow (<2 m/s or <0.05 m/s)  
- **High calm frequency**: Often stagnant → potential for thermal stress
- **Low calm frequency**: Always flowing → low thermal stress

---

## 6️⃣ LARVAL TRANSPORT METRICS

### `wind_transport_potential` / `current_transport_potential` (m/s)
**Formula:** directional_constancy × mean_speed  
- **High**: Strong, consistent flow → long-distance directional dispersal
- **Low**: Weak or variable flow → limited, local dispersal

### `wind_dispersal_symmetry` / `current_dispersal_symmetry` (0-1)
**Formula:** 1 - |mean_u / mean_speed|  
- **1.0**: Omnidirectional (perfectly symmetric)
- **0.0**: Strongly asymmetric (unidirectional)

### `wind_retention_potential` / `current_retention_potential`
**Formula:** (1 - constancy) × (1 / mean_speed)
- **High**: Variable, weak flow → larvae retained locally → self-recruitment
- **Low**: Strong, consistent flow → larvae exported → connectivity to distant sites

---

## 7️⃣ CATEGORICAL VARIABLES

### `wind_exposure_category` / `current_exposure_category`
**Categories:** Sheltered / Moderate / Exposed (or Low_Flow / Moderate_Flow / High_Flow)  

### `directional_consistency_category`
**Categories:** Variable / Moderate / Consistent  

### `seasonal_variability_category`
**Categories:** Low / Moderate / High  

### `wind_regime_type` / `current_regime_type`
**Categories:** Steady_Strong / Variable_Strong / Steady_Weak / Variable_Weak  

---

## 📊 QUICK REFERENCE TABLE

| Metric | Range | Low = | High = | Key For |
|--------|-------|-------|--------|---------|
| `mean_wind_speed` | 2-8 m/s | Sheltered | Exposed | Physical stress, morphology |
| `directional_constancy` | 0-1 | Variable | Consistent | Connectivity, IBD |
| `transport_potential` | 0-10 | Local retention | Long dispersal | Gene flow magnitude |
| `seasonal_range` | 0-4 m/s | Stable | Variable | Plasticity selection |
| `trade_wind_component` | -5 to +5 | Sheltered | Exposed to trades | GBR-specific |
| `retention_potential` | 0-2 | Export | Self-recruitment | Population structure |

---
