# bclake

Hydrological simulations based on SWAM cellular automata model

Ver. 2 20181231

- Version with CRU 10m climatology
- Extract CRU, with rhum + tmin and tmax (from tmp + dtr)
  - Convert to daily using splash code
  - Run splashf including sunlight hours
  - Run CRWE to get open water evaporation
- Extract/resample to model DEM

To do :

- ------- UPDATE SWAM CODE --------
- Check and modify dewpoint temperature calculations

## Methods

  1. Run/modify 'getClimDailyGrid.R' to extract cru data and calculate daily forcings PRE, PET, EVP and CN. [Condensation is not yet used]
  2. Run/modify 'reGrid.R' to resample forcing data onto the DEM used for the hydrological model
  3. Run/modify 'runSwam*.R' to calculate the hydrological flow