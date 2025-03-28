# **Apex scavenger declines have cascading effects on soil biogeochemistry**

### **Corresponding Author:** Torrey Stephenson (torreys@uidaho.edu), University of Idaho
### **Co-Authors:** Savannah Bartel, Ernest Osburn, Michael Strickland, David W. Crowder, Menna Jones, Kawinwit Kittipalawattanapol, Calum X. Cunningham, Tara Hudiburg, Andrew Storfer, Julia Piaskowski, Laurel M. Lynch

### **Abstract**
Decomposition of animal carcasses provides unique and ephemeral resources that can cycle through ecosystems more efficiently than plant litter. Despite previous work demonstrating how the loss of apex consumers can affect carcass decomposition rates across landscapes and ecosystem productivity, the consequences of trophic downgrading for terrestrial biogeochemistry remain understudied. Using manipulative field experiments conducted across a natural gradient in Tasmanian devil density, we assessed how declines in apex scavenger populations influenced carcass persistence, soil biogeochemistry, and microbial ecology in summer and winter. Devil declines extended carcass persistence, concentrating carcass-derived nutrient input into soil, reducing soil microbiome diversity, and selecting for faster-growing taxa and pathogens. In the winter, lower invertebrate and microbial activity reduced carcass consumption regardless of devil populations, dampening soil biogeochemical responses. Further, while less efficient scavenger guilds clearly facilitated carcass consumption, they did not fill the functional role of apex scavengers. Our study highlights the unique role apex scavengers play in maintaining ecosystem function by linking decomposition dynamics with broader ecological consequences including plant productivity, carbon sequestration, and pathogen risk.


### **Data Description:**
Soil chemistry and microbial community data is from the top 7 cm of mineral soil collected from field sites across Tasmania in August-September 2022 ("Winter") or February-March 2023 ("Summer"). Complete experimental design, sampling, and analysis methods are described in the published manuscript.

Names and units for soil chemistry data (SummerCarcassData.xlsx) are as follows:

Most\* abbreviations/units used in this file are consistent in other data files, so only this file has column-specific detail

A: year of sampling. 2023 = "summer", 2022 = "winter"
B: sample number.
C: trueID. Internal lab shorthand for keeping track of samples
D: site. G = Low Devils, A = Medium Devils, T = High Devils, B = No Devils
E: block. experimental block/replicate
F: trt1. "Carcass Access" - whether carcass was in an exclusion cage (E) or staked with no cage (S)
G: trt2. "Location" - whether the soil sample was taken from underneath the carcass treatment (T), or the paired control (C)
H: time. Days after carcass placement that the sample was taken. 0 = day of placement, 1 = +5 days, 2 = +10 days, 3 = +15 days, 4 = +30 days.
I: dark. Indicates whether visible color was present in the soil extract; relevant for soil PO4 determination.
J: day_remove. Day of experiment on which all soft tissue had been removed from the carcass.
K: pct.consumed. The percent of the carcass which was consumed at each sampling point/location. Data are incomplete in this file; this is rectified by the presence of "day removed both yrs.xlsx", which contains complete data.
L: Wm. Soil moisture (%), expressed as a decimal.
M-Q: Ammonium, Nitrate, Phosphate, Dissolved Organic Carbon, Total Dissolved Nitrogen in soil samples, expressed as micrograms per gram of dry soil.
R: EC. Soil electrical conductivity, expressed as microsiemens per cm.
S: pH. Soil pH, expressed as pH units.
T-AA: Relative abundance of various bacterial phyla. "Other" is the sum of all rare phyla. See published methods for details.
AB: rK_ratio. The ratio of "r-selected" to "K-selected" phyla.
AC: 16S_shandiv. Shannon Diversity (H) of the soil bacterial community.
AD-AL: Relative abundance of various fungal phyla. "Other" is the sum of all rare phyla. See published methods for details.
AM: ITS_shandiv. Shannon Diversity (H) of the soil fungal community.

\* "time" column varies from summer to winter. In winter, only 2 sampling periods occurred so the final sampling point, equivalent to '4' in the summer data is labeled as '1'. Only relevant for "timefinal_bothyears.xlsx" file.
