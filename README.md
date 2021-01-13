# Field Enclosure 2017 [FEn17]
Scripts associated with an enclosure experiment looking at food-web changes based on mussel presence/community.
50 enclosures placed in the Kiamichi River: 10 Amblema plicata dominated alive, 10 Actinonaias ligamentina dominated alive, 10 Amblema plicata dominated shams, 10 Actinaias ligamentina dominated shams, and 10 true controls. 
X species dominated communities had 7 dominants and 3 sub-ordinant individuals.

##Order to run the scripts

###Slurry Analysis.R
Script that loads basic data, calculates Chlorophyl density and AFDM at each enclosure, results in a master matrix called basalres that contains chlorophyl and AFDM measures for each enclosure at each week.

###spatialanlaysis.R
Script that creates a raster that shows an approximation of the relative position of each enclosure. Also creates a raster map of different physiochemical characteristics and primary production.

###Inverts.R
Takes raw data and produces counts, sizes, and biomass densities. Need to use densities because sampling was not consistent (basket recovery was not perfect)

###NDS_analysis.R
Based on code by Dr. Thomas Parr, used to analyze the nutrient diffusing substrata placed within the river during our experiment.

###Stoich.R
Script to evaluate potential stoichiometry differences in the periphyton samples on tiles in each enclosure. 

###Model Script.R

##Products from this project
Hopper, G. W., DuBose, T. P., Gido, K. B., & Vaughn, C. C. (2019). Freshwater mussels alter fish distributions through habitat modifications at fine spatial scales. Freshwater Science, 38(4), 702-712.