### Field Enclosure 2017 [FEn17]
Scripts associated with an enclosure experiment looking at food-web changes based on mussel presence/community.
50 enclosures placed in the Kiamichi River: 10 Amblema plicata dominated alive, 10 Actinonaias ligamentina dominated alive, 10 Amblema plicata dominated shams, 10 Actinaias ligamentina dominated shams, and 10 true controls. 
X species dominated communities had 7 dominants and 3 sub-ordinant individuals.

##Order to run the scripts
#Slurry Analysis
Script that loads basic data, calculates Chlorophyl density and AFDM at each enclosure, results in a master matrix called basalres that contains chlorophyl and AFDM measures for each enclosure at each week

#spatialanlaysis
Script that creates a raster that shows an approximation of the relative position of each enclosure. 
Obvious AFDM spatial autocorrelation
TODO: bring in basalres to map these variables as well

#ndsanalysis
Based on code by Dr. Thomas Parr, used to analyze the nutrient diffusing substrata placed within the river during our experiment.

#Inverts
Takes raw data and produces counts, sizes, and biomass densities. Need to use densities because sampling was not consistent (basket recovery was not perfect)

#InvFunction
preliminary script to do functional anlaysis for the invertebrates 
Traci is working on this