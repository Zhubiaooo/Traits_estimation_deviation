The following files allow one to reproduce analyses in the manuscript entitled "Traits estimated when grown alone may underestimate competitive advantage and invasiveness of exotic species" by Biao Zhu, Chunqiang Wei, Hao Zhou and Xinmin Lu.

DATA & FILE OVERVIEW

In *Data* folder
1)  All_species.newick
2)  Total_species_list.xlsx
3)  Common_species_list.xlsx
4)  Pot_traits_database.xlsx
5)  Pot_traits_mean.xlsx
6)  Field_traits_database(SLA).xlsx
7)  Field_traits_database(Hmax&AGB).xlsx
8)  Imputation_data_in_field.xlsx
9)  Field_traits_mean.xlsx
10) Survival_database_in_field.xlsx
11) Field_composition_database.xlsx

In *Code* folder
1)  Fig 1.R
2)  Fig 2 & Fig 4 & Fig S2.R
3)  Fig 3 & Fig S1 & Fig S3 & Fig S4.R
4)  Imputation missing valus in field experiment


#########################################################################

DATA-SPECIFIC INFORMATION FOR: All_species.newick

Description: Phylogenetic relationship of all experimental species in pot experiment and field experiment (N = 172)

#########################################################################

DATA-SPECIFIC INFORMATION FOR: Total_species_list.xlsx

1. Number of variables: 4
2. Number of cases/rows: 173
3. Variable List:  
    * Species: Latin name of species
    * Species in KEW: Latin name of species based on Plants of the World Online 
    * Origin: Geographical origin of species
    * Lifeform: Lifeform of species

#########################################################################

DATA-SPECIFIC INFORMATION FOR: Common_species_list.xlsx

1. Number of variables: 7
2. Number of cases/rows: 65
3. Variable List:
    * Species: Latin name of species
    * Origin: Geographical origin of species
    * Genus: Taxonomic level generic names
    * Family: Taxonomic level family name
    * Lifeform: Lifeform of species
    * Species_AGB: The list of species used in comparing Hmax and AGB under different types of experiments
    * Species_SLA: The list of species used in comparing SLA under different types of experiments

#########################################################################

DATA-SPECIFIC INFORMATION FOR: Pot_traits_database.xlsx

1. Number of variables: 9
2. Number of cases/rows: 321
3. Variable List:
    * Pots_num: Number of experimental pot
    * Sp_code: Species number
    * Repeats: Repetition number of experimental species
    * Species: Latin name of species
    * Origin:  Geographical origin of species
    * Lifeform: Lifeform of species
    * SLA: Specific leaf area (cm2/g)
    * Hmax:  Maximum height (cm)
    * AGB: Aboveground biomass (g)

#########################################################################

DATA-SPECIFIC INFORMATION FOR: Pot_traits_mean.xlsx

Several metrics obtained from the community information at each unique sampling event.

1. Number of variables: 4
2. Number of cases/rows: 65
3. Variable List:
    * Species: Latin name of species
    * SLA: Specific leaf area (cm2/g)
    * Hmax:  Maximum height (cm)
    * AGB: Aboveground biomass (g)

#########################################################################

DATA-SPECIFIC INFORMATION FOR: Field_traits_database(SLA).xlsx

1. Number of variables: 4
2. Number of cases/rows: 254
3. Variable List:
    * Species: Latin name of species
    * SLA: Specific leaf area (cm2/g)
    * Origin:  Geographical origin of species
    * Lifeform: Lifeform of species

#########################################################################

DATA-SPECIFIC INFORMATION FOR: Field_traits_database(Hmax&AGB).xlsx

1. Number of variables: 8
2. Number of cases/rows: 405
3. Variable List:
    * Plot_num: Number of experimental plot
    * Block: Number of experimental block
    * Sp_num:  Species number
    * Species: Latin name of species
    * Field_Hmax:  Maximum height (cm)
    * Field_AGB: Aboveground biomass (g)
    * Origin:  Geographical origin of species
    * Lifeform: Lifeform of species

#########################################################################

DATA-SPECIFIC INFORMATION FOR: Imputation_data_in_field.xlsx

1. Number of variables: 4
2. Number of cases/rows: 89
3. Variable List:
    * Species: Latin name of species
    * SLA: Specific leaf area (cm2/g)
    * Hmax:  Maximum height (cm)
    * AGB: Aboveground biomass (g)

#########################################################################

DATA-SPECIFIC INFORMATION FOR: Field_traits_mean.xlsx

1. Number of variables: 5
2. Number of cases/rows: 65
3. Variable List:
    * Species: Latin name of species
    * SLA: Specific leaf area (cm2/g)
    * SLA_imp: Specific leaf area (cm2/g) inferred by missForest
    * Hmax:  Maximum height (cm)
    * AGB: Aboveground biomass (g)

#########################################################################

DATA-SPECIFIC INFORMATION FOR: Survival_database_in_field.xlsx

1. Number of variables: 4
2. Number of cases/rows: 123
3. Variable List:
    * Species: Latin name of species
    * Origin:  Geographical origin of species
    * Lifeform: Lifeform of species
    * Survival: Survival rate of species in the first year in field experiments

#########################################################################

DATA-SPECIFIC INFORMATION FOR: Field_composition_database.xlsx

1. Number of variables: 16
2. Number of cases/rows: 405
3. Variable List:
    * Plot_num: Number of experimental plot
    * Block: Number of experimental block
    * Sp_num:  Species number
    * Species: Latin name of species
    * Field_SLA: Specific leaf area (cm2/g) measured in the field
    * Field_SLA_imp: Specific leaf area (cm2/g) inferred by missForest measured in the field
    * Field_Hmax:  Maximum height (cm) measured in the field
    * Field_AGB: Aboveground biomass (g) measured in the field
    * 2019_AGB: Total biomass of species in each sample plot in the first year of field experiment
    * 2019_total: Total biomass of each plot in the first year of field experiment
    * rebio2019: Relative biomass of each species in each plot in the first year of the field experiment
    * 2020_AGB: Total biomass of species in each sample plot in the second year of field experiment
    * 2020_total: Total biomass of each plot in the second  year of field experiment
    * rebio2020: Relative biomass of each species in each plot in the second year of the field experiment
    * Origin:  Geographical origin of species
    * Lifeform: Lifeform of species

#########################################################################


CODE INFORMATION

1)  Fig 1.R

This script contains all the data analysis and visualization code of Fig.1.

2)  Fig 2 & Fig 4 & Fig S2.R

This script contains all the data analysis and visualization code of Fig 2 & Fig 4 & Fig S2.

3)  Fig 3 & Fig S1 & Fig S3 & Fig S4.R

This script contains all the data analysis and visualization code of Fig 3 & Fig S1 & Fig S3 & Fig S4.

4)  Imputation missing valus in field experiment.R

This script is used to interpolate and infer missing specific leaf area data from field experiments.
