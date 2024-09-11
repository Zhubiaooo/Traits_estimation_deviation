The following files allow one to reproduce analyses in the manuscript entitled "Traits estimated when grown alone may underestimate competitive advantage and invasiveness of exotic species" by Biao Zhu & Wei Chen.

DATA & FILE OVERVIEW

In *Data* folder
1)  iq_tree.treefile
2)  Common_species_list.xlsx
3)  Pot_traits_database0831.xlsx
4)  Pot_traits_mean0831.xlsx
5)  Field_traits_database0831_2.xlsx
6)  Field_traits_mean0831_2.xlsx
7)  all_row_data0829.xlsx
8)  Imputation_data0831_2.xlsx

In *Code* folder
1)  Fig 1-0831.R
2)  Fig 2-0831.R
3)  Fig 3 & Fig S9 & Fig S11-0831.R
4)  Fig 4-0831.R
5)  Fig 5-0831 & Fig S3-0831.R
6)  Fig S1-0831.R
7)  Fig S2-0831.R
8)  Fig S4-0831.R
9)  Fig S5-0831.R
10)  Fig S6-0831.R
11)  Fig S7-0831.R
12)  Fig S8-0831.R
13)  Fig S10-0831.R
14)  Traits imputation.R

#########################################################################

DATA-SPECIFIC INFORMATION FOR: iq_tree.treefile

Description: Phylogenetic relationship of all experimental species in pot experiment and field experiment

#########################################################################

DATA-SPECIFIC INFORMATION FOR: Common_species_list.xlsx

1. Number of variables: 8
2. Number of cases/rows: 65
3. Variable List:
    * Species: Latin name of species
    * Origin: Geographical origin of species
    * Genus: Taxonomic level generic names
    * Family: Taxonomic level family name
    * Lifeform: Lifeform of species
    * Species_AGB: The list of species used in comparing Hmax and AGB under different types of experiments
    * Species_SLA: The list of species used in comparing SLA under different types of experiments
    * Invasive: invasive vs. not invasive

#########################################################################

DATA-SPECIFIC INFORMATION FOR: Pot_traits_database0831.xlsx

1. Number of variables: 13
2. Number of cases/rows: 320
3. Variable List:
    * species_num: Species number 
    * Species: Latin name of species
    * Repeats: Repetition number of experimental species
    * Pots_num: Number of experimental pot
    * SLA: Specific leaf area (cm2/g)
    * Hmax_time1:  First measurement of plant height (cm)
    * Hmax_time2:  Second measurement of plant height (cm)
    * Hmax:  Maximum height (cm)
    * AGB: Aboveground biomass (g)
    * Origin:  Geographical origin of species
    * Genus: Taxonomic level generic names
    * Family: Taxonomic level family name
    * Lifeform: Lifeform of species

#########################################################################

DATA-SPECIFIC INFORMATION FOR: Pot_traits_mean0831.xlsx

Several metrics obtained from the community information at each unique sampling event.

1. Number of variables: 6
2. Number of cases/rows: 65
3. Variable List:
    * Species: Latin name of species
    * SLA: Specific leaf area (cm2/g)
    * Hmax_time1:  First measurement of plant height (cm)
    * Hmax_time2:  Second measurement of plant height (cm)
    * Hmax:  Maximum height (cm)
    * AGB: Aboveground biomass (g)

#########################################################################

DATA-SPECIFIC INFORMATION FOR: Field_traits_database0831_2.xlsx

1. Number of variables: 7
2. Number of cases/rows: 439
3. Variable List:
    * Block: Plot code of experimental 
    * Species: Latin name of species
    * Origin:  Geographical origin of species
    * Lifeform: Lifeform of species
    * SLA: Specific leaf area (cm2/g)
    * Hmax:  Maximum height (cm)
    * AGB: Aboveground biomass (g)

#########################################################################

DATA-SPECIFIC INFORMATION FOR: Field_traits_mean0831_2.xlsx

1. Number of variables: 6
2. Number of cases/rows: 65
3. Variable List:
    * Species: Latin name of species
    * SLA: Specific leaf area (cm2/g)
    * SLA_imp: Specific leaf area (cm2/g) inferred by missForest
    * Hmax:  Maximum height (cm)
    * AGB: Aboveground biomass (g)
    * Origin:  Geographical origin of species

#########################################################################

DATA-SPECIFIC INFORMATION FOR: Imputation_data0831_2.xlsx

1. Number of variables: 4
2. Number of cases/rows: 89
3. Variable List:
    * Species: Latin name of species
    * SLA: Specific leaf area (cm2/g)
    * Hmax:  Maximum height (cm)
    * AGB: Aboveground biomass (g)

#########################################################################

DATA-SPECIFIC INFORMATION FOR: all_row_data0829.xlsx

1. Number of variables: 16
2. Number of cases/rows: 405
3. Variable List:
    * Plot_num: Subplot code of field exp.
    * Block: Plot code of field exp.
    * Species_num:  Species number
    * Species: Latin name of species
    * Seed_source: Species population of field exp.
    * Origin:  Geographical origin of species
    * Lifeform: Lifeform of species
    * rebio2021: Relative biomass of each species in each plot in the first year of the field experiment
    * rebio2022: Relative biomass of each species in each plot in the second year of the field experiment

#########################################################################


CODE INFORMATION

1)  Fig 1-0831.R

This script contains all the data analysis and visualization code of Fig 1.

2)  Fig 2-0831.R

This script contains all the data analysis and visualization code of Fig 2.

3)  Fig 3 & Fig S9 & Fig S11-0831.R

This script contains all the data analysis and visualization code of Fig 3 & Fig S9 & Fig S11.

4)  Fig 4-0831.R

This script contains all the data analysis and visualization code of Fig 3 & Fig S9 & Fig S11.

5)  Fig 5-0831 & Fig S3-0831.R

This script contains all the data analysis and visualization code of Fig 5 $ Fig S3.

6)  Fig S1-0831.R

This script contains all the data analysis and visualization code of Fig S1.

7)  Fig S2-0831.R

This script contains all the data analysis and visualization code of Fig S2.

8)  Fig S4-0831.R

This script contains all the data analysis and visualization code of Fig S4.

9)  Fig S5-0831.R

This script contains all the data analysis and visualization code of Fig S5.

10)  Fig S6-0831.R

This script contains all the data analysis and visualization code of Fig S6.

11)  Fig S7-0831.R

This script contains all the data analysis and visualization code of Fig S7.

12)  Fig S8-0831.R

This script contains all the data analysis and visualization code of Fig S8.

13)  Fig S10-0831.R

This script contains all the data analysis and visualization code of Fig S10.

14)  Traits imputation.R

This script is used to interpolate and infer missing specific leaf area data from field experiments.
