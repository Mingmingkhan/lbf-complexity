# lbf-complexity
Code for Khan et al. (2025) Metacommunity structural changes of Antarctic benthic invertebrates over the late Maastrichtian 

# Abstract
Seymour (Marambio) Island, Antarctica has one of the most expanded onshore Cretaceous-Paleogene sedimentary successions in the world. The deposition of the López de Bertodano Formation (~70-65.6 Ma) covered a time of fluctuating sea temperatures, including cold snaps, and warming linked to Deccan Traps volcanism. Here, we study community dynamics of latest Cretaceous (Maastrichtian) Antarctic invertebrates using fossils from the Zinsmeister Collection, Paleontological Research Institution, USA, in order to assess ecological complexity prior to the Cretaceous-Paleogene (K-Pg) mass extinction. Our data set consists of 7400 fossils from 85 genera across bivalves, gastropods, cephalopods, echinoderms, brachiopods, scaphopods, polychaetes and octocorals, from 324 localities spanning six informal sub-units, KLBs 5-9. Due to positional uncertainty of the KLB boundaries, we performed sensitivity analyses to ensure robust results. We found that the number of significantly non-random taxonomic co-occurrences and complexity increased throughout this period. To investigate metacommunity structure that may arise from taxa interactions or environmental filtering, we used the Elements of Metacommunity Structure framework, where we found that taxa replacement, rather than nestedness, increased through time, also highlighting complexity. However, our sensitivity analyses found that our metacommunity results could not be distinguished from sampling biases in the most conservative sensitivity test. Thus, whilst we can be confident that ecological complexity increased throughout the Maastrichtian, the detailed community mechanisms behind this increase cannot be firmly established; nonetheless, this result reinforces the presence of a single, rather than two-fold, K-Pg extinction in the southern high latitudes.

# Using the repository 

1) Site-by-taxa occurrences of the LBF fossils can be downloaded from the UK Polar Data Centre: Khan, T.M., Whittle, R.J., Witts, J.D., Griffiths, H.J., Manica, A., & Mitchell, E.G. (2025). Invertebrate fossil occurrences from the Cretaceous Lopez de Bertodano Formation (Maastrichtian), Seymour Island, housed in the Paleontological Research Institution (PRI), Ithaca, NY (Version 1.0) [Data set]. NERC EDS UK Polar Data Centre. https://doi.org/10.5285/c1252fb5-8075-431c-a78b-1c1d8da37c5e
   
2) The above file is also present in this repository.

3) R scripts present:
      01-Analyses_KLB7.R - creates the occurrence matrices for KLB 7 under sensitivity treatments and does co-occurrence and metacommunity analyses.
      02-Analyses_KLB8.R - creates the occurrence matrices for KLB 8 under sensitivity treatments and does co-occurrence and metacommunity analyses.
      03-Analyses_KLB9.R - creates the occurrence matrices for KLB 9 under sensitivity treatments and does co-occurrence and metacommunity analyses.
      04-Ordination.R - visualization of the reciprocally averaged occurrence data for use in metacommunity analyses, and a comparison with NMDS analyses.
      05-all_cooccur.R - boxplots of the significant co-occurrences based on sensitivity treated data.
      06-all_metacom_plots.R - plots the results of the metacommunity analyses.
      07-coverage_thinning.R - assesses rarefaction and coverage based rarefaction; spatially explicit subsampling; plots results.

Note that all figures have been modified in Inkscape. 

   
   



