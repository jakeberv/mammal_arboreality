# mammal_arboreality
R code for reconstructing the evolution of mammal arboreality

Ecological selectivity and the evolution of mammalian arboreality
across the K–Pg mass extinction

Running title: Arboreal mammals and the K–Pg extinction

Jonathan J. Hughes(1), Jacob S. Berv(1,2,3), Stephen G. B. Chester(4,5,6), Eric J. Sargis(7,8), Daniel J. Field(9)
1. Department of Ecology & Evolutionary Biology, Cornell University, Ithaca, NY, USA
2. Department of Ecology & Evolutionary Biology, University of Michigan, Ann Arbor, MI, USA 
3. University of Michigan Museum of Paleontology, University of Michigan, Ann Arbor, MI, USA 
4. Department of Anthropology, Brooklyn College, City University of New York, Brooklyn, NY, USA
5. Department of Anthropology, The Graduate Center, City University of New York, New York, NY, USA
6. New York Consortium in Evolutionary Primatology, New York, NY, USA
7. Department of Anthropology, Yale University, New Haven, CT, USA
8. Divisions of Vertebrate Paleontology and Vertebrate Zoology, Yale Peabody Museum of Natural History, New Haven, CT, USA
9. Department of Earth Sciences, University of Cambridge, Cambridge, UK


Meredith_2011.R and Upham_2019.R are respectivley, R scripts for producing analyses of ancestral ecologies using the Meredith et al 2011 and Upham et al 2019 consensus topologies. In the case of Upham et al 2019, additional analyses are performed on a sample of 1000 posterior trees. These primary scripts rely on and source functions contained within rate_through_time.R, helpers.R, and simmap_parallel.R R documents. 

mammal_data_v2 is the dataset used for reconstructions on the Meredith_2011 topology.
mammal_data_v2_upham is the dataset used for reconstructions on the Upham_2019 topology.

Analysis output files are deposited inside analyses_out