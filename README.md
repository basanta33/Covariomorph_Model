# A covarion model for phylogenetic estimation using discrete morphological datasets
authors: Basanta Khakurel and Sebastian Höhna

# Summary of the study
In this manuscript, we present and original study exploring the usage of covarion model in discrete morphological datasets for phylogenetic tree estimation.

This compressed file archive contains all the data and scripts used for the simulation-related and empirical phylogenetic analyses in RevBayes.


## Funding Sources
This work was supported by the European Union (ERC, MacDrive, GA 101043187).
Views and opinions expressed are however those of the authors only and do not necessarily reflect those of the European Union or the European Research Council Executive Agency. Neither the European Union nor the granting authority can be held responsible for them.


# Description of data files and file structure. Folders are indicated by an asterisk.

`Khakurel_Höhna_electronic_supp_mat.pdf`: this file includes:
- Introductory figures and sections for the Covariomorph model.
- References to the original studies from which empirical datasets are obtained
- Results not included in the main manuscript.
The following supplementary figures and tables can be found in the electronic supplementary material:
	Figure S1: Graphical representation of the Covariomorph model
	Figure S2: Results of the validation analysis for the model implementation.
	Figure S3: Expansion of virtual states from observed states
	Figure S4: Rate categories from a discretized lognormal distribution
	Figure S5: Size and parameter estimates of the empirical datasets
	Figure S6: Marginal likelihood estimates of Rays and Sharks dataset with varying number of covariomorph rate categories
	Figure S7: Robinson-Foulds difference between the Maximum a posteriori tree from the Mk model and covariomorph model with varying rate categories
	Figure S8: Differences between trees sampled by different models for the Rays dataset
	Figure S9: Differences between trees sampled by different models for the Sharks dataset

	Table S1: Bayes factor for covariomorph model varying number of rate categories.

* `empirical_data` folder: contains 164 empirical matrices that are used to test the covariomorph model in this study. The source studies are cited in the electronic supplementary materials.

* `scripts` folder: contains the scripts used for phylogenetic analyses using the soured data files described above.
Please refer to our analysis pipeline for proper use of the following scripts:

	* `RevScripts` folder: contains scripts for use with RevBayes
		- `Marginal_Likelihood.Rev`: script to run the marginal likelihood estimation in RevBayes.
		- `mcmc.Rev`: script containing the analysis settings and Markov chain Monte Carlo (MCMC) for covariomorph model in RevBayes
		- `model_Covariomorph.Rev`: general model settings for the covariomorph model
		- `model_Tree.Rev`: script containing the setting for prior on the topology
		- `simulate_Covariomorph.Rev`: script to simulate datasets under the covariomorph model

	* `RScripts` folder: contains the scripts to summarize the output from RevBayes and scripts to plot the results.
		- `convert_state_labels.r`: script to convert the state labels for the simulated data from virtual states to observed states
		- `figure3_simulations_sd_vs_sr.r`: script to generate Figure 3 from the manuscript. This script plots the standard deviation against the switching rate for simulation results.
		- `figure4_empirical_sd_vs_sr.r`: script to generate Figure 4 from the manuscript. This script plots the standard deviation against the switching rate for empirical results.
		- `figure5_posterior_rays_sharks.r`: script to generate Figure 5 from the manuscript. This script plots the posterior distribution of the analyzed datasets.
		- `figure6_cpp_rays_sharks.r`: script to generate Figure 6 from the manuscript. This script calculates the posterior probabilities of the clades and plots them.
		- `summarize_empirical_covariomorph.r`: script to obtain results from the empirical analysis using the covariomorph model
		- `summarize_simulations.r`: script to obtain results from the simulated analyses.

# Analysis pipeline

We assume you have RevBayes installed and it can be used with `rb` in command line.
Note that you can use the MPI version of RevBayes instead.

To simulate datasets using the covariomorph model:
 > rb --file scripts/RevScripts/simulate_Covariomorph.Rev --args ${number_of_datasets} ${number_of_characters} ${output_directory} ${switching_rate} ${number_of_character_states} ${number_of_rate_categories}

 After you simulate, you need to convert the state labels from virtual expanded states to observed states. You can do so be using:
 > Rscript scripts/RScripts/convert_state_labels.r -d "output_directory" -n "number_of_datasets" -s "number_of_character_states" -c "number_of_rate_categories"


To use the covariomorph model with a dataset:
 > rb --file scripts/RevScripts/mcmc.Rev --args ${dataset_name} ${data_directory} Covariomorph ${output_directory} ${number_of_mcmc_chains} ${number_of_rate_categories}


To do a marginal likelihood estimation analysis:
 >rb --file scripts/RevScripts/Marginal_Likelihood.Rev --args ${dataset_name} ${data_directory} Covariomorph ${output_directory} ${number_of_mcmc_chains} ${number_of_rate_categories}


To summarize results from the posterior of the empirical analyses:
(this code assumes that you have a 4 rate category Covariomorph model and a list of data files that are in `completed_list.txt`)
 > Rscript scripts/RScripts/summarize_empirical_covariomorph.r -f completed_list.txt -d "output_directory_from_mcmc" -m Covariomorph


To summarize results from the posterior of the simulated analyses:
> Rscript scripts/RScripts/summarize_simulations_median.r -d "output_directory_from_mcmc" -n "number_of_datasets_used"


To generate figures:
 > Rscript scripts/RScripts/figure3_simulations_sd_vs_sr.r
 > Rscript scripts/RScripts/figure4_empirical_sd_vs_sr.r
 > Rscript scripts/RScripts/figure5_posterior_sharks_rays.r
 > Rscript scripts/RScripts/figure6_cpp_sharks_rays.r
