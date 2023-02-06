# Supplementary to 
### *In Silico Engineering of Ion-Exchanged Zeolites for High-Performance Carbon Capture in Psa Processes*

*Zijun Deng, Arun Gopalan, Lev Sarkisov*


This repository contains some of the pre and post-processing scripts and simulation inputs to the article *In Silico Engineering of Ion-Exchanged Zeolites for High-Performance Carbon Capture in Psa Processes* by Zijun Deng, Arun Gopalan, Lev Sarkisov.2023. 

You can find the publication here: 

[In Silico Engineering of Ion-Exchanged Zeolites for High-Performance Carbon Capture in PSA Processes](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4306028)


Note:
-----
The scripts are provided as-is, and are not guaranteed to work. If you have any questions, please contact the corresponding author of the article.

##	Generate random Al distribution
Requirement: Jupyter Notebook, Python
1.	Copy the file “LTA_my.cif” to the folder “C:\\LTA\\”
2.	Change the variable i (to determine the number of Al to be replaced by Si script) in the script “Generate random Al distribution.ipynb”.
3.	Run the script “Generate random Al distribution.ipynb” using Jupyter Notebook, pick the Al distribution with the closest Warren-Cawley parameter to the average.
[Generate Al distribution notebook](./generate_random_Al_dist.ipynb)
4.	Override the atom distribution content in “LTA_my.cif” with the Al distribution content picked.

##	Equilibrium cation distribution using parallel tempering
Requirement: RASPA
1.	Copy the file “LTA_my.cif” in 1.4 to the folder “Restart”
2.	Change the value of “CreateNumberOfMolecules” (to determine the number of Na and K in a unit cell) in the script “Restart\\ simulation.input”.
3.	Run the simulation using RASPA,

##	Identify cation distribution

Please refer to the [Cation distribution notebook](./cation_distribution.ipynb)



##	Predict adsorption uptake using GCMC simulation
Requirement: RASPA
1.	Copy the file “LTA_my.cif” in 1.4 to the folders “co2” and “n2”.
2.	Copy the file “Restart\\Restart\\System_0\\ restart_LTA_my_1.1.1_303.000000_0” to the folders “co2\\RestartInitial\\System_0\\” and “n2\\RestartInitial\\System_0\\”; change the name of the file in both folders to “restart_LTA_my_1.1.1_298.150000_100000”.
3.	Run the simulations using RASPA,

##	Obtain DSL model
Requirement: Origin LAB
1.	Collect the adsorption uptakes predicted in 4.3, and generate adsorption isotherms at 263.15K, 298.15K, 308.15K, and 323.15K.
2.	Fit these data points in a DSL model using Origin LAB.

6.	Process simulation and optimization
Requirement: MATLAB
1.	Override the DSL parameters the row 16 of the table “PSA\\Params.mat\\IsothermPar“ with your DSL model obtained in 5.2. 
2.	Run the script “PSA\\run_FullOptimization.sh” using MATLAB,

