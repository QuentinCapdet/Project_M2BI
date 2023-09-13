# Project_M2BI


## Setup your environment

Install [miniconda](https://docs.conda.io/en/latest/miniconda.html).

Clone the repository

```
git clone https://github.com/QuentinCapdet/Project_M2BI.git
```

Move to the new directory

```
cd Project_M2BI
```

Create the Project_M2BI conda environment:

```
conda env create -f py35.yml
```

Load the `Project_M2BI` conda environment:

```
conda activate Project_M2BI
```
To deactivate an conda active environment, use

```
conda deactivate
```
1) Run the file Projet_access_solvant.py to have all the results on the terminal

```
python Projet_access_solvant.py --pdb /YOUR/PATH/Project_M2BI/PDB/7kh5.pdb
````
2) The file Projet_access_solvant_results.py alows you to have the results in a text file found in the OUTPUTS/ directory.

3) The file RESULTS-PROJECT.ipynb contains grpahs produced in relation to the results obtained (lack of time, no commments to properly eplain the results apart from the last figure which is explained in the written report  !)
