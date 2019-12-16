# PharMD - extraction of pharmacophores from MD trajectories and virtual screening

PharMD is a tool to retrieve pharmacophore models from MD trajectories of protein-ligand complexes, identification of redundant pharmacophores and virtual screening with multiple pharmacophore models using different scoring schemes. 

## Dependency

`mdtraj >= 1.9.3`  
`plip >= 1.4.2`  
`pmapper >= 0.3.1`  
`psearch >= 0.0.2`

## Installation
```text
pip install pharmd
```

## Usage

### Retrieve pharmacophores from an MD trajectory

To retrieve individual snapshots of MD trajectory `mdtraj` package is used. 
Therefore the `md2pharm` utility takes the same arguments as `mdconvert` utility from `mdtraj`. 
Thus you may extract only specified frames not all of them. 
You have to specify ligand code as it is given in PDB topology file.
Individual frames will be stored in a single PDB file without solvent molecules.
Pharmacophore models for each frame in xyz-format will be stored in the same directory as output pdb-file. 

```bash
md2pharm -i md.xtc -t md.pdb -s 10 -g LIG -o pharmacophores/frames.pdb
```

### Retrieve non-redundant pharmacophores

Similar pharmacophores are recognized by identical 3D pharmacophore hashes. 
It is expected that pharmacophores with identical hashes would have RMSD less than the specified binning step.
By default binning step equals to 1A.
Pharmacophores with distinct hashes are stored in a specified directory. Optionally one may provide a path where to store hashes for al pharmacophores.   

```bash
get_distinct -i pharmacophores/ -o distinct_pharmacophores/
```

### Perform virtual screening using multiple non-redundant pharmacophores

`screen_db` utility from `psearch` package is used for this purpose.
Therefore you have to generate database of compound conformers and their pharmacophore representations using utilities from `psearch` package. 
At this step you may specify a desired binning step value which will be used further in screening (default is 1).

```bash
prepare_db -i input.smi -o compounds.db -c 2 -v 
```

If you would like to calculate scoring based on Conformer Coverage Approach you have to specify `--conf` argument for `screen_db`. 
Then all conformers of a compound matching pharmacophore models will be retrieved as hits (may be slower). 
Otherwise only the first matching conformer will be returned.

It is recommended to restrict screening to complex pharmacophores having at least four features, because less complex models would retrieve many irrelevant compounds.

```bash
screen_db -i compounds.db -q distinct_pharmacophores/ -o screen/ --conf -c 2 -f 4
```

Multiple txt-files will be created in the output directory containing hit lists retrieved by individual pharmacophore models.

### Calculate compound scores based on multiple hit lists

The advantage of ensemble scoring is that you do not need validate individual models and select best performing ones.
Ensemble scoring is calculated by:   
1. Conformer Coverage Approach (CCA) - the score is equal to the percentage of conformers matching at least one of supplied pharmacophore models.
2. Common HIts Approach (CHA) - the score is equal to the percentage of models matched at least one conformer of a compound.

In the case of CCA scoring you have to supply the database of screened compounds as an additional parameter.
```bash
get_scores -i screen/ -o cca_scores.txt -s cca -d compounds.db
```
 
## Documentation
All utilities have `-h` option to get help pages with descriptions of all available arguments. 


## Citation
Virtual Screening Using Pharmacophore Models Retrieved from Molecular Dynamic Simulations  
Pavel Polishchuk, Alina Kutlushina, Dayana Bashirova, Olena Mokshyna, Timur Madzhidov  
*Int. J. Mol. Sci.* **2019**, 20(23), 5834  
https://doi.org/10.3390/ijms20235834


## Issues
Currently there is an issue with installation of dependencies. `plip` requires `openbabel` which causes an error during installation via `pip`. Therefore it is recommended to solve dependencies manually and use `pip install -U --no-deps pharmd` to install `pharmd` ignoring dependencies.

## License
BSD-3 clause
