conda activate instanexus
conda env update -f environment.yml
conda list
conda info --envs
conda deactivate
conda env remove --name instanexus








eval "$(/home/psq/miniconda3/bin/conda shell.zsh hook)"
conda env create -f environment.linux.yml
conda activate instanexus

navigate to: https://github.com/Multiomics-Analytics-Group/InstaNexus/wiki/Reproducibility-test

# we added biopython to yml file
conda env update --file environment.linux.yml --prune

remind user to cd into src
  - better: do not use relative paths, somehow
