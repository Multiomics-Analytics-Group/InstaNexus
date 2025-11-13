# Notes reproducibility

Remind user to cd into src
- better: do not use relative paths, somehow

## basic commands for activate and deactivate an environment
```bash
conda activate instanexus
conda env update -f environment.yml
conda list
conda info --envs
conda deactivate
conda env remove --name instanexus
```

```bash
eval "$(/home/psq/miniconda3/bin/conda shell.zsh hook)"
conda env create -f environment.linux.yml
conda activate instanexus
```

navigate to: https://github.com/Multiomics-Analytics-Group/InstaNexus/wiki/Reproducibility-test


## Update your environment after updating the file with new libraries/modules

we added biopython to yml file

```bash
conda env update --file environment.linux.yml --prune
```

## Update the main from a different branch anche check what has been changed

```bash
git branch -r
git log origin/tests-pasdom --oneline
git checkout main
git fetch origin
git merge origin/tests-pasdom
git push origin main
```
1. You first listed all remote branches and inspected the commit history of the remote branch `origin/tests-pasdom` using `git branch -r` and `git log`.
2. Then, you switched to your local `main` branch, fetched the latest updates from the remote repository, and merged the changes from `origin/tests-pasdom` into `main`.
3. Finally, you pushed the updated `main` branch (now including the merged changes) to the remote repository.
