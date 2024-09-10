# PACKAGENAME

Variables to replace across repository:

- PACKAGENAME
- PARSERNAME
- AUTHORNAME
- REPONAME

## Installation

```bash
conda create -n PACKAGENAME python=3.10
conda activate PACKAGENAME
git clone --shared https://github.com/AUTHORNAME/REPONAME.git
cd REPONAME
pip install -e . 
```

## Usage

```bash
python -m PACKAGENAME --help
```
