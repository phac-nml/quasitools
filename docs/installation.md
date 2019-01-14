# Installation

## Conda

Conda is an open source package management system and environment management system that allows users to install software quickly and easily on a variety of operating systems.

The simpliest way to install quasitools is with Conda. Install quasitools from [Bioconda](https://bioconda.github.io/) with [Conda](https://conda.io/docs/) ([installation instructions](https://bioconda.github.io/#install-conda)).

Adding Conda channels for the first time:

```
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
```

Installing quasitools with Conda:

```
conda install quasitools
```

If everything was successful, you should be able to verify your installation by running ```quasitools --help``` and see a quasitool's help message.

## Galaxy

Install quasitools from the [main Galaxy toolshed](https://toolshed.g2.bx.psu.edu/view/nml/quasitools/9fb9fed71486).
