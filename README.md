
Stoma Simulator
==

Biomechanical simulations of plant stomata.

The following instructions have been tested on macOS 10.12.5.

### Python environment

To run the StomaSimulator software you will need Python 2.7.
You will also need to install the project's dependencies.

The following instructions use [Miniconda](https://conda.io/miniconda.html) although a 
Python "virtual environment" could be used instead.
Install Miniconda by downloading and running the installer.
To setup the environment issue the following commands:

```bash
conda create -n stomaenv --yes python=2.7       # create Python 2.7 environment
source activate stomaenv                        # activate env.
```

### Set up `stomasimulator`

Clone this repository (if you haven't already) and 'cd' into its directory
Run the following to setup `stomasimulator`:
```bash
conda install --file requirements.txt --yes     # install dependencies
python setup.py install                         # to install stomasimulator
```

To uninstall `stomasimulator` remove the miniconda environment using
```bash
conda remove -n stomaenv --all --yes
```

### Finite element software

This step is only required if you are going to rerun the simulations.

Download and install the [FEBio](http://www.febio.org) finite element software (v2.3.0 was used).

Update the path to the executable in the `stomasimulator.yml` file at the top level of the `stomasimulator` package.
You will have to run `python setup.py install` again. 

See [Maas *et al.* (2012)](http://dx.doi.org/10.1115/1.4005694) for details of FEBio.
 

### Running `stomasimulator`

Activate the environment. Typing `stomasimulator -h` will show:
```bash
$ stomasimulator -h
************************************************************************************************************************
Reading configuration file '<path/to/repo>/stomasimulator/stomasimulator/stomasimulator.yml'...
Finished reading the configuration file

Configuration settings for stomasimulator...
FEBio executable     : <path/to/febio>/FEBio2
Model data directory : <path/to/repo>/stomasimulator/stomasimulator/.cfg
************************************************************************************************************************
usage: stomasimulator <command> [<args>]

The commands are:
  cfg          Output the configuration settings and model labels used by stomasimulator
  feb          Output a FEB file for the specified stoma model (the simulation can be run)
  opt          Optimise the material parameters for a given stoma model
  xplt         Process an XPLT file

Execute functions within the stomata biomechanics suite

positional arguments:
  command     Subcommand to run

optional arguments:
  -h, --help  show this help message and exit
```
where `<path/to/repo>` will show your path.
These commands each have the `-h` option to show help.

### Display the graphs

The graph plotting assumes you have [gnuplot](http://www.gnuplot.info/) installed.

In a terminal, go to the `figures/plotting` directory.

Start `gnuplot` and run the scripts using, for example, `call 'plot-fig1.gp'`. 

### Rerun simulations

Rerunning the simulations will take several hours on one machine.

To rerun the simulations, go to the `figures/simulations` directory.

Run `./create-figures.sh` to recreate all of the simulations.

 
 
 
 
 
 
 

