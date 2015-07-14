# simple analysis - A package to read and analyze data saved by LabRad

This is mainly meant to navigate in the forest of directories where LabRad saves the data. 
It has some data fitting capabilities


## Installation:


Required packages: treedict, pylab, numpy

In order to find the files it is best to add an environment variable DV_PATH.
The easiest way to do this under linux is to add following line into your .bashrc

```bash
export DV_PATH="/path/to/your/datavault/"
```

It is also convenient to add simple analysis to your PYTHONPATH

```bash
export PYTHONPATH="/path-to-simpleanalysis/"
```

## Basic Usage:


```python

# Import the get_data module:

>>>> import simple_analysis.get_data as gd

# Analyze a day full of measurements

>>>> md = gd.MeasDay('2015Mar17')

# Let's see what data files we have:

>>>> print md.file_list

# Load a single file:

>>>> time_str = '2107_26'
>>>> md.read_file(time_str)

# Plot a dataset

>>>> md.plot_file(time_str)


# Have a look at the saved parameters

>>>> md.param_dict[time_str].keys()

```

## Data fitting

There are a couple of old fitting routines in the module fitdata. To get an overview type:

```python
import simple_analysis.fitdata as fd
help(fd)
```

The module fit_heating is based on fit_data and fits a heating rate to multiple 
spectra of the blue and red sideband.

A newer fitting class is called datafit. It is based on the lmfit class. 
See the file heatingraby.py for fitting heating rates to Rabi flops.
