# Copula-Based Multi-User Performance Bounds

This repository contains supplementary material for the papers "Copula-Based
Bounds for Multi-User Communications - Part I: Average Performance" and "Part
II: Outage Performance" (Karl-L. Besser, Eduard Jorswieck, 2020, [doi:XXX]()
and [doi:XXX]()).

The idea is to give an interactive version of the calculations and presented
concepts to the reader. One can also change different parameters and explore
different behaviors on their own.

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gl/klb2%2Fcopula-mu-performance-bounds/master)


## File List
The following files are provided in this repository:

### Part I: Average Performance
* [Motivation.ipynb](https://mybinder.org/v2/gl/klb2%2Fcopula-mu-performance-bounds/master?filepath=Motivation.ipynb):
  Jupyter notebook that contains a simple motivating example about channels
  with negative correlation.
* `Bounds-Expectation.nb`: Mathematica notebook that shows the bounds on the
  expected values of sum rate, MAC rate, and SINR.
* [Expected Value Bounds.ipynb](https://mybinder.org/v2/gl/klb2%2Fcopula-mu-performance-bounds/master?filepath=Expected%20Value%20Bounds.ipynb): Jupyter notebook that
  shows the bounds on the expected values of sum rate, MAC rate, and SINR.
* `expectation_bounds.py`: Calculations of the bounds for the expected values
  of the sum rate, MAC rate and SINR.

### Part II: Outage Performance
* [Basics - Copula.ipynb](https://mybinder.org/v2/gl/klb2%2Fcopula-mu-performance-bounds/master?filepath=Basics%20-%20Copula.ipynb): Jupyter notebook that contains plots about basic copula theory.
* [Probability Bounds.ipynb](https://mybinder.org/v2/gl/klb2%2Fcopula-mu-performance-bounds/master?filepath=Probability%20Bounds.ipynb): Jupyter notebook that
  contains plots and calculations about the bounds on the probability of sum
  and product of uniform distributions.
* [Probability Bounds - Rayleigh Fading.ipynb](https://mybinder.org/v2/gl/klb2%2Fcopula-mu-performance-bounds/master?filepath=Probability%20Bounds%20-%20Rayleigh%20Fading.ipynb): Jupyter
  notebook that contains the probability bounds for the product of two Rayleigh
  fading channel gains, similar to the model for RIS channels.
* `copulas.py`: Python module that provides classes for different copulas.
* `probability_bounds.py`: Calculations of the bounds for the probability of
  functions of random variables.
* `comparison_correlation_model.py`: Python module that provides a comparison
  of the copula-based bounds with a (linear) correlation model from literature.


## Usage
### Running it online
The easiest way is to use services like [Binder](https://mybinder.org/) to run
the notebook online. Simply navigate to
[https://mybinder.org/v2/gl/klb2%2Fcopula-mu-performance-bounds/master](https://mybinder.org/v2/gl/klb2%2Fcopula-mu-performance-bounds/master)
to run the notebooks in your browser without setting everything up locally.

### Local Installation
If you want to run it locally on your machine, Python3 and Jupyter are needed.
The present code was developed and tested with the following versions:
- Python 3.8
- Jupyter 1.0
- numpy 1.18
- scipy 1.4

Make sure you have [Python3](https://www.python.org/downloads/) installed on
your computer.
You can then install the requires packages (including Jupyter) by running
```bash
pip3 install -r requirements.txt
jupyter nbextension enable --py widgetsnbextension
```
This will install all the needed packages which are listed in the requirements 
file. The second line enables the interactive controls in the Jupyter
notebooks.

Finally, you can run the Jupyter notebooks with
```bash
jupyter notebook 'Basics - Copula.ipynb'
```

## More Applications
Some applications, that are mentioned in the papers, are discussed in more
detail in the following publications:

* Ergodic Secret-Key Capacity:
  [https://gitlab.com/klb2/bounds-secret-key-capacity](https://gitlab.com/klb2/bounds-secret-key-capacity)
* Secrecy Outage Probability:
  [https://gitlab.com/klb2/bounds-secrecy-outage](https://gitlab.com/klb2/bounds-secrecy-outage)
  ([arXiv:2004.06644](https://arxiv.org/abs/2004.06644))


## Acknowledgements
This research was supported in part by the Deutsche Forschungsgemeinschaft
(DFG) under grant JO 801/23-1.


## License and Referencing
This program is licensed under the GPLv3 license. If you in any way use this
code for research that results in publications, please cite our original
article listed above.
