# Copula-Based Multi-User Performance Bounds

This repository contains supplementary material for the papers "Copula-Based
Multi-User Performance Bounds -- Part I: Theory" and "Part II: Applications"
(Karl-L. Besser, Eduard Jorswieck, 2020, [https://doi.org/XXX](doi:XXX) and
[https://doi.org/XXX](doi:XXX)).

The idea is to give an interactive version of the calculations and presented
concepts to the reader. One can also change different parameters and explore
different behaviors on their own.

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gl/)


## File List
The following files are provided in this repository:

### Part I: Theory
* [Basics - Copula.ipynb](https://mybinder.org/v2/gl/): Jupyter notebook that
  contains plots about basic copula theory.
* [Probability Bounds.ipynb](https://mybinder.org/): Jupyter notebook that
  contains plots and calculations about the bounds on the probability of sum
  and product of uniform distributions.
* `copulas.py`: Python module that provides classes for different copulas.

### Part II: Applications
* `Bounds-Expectation.nb`: Mathematica notebook that shows the bounds on the
  expected values of sum rate, MAC rate, and SINR.
* [Expected Value Bounds.ipynb](https://mybinder.org/): Jupyter notebook that
  shows the bounds on the expected values of sum rate, MAC rate, and SINR.
* [Probability Bounds - Rayleigh Fading.ipynb](https://mybinder.org/): Jupyter
  notebook that contains the probability bounds for the product of two Rayleigh
  fading channel gains, similar to the model for RIS channels.
* `expectation_bounds.py`: Calculations of the bounds for the expected values
  of the sum rate, MAC rate and SINR.


## Usage
### Running it online
The easiest way is to use services like [Binder](https://mybinder.org/) to run
the notebook online. Simply navigate to [https://mybinder.org/v2/gl/](TODO) to
run the notebook in your browser without setting everything up locally.

### Local Installation
If you want to run it locally on your machine, Python3 and Jupyter are needed.

Make sure you have [Python3](https://www.python.org/downloads/) installed on
your computer.
You can then install the requires packages (including Jupyter) by running
```
pip3 install -r requirements.txt
jupyter nbextension enable --py widgetsnbextension
```
This will install all the needed packages which are listed in the requirements 
file. The second line enables the interactive controls in the Jupyter
notebooks.

Finally, you can run the Jupyter notebooks with
```
jupyter notebook 'Basics - Copula.ipynb'
```

## More Applications (also presented in Part II)
[https://doi.org/XXX](Part II of the paper) shows some applications of the
theory presented in [https://doi.org/XXX](Part I).

`TODO:`

* Ergodic Secret-Key Capacity:
  [https://gitlab.com/klb2/bounds-secret-key-capacity](https://gitlab.com/klb2/bounds-secret-key-capacity)
  ([https://doi.org/XXX](doi:XXX))
* Secrecy Outage Probability:
  [https://gitlab.com/klb2/bounds-secrecy-outage](https://gitlab.com/klb2/bounds-secrecy-outage)
  ([https://doi.org/XXX](doi:XXX), [](arXiv:XXX))


## Acknowledgements
This research was supported in part by the Deutsche Forschungsgemeinschaft
(DFG) under grant JO 801/23-1.


## License and Referencing
This program is licensed under the GPLv3 license. If you in any way use this
code for research that results in publications, please cite our original
article listed above.
