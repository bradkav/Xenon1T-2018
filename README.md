# Xenon1T-2018

[![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/bradkav/Xenon1T-2018/master?filepath=Xenon1T-limit.ipynb) [![MIT Licence](https://badges.frapsoft.com/os/mit/mit.svg?v=103)](https://opensource.org/licenses/mit-license.php)

This code is a first attempt at reproducing (at least approximately) the Xenon1T-2018 upper limits on the DM-proton scattering cross section, announced by the Xenon1T collaboration on 28th May. The Xenon1T-2018 paper can be found here: [arXiv:1805.12562](https://arxiv.org/abs/1805.12562).

The python code is in the form of a jupyter notebook: [Xenon1T-limit.ipynb](Xenon1T-limit.ipynb). **You can run the notebook in your browser by clicking this badge: [![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/bradkav/Xenon1T-2018/master?filepath=Xenon1T-limit.ipynb)**. It's also available as an executable python script: [Xenon1T-limit.py](Xenon1T-limit.py).

Feel free to use this code however the hell you want. For questions, comments or bug reports, please contact Bradley J Kavanagh (bradkav@gmail.com).

### Notes & Caveats

* The code roughly reproduces the Xenon1T-2018 limit for standard spin-independent interactions (within a factor of about 2). However, it is straightforward to use it to recast results for other types of interactions (and I might add/calculate those limits myself at some point).
* At the moment, the code does a very rough Poisson upper limit, including no background uncertainties and using no spectral information from the observed events themselves. This should be improved as more information becomes available.
* With the release of the paper on arXiv, things have changed a little bit. Because we don't use event-by-event information, our best bet at getting a reasonable limit is to minimise the number of events (because not all the events look like WIMP recoils). So we use the 0.65t exposure in the signal reference region when we reproduce the limit.


### Requirements

* python3
* Standard python libraries like [`numpy`](http://www.numpy.org) and [`scipy`](https://www.scipy.org)
* [`WIMpy`](https://github.com/bradkav/WIMpy_NREFT), for calculating signal spectra (or your favourite DM signal calculator)


### Version history

**1st June 2018:** Updated to include the correct efficiency and updated numbers from the paper - [arXiv:1805.12562](https://arxiv.org/abs/1805.12562). Now uses the 0.65t information to set the limit.  
**28th May 2018 (2):** Updated to include Brazilian bands (updated folder structure a bit)  
**28th May 2018 (1):** Initial version, using only a simple Poisson upper limit, from crappy screenshots of the Xenon1T webcast.

### License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
