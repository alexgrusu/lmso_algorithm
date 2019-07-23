# lmso_algorithm

The least-mean-square (LMS) and the normalized least-mean-square (NLMS) algorithms require a trade-off between fast convergence 
and low misadjustment, obtained by choosing the control parameters. In general, time variable parameters are proposed 
according to different rules. Many studies on the optimization of the NLMS algorithm imply time variable control parameters 
according some specific criteria.

The optimized LMS (LMSO) algorithm [1] for system identification is developed in the context of a state variable model, assuming 
that the unknown system acts as a time-varying system, following a first-order Markov model [2]. 

The proposed algorithm follows an optimization problem and introduces a variable step-size in order to minimize the system misalignment.

Installation
============

This poject is on pypi <https://pypi.org/project/lmso-algorithm/>

To install from pypi: 

	pip3 install lmso_algorithm
  
From this repository::

	git clone https://github.com/alexgrusu/lmso_algorithm.git
	python3 setup.py develop


Usage
=====

Supported features:

- The Optimized Least-Mean-Square (LMSO) Algorithm;
- WGN and AR(1) generators;
- NMSE computation;
- Echo path processing.

Example use:

    <example/test.py>.


[1] A. G. Rusu, S. Ciochină, and C. Paleologu, “On the step-size optimization of the LMS algorithm,” in Proc. IEEE TSP, 2019, 6 pages.

[2] G. Enzner, H. Buchner, A. Favrot, and F. Kuech, “Acoustic echo control,” in Academic Press Library in Signal Processing, 
vol. 4, ch. 30, pp. 807–877, Academic Press 2014.

