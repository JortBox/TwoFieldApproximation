# TwoFieldApproximation
GitHub page with the purpose to present the code accompanying thesis. The code as presented here can not be considered a complete usable module but as part of the final thesis product.

## Notes

- The PyTransport directory within this repository is written by *David J. Mulryne and John W. Ronayne. PyTransport: A Python package for the calculation of inflationary correlation functions. arXiv e-prints, art. arXiv:1609.00381, September 2016*. See also: https://transportmethod.com/pytransport/. It is added to this repository for the sake of convenience
- The code is still heavily under development and not very user friendly.



```
usage: trajectory_analyser.py [-h] [-i I] [-f F] [-steps STEPS] [-twopt TWOPT]
                              [-bispec BISPEC] [-f1 F1] [-f2 F2] [-V V]
                              [-lab LAB] [-beta BETA] [-alp ALP] [-M M] [-a A]
                              potential c

positional arguments:
  potential       potential name
  c               (float) c parameter supergravity Inflation. Default: 1

optional arguments:
  -h, --help      show this help message and exit
  -i I            (float) initial e-folding. Default: 0
  -f F            (float) final e-folding. Default: 60
  -steps STEPS    (int) number of steps. Default: 8000
  -twopt TWOPT    (bool) whether to perform twopt function. Default: False
  -bispec BISPEC  (bool) whether to perform bispec function. Default: False
  -f1 F1          (float) initial value field 1. Default: 1
  -f2 F2          (float) initial value field 2. Default: 10
  -V V            (float) V parameter Orbital Inflation. Default: 1
  -lab LAB        (float) lambda parameter Orbital Inflation. Default: 0
  -beta BETA      (float) beta parameter Orbital Inflation. Default: 0
  -alp ALP        (float) alpha parameter supergravity Inflation. Default: 1
  -M M            (float) M parameter supergravity Inflation. Default: 10
  -a A            (float) a parameter supergravity Inflation.Default: 0
```
