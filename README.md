# TwoFieldApproximation
GitHub page with the purpose to present the code accompanying thesis. The code as presented here can not be considered a complete usable module but as part of the final thesis product.

## Notes

- The PyTransport directory within this repository is written by *David J. Mulryne and John W. Ronayne. PyTransport: A Python package for the calculation of inflationary correlation functions. arXiv e-prints, art. arXiv:1609.00381, September 2016*. See also: https://transportmethod.com/pytransport/. It is added to this repository for the sake of convenience
- The code is still heavily under development and not very user friendly.


The text below describes how to use '''trajectory_analyser.py'''. The main piece of code that calculates trajectory properties after an evolution by PyTransport.

```
usage: trajectory_analyser.py [-h] [-i I] [-f F] [-steps STEPS] [-twopt TWOPT]
                              [-bispec BISPEC] [-f1 F1] [-f2 F2] [-V V]
                              [-lab LAB] [-beta BETA] [-alp ALP] [-M M] [-a A] [-c C]
                              potential

positional arguments:
  potential       potential name (str)
  

optional arguments:
  -h, --help      show this help message and exit
  -i I            (float) initial e-folding. Default: 0 (float)
  -f F            (float) final e-folding. Default: 60 (float)
  -steps STEPS    (int) number of steps. Default: 8000 (int)
  -twopt TWOPT    (bool) whether to perform twopt function. Default: False (bool)
  -bispec BISPEC  (bool) whether to perform bispec function. Default: False (bool)
  -f1 F1          (float) initial value field 1. Default: 1 (float)
  -f2 F2          (float) initial value field 2. Default: 10 (float)
  -V V            (float) V parameter Orbital Inflation. Default: 1 (float)
  -lab LAB        (float) lambda parameter Orbital Inflation. Default: 0 (float)
  -beta BETA      (float) beta parameter Orbital Inflation. Default: 0 (float)
  -alp ALP        (float) alpha parameter supergravity Inflation. Default: 1 (float)
  -M M            (float) M parameter supergravity Inflation. Default: 10 (float)
  -a A            (float) a parameter supergravity Inflation.Default: 0 (float) 
  -c C            (float) c parameter supergravity Inflation. Default: 1 (float)
```

The preliminary script to patch together multiple evolutions of orbital inflation based on the EGNO model for a set of parameters can be found in '''patching.py'''.
