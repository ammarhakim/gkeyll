giADAS data
---------

ADAS data is not included in gkylzero and must be downloaded by user.

To download ADAS data and convert to format accessible by gkylzero,
run machines/mkdeps.<machine_name>.sh with flag --download-adas=yes

Radiation fits were originally are produced with [code maintained by J. Roeltgen](https://github.com/jRoeltgen/radiation_operator.git), requiring Matlab. Check this branch xxx for that version.

Now, we ported the matlab implemenation to python and the fits can be produced natively. More updates here soon, but check the main repository here: https://github.com/jaimix4/kin-rad-gkeyll.git
