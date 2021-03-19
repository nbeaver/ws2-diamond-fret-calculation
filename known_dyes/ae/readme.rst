This says:

    OBS: Make sure your absorption spectra are in molar absorptivities (units 1/(M cm)).

    Overlap integral between 'ATTO_425_absorption' and 'ATTO_390_emission':
    J = 7.529e+09 nm^4/(M*cm)

which is the same as the FPbase value
once we multiply by the absorption coefficient 45,000/(M cm).

::

    You have: 7.529e+09 nm^4 45e3 / ((mol/L) cm)
    You want: 1e15 nm^4 / ((mol/L) cm)
    	* 0.338805
    	/ 2.9515503
