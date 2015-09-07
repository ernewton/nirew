FUNCTION waterindex, lambda, spectra
    sp = [[2.070,2.090],[2.235,2.255],[2.360,2.380]]

    FOR i=0,2 DO BEGIN
        sp[*,i] = MEDIAN(spectra[where(lambda ge sp[0,i] and lambda le sp[1,i])])
    ENDFOR
    sp = sp[0,*]
    index = sp[0] * sp[2] / sp[1]^2.

    return, index
END
