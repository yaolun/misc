def radial_chisq(modelrange, filename_ext, indir, plotpath):
    import numpy as np
    import matplotlib.pyplot as plt
    from astropy.io import ascii

    model_start, model_end = modelrange
    
