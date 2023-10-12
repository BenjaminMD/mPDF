class AttrDict:
    def __init__(self, d):
        for k, v in d.items():
            if isinstance(v, dict):
                setattr(self, k, AttrDict(v))
            else:
                setattr(self, k, v)

    def __repr__(self):
        return str(self.__dict__)
    

def r_gr_mgr(fit):
    prof = fit.recipe._contributions['PDF'].profile

    r = prof.x
    gobs = prof.y
    gcalc = fit.recipe._contributions['PDF'].evaluate()
    magcalc = mfit.mpdf(fit.recipe.ordscale.value,fit.recipe.parascale.value,fit.recipe.xi.value) # just mPDF
    baseline = 1.10 * gobs.min()
    gdiff = gobs - gcalc
    
    return r, gobs, gcalc, magcalc, baseline, gdiff