

from ecolopy_dev import Community
from ecolopy_dev.utils import draw_shannon_distrib

abd = Community('Carthamus_tinctorius_1kfrac.txt')
print abd

abd.fit_model('ewens')
abd.set_current_model('ewens')
ewens_model = abd.get_model('ewens')
print ewens_model

abd.fit_model('lognormal')
abd.set_current_model('lognormal')
lognormal_model = abd.get_model('lognormal')
print lognormal_model

abd.fit_model('etienne')
abd.set_current_model('etienne')
etienne_model = abd.get_model('etienne')
print etienne_model

lrt = abd.lrt('ewens', 'etienne')
best = 'ewens' if lrt > 0.05 else 'etienne'
print best

abd.generate_random_neutral_distribution(model=best)

pval, neut_h = abd.test_neutrality (model=best,
                                    gens=10000, full=True)
draw_shannon_distrib(neut_h, abd.shannon)
print pval
