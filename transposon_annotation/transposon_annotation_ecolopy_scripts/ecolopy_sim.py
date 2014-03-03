import matplotlib     # These two lines is because you're in a cluster
matplotlib.use('Agg') # and there may be no "X".
from ecolopy_dev import Community
from ecolopy_dev.utils import draw_shannon_distrib

com = Community('Fulcaldea_stuessyi_1kfrac.txt')
print com

com.fit_model('ewens')
com.set_current_model('ewens')
ewens_model = com.get_model('ewens')
print ewens_model

com.fit_model('lognormal')
com.set_current_model('lognormal')
lognormal_model = com.get_model('lognormal')
print lognormal_model

com.fit_model('etienne')
com.set_current_model('etienne')
etienne_model = com.get_model('etienne')
print etienne_model

tmp = {}
for met in ['fmin', 'slsqp', 'l_bfgs_b', 'tnc']:
    print 'Optimizing with %s...' % met
    try:
        com.fit_model(name='etienne', method=met, verbose=False)
        model = com.get_model('etienne')
        tmp[met] ={}
        tmp[met]['model'] = model
        tmp[met]['theta'] = model.theta
        tmp[met]['I']     = model.I
        tmp[met]['m']     = model.m
        tmp[met]['lnL']   = model.lnL
    except Exception as e:
        print '    optimization failed: ' + e.args[0]

# in case optimization by fmin failed to found correct values for theta and m:
if not (1 <= tmp['fmin']['theta'] < com.S and \
        1e-50 <= tmp['fmin']['m'] < 1-1e-50):
    del (tmp['fmin'])

# find the model with the higher likelihood:
met = min(tmp, key=lambda x: tmp[x]['lnL'])

# load it as 'etienne' model
com.set_model(tmp[met]['model'])

lrt = abd.lrt('ewens', 'etienne')
best = 'ewens' if lrt > 0.05 else 'etienne'
print 'Best model by LRT was: ' + best

com.generate_random_neutral_distribution(model=best)

pval, neut_h = com.test_neutrality (model=best, gens=10000, full=True)
#draw_shannon_distrib(neut_h, abd.shannon)
draw_shannon_distrib(neut_h, abd.shannon, outfile='Fulcaldea_stuessyi_1kfrac_shannon_dist.pdf', filetype='pdf')
print 'P-value for neutrality test was: ' + pval

out = open('Fulcaldea_stuessyi_1kfrac_shannon_neutral_data.tsv', 'w')
out.write('# shannon:' + str(shannon) + '\n')
out.write('\n'.join([str(s for s in neut_h)]) + '\n')
out.close()

com.dump_abundance('Fulcaldea_stuessyi_1kfrac_ecolopy.pik')
