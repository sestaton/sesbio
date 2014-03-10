import matplotlib     
matplotlib.use('Agg') 
from ecolopy_dev import Community
from ecolopy_dev.utils import draw_shannon_distrib

##test_abund.txt would be the genome abundance / the diploid number of chromosomes
## j_tot is obtained by taking the total number of individuals
com = Community('test_log_abund.txt', j_tot=2150924886)
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
likelihoods = []
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
        # in case you reach two times the same likelyhood it may not be necessary
        # to go on with other optimization strategies... 
        # of course if time is not limiting it is not worth to check :)
        if round(model.lnL,1) in likelihoods:
            break
        likelihoods.append(round(model.lnL, 1))
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

lrt = com.lrt('ewens', 'etienne')
best = 'ewens' if lrt > 0.05 else 'etienne'
print 'Best model by LRT was: ' + best

com.generate_random_neutral_distribution(model=best)

## NB: There appears to be a bug with the method='loglike' function,
##     so it is only possible to compare Shannon's entropy
pval, neut_h = com.test_neutrality (model=best, gens=10000, full=True)
draw_shannon_distrib(neut_h, com.shannon, outfile='test_log_shannon_dist.pdf', filetype='pdf')
print 'P-value for neutrality test was: ', pval

out = open('test_log_shannon_neutral_data.tsv', 'w')
out.write('# shannon:' + str(com.shannon) + '\n')
out.write('\n'.join([str(s) for s in neut_h]) + '\n')
out.close()

com.dump_community('test_log_ecolopy.pik')
