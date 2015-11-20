from __future__ import print_function
from scipy import stats
import argparse
import os

def get_maxl_params(lines):
    maxL_paramline = lines.index('Done!') - 2
    params = lines[maxL_paramline].split()
    
    return params[2:]

parser = argparse.ArgumentParser(description='Gather parameters from dfe_ewk output and calculates the p-values')
parser.add_argument('-i', '--in_neutral', required=True, dest='infile_1', help = 'The input ile is the output from dfe_ewk program run with neutral')
parser.add_argument('-I', '--in_none', required=True, dest='infile_2', help = 'The input ile is the output from dfe_ewk program run with none')
parser.add_argument('-o', '--out', required=True, dest='outfile', help = 'Output fwith parameters and p-values for the models')
args = parser.parse_args()

rep = args.infile_1.split('.')[-3]

with open(args.infile_1) as infile:
    lines = [i.rstrip() for i in infile]
    neutral_params = get_maxl_params(lines)
    ln0 = float(neutral_params[-1])

with open(args.infile_2) as infile:
    lines = [i.rstrip() for i in infile]
    none_params = get_maxl_params(lines)
    ln1 = float(none_params[-1])

params = neutral_params + none_params
                                                                    
diff = 2*(ln1 - ln0)
df = 1
pval = str(stats.chi2.sf(diff, df))


params.append(str(diff))
params.append(pval)
params.append(rep)
out_string = '\t'.join(params)


with open(args.outfile, 'a') as outfile:
    if os.stat(args.outfile).st_size == 0: # print header if file is empty
        print('theta_neutral', 'lnL_neutral', 'theta_none', 'gamma_none', 'lnL_none', '2deltaL', 'p_value', 'replicate', file=outfile, sep='\t')
    print(out_string, file= outfile, sep='\t')
