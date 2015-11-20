from __future__ import print_function
import sys
import numpy as np
#from scipy import integrate
from math import exp, log
import argparse
sys.path.append('/home/bo1pgc/local/python')
from pygsl import integrate

def log_factorial(n):
    log_sum = 0
    for i in xrange(1, n+1):
        log_sum += log(i)
    return log_sum


def n_choose_i(n, i):

    """n choose i for use in  Equation 2 Bustamante et al. 2001 Genetics 159:1779-1788"""

    if i == 0 or i == n:
        nci = 1
    else:
        nci = exp(log_factorial(n) - log_factorial(i) - log_factorial(n - i))

    return nci


def integrand_eq2(q, args):
    """Equation 2 Bustamante et al. 2001 Genetics 159:1779-1788
    The q(1-q) in the denominator term has been canceled out using the binomial probability
    density
    
    Kai's program uses is twice the gamma in the Bustamante paper, hence the 2*gamma is replaced
    by gamma
    """

    gamma = args[0]
    nci = args[1]
    n = args[2]
    i = args[3]

    return ((1 - exp(-1 * gamma * (1 - q))) / (1 - exp(-1 * gamma))) * (nci * (q ** (i - 1)) * ((1 - q) ** (n - i - 1)))

def integrand_eq2_neg(q, args):
    """Equation 2 Bustamante et al. 2001 Genetics 159:1779-1788                                                                                                                  
    The right hand side of Eqaution 1 is modified in this function to 
    
    (exp(gamma) - exp(gamma*q))/(exp(gamma) - 1)
    
    This is used for cased when gamma is large and negative and may 
    cause overflow when in the exponential term
    """

    gamma = args[0]
    nci = args[1]
    n = args[2]
    i = args[3]

    #print(gamma, nci, n, i)
    
    
    return ((exp(gamma) - exp(gamma * q)) / (exp(gamma) - 1))  * (nci * (q ** (i - 1)) * ((1 - q) ** (n - i - 1)))

    #return ((exp(gamma) - exp(gamma * q)) / (exp(gamma) - exp(-gamma))) * (nci * (q ** (i - 1)) * ((1 - q) ** (n - i - 1)))


def neutral_integrand_eq2(q, args):

    """Equation 2 Bustamante et al. 2001 Genetics 159:1779-1788 for the neutral case where gamma is zero
    Equation 1 reduced to 2/ q when gamma is zero, the q term in the denominator for the
    neutral version of Equation 2 has been canceled out using the binomial probability
    density"""

    nci = args[0] 
    n = args[1]
    i = args[2]

    return nci * (q ** (i - 1)) * ((1 - q) ** (n - i))


def get_sfs(integral, theta_ls):

    """
    Generates the sfs using Equation 3 from Bustamanted et al. 2001
    by randomly sampling from poisson distribution where lambda = theta_ls*integral
    """

    x = np.random.poisson(theta_ls*integral)

    return x


def main():

    parser = argparse.ArgumentParser(description='Simulation of selected and neutral site frequency spectra using the Poisson Random Field Model (PRF)')
    parser.add_argument('-n', '--nsam', type=int, required=True, dest='sample', help='The sample size')
    parser.add_argument('-t', '--theta', type=float, required=True, dest='theta', help='The scaled mutation parameter per site')
    parser.add_argument('-g', '--gamma', type=float, required=True, dest='gamma', help='The scaled selection coefficient (Ns). Can be positive or negative')
    parser.add_argument('-N', '--Nsites', required=True, type=int, dest='length_S', help = 'Size of region for neutral sites in base pairs to simulate using the PRF [default: 500] ')
    parser.add_argument('-S', '--Ssites', required=True, type=int, dest='length_N', help = 'Size of region for selected sites in base pairs to simulate using the PRF [default: 500] ')
    parser.add_argument('-o', '--out', required=True, dest='outfile', help = 'Output file for simulated Neutral and Selected sites site frequency spectra')
    parser.add_argument('-s', '--sfs', default='folded', required=False, dest='fold', help = 'Specify folded or unfolded [ default: folded]')
    args = parser.parse_args()

    n = args.sample
    
    Ls_N = args.length_N
    theta_ls_N = args.theta * Ls_N
    
    Ls_S = args.length_S
    theta_ls_S = args.theta * Ls_S # specified theta per site at input, so multiply by length

    gamma = args.gamma # gamma in Kai's program is twice that in the Bustamante paper so multiply by a half


    if args.fold == 'folded':
        fold = True
    elif args.fold == 'unfolded':
        fold = False
    else:
        print('Invalid option for -s, it should be folded or unfolded')
        sys.exit(1)
    
    with open(args.outfile, 'w') as outfile:
    
        nl = 1 # consider adding option to have multiple replicate of varying size output
        sfs = []
        workspace = integrate.workspace(1000000)
        if fold:
            print('S', file=outfile)
            print(nl, file=outfile)
            for i in xrange(1, (n + 2) / 2):
                if gamma == 0: # still print "Selected" sites in output when gamma is zero
                    if i < (n - i):
                        nci = n_choose_i(n, i)
                        neutral_integrand_eq2_func = integrate.gsl_function(neutral_integrand_eq2, [nci, n, i])
                        F_i_gamma = integrate.qag(neutral_integrand_eq2_func, 0, 1, 1e-8, 1e-8, 1000, integrate.GAUSS61, workspace)[1]
                        
                        nci = n_choose_i(n, n - i)
                        neutral_integrand_eq2_func = integrate.gsl_function(neutral_integrand_eq2, [nci, n, n - i])
                        F_n_minus_i_gamma = integrate.qag(neutral_integrand_eq2_func, 0, 1, 1e-8, 1e-8, 1000, integrate.GAUSS61, workspace)[1]
    
                        G_i_gamma = F_i_gamma + F_n_minus_i_gamma # equation 6 condition i < n - i
                        
                        xN = get_sfs(G_i_gamma, theta_ls_N)
                        print(i, G_i_gamma*theta_ls_N)
                        sfs.append(str(xN))
    
                    elif i == n / 2:
                        nci = n_choose_i(n, n / 2)
                        neutral_integrand_eq2_func = integrate.gsl_function(neutral_integrand_eq2, [nci, n, n / 2])
                        F_n_over_two_gamma = integrate.qag(neutral_integrand_eq2_func, 0, 1, 1e-8, 1e-8, 1000, integrate.GAUSS61, workspace)[1]
                        xN = get_sfs(F_n_over_two_gamma, theta_ls_N)
                        print(i, F_n_over_two_gamma * theta_ls_N)
                        sfs.append(str(xN))
                        
                elif gamma < 0:
                    if i < (n - i):
                        nci = n_choose_i(n, i)
                        integrand_eq2_func = integrate.gsl_function(integrand_eq2_neg, [gamma, nci, n, i])
                        F_i_gamma = integrate.qag(integrand_eq2_func, 0, 1, 1e-8, 1e-8, 1000, integrate.GAUSS61, workspace)[1]
                        
                        
                        nci = n_choose_i(n, n - i)
                        integrand_eq2_func = integrate.gsl_function(integrand_eq2_neg, [gamma, nci, n, n - i])
                        F_n_minus_i_gamma = integrate.qag(integrand_eq2_func, 0, 1, 1e-8, 1e-8, 1000, integrate.GAUSS61, workspace)[1]

                        G_i_gamma = F_i_gamma + F_n_minus_i_gamma # equation 6 condition i < n - i                                                                               
                        
                        #print(G_i_gamma, F_i_gamma, F_n_minus_i_gamma)
                        xS = get_sfs(G_i_gamma, theta_ls_S)
                        print(i, G_i_gamma * theta_ls_S)
                        sfs.append(str(xS))
                    
                    elif i == n / 2: # floor division in python 2
                        nci = n_choose_i(n, n / 2)
                        integrand_eq2_func = integrate.gsl_function(integrand_eq2_neg, [gamma, nci, n, n / 2])
                        F_n_over_two_gamma = integrate.qag(integrand_eq2_func, 0, 1, 1e-10, 1e-10, 1000, integrate.GAUSS61, workspace)[1]
                        xS = get_sfs(F_n_over_two_gamma , theta_ls_S)
                        #print(F_n_over_two_gamma)
                        print(i, F_n_over_two_gamma * theta_ls_S)
                        sfs.append(str(xS))
                    
    
                else:
                    if i < (n - i):
                        nci = n_choose_i(n, i)
                        integrand_eq2_func = integrate.gsl_function(integrand_eq2, [gamma, nci, n, i])
                        F_i_gamma = integrate.qag(integrand_eq2_func, 0, 1, 1e-8, 1e-8, 1000, integrate.GAUSS61, workspace)[1]
    
                        nci = n_choose_i(n, n - i)
                        integrand_eq2_func = integrate.gsl_function(integrand_eq2, [gamma, nci, n, n - i])
                        F_n_minus_i_gamma = integrate.qag(integrand_eq2_func, 0, 1, 1e-8, 1e-8, 1000, integrate.GAUSS61, workspace)[1]
    
                        G_i_gamma = F_i_gamma + F_n_minus_i_gamma # equation 6 condition i < n - i
                        xS = get_sfs(G_i_gamma, theta_ls_S)
                        print(i, G_i_gamma * theta_ls_S)
                        sfs.append(str(xS))
    
                    elif i == n / 2: # floor division in python 2
                        nci = n_choose_i(n, n / 2)
                        integrand_eq2_func = integrate.gsl_function(integrand_eq2, [gamma, nci, n, n / 2])
                        F_n_over_two_gamma = integrate.qag(integrand_eq2_func, 0, 1, 1e-8, 1e-8, 1000, integrate.GAUSS61, workspace)[1]
                        xS = get_sfs(F_n_over_two_gamma , theta_ls_S)
                        print(i, F_n_over_two_gamma * theta_ls_S)
                        sfs.append(str(xS))
    
            sfs_string = ', '.join(sfs)
            print(n, Ls_S, sfs_string, sep=', ', file=outfile)
            del sfs[:]
    
    
            # Neutral sites
            print('N', file=outfile)
            print(nl, file=outfile)
            for i in xrange(1, (n + 2) / 2):
                if i < (n - i):
                    nci = n_choose_i(n, i)
                    neutral_integrand_eq2_func = integrate.gsl_function(neutral_integrand_eq2, [nci, n, i])
                    F_i_gamma = integrate.qag(neutral_integrand_eq2_func, 0, 1, 1e-8, 1e-8, 1000, integrate.GAUSS61, workspace)[1]
    
                    nci = n_choose_i(n, n - i)
                    neutral_integrand_eq2_func = integrate.gsl_function(neutral_integrand_eq2, [nci, n, n - i])
                    F_n_minus_i_gamma = integrate.qag(neutral_integrand_eq2_func, 0, 1, 1e-8, 1e-8, 1000, integrate.GAUSS61, workspace)[1]
    
                    G_i_gamma = F_i_gamma + F_n_minus_i_gamma # equation 6 condition i < n - i
                    xN = get_sfs(G_i_gamma, theta_ls_N)
                    #print(i, G_i_gamma * theta_ls_N)
                    sfs.append(str(xN))
    
                elif i == n / 2:
                    nci = n_choose_i(n, n / 2)
                    neutral_integrand_eq2_func = integrate.gsl_function(neutral_integrand_eq2, [nci, n, n / 2])
                    F_n_over_two_gamma = integrate.qag(neutral_integrand_eq2_func, 0, 1, 1e-8, 1e-8, 1000, integrate.GAUSS61, workspace)[1]
                    xN = get_sfs(F_n_over_two_gamma, theta_ls_N)
                    #print(i, F_n_over_two_gamma * theta_ls_N)
                    sfs.append(str(xN))
    
            sfs_string = ', '.join(sfs)
            print(n, Ls_N, sfs_string, sep=', ', file=outfile)
            del sfs[:]
    
        else: # produce unfolded sfs
            print('S', file=outfile)
            print(nl, file=outfile)
            for i in xrange(1, n):
                if gamma == 0:
                    nci = n_choose_i(n , i)
                    neutral_integrand_eq2_func = integrate.gsl_function(neutral_integrand_eq2, [nci, n, i])
                    eq2_eval_neut = integrate.qag(neutral_integrand_eq2_func, 0, 1,  1e-8, 1e-8, 1000, integrate.GAUSS61, workspace)[1]
                    xN = get_sfs(eq2_eval_neut, theta_ls_N)
                    # print(i, xN)
                    sfs.append(str(xN))
                    
                elif gamma < 0:
                    nci = n_choose_i(n, i)
                    integrand_eq2_func = integrate.gsl_function(integrand_eq2_neg, [gamma, nci, n, i])
                    eq2_eval_sel = integrate.qag(integrand_eq2_func, 0, 1,  1e-8, 1e-8, 1000, integrate.GAUSS61, workspace)[1]
                    xS = get_sfs(eq2_eval_sel, theta_ls_S)
                    # print(i, eq2_eval_sel)                                                                                                                                     
                    sfs.append(str(xS))
                    
                else:
                    nci = n_choose_i(n, i)
                    integrand_eq2_func = integrate.gsl_function(integrand_eq2, [gamma, nci, n, i])
                    eq2_eval_sel = integrate.qag(integrand_eq2_func, 0, 1,  1e-8, 1e-8, 1000, integrate.GAUSS61, workspace)[1]
                    xS = get_sfs(eq2_eval_sel, theta_ls_S)
                    # print(i, eq2_eval_sel)
                    sfs.append(str(xS))
    
    
            sfs_string = ', '.join(sfs)
            print(n, Ls_S, sfs_string, sep=', ', file=outfile)
            del sfs[:]
    
            print('N', file=outfile)
            print(nl, file=outfile)
            for i in xrange(1, n):
                nci = n_choose_i(n , i)
                neutral_integrand_eq2_func = integrate.gsl_function(neutral_integrand_eq2, [nci, n, i])
                eq2_eval_neut = integrate.qag(neutral_integrand_eq2_func, 0, 1, 1e-8, 1e-8, 1000, integrate.GAUSS61, workspace)[1]
                xN = get_sfs(eq2_eval_neut, theta_ls_N)
                # print(i, xN)
                sfs.append(str(xN))
    
            sfs_string = ', '.join(sfs)
            print(n, Ls_N, sfs_string, sep=', ', file=outfile)
            del sfs[:]
    
if __name__ == '__main__':
    main()
