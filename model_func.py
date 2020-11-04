#!/usr/bin/env python3
# PYTHON_PREAMBLE_START_STANDARD:{{{

# Christopher David Cotton (c)
# http://www.cdcotton.com

# modules needed for preamble
import importlib
import os
from pathlib import Path
import sys

# Get full real filename
__fullrealfile__ = os.path.abspath(__file__)

# Function to get git directory containing this file
def getprojectdir(filename):
    curlevel = filename
    while curlevel is not '/':
        curlevel = os.path.dirname(curlevel)
        if os.path.exists(curlevel + '/.git/'):
            return(curlevel + '/')
    return(None)

# Directory of project
__projectdir__ = Path(getprojectdir(__fullrealfile__))

# Function to call functions from files by their absolute path.
# Imports modules if they've not already been imported
# First argument is filename, second is function name, third is dictionary containing loaded modules.
modulesdict = {}
def importattr(modulefilename, func, modulesdict = modulesdict):
    # get modulefilename as string to prevent problems in <= python3.5 with pathlib -> os
    modulefilename = str(modulefilename)
    # if function in this file
    if modulefilename == __fullrealfile__:
        return(eval(func))
    else:
        # add file to moduledict if not there already
        if modulefilename not in modulesdict:
            # check filename exists
            if not os.path.isfile(modulefilename):
                raise Exception('Module not exists: ' + modulefilename + '. Function: ' + func + '. Filename called from: ' + __fullrealfile__ + '.')
            # add directory to path
            sys.path.append(os.path.dirname(modulefilename))
            # actually add module to moduledict
            modulesdict[modulefilename] = importlib.import_module(''.join(os.path.basename(modulefilename).split('.')[: -1]))

        # get the actual function from the file and return it
        return(getattr(modulesdict[modulefilename], func))

# PYTHON_PREAMBLE_END:}}}

def getinputdict(loglineareqs = True):
    inputdict = {}

    inputdict['paramssdict'] = {'GAMMA': 1, 'BETA': 0.96 ** (1/4), 'ETA': 2, 'ALPHA': 0.3, 'RHO_A': 0.95, 'Abar': 1, 'Pistar': 1.02 ** (1/4), 'PHIpi': 1.5, 'LAMBDA': 0.3, 'SIGMA': 8}
    inputdict['states'] = ['nu_tm1', 'A']
    inputdict['controls'] = ['C', 'Rp', 'W', 'L', 'Y', 'MC', 'Omega', 'I', 'Pi', 'U', 'V', 'PstaroverP']
    inputdict['shocks'] = ['epsilon_I']

    # equations:{{{
    inputdict['equations'] = []

    # household
    if loglineareqs is True:
        inputdict['equations'].append('-GAMMA * C = Rp - GAMMA * C_p')
    else:
        inputdict['equations'].append('C^(-GAMMA) = BETA*Rp*C_p^(-GAMMA)')
    if loglineareqs is True:
        inputdict['equations'].append('I = Rp + Pi_p')
    else:
        inputdict['equations'].append('I = Rp * Pi_p')
    if loglineareqs is True:
        inputdict['equations'].append('W - GAMMA * C = ETA * L')
    else:
        inputdict['equations'].append('W*C^(-GAMMA) = L^(ETA)')

    # firm production
    if loglineareqs is True:
        inputdict['equations'].append('MC = W - A')
    else:
        inputdict['equations'].append('MC = W / A')
    if loglineareqs is True:
        inputdict['equations'].append('nu_tm1_p + Y = A + L')
    else:
        inputdict['equations'].append('nu_tm1_p * Y = A * L')
    if loglineareqs is True:
        inputdict['equations'].append('Omega - Y = - (MC_ss * nu_tm1_ss) / (1 - MC_ss * nu_tm1_ss) * (MC + nu_tm1_p)')
    else:
        inputdict['equations'].append('Omega / Y = 1 - MC * nu_tm1_p')

    # firm pricing
    if loglineareqs is True:
        inputdict['equations'].append('(1 - LAMBDA) * Pi_ss ** (SIGMA - 1) * Pi = LAMBDA * PstaroverP_ss**(1-SIGMA) * PstaroverP')
    else:
        inputdict['equations'].append('1 = LAMBDA * PstaroverP ** (1 - SIGMA) + (1 - LAMBDA) * Pi ** (SIGMA - 1)')
    if loglineareqs is True:
        inputdict['equations'].append('U + PstaroverP = V')
    else:
        inputdict['equations'].append('U * PstaroverP = V')
    if loglineareqs is True:
        inputdict['equations'].append('U_ss * U = Y_ss * Y + Pi_ss**(SIGMA-1) * (1-LAMBDA) * 1 / Rp_ss * U_ss * ( (SIGMA-1) * Pi_p - Rp + U_p )')
    else:
        inputdict['equations'].append('U = Y + Pi_p**(SIGMA-1) * (1 - LAMBDA) * 1 / Rp * U_p')
    if loglineareqs is True:
        inputdict['equations'].append('V_ss * V = Y_ss * SIGMA / (SIGMA - 1) * MC_ss * (Y + MC) + Pi_ss**SIGMA * (1 - LAMBDA) * 1 / Rp_ss * V_ss * ( SIGMA * Pi_p - Rp + V_p )')
    else:
        inputdict['equations'].append('V = Y * SIGMA / (SIGMA - 1) * MC + Pi_p ** SIGMA * (1 - LAMBDA) * 1 / Rp * V_p')
    if loglineareqs is True:
        inputdict['equations'].append('nu_tm1_ss * nu_tm1_p = (1 - LAMBDA) * nu_tm1_ss * Pi_ss**SIGMA * (nu_tm1 + SIGMA * Pi) - SIGMA * LAMBDA * PstaroverP_ss**(-SIGMA) * PstaroverP')
    else:
        inputdict['equations'].append('nu_tm1_p = (1 - LAMBDA) * nu_tm1 * Pi ** SIGMA + LAMBDA * PstaroverP ** (-SIGMA)')
    
    # exogenous process
    if loglineareqs is True:
        inputdict['equations'].append('A_p = RHO_A * A')
    else:
        inputdict['equations'].append('log(A_p) = RHO_A*log(A) + (1 - RHO_A) * log(Abar)')

    # monetary policy
    if loglineareqs is True:
        inputdict['equations'].append('I = Rp + PHIpi * Pi + epsilon_I')
    else:
        inputdict['equations'].append('I = Rp * Pistar * (Pi / Pistar) ** PHIpi * exp(epsilon_I)')


    # resource
    if loglineareqs is True:
        inputdict['equations'].append('C = Y')
    else:
        inputdict['equations'].append('C = Y')
        
    # equations:}}}

    p = inputdict['paramssdict']
    p['Pi'] = p['Pistar']
    p['PstaroverP'] = ((1 - (1 - p['LAMBDA']) / p['Pi'] ** (1 - p['SIGMA'])) / p['LAMBDA']) ** (1 / (1 - p['SIGMA']))
    p['MC'] = (p['SIGMA'] - 1) / p['SIGMA'] * (1 - (1 - p['LAMBDA']) * p['BETA'] * p['Pi']**p['SIGMA']) / (1 - (1 - p['LAMBDA']) * p['BETA'] * p['Pi'] ** (p['SIGMA'] - 1)) * p['PstaroverP']
    p['nu_tm1'] = p['LAMBDA'] * p['PstaroverP'] ** (-p['SIGMA']) / (1 - (1 - p['LAMBDA']) * p['Pi'] ** p['SIGMA'])

    p['A'] = p['Abar']
    p['W'] = p['A'] * p['MC']
    p['Rp'] = 1/p['BETA']
    p['I'] = p['Rp'] * p['Pi']

    p['L'] = (p['W'] * (p['A'] / p['nu_tm1']) ** (-p['GAMMA'])) ** (1 / (p['ETA'] + p['GAMMA']))
    p['Y'] = p['A'] * p['L'] / p['nu_tm1']
    p['C'] = p['Y']

    p['Omega'] = p['Y'] * (1 - p['MC'] * p['nu_tm1'])

    p['U'] = p['Y'] / (1 - p['Pi'] ** (p['SIGMA'] - 1) * (1 - p['LAMBDA']) / p['Rp'])
    p['V'] = p['Y'] * p['SIGMA'] / (p['SIGMA'] - 1) * p['MC'] / (1 - p['Pi'] ** p['SIGMA'] * (1 - p['LAMBDA']) / p['Rp'])

    if loglineareqs is True:
        inputdict['loglineareqs'] = True
    else:
        inputdict['logvars'] = inputdict['states'] + inputdict['controls']
    inputdict['irfshocks'] = ['A', 'epsilon_I']

    # save stuff
    inputdict['savefolder'] = __projectdir__ / Path('temp/')

    # main vars
    inputdict['mainvars'] = ['C', 'Rp', 'Pi', 'I', 'L', 'Omega']

    return(inputdict)


def check():
    inputdict_loglin = getinputdict(loglineareqs = True)
    inputdict_log = getinputdict(loglineareqs = False)
    importattr(__projectdir__ / Path('submodules/dsge-perturbation/dsgediff_func.py'), 'checksame_inputdict')(inputdict_loglin, inputdict_log)
    

def dsgefull():
    inputdict = getinputdict()
    # importattr(__projectdir__ / Path('submodules/dsge-perturbation/dsge_bkdiscrete_func.py'), 'discretelineardsgefull')(inputdict)
    importattr(__projectdir__ / Path('submodules/dsge-perturbation/dsge_bkdiscrete_func.py'), 'discretelineardsgefull')(inputdict)


# Run:{{{1
check()
dsgefull()
