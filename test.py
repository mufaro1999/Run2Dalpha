from time import time
from TwoDAlphabet import plot
from TwoDAlphabet.twoDalphabet import MakeCard, TwoDAlphabet
from TwoDAlphabet.alphawrap import BinnedDistribution, ParametricFunction
from TwoDAlphabet.helpers import make_env_tarball, cd, execute_cmd
from TwoDAlphabet.ftest import FstatCalc
import os
import numpy as np

# Helper function to get region names
def _get_other_region_names(Pass_reg_name):
    '''
    If Passing e.g. "Fail", will return ("Fail", "Pass")
    '''
    return Pass_reg_name, Pass_reg_name.replace('Fail','Loose'),Pass_reg_name.replace('Fail','Pass')
    #return Pass_reg_name, Pass_reg_name.replace('Fail','Pass')

# Helper function to generate constraints for parametric Transfer Functions
# Change values as you see fit
def _generate_constraints(nparams):
    out = {}
    for i in range(nparams):
        if i == 0:
            out[i] = {"MIN":-50,"MAX":50}
        else:
            out[i] = {"MIN":-50,"MAX":50}
    return out

# Dict to store transfer function forms and constraints
_rpf_options = {
    '0x0': {
        'form': '(@0)',
        'constraints': _generate_constraints(1)
    },
    '1x0': {
        'form': '(@0+@1*x)',
        'constraints': _generate_constraints(2)
    },
    '0x1': {
        'form': '(@0+@1*y)',
        'constraints': _generate_constraints(2)
    },
    '1x1': {
        'form': '(@0+@1*x)*(1+@2*y)',
        'constraints': _generate_constraints(3)
    },
    '2x0': {
        'form': '(@0+@1*x+@2*x**2)*(@3)',
        'constraints': _generate_constraints(4)
    },
    '2x1': {
        'form': '(@0+@1*x+@2*x**2)*(1+@3*y)',
        'constraints': _generate_constraints(4)
    },
    '2x2': {
        'form': '(@0+@1*x+@2*x**2)*(1+@3*y+@4*y**2)',
        'constraints': _generate_constraints(4)
    }
}

# Helper function for selecting the signal from the ledger
def _select_signal(row, args):
    signame = args[0]
    poly_order = args[1]
    if row.process_type == 'SIGNAL':
        if signame in row.process:
            return True
        else:
            return False
    elif 'Background' in row.process:
        if row.process == 'Background_'+poly_order:
            return True
        elif row.process == 'Background':
            return True
        else:
            return False
    else:
        return True


# Make the workspace
def make_workspace():
    # Create the workspace directory, using info from the specified JSON file 'analysis_note_CR.json'
    twoD = TwoDAlphabet('XHY_fits', 'analysis_note_CR.json', loadPrevious=False)

    # 2DAlphabet wasn't intended for an analysis like this, so the default function
    # for Looping over all regions and for a given region's data histogram, subtracting
    # the list of background histograms, and returning a data-bkgList is called initQCDHists.
    # This is b/c QCD multijet is the main background we usually estimate via 2DAlphabet
    bkg_hists = twoD.InitQCDHists()
    print('bkg_hists = {}'.format(bkg_hists))

    # Now, we loop over "Pass" and "Fail" regions and get the binnings
    for f,l,p in [_get_other_region_names(r) for r in twoD.ledger.GetRegions() if 'Fail' in r]:
        binning_f, _ = twoD.GetBinningFor(f)
    Fail_name = 'Background_'+f
    bkg_f = BinnedDistribution(
		Fail_name, bkg_hists[f],
		binning_f, constant=False
	)
    twoD.AddAlphaObj('Background',f,bkg_f)

    # 1x0 TF to go between F->L
    qcd_rpfL = ParametricFunction(
        Fail_name.replace('Fail','rpfL'),
        binning_f, _rpf_options['2x1']['form'],
        constraints=_rpf_options['2x1']['constraints']
    )
    # 2x1 TF to go between L->T
    qcd_rpfT = ParametricFunction(
        Fail_name.replace('Fail','rpfT'),
        binning_f, _rpf_options['1x0']['form'],
        constraints=_rpf_options['1x0']['constraints']
    )

    # qcd_l = qcd_f*rpfL
    bkg_l = bkg_f.Multiply(Fail_name.replace('Fail','Loose'),qcd_rpfL)
    # qcd_p = qcd_l*rpfT
    bkg_p = bkg_l.Multiply(Fail_name.replace('Fail','Pass'),qcd_rpfT)

    twoD.AddAlphaObj('Background',l,bkg_l,title='Background')
    twoD.AddAlphaObj('Background',p,bkg_p,title='Background')
    '''for f, p in [_get_other_region_names(r) for r in twoD.ledger.GetRegions() if 'Fail' in r]:
        print(f, p)
        # get the binning for the Fail region
        binning_f, _ = twoD.GetBinningFor(f)
        # you can change the name as you see fit
        Fail_name = 'Background_'+f
        # this is the actual binned distribution of the Fail
        bkg_f = BinnedDistribution(Fail_name, bkg_hists[f], binning_f, constant=False)
        # now we add it to the 2DAlphabet ledger
        twoD.AddAlphaObj('Background',f, bkg_f)

    # Now let's make the parametric transfer function to go from Fail -> Pass
    bkg_rpf = ParametricFunction(
        Fail_name.replace('Fail','rpf'),        # this is our Pass/Fail ratio
        binning_f,                              # we use the binning from Fail
        _rpf_options['0x0']['form'],          # let's make it constant in Ias (and ProbQ obviously)
        _rpf_options['0x0']['constraints']    # use the default constraints [0,5]
    )

    # now define the bkg in Pass as the bkg in Fail multiplied by the transfer function (bkg_rpf)
    bkg_p = bkg_f.Multiply(Fail_name.replace('Fail','Pass'), bkg_rpf)

    # then add this to the 2DAlphabet ledger
    twoD.AddAlphaObj('Background',p,bkg_p,title='Background')

    # and save it out

    '''
    twoD.Save()


# function for perfomring the fit
def perform_fit(signal='MX3000_MY300', extra='--robustFit 1'):
    '''
        extra (str) = any extra flags to Pass to Combine when running the ML fit
    '''
    # this is the name of the directory created in the workspace function
    working_area = 'XHY_fits'
    # we reuse the workspace from the last step.
    # The runConfig.json is copied from the origin JSON config file,
    # and we must specify that we want to load the previous workspace
    twoD = TwoDAlphabet(working_area, '{}/runConfig.json'.format(working_area), loadPrevious=True)

    # Now we create a ledger and make a new area for it with a Combine card
    # this select() method uses lambda functions. Will explain later
    subset = twoD.ledger.select(_select_signal, "{}".format("MX3000_MY300"), '')
    twoD.MakeCard(subset, '{}_area'.format(signal))

    # perform fit
    twoD.MLfit('{}_area'.format(signal), rMin=-5, rMax=20, verbosity=1, extra=extra)

def plot_fit(signal='MX3000_MY300'):
    working_area = 'XHY_fits'
    twoD = TwoDAlphabet(working_area, '{}/runConfig.json'.format(working_area), loadPrevious=True)
    subset = twoD.ledger.select(_select_signal, '{}'.format("MX3000_MY300"), '')
    twoD.StdPlots('{}_area'.format(signal), subset)

def GOF(signal='MX3000_MY300',condor=True, extra=''):
    working_area = 'XHY_fits'
    twoD = TwoDAlphabet(working_area, '{}/runConfig.json'.format(working_area), loadPrevious=True)
    if not os.path.exists(twoD.tag+'/{}_area/card.txt'.format(signal)):
        print(twoD.tag+'/'+signal+'_area/card.txt does not exist, making Combine card')
        subset = twoD.ledger.select(_select_signal, signal, '')
        twoD.MakeCard(subset, '{}_area'.format(signal))
    if condor:
        twoD.GoodnessOfFit(
            '{}_area'.format(signal),      # tag
            ntoys=500,          # number of toys to generate
            freezeSignal=0,     # whether to freeze signal
            extra=extra,        # any extra commands to Pass
            condor=True,        # ship GOF toys to condor
            njobs=10            # N jobs will run the ntoys
        )
    else:
        twoD.GoodnessOfFit(
            '{}_area'.format(signal),
            ntoys=100,          # this will take ages to run without condor
            freezeSignal=0,
            extra=extra,
            condor=False
        )

def plot_GOF(signal='MX3000_MY300',condor=True):
    working_area = 'XHY_fits'
    plot.plot_gof('{}'.format(working_area), '{}_area'.format(signal), condor=condor)

def load_RPF(twoD):
    '''
        loads the rpf parameter values for use in toy generation
    '''
    params_to_set = twoD.GetParamsOnMatch('rpf.*', 'Signal', 'b')
    return {k:v['val'] for k,v in params_to_set.items()}

def SignalInjection(r, signal='MX3000_MY300',condor=True):
    working_area = 'XHY_fits'
    twoD = TwoDAlphabet(working_area, '{}/runConfig.json'.format(working_area), loadPrevious=True)
    #params = load_RPF(twoD)
    twoD.SignalInjection(
        '{}_area'.format(signal),
        injectAmount = r,       # injected signal xsec (r=0 : bias test)
        ntoys=500,              # will take forever if not on condor
        blindData = True,       # make sure you're blinding if working with data
        #setParams = params,     # give the toys the same RPF params
        verbosity = 0,          # you can change this if you need
        njobs=10,
        condor = condor
    )

def plot_SignalInjection(r, signal='MX3000_MY300', condor=False):
    working_area = 'XHY_fits'
    plot.plot_signalInjection(working_area, '{}_area'.format(signal), injectedAmount=r, condor=condor)

def Impacts(signal='MX3000_MY300'):
    working_area = 'XHY_fits'
    twoD = TwoDAlphabet(working_area, '{}/runConfig.json'.format(working_area), loadPrevious=True)
    twoD.Impacts('{}_area'.format(signal), cardOrW='card.txt', extra='-t 1')

def run_limits(signal='MX3000_MY300'):
    working_area = 'XHY_fits'
    twoD = TwoDAlphabet(working_area, '{}/runConfig.json'.format(working_area), loadPrevious=True)
    twoD.Limit(
        subtag='{}_area'.format(signal),
        blindData=False,        # BE SURE TO CHANGE THIS IF YOU NEED TO BLIND YOUR DATA
        verbosity=1,
        condor=False
    )

if __name__ == "__main__":
    make_workspace()
    perform_fit(signal='MX3000_MY300', extra='--robustHesse 1')
    plot_fit(signal='MX3000_MY300')
    GOF(signal='MX3000_MY300',condor=False)
    plot_GOF(signal='MX3000_MY300',condor=False)
    
