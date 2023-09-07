from TwoDAlphabet.twoDalphabet import TwoDAlphabet
from TwoDAlphabet.alphawrap import BinnedDistribution, ParametricFunction
from TwoDAlphabet.helpers import open_json
import json
import ROOT
from collections import OrderedDict
import os

def _generate_constraints(nparams):
    out = {}
    for i in range(nparams):
        if i == 0:
            out[i] = {"MIN":0,"MAX":1}
        else:
            out[i] = {"MIN":-5,"MAX":5}
    return out

def _get_other_region_names(pass_reg_name):
    return pass_reg_name, pass_reg_name.replace('Fail','Loose'),pass_reg_name.replace('Fail','Pass')
_rpf_options = {
    '0x0': {
        'form': '0.1*(@0)',
        'constraints': _generate_constraints(1)
    },
    '1x0': {
        'form': '0.1*(@0+@1*x)',
        'constraints': _generate_constraints(2)
    },
    '0x1': {
        'form': '0.1*(@0+@1*y)',
        'constraints': _generate_constraints(2)
    },
    '1x1': {
        'form': '0.1*(@0+@1*x)*(1+@2*y)',
        'constraints': _generate_constraints(3)
    },
    '1x2': {
        'form': '0.1*(@0+@1*x)*(1+@2*y+@3*y*y)',
        'constraints': _generate_constraints(4)
    },
    '2x1': {
        'form': '0.1*(@0+@1*x+@2*x**2)*(1+@3*y)',
        'constraints': _generate_constraints(4)
    },
    '2x2': {
        'form': '0.1*(@0+@1*x+@2*x**2)*(1+@3*y+@4*y**2)',
        'constraints': _generate_constraints(5)
    },
    '2x3': {
        'form': '0.1*(@0+@1*x+@2*x*x)*(1+@3*y+@4*y*y+@5*y*y*y)',
        'constraints': _generate_constraints(6)
    },
    '3x2': {
        'form': '0.1*(@0+@1*x+@2*x*x+@3*x*x*x)*(1+@4*y+@5*y*y)',
        'constraints': _generate_constraints(6)
    }
}

def _getParams(line):
    '''
    helper function for makeRPF(), acts on each line of the rpf_params .txt 
    file to get the final parameters. This way, whole file isn't saved in memory.
    param format:
            Background_CR_rpfT_par0: 0.999999999868 +/- 0.148048177441
    '''
    # get param name (everything before colon)
    name =  line.split(':')[0]
    # get actual param value
    param = line.split(':')[1].split('+/-')[0].strip()

    return (name, float(param))

def makeRPF(fitDir, tmpDir, paramFile, binnings, poly):
    '''
    fitDir [str] = '/path/to/2DAlphabet/workspace/dir/'
    tmpDir [str] = name of temporary directory for 2DAlphabet workspace that gets made
    paramFile [str] = '/path/to/rpf_params_fitb.txt'
    binnings [dict] = {'X':{MIN,MAX,NBINS},'Y':{MIN,MAX,NBINS}}
    poly [str] = RPF polynomial order, e.g. '2x2'

    * NOTE * assumes that the 2DAlphabet fit has already been performed 
    and the workspace exists.

    Gets the binning from the 2DAlphabet-generated workspace and the 
    RPF parameters from the 2DAlphabet-generated .txt file, then 
    constructs the ParametricFunction associated with those parameters.
    '''
    # hard-code the region detection as well. L = 'rpfL', T = 'rpfT'
    rpfName = 'rpfL' if 'rpfL' in paramFile else 'rpfT'

    # first, we need to ensure the Ratio shapes have the same binning as the TIMBER selection histograms
    print('Performing rebinning of X and Y axis for ratio shapes:')
    c = open_json('/uscms/home/mchitoto/nobackup/XtoHY/CMSSW_10_6_14/src/test4.json')
    print(binnings.items())
    for axis, vals in binnings.items():
        for val in vals.keys():
            print('val = {}').format(val)
            print('axis = {}').format(axis)
            print('c = {}').format(c)
            #print('Replacing {} {} : {} -> {}'.format(axis,val,c['BINNING']['default'][axis][val],binnings[axis][val]))
            c['BINNING']['default'][axis][val] = binnings[axis][val]
    # now, write this to a new json file.
    with open('temp_config.json','w') as json_file:
        json.dump(c, json_file, indent=4)

    # now, use this new json file when making the Ratio shapes
    # we have to create a new working directory, since the existing fit directories get the binning from binnings.p
    # If we specify loadPrevious=True, then it will simply use the old binnings. If we use loadPrevious=False, it'll overwrite the old FLT fit dir which is undesirable
    # so, we make a new 2DAlphabet workspace with the new binnings created in the previous step, and go from there.
    workspaceDir = tmpDir+'_'+rpfName       # the workspace created with the new binnings
    twoD = TwoDAlphabet(workspaceDir, 'temp_config.json', loadPrevious=False) 
    # for now we'll just hard-code this. We're looking in the CR, so we'll get the binning from 'CR_fail' region
    for f,l,p in [_get_other_region_names(r) for r in twoD.ledger.GetRegions() if 'Fail' in r]:
        print('regions detected in workspace: {}, {}, {}\n'.format(f,l,p))
        binning_f, _ = twoD.GetBinningFor(f)
    fail_name = 'Background_{}'.format(f)
    
    # construct the actual ParametricFunction
    rpf_func = ParametricFunction(
        fail_name.replace('Fail',rpfName),
        binning_f, _rpf_options[poly]['form'],
        constraints=_rpf_options[poly]['constraints']
    )

    # get rpf params from the .txt file
    params = OrderedDict()
    with open(paramFile) as f:
        for line in f:
            # output of _getParams() is a tuple (paramName, paramVal)
            paramName, paramVal = _getParams(line)
            params[paramName] = paramVal
    
    # quick check to ensure that the number of params matches the polynomial order
    nParams = _rpf_options[poly]['form'].count('@')
    assert(len(params.keys()) == nParams)

    # status update
    print('{} polynomial : {}\n'.format(rpfName, _rpf_options[poly]['form']))
    # update user and set parameter values in the rpf_func ParametricFunction object
    count = 0
    for iparam, param_val in params.items():
        print('Parameter @{} = {} = {}'.format(count, iparam, param_val))
        rpf_func.setFuncParam(iparam, param_val)
        count += 1

    # save the histogram in a root file
    rpfFile = ROOT.TFile.Open('{}.root'.format(rpfName), 'RECREATE')
    rpfFile.cd()

    # create and fill RPF shape histogram
    out_hist = rpf_func.binning.CreateHist(rpfName,cat='')
    print('{} NbinsX: {}'.format(rpfName,out_hist.GetNbinsX()))
    print('{} NbinsY: {}'.format(rpfName,out_hist.GetNbinsY()))
    for ix in range(1, out_hist.GetNbinsX()+1):
        for iy in range(1, out_hist.GetNbinsY()+1):
            out_hist.SetBinContent(ix, iy, rpf_func.getBinVal(ix, iy))

    out_hist.SetDirectory(0)
    out_hist.Write()
    rpfFile.Close()

    return out_hist