import re
import ROOT
import sys
import os
import random
from makeRPF import makeRPF


# for constructing data-like toys in the SR from the RPF shapes generated in the CR FLT fit


# global variables needed for makeRPF()
fitDir = '/uscms/home/mchitoto/nobackup/XtoHY/CMSSW_10_6_14/src/XHY_fits'
rpfL_params = '{}/MX3000_MY300_area/rpf_params_Background_rpfL_fitb.txt'.format(fitDir)
rpfT_params = '{}/MX3000_MY300_area/rpf_params_Background_rpfT_fitb.txt'.format(fitDir)
# binnings for the ratio shapes in makeRPF() - {'X':{MIN,MAX,NBINS},'Y':{MIN,MAX,NBINS}}
binnings = {
    'X': {'MIN':0,'MAX':500,'NBINS':50},
    'Y': {'MIN':1000,'MAX':4000,'NBINS':50}
}

def constructBkg(bkg, region, HT, tagger='particleNet'):
    '''
    if region == 'Tight':
        region = 'Pass'
    
    bkg [str] = 'ttbar', 'WJets', 'ZJets'
    region [str] = 'Loose', 'Tight'
    tagger [str] = 'particleNet', 'deepTag'

    sums up the total for a given background in a given region, over all years.
    Ideally, this should be used in subtractBkg(), but I created this function for generatePDF() after having written subtractBkg()
    might add that later.

    returns TH2D of total bavkground for all years in the given region
    '''
    histName = region   # SR DATA HISTO
    rootfile_dir = '/uscms/home/mchitoto/nobackup/XtoHY/CMSSW_10_6_14/src/'
    # loop through all years
    print('Generating total {} in SR {}'.format(bkg,region))
    totalBkg = None
    for year in ['16APV']:
        print('Acquring {} {} - {}'.format(bkg,year,histName))
        fileName = '2dhist_CorrectSRTTToHadronic_rescaled.root'
        f = ROOT.TFile.Open(fileName)
        bkgTemp = f.Get(histName)
        if not totalBkg:
            # the loop has just started, set up this year's bkg in the given region as our current total bkg
            print(region)
            totalBkg = bkgTemp.Clone(region)
            totalBkg.SetDirectory(0)
        else:
            # totalBkg already exists, simply append to it
            totalBkg.Add(bkgTemp)
        totalBkg.SetDirectory(0)
        f.Close()
    totalBkg.SetDirectory(0)
    return totalBkg

def subtractBkg(HT):
    '''
    takes data from SR fail and subtracts ttbar and V+Jets bkgs from it

    requires the HT cut that was performed when making the selection

    returns a TH2D of the data minus background (giving QCD est) in SR fail
    '''
    tagger = 'particleNet'
    region = 'Fail'
    rootfile_dir = '/uscms/home/ammitra/nobackup/TopPhi_analysis/CMSSW_11_1_4/src/TopHBoostedAllHad/rootfiles'
    dataFile = 'THselection_HT{}_Data_Run2.root'.format(HT)
    
    print('Subtracting backgrounds from data in {}...'.format(region))

    # for this analysis we are considering ttbar (combined semilep+allhad) and V+jets
    # different analyses could modify this dict and the following loop logic
    bkg_SR_fail = {'ttbar':None, 'Data_Run2':None}
    # loop through all backgrounds to be subtracted
    for bkg in bkg_SR_fail.keys():
        # loop through all years
        for year in ['16APV']:
            # get m_t vs m_phi(H) SR fail histogram (nominal) - same hist name for all files
            histName = 'Fail'
            if ('Data' in bkg):
                # data has been added to this bkg dict, just for debug purposes later
                fileName = '2dhist_CorrectSRJetHT2016.root'
                print('Opening Run2 Data file {}'.format(fileName))
                f = ROOT.TFile.Open(fileName)
                print('Acquiring histogram {}'.format(histName))
                bkg_temp = f.Get('Fail')
                bkg_SR_fail[bkg] = bkg_temp.Clone('Fail')
                bkg_SR_fail[bkg].Reset()
                bkg_SR_fail[bkg].SetDirectory(0)
                bkg_SR_fail[bkg].Add(bkg_temp)
                f.Close()
            else:
                # open the specific background file (nominal) for given year.
                # file will look like, e.g. ..../THselection_ttbar_16APV.root 
                fileName = '2dhist_CorrectSRTTToHadronic_rescaled.root'
                print('Opening {}'.format(fileName))
                f = ROOT.TFile.Open(fileName)
                print('Acquiring histogram {}'.format(histName))
                bkg_temp = f.Get(histName)
                # check if bkg histogram exists yet in our dict
                if not bkg_SR_fail[bkg]:
                    bkg_SR_fail[bkg] = bkg_temp.Clone(region)
                    bkg_SR_fail[bkg].Reset()
                    bkg_SR_fail[bkg].SetDirectory(0)
                bkg_SR_fail[bkg].Add(bkg_temp)
                f.Close()

    # at this point the dict should contain the bkgs and data in SR_fail. Now just subtract
    dataMinusBkg = bkg_SR_fail['Data_Run2'].Clone(region)
    print('Performing bkg subtraction from Run2 Data')
    for bkg, hist in bkg_SR_fail.items():
        if 'Data' in bkg:
            pass
        else:   # perform subtraction via Add() TH2D method
            print('Subtracting {} from Run2 data'.format(bkg))
            dataMinusBkg.Add(hist, -1.)
            # test bkgs subtracted individually, just for debug
            dataMinusIndividualBkg = bkg_SR_fail['Data_Run2'].Clone('dataMinus{}_only'.format(bkg))
            dataMinusIndividualBkg.Add(hist, -1.)
            bkg_SR_fail['DataMinus{}'.format(bkg)] = dataMinusIndividualBkg

    # add data minus bkgs hist to dict for debug
    bkg_SR_fail['Data_minus_background'] = dataMinusBkg

    # debug file to check that results make sense
    test = ROOT.TFile.Open('bkgs_SR_fail.root','RECREATE')
    test.cd()
    for bkg,hist in bkg_SR_fail.items():
        cloned = hist.Clone('{}_Fail'.format(bkg))
        cloned.SetDirectory(0)
        cloned.Write()
    test.Close()
    
    # return TH2D corresponding to QCD estimate in SR fail
    return bkg_SR_fail['Data_minus_background']

def Multiply(h1, h2, name=''):
    '''
    h1, h2 = TH2Ds with different binning (nBins_h1 > nBins_h2)
    namne [str] = name of output hist
    Loops over bins in h1, get the corresponding value at that coordinate in h2.
    Then, multiply the value at h2 coord by the bin value in h1
    returns TH2D corresponding to h1 * h2
    '''
    print('Multiplying {} x {}'.format(h1.GetName(),h2.GetName()))
    nh1 = (h1.GetNbinsX(),h1.GetNbinsY())
    nh2 = (h2.GetNbinsX(),h2.GetNbinsY())
    print('{} binning : {}\n{} binning : {}'.format(h1,nh1,h2,nh2))
    finalHist = h1.Clone(name)
    # loop over bins in h1
    for i in range(0, h1.GetNbinsX()+1):
        for j in range(0, h1.GetNbinsY()+1):
            h1Val = h1.GetBinContent(i,j)
            if h1Val < 0:     # don't allow negative yields
                h1Val = 0
            xVal = h1.GetXaxis().GetBinCenter(i)
            yVal = h1.GetYaxis().GetBinCenter(j)
            ih2 = h2.GetXaxis().FindBin(xVal)
            jh2 = h2.GetYaxis().FindBin(yVal)
            h2Val = h2.GetBinContent(ih2,jh2)
            finalHist.SetBinContent(i,j,h1Val*h2Val)
    finalHist.SetDirectory(0)
    return finalHist

def getRatioHist(fileName,ratioName):
    '''
    gets the ratio histogram from the file, if it's already been created by makeRPF()
    fileName [str] = name of ratio file
    ratioName [str] = name of histogram in ratio file
    '''
    f = ROOT.TFile.Open(fileName)
    h = f.Get(ratioName)
    h.SetDirectory(0)
    f.Close()
    return h

def constructQCD(HT):
    '''
    construct QCD estimate in SR loose and pass regions from data-bkg (QCD) in SR fail
    requires the HT cut made on the selection
    '''
    # first, generate R_LF and R_TL ratios and get histograms (will also store in root file)
    # this step requires making a new 2DAlphabet workspace with the proper binning (see makeRPF.py), which is time consuming
    # therefore, we do a check to see if the temporary workspace directories (made during makeRPF()) exist yet, and if so we just skip this step.
    isDir_rpfL = os.path.isdir('tmpDir_rpfL')
    isDir_rpfT = os.path.isdir('tmpDir_rpfT')
    print(isDir_rpfL)
    if isDir_rpfL:
        rpfL = getRatioHist('rpfL.root','rpfL')
    else:
        print('rpfT not created yet, creating ratio shape now...')
        rpfL = makeRPF(fitDir,'tmpDir',rpfL_params,binnings,'1x1')
    if isDir_rpfT:
        rpfT = getRatioHist('rpfT.root','rpfT')
    else:
        print('rpfT not created yet, creating ratio shape now...')
        rpfT = makeRPF(fitDir,'tmpDir',rpfT_params,binnings,'1x0')

    # next, perform background subtraction (will also generate data minus bkg root file)
    QCD_SR_Fail = subtractBkg(HT)

    # book an output root file
    outFile = ROOT.TFile.Open('QCD_SR_distributions.root','RECREATE')
    outFile.cd()
    # get QCD distribution in SR Loose and Tight, write histograms to outFile
    qcdFail = QCD_SR_Fail.Clone('QCD_Fail')
    qcdFail.SetDirectory(0)
    qcdFail.Write()
    qcdLoose = qcdFail.Clone('QCD_Loose')
    qcdLoose = Multiply(qcdFail, rpfL, 'QCD_Loose')    
    #qcdLoose.Multiply(rpfL)    # I'll use my own multiply() method instead of built-in one so that we can zero negative bins
    qcdLoose.SetDirectory(0)
    qcdLoose.Write()
    qcdL_tmp = qcdLoose.Clone('qcdL_tmp')
    qcdTight = Multiply(qcdL_tmp, rpfT, 'QCD_Tight')
    #qcdTight.Multiply(rpfT)    # again, Multiply() TH2D method would work for same-binned histos, but want to zero negative bins so use my own implementation
    qcdTight.SetDirectory(0)
    qcdTight.Write()
    outFile.Close()

    # return qcd estimate in all three regions, just in case you want to use this function later
    return (qcdFail, qcdLoose, qcdTight)

def getCumulativePDF(pdf, pdf_name):
    print('Creating cumulative PDF : {}'.format(pdf_name))
    nx = pdf.GetNbinsX()
    ny = pdf.GetNbinsY()
    cPDF = ROOT.TH1F(pdf_name,"",nx*ny,0,nx*ny)
    print('input:  {} - {}\noutput: {} - {}'.format(type(pdf),pdf.GetName(),type(cPDF),pdf_name))
    cumulativeBin = 0
    for i in range(1,nx+1):
        for j in range(1,ny+1):
            cumulativeBin += 1
            pdfVal = pdf.GetBinContent(i,j)+cPDF.GetBinContent(cumulativeBin-1)
            cPDF.SetBinContent(cumulativeBin,pdfVal)
    return cPDF

def generatePDF(region, HT):
    '''
    generates PDFs from ttbar/V+Jet simulation and QCD from ratios in given region
    region [str] = 'Loose', 'Tight'
    HT (int) = HT used to make the selection
    returns: scaled PDF in region (TH2D), cumulativePDF in region (TH2D), number of events from estimate
    '''
    assert(region=='Loose' or region=='Tight')
    # get estimated QCD Loose, Tight distributions in SR (Fail is included, but not necessary - hence '_')
    _, qcdLoose, qcdTight = constructQCD(HT)
    # get the distributions of bkgs for Loose, Tight regions and store in dict
    PDFs = {
        region:{'ttbar':None,'QCD_scaled':None,'cumulativePDF':None,'nRegion':None}
        }
    if (region=='Loose'):
        PDFs[region]['QCD']=qcdLoose.Clone('QCD_Loose')
    else:
        PDFs[region]['QCD'] = qcdTight.Clone('QCD_Tight')
    





    for region, backgrounds in PDFs.items():

        print('Generating PDF in {}'.format(region))
        for bkg in backgrounds:
            if ('QCD' in bkg) or ('cumulative' in bkg) or ('nRegion' in bkg):
                continue
            else:
                # get the cumulative (non-QCD) background from simulation in the given region
                totalBkg = constructBkg(bkg,region,HT)
                totalBkg.SetDirectory(0)
                PDFs[region][bkg] = totalBkg
                PDFs[region][bkg].SetDirectory(0)
                # add it to the QCD in the given region
                PDFs[region]['QCD'].Add(totalBkg)
                PDFs[region]['QCD'].SetDirectory(0)

    for region in PDFs.keys():
        # grab the QCD distribution for this region
        pdf = PDFs[region]['QCD'].Clone('pdf_{}'.format(region))
        pdf.SetDirectory(0)
        nRegion = pdf.Integral(1,pdf.GetNbinsX()+1,1,pdf.GetNbinsY()+1)
        PDFs[region]['nRegion'] = nRegion
        print('nEvents in {} from pdf: {}'.format(region,nRegion))
        print('Scaling PDF - {}'.format(region))
        pdf.Scale(1/nRegion)
        pdf.SetDirectory(0)
        PDFs[region]['QCD_scaled'] = pdf
        # get cumulative PDF
        cumulativePDF = getCumulativePDF(pdf,'cumulative_pdf_{}'.format(region))
        cumulativePDF.SetDirectory(0)
        PDFs[region]['cumulativePDF'] = cumulativePDF

    # debug - send to output file
    oFile = ROOT.TFile.Open('PDF_info_{}.root'.format(region),'RECREATE')
    for region, backgrounds in PDFs.items():
        for bkg in backgrounds:
            if bkg == 'QCD':
                oHist = PDFs[region][bkg].Clone('PDF_{}'.format(region))
            elif 'scaled' in bkg:
                oHist = PDFs[region][bkg].Clone('PDF_{}_scaled'.format(region))
            elif 'cumulative' in bkg:
                oHist = PDFs[region][bkg].Clone('PDF_{}_cumulative'.format(region))
            elif bkg == 'nRegion':
                continue
            else:
                oHist = PDFs[region][bkg].Clone('{}_fromSimulation_{}'.format(bkg,region))
            oHist.SetDirectory(0)
            oHist.Write()
    oFile.Close()

    return PDFs[region]['QCD_scaled'], PDFs[region]['cumulativePDF'], PDFs[region]['nRegion']


def findPDFintersection(rand,cumulativePDF):
    '''
    rand [float] - random float [0, 1)
    cumulativePDF [TH1F]
    returns global bin for PDF intersection
    '''
    intFound = False
    intBin = None
    for i in range(1,cumulativePDF.GetNbinsX()+1):
        pdfVal = cumulativePDF.GetBinContent(i)
        if (pdfVal >= rand):
            #return i
            intFound = True
            intBin = i
            break
    if intFound:
        return intBin
    else:
        print("Intersection with PDF not found, something is wrong")
        return intBin
    
def globalBinTo2D(pdf,globalBin):
    '''
    globalBins start from 1
    globalBin 1 = (1,1), 2 = (1,2) and so on
    '''
    globalBin = globalBin - 1
    NX = pdf.GetNbinsX()
    NY = pdf.GetNbinsY()
    nx = int(globalBin)/int(NY)+1
    ny = globalBin%NY+1
    return nx,ny

def generateToy(pdf, cumulativePDF, nEvents, name):
    '''
    generates toy data from the cumulative PDF
    '''
    toy = pdf.Clone(name)
    toy.Reset()
    toy.SetDirectory(0)
    maxCDF =cumulativePDF.GetMaximum()
    print('max cdf is {}'.format(maxCDF))
    for i in range(int(nEvents)):
        rand = random.uniform(0,maxCDF)
        globalBin = findPDFintersection(rand,cumulativePDF)
        nx, ny = globalBinTo2D(pdf,globalBin)
        nx = int(nx)
        ny = int(ny)
        toy.SetBinContent(nx,ny,toy.GetBinContent(nx,ny)+1)
    return toy

def getFailTemplate(rFile):
    '''
    Gets the actual data distribution in the SR fail region, so that it can be used for 2DAlphabet when generating the fits
    rFile [str] = /path/to/TIMBER/selection/DataRun2/file
    '''
    f = ROOT.TFile.Open(rFile)
    h = f.Get('Fail')
    h.SetDirectory(0)
    f.Close()
    return h	

# -----------------------------------------------------------------
# MAIN LOOP
# -----------------------------------------------------------------
if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('--HT', type=str, dest='HT',
                        action='store', required=True,
                        help='Value of HT to cut on')
    args = parser.parse_args()

    out = {
	'Loose': {'pdf':None, 'cumulative_pdf':None, 'nEvents':None, 'toy':None},
        'Tight': {'pdf':None, 'cumulative_pdf':None, 'nEvents':None, 'toy':None}
    }
    for region in out.keys():
        pdf, cumulativePDF, nEvents = generatePDF(region, args.HT)
	# we want the histogram to have the same name as the actual TIMBER data selection files, for ease of use when running 2DAlpabet.
        toy = generateToy(pdf,cumulativePDF,nEvents,region)
        out[region]['pdf'] = pdf
        out[region]['cumulative_pdf'] = cumulativePDF
        out[region]['nEvents'] = nEvents
        out[region]['toy'] = toy

    # get the SR fail data distribution
    SR_fail = getFailTemplate('/uscms/home/mchitoto/nobackup/XtoHY/CMSSW_10_6_14/src/2dhist_CorrectSRJetHT2016.root')

    f = ROOT.TFile.Open('2dhist_30ToySRData2016.root','RECREATE')
    f.cd()
    for region, elements in out.items():
        for element in elements:
            if element == 'nEvents':
                continue
            if ((element =='toy') & (region =='Tight')):
                out[region][element].SetDirectory(0)
                out[region][element].Write('Pass')
                continue
            out[region][element].SetDirectory(0)
            out[region][element].Write()
    SR_fail.Write()
    f.Close()



    