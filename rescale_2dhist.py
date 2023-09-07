import ROOT

# open files
QCD = ROOT.TFile.Open('correctHistos/2dhist_SRTToHadronic.root','READ')
#QCD = ROOT.TFile.Open('correctHistos/2dhist_TightCRMX3000_MY300.root','READ')

# get bkg distributions 
QCDpass = QCD.Get('Pass')
QCDfail = QCD.Get('Fail')
QCDloose = QCD.Get('Loose')

# TH2::SetDirectory(0) if you don't want a histogram to be added to any directory, so that
# you can open and close root files without garbage collecting the histogram
QCDpass.SetDirectory(0)
QCDloose.SetDirectory(0)
QCDfail.SetDirectory(0)


# close the files so they won't get written to accidentally
QCD.Close()
fghruns_2016 = 19.528/137.6
#scale = ((3.57606*(10**6))/65337))*fghruns_2016 #1500 to 2000 QCD
#scale =((1.87317*(10**7))/100484)*fghruns_2016 #1000 to 1500 QCD
scale =  (145268./76035)*fghruns_2016 #TTToHadronic Scaling
#scale = (195./91141)*fghruns_2016 #MX3000_MY300 scaling
# create output file 
#out = ROOT.TFile.Open('2dhist_TightCRMX3000_MY300_rescaled.root','RECREATE')
out = ROOT.TFile.Open('2dhist_CorrectSRTTToHadronic_rescaled.root','RECREATE')

# create new Pass and Fail histos
newFail = QCDfail.Clone('Fail')	# clone so that we get same binning and title
newFail.Scale(scale)

newLoose = QCDloose.Clone('Loose')	# clone so that we get same binning and title
newLoose.Scale(scale)

newPass = QCDpass.Clone('Pass')
newPass.Scale(scale)
#newPass.Reset()
#newPass.SetDirectory(0)
#newPass.Add(QCDpass)

# cd to output file to write histos
out.cd()
newPass.Write()
newFail.Write()
newLoose.Write()

out.Close()
