import ROOT

# open files
QCD = ROOT.TFile.Open('2dhist_SRQCD_unscaled.root','READ')
tt = ROOT.TFile.Open('2dhist_SRTTToHadronic.root','READ')

# get bkg distributions 
QCDpass = QCD.Get('Pass')
QCDfail = QCD.Get('Fail')
QCDloose = QCD.Get('Loose')
ttpass = tt.Get('Pass')
ttfail = tt.Get('Fail')
ttloose = tt.Get('Loose')

# TH2::SetDirectory(0) if you don't want a histogram to be added to any directory, so that
# you can open and close root files without garbage collecting the histogram
QCDpass.SetDirectory(0)
QCDfail.SetDirectory(0)
QCDloose.SetDirectory(0)
ttpass.SetDirectory(0)
ttfail.SetDirectory(0)
ttloose.SetDirectory(0) 

# close the files so they won't get written to accidentally
QCD.Close()
tt.Close()

# create output file 
out = ROOT.TFile.Open('2dhist_data_unscaled.root','RECREATE')

# create new Pass and Fail histos
newFail = QCDfail.Clone('Fail')	# clone so that we get same binning and title
newFail.Reset()	# empty the histo
newFail.SetDirectory(0)
newFail.Add(QCDfail)	# add QCD fail distribution
newFail.Add(ttfail)	# add ttbar fail distribution
#newFail.Scale(0.2)

# create new Pass and Fail histos
newLoose = QCDloose.Clone('Loose')	# clone so that we get same binning and title
newLoose.Reset()	# empty the histo
newLoose.SetDirectory(0)
newLoose.Add(QCDloose)	# add QCD fail distribution
newLoose.Add(ttloose)	# add ttbar fail distribution
#newFail.Scale(0.2)

newPass = QCDpass.Clone('Pass')
newPass.Reset()
newPass.SetDirectory(0)
newPass.Add(QCDpass)
newPass.Add(ttpass)
#newPass.Scale(0.2)
# cd to output file to write histos
out.cd()
newPass.Write()
newLoose.Write()
newFail.Write()


out.Close()