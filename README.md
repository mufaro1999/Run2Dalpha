# Run2Dalpha
  It's pretty straightforward to run 2dalphabet. You write a python script that takes in a json file with the locations of all the input histograms
  that you need. The process of making the nominal and systematics histograms for this analysis is detailed in my other repo: https://github.com/mufaro1999/TagNTrain/tree/master/plotting .
  Here I ran 2dAlphabet in different ways. Since I split my histos into 3 separate fail, loose and pass regions, I needed to get 2 deifferent transfer functions to go from fail to loose and go from loose to pass
  in the background estimate. I ran TH_fail_loose.py which took in AN_failtoloose.json and ran this with multiple transfer functions while doing Ftests and GoF to determine which tranfer function to pick. It spit out the results to XHY_CRfl.. After getting
  that transfer function I then do this process with TH_loose_pass.py and AN_loosetopass.json to get the best transfer function to use here.It spit out the results to XHY_CRlp. 
  After getting these I then run the full fail to loose to pass fit which is test.py with the analysis_note_CR.json. It spit out the result to XHY_fits.
