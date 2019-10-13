from ROOT import *

class Histo:
        def __init__(self, histo_name, N_bins, X_min, X_max,  color):
                self.histo=TH1D(histo_name,histo_name, N_bins, X_min, X_max)
                self.histo.SetFillColor(color)
                self.histo_name = histo_name
        def reset(self):
                self.histo=TH1D("", "", 0, 0, 0)
                self.histo.SetFillColor(0)
                self.histo_name = "None"
                
c_WJetsToLNu = Histo('WJets #rightarrow l #nu', Nbins, Xmin, Xmax, 2) #42
c_TTJets_Dilept = Histo('TTJets dilept', Nbins, Xmin, Xmax, 3) #kOrange
c_TTJets_SingleLeptFromT =  Histo('TTJets singleLeptFromT', Nbins, Xmin, Xmax, 4) #kRed
c_TTJets_SingleLeptFromTbar =  Histo('TTJets singleLeptFromTbar', Nbins, Xmin, Xmax, 5) #kBlue
c_DYJetsToLL_M_50_FXFX = Histo('DYJetsToLL_M-50_FXFX', Nbins, Xmin, Xmax, 6) #kGreen
c_DYJetsToLL_M_10_50_MLM = Histo('DYJetsToLL_M_10-50', Nbins, Xmin, Xmax, 7) #kCyan
c_diBoson = Histo('di-boson', Nbins, Xmin, Xmax, kOrange) #kMagenta
c_singleTop = Histo('single top', Nbins, Xmin, Xmax, 9) #2

files={
        'WJets_new' : ['WJets_total_2017_check/merged.root', c_WJetsToLNu, float(52940.0), float(74799469), 'mu'], #sum Files 74981747
        'TTJets_Dilept' : ['TTJets_DiLept_corrected_check/merged.root', c_TTJets_Dilept, float(54.23), float(28380110), 'mu'], #26697032
        'TTJets_SingleLeptFromT' : ['TTjets_SingleLeptFromT_new_check/merged.root', c_TTJets_SingleLeptFromT, float(109.6), float(58019135), 'mu'],#60030172#61077143
        'TTJets_SingleLeptFromTbar' : ['TTjets_SingleLeptFromTbar_new_check/merged.root', c_TTJets_SingleLeptFromTbar, float(109.6), float(44836885), 'mu'],#52680113#55108932#
        'DYJetsToLL_M-10-50_MLM' : ['DYJetsToLL_M_10to50_pu_trig_check/merged.root', c_DYJetsToLL_M_10_50_MLM, float(18610.0), float(39520441), 'mu'],
        'DY_new' : ['DYJest_M-50_total_2017_check/merged.root', c_DYJetsToLL_M_50_FXFX, float(6225.42), float(208046450), 'mu'], #sum files 207517695
        
        #di-boson
        
        'WZTo1L3Nu_13TeV' : ['WZTo1L3Nu_13TeV_pileup_trig_new_check/merged.root', c_diBoson, float(3.033), float(4994395), 'mu'],
        'ZGTo2MuG_MMuMu' : ['ZGTo2MuG_MMuMu_pileup_trig_new_check/merged.root', c_diBoson, float(117.864), float(137800), 'mu'],
        'ZZTo4L_13TeV' : ['ZZTo4L_13TeV_pileup_trig_new_check/merged.root', c_diBoson, float(1.256), float(15101000), 'mu'],
        'WGToLNuG' : ['WGToLNuG_pileup_trig_new_check/merged.root', c_diBoson, float(489), float(6283083), 'mu'],
        'WWTo2L2Nu' : ['WWTo2L2Nu_pileup_trig_new_check/merged.root', c_diBoson, float(11.08), float(16000), 'mu'],

        ##                'ZZTo2Q2Nu' : ['ZZTo2Q2Nu_corrected_check/merged.root', c_diBoson, float(4.04), float(49262784), 'mu'],
                
        'WWToLNuQQ' : ['WWToLNuQQ_total_new_check/merged.root', c_diBoson, float(45.99), float(18674262), 'mu'],
        'WZTo1L1Nu2Q' : ['WZTo1L1Nu2Q_pileup_trig_new_check/merged.root', c_diBoson, float(11.66), float(19086373), 'mu'],
        'WZTo3LNu' : ['WZTo3LNu_pileup_trig_new_check/merged.root', c_diBoson, float(5.052), float(10881711), 'mu'],
        #'ZZTo2L2Q' : ['ZZTo2L2Q_pileup_trig_new_check/merged.root', c_diBoson, float(3.688), float(17047050), 'mu'],
        'WZTo2L2Q' : ['WZTo2L2Q_pileup_trig_new_check/merged.root', c_diBoson, float(6.331), float(27503861), 'mu'],
                
        #single top
                
        'ST_tW_top_5f_inclusiveDecays' : ['ST_tW_top_5f_inclusiveDecays_pu_trig_new_check/merged.root', c_singleTop, float(34.91), float(7746324), 'mu'],
        'ST_tW_antitop_5f_inclusiveDecays' : ['ST_tW_antitop_5f_inclusiveDecays_pu_trig_new_check/merged.root', c_singleTop, float(34.97), float(7038570), 'mu'],
        'ST_s-channel_4f_leptonDecays' : ['ST_s-channel_4f_leptonDecays_check/merged.root', c_singleTop, float(3.74), float(7952142), 'mu']
}


