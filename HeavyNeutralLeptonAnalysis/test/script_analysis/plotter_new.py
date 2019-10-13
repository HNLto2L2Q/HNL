from ROOT import *
from bkg_all import *
from math import *
from numpy import *

ratio = True

histo= 'histos_name'
N_bins= Nbins
X_min= Xmin
X_max= Xmax

data = TH1D("", "", N_bins, X_min, X_max)
histo_bkg_all = TH1D("", "", N_bins, X_min, X_max)
frame = TH1D("", "", frame_bins, frame_min, frame_max)


for name, sample_name in files.items():
        
    inputFile=TFile(sample_name[0])
    mc=inputFile.Get(histo)#DiMuMass           
    ciola = 41557*sample_name[2]/sample_name[3]#23765
        
    sample_name[1].histo.Add(mc, ciola)



ciao = TCanvas("canvas", "canvas", 900,900,900,900)
ciao.cd()
ciao.SetLogy()

stack = THStack('stack','stack')

stack.Add(c_diBoson.histo)
stack.Add(c_singleTop.histo)
stack.Add(c_WJetsToLNu.histo)
stack.Add(c_TTJets_Dilept.histo)
stack.Add(c_TTJets_SingleLeptFromT.histo)
stack.Add(c_TTJets_SingleLeptFromTbar.histo)
stack.Add(c_DYJetsToLL_M_10_50_MLM.histo)
stack.Add(c_DYJetsToLL_M_50_FXFX.histo)

c_diBoson.histo.Rebin(num_rebin)
c_singleTop.histo.Rebin(num_rebin)
c_WJetsToLNu.histo.Rebin(num_rebin)
c_TTJets_Dilept.histo.Rebin(num_rebin)
c_TTJets_SingleLeptFromT.histo.Rebin(num_rebin)
c_TTJets_SingleLeptFromTbar.histo.Rebin(num_rebin)
c_DYJetsToLL_M_10_50_MLM.histo.Rebin(num_rebin)
c_DYJetsToLL_M_50_FXFX.histo.Rebin(num_rebin)

    
inputFile=TFile('new_data_all_check/merged.root')
data=inputFile.Get(histo)#DiMuMass
data.Rebin(num_rebin)

lista = TList()
lista.Add(c_diBoson.histo)
lista.Add(c_singleTop.histo)
lista.Add(c_WJetsToLNu.histo)
lista.Add(c_TTJets_Dilept.histo)
lista.Add(c_TTJets_SingleLeptFromT.histo)
lista.Add(c_TTJets_SingleLeptFromTbar.histo)
lista.Add(c_DYJetsToLL_M_10_50_MLM.histo)
lista.Add(c_DYJetsToLL_M_50_FXFX.histo)


histo_bkg_all.Merge(lista)

data.SetMarkerStyle(8);
data.SetMarkerColor(kBlack);
data.SetLineColor(kBlack);

frame.Draw()
frame.GetXaxis().SetTitle("labelX")
frame.GetYaxis().SetTitle("labelY")
frame.SetMaximum(10e+9)
frame.SetMinimum(10)
#frame.GetYaxis().SetRange(100, 10000000)
frame.SetStats(0)
stack.Draw("same HISTO")

data.Draw("EP same")

gPad.Update()
gPad.RedrawAxis()


legendmass4l = TLegend(0.47,0.68,0.67,0.89)
legendmass4l.SetTextSize(0.030)
legendmass4l.SetLineColor(0)
legendmass4l.SetLineWidth(1)
legendmass4l.SetFillColor(kWhite)
legendmass4l.SetBorderSize(0)
    
legendmass4l.AddEntry(c_singleTop.histo, c_singleTop.histo_name, "f")
legendmass4l.AddEntry(c_diBoson.histo, c_diBoson.histo_name, "f")
legendmass4l.AddEntry(c_WJetsToLNu.histo, c_WJetsToLNu.histo_name, "f")
legendmass4l.AddEntry(c_TTJets_Dilept.histo, c_TTJets_Dilept.histo_name, "f")
legendmass4l.AddEntry(c_TTJets_SingleLeptFromT.histo, c_TTJets_SingleLeptFromT.histo_name, "f")
legendmass4l.AddEntry(c_TTJets_SingleLeptFromTbar.histo, c_TTJets_SingleLeptFromTbar.histo_name, "f")
legendmass4l.AddEntry(c_DYJetsToLL_M_10_50_MLM.histo, c_DYJetsToLL_M_10_50_MLM.histo_name, "f")
legendmass4l.AddEntry(c_DYJetsToLL_M_50_FXFX.histo, c_DYJetsToLL_M_50_FXFX.histo_name, "f")


legendmass4l.AddEntry(data,"Data","lep")

legendmass4l.Draw("same")


if(ratio == True):
    frame2 = TH1D("","", frame_bins, frame_min, frame_max)
    frame2.SetMaximum(2)
    frame2.SetMinimum(0)
    frame2.SetStats(0)
    frame2.GetYaxis().SetTitle("Data / mc")
    frame2.GetYaxis().SetNdivisions(508)
    h3 = data.Clone("h3")
    h3.Sumw2()
    h3.Divide(histo_bkg_all)
    h3.SetMarkerStyle(8)
    h3.SetMarkerColor(kBlack)
    h3.SetLineColor(kBlack)
    canvasratio = 0.3
    ciao.SetBottomMargin(canvasratio + (1-canvasratio)*ciao.GetBottomMargin()-canvasratio*ciao.GetTopMargin())
    canvasratio = 0.16
    ratioPad =  TPad("BottomPad","",0,0,1,1)
    ratioPad.SetTopMargin((1-canvasratio) - (1-canvasratio)*ratioPad.GetBottomMargin()+canvasratio*ratioPad.GetTopMargin())
    ratioPad.SetFillStyle(4000)
    ratioPad.SetFillColor(4000)
    ratioPad.SetFrameFillColor(4000)
    ratioPad.SetFrameFillStyle(4000)
    ratioPad.SetFrameBorderMode(0)
    ratioPad.SetGrid(1,1)
    ratioPad.Draw("HISTO")
    #frame2.Draw("same HISTO")
    ratioPad.cd()
    frame2.Draw("same HISTO")
    h3.Draw("EPsame HISTO")
        
print "muori"

ciao.SaveAs('plots/name_canvas_try.png')


