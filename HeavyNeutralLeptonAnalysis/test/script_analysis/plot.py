#!/usr/bin/env python
import yaml
import os, re

#input_data='data_2017_4mu.txt'
#input_bkg='bkg_4mu_2017_new.txt'

with open('plotter.yml','r') as f:
    doc = yaml.load(f)
    histos = doc['histos'].keys()
       
for histo in histos:

    input_histo = doc['histos'][histo]['name_histo']
    input_Nbins = doc['histos'][histo]['Nbins']
    input_Xmin = doc['histos'][histo]['Xmin']
    input_Xmax = doc['histos'][histo]['Xmax']
    input_name_canvas = histo
    input_labelX = doc['histos'][histo]['labelX']
    input_labelY = doc['histos'][histo]['labelY']
    input_frame_Nbins = doc['histos'][histo]['frame_Nbins']
    input_frame_Xmin = doc['histos'][histo]['frame_Xmin']
    input_frame_Xmax = doc['histos'][histo]['frame_Xmax']
    input_rebin = doc['histos'][histo]['rebin']


    os.system("cat bkg_all_new.py | sed 's?Nbins?"+str(input_Nbins)+"?g' | sed 's?Xmin?"+str(input_Xmin)+"?g' | sed 's?Xmax?"+str(input_Xmax)+"?g' > bkg_confing.py")
    os.system("cat plotter_new.py | sed 's?bkg_all?bkg_confing?g' | sed 's?histos_name?"+input_histo+"?g' | sed 's?Nbins?"+str(input_Nbins)+"?g' | sed 's?Xmin?"+str(input_Xmin)+"?g' | sed 's?Xmax?"+str(input_Xmax)+"?g' | sed 's?frame_bins?"+str(input_frame_Nbins)+"?g' | sed 's?frame_min?"+str(input_frame_Xmin)+"?g' | sed 's?frame_max?"+str(input_frame_Xmax)+"?g' | sed 's?name_canvas?"+input_name_canvas+"?g' | sed 's?labelX?"+input_labelX+"?g' | sed 's?labelY?"+input_labelY+"?g' | sed 's?num_rebin?"+str(input_rebin)+"?g' > plotter_exe.py")
    os.system("python -i plotter_exe.py")

