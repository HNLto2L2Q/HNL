#!/usr/bin/env python
import os, re
import commands
import math, time
import sys
from  optparse  import OptionParser

print 
print 'START'
print 
#M_ele_M2_31.15 474
#M_ele_M2_40.54 335
#M_ele_M3_32.43 459
#M_ele_M4_16.14 398
#M_ele_M5_18.37 478
#M_ele_M6_12.48 118
#M_ele_M8_5.45 350
parser = OptionParser()
parser.add_option("-f","--file",dest="file",help="list of signals",action="store",type="string")
(options, args) = parser.parse_args()
########   customization  area #########
#NumberOfJobs= 62 # number of jobs to be submitted please check the number of jobs is correct
interval = 1
OutputFileNames = "Analysis_output" 
ScriptName = "Reco_1.py" 
OutputDir = options.file #"M_mu_M8_0.42"
FileList = "list_"+OutputDir+".list" 
count = 0
with open(FileList, 'r') as f:
    for line in f:
        count += 1
NumberOfJobs= count
queue = "localgrid "
########   customization end   #########

path = os.getcwd()
print
print 'do not worry about folder creation:'
os.system("rm -r tmp")
os.system("mkdir tmp")
os.system("mkdir "+OutputDir) # change the name of directroy and don't forget change it below 
print

##### loop for creating and sending jobs #####
for x in range(1, int(NumberOfJobs)+1):
   ##### creates directory and file list for job #######
   os.system("mkdir tmp/"+str(x))
   os.chdir("tmp/"+str(x))
   os.system("sed '"+str(1+interval*(x-1))+","+str(interval*x)+"!d' ../../"+FileList+" > list.txt ")
   
   ##### creates jobs #######
   with open('job.sh', 'w') as fout:
      fout.write("#!/bin/sh\n")
      fout.write("echo\n")
      fout.write("echo\n")
      fout.write("echo 'START---------------'\n")
      fout.write("echo 'WORKDIR ' ${PWD}\n")
      #fout.write("source /afs/cern.ch/cms/cmsset_default.sh\n")
      fout.write("cd "+str(path)+"\n")
      #fout.write("source $VO_CMS_SW_DIR/cmsset_default.sh")
      fout.write("eval `scram runtime -sh`")
      #fout.write("cmsenv\n")
      fout.write("export X509_USER_PROXY=/user/$USER/x509up_u23054\n")
      fout.write("cmsRun "+ScriptName+" outputFiles='"+OutputDir+"/"+OutputFileNames+"_"+str(x)+".root'  inputFile='tmp/"+str(x)+"/list.txt'\n")
      fout.write("echo 'STOP---------------'\n")
      fout.write("echo\n")
      fout.write("echo\n")
   os.system("chmod 755 job.sh")
   
   ###### sends bjobs ######
   os.system("qsub -q "+queue+" -o logs job.sh")
   print "job nr " + str(x) + " submitted"
   
   os.chdir("../..")
   
print
print "your jobs:"
os.system("qstat -u $USER")
print
print 'END'
print
