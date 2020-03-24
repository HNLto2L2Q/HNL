#!/usr/bin/env python
import os, re
import commands
import math, time
import sys
import subprocess
from sig_list import *

queque="cms"

for items, infos in sig_files.items():
   txt_files=infos[0]
   directory_output=infos[1]

   number = subprocess.check_output("cat /user/atalierc/CMSSW_9_4_13_patch4/src/HNL/HeavyNeutralLeptonAnalysis/test/skimmed_signal/"+txt_files+" | wc -l", shell = True)
   interval = 100
   NumberOfJobs = (int(number.strip())/interval) + 1
   #NumberOfJobs = int(number.strip())

   File='/user/atalierc/CMSSW_9_4_13_patch4/src/HNL/HeavyNeutralLeptonAnalysis/test/skimmed_signal/'+txt_files
   Name=directory_output

   print File
   print Name

   for x in range(1, int(NumberOfJobs)+1):

      print x
      os.system("sed '"+str(1+interval*(x-1))+","+str(interval*x)+"!d' "+File+" > list_"+Name+str(x)+".txt ")
      #   print x
      os.system("cat CloneTree_new.C | sed 's?sample?"+Name+"?g' > CloneTree_new"+Name+str(x)+".C")
      #   print x
      os.system("cat read_condor.sh | sed 's?CloneTree_new?CloneTree_new"+Name+str(x)+"?g' | sed 's?sample?"+Name+"?g' |  sed 's?input_root?list_"+Name+str(x)+".txt?g' | sed 's?dir_skimmed?SKIMMER_"+Name+"?g' | sed 's?dir_jobs?jobs_new_"+Name+"?g' | sed 's?dir_list?list_new_"+Name+"?g' > skimmer_"+Name+str(x)+".sh")
      #   print x
      #os.system("cat condor_ex.cfg | sed 's?input_sh?skimmer_"+Name+str(x)+".sh?g' | sed 's?input.txt?"+File+"?g'> condor_ex_new_"+Name+str(x)+".cfg")
      #   print x
      #os.system("condor_submit -name ettore condor_ex_new_"+Name+str(x)+".cfg")
      os.system("qsub -q "+queque+" -o logs.stdout -e logs.stderr "+"skimmer_"+Name+str(x)+".sh")
      #print "job nr " + str(x) + " submitted"
      print "qsub -q "+queque+" -o logs.stdout -e logs.stderr job.sh"