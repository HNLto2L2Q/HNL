#!/usr/bin/env python
import os, re
import commands
import math, time
import sys
import subprocess


number = subprocess.check_output("cat /lustre/home/taliercio/SL7/CMSSW_9_4_13_patch4/src/condor/prova.txt | wc -l", shell = True)
interval = 2
NumberOfJobs = (int(number.strip())/interval) + 1


File='/lustre/home/taliercio/SL7/CMSSW_9_4_13_patch4/src/condor/prova.txt'
Name='prova'

for x in range(1, int(NumberOfJobs)+1):

   print x
   os.system("sed '"+str(1+interval*(x-1))+","+str(interval*x)+"!d' "+File+" > list"+str(x)+".txt ")

   os.system("cat read_condor.sh | sed 's?input_root?list"+str(x)+".txt?g' | sed 's?dir_skimmed?SKIMMER_"+Name+"?g' | sed 's?dir_jobs?jobs_new_"+Name+"?g' | sed 's?dir_list?list_new_"+Name+"?g' > skimmer"+str(x)+".sh")
   os.system("cat condor_ex.cfg | sed 's?input_sh?skimmer"+str(x)+".sh?g' | sed 's?input.txt?"+File+"?g'> condor_ex_new_"+str(x)+".cfg")

   os.system("condor_submit -name ettore condor_ex_new_"+str(x)+".cfg")
