import os
import yaml
import glob
import re
from datetime import datetime

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--ymls", type=str, nargs='*', required=True, help="yml files" )
parser.add_argument("--outputPath", type=str, required=True, help="output path" )
parser.add_argument("--cfgfile", type=str, required=True, help="path to cfgfile to be used" )
parser.add_argument("--masses", type=int, choices=[1,2,3,4,5,6,8,10,15,20], nargs='*',
                                help="HNL mass value [1,2,3,4,5,6,8,10,15,20]" )
parser.add_argument("--newIVF", action='store_true', default=False, help="Process events with the modify IVF")


now = datetime.now()
time = now.strftime("%y%m%d")

args = parser.parse_args()
ymls = args.ymls
masses = args.masses
cfgfile = os.path.abspath(args.cfgfile)
newIVF = args.newIVF

tag = "newIVF_{}".format(time) if newIVF else "usualIVF_{}".format(time)


print("Running on cfg file:{cfg}".format(cfg=cfgfile))

def createFolder(path):
    if os.path.isdir(path) == False:
        if path.find("pnfs") > -1:
            print("Creating folder on pnfs....")
            ret = os.system('ssh $USER@rw.iihe.ac.be \'mkdir -p {}\''.format(path))
            if ret:
                print("..Creation failed!!")
                exit(69)
            else:
                print("...Creation done!!")
        else:
            os.system('mkdir -p ' + path)



inputs_folder  = "FARM/inputs/"
logs_folder     = "FARM/logs/"
outputs_folder = os.path.abspath(args.outputPath)
cmssw_base = os.getenv('CMSSW_BASE')
if not cmssw_base:
    print("Please, set your CMSSW environment before")
    exit(12)

cmssw_src_folder = cmssw_base + "/src/"

bash_script_tmp = """#!/bin/bash
set -e

source $VO_CMS_SW_DIR/cmsset_default.sh
pushd {cmssw_src}
eval `scramv1 runtime -sh`
popd

{command}

{copy_cmd} {outFile} {output_path}
"""
for yml in ymls:
    with open(yml,'r') as f:

        pwd = os.getcwd()

        yml_name = yml.split('.')[0]

        process_folder_name = inputs_folder + yml_name + "/" + tag
        logs_folder_name    = logs_folder   + yml_name + "/" + tag
        createFolder(process_folder_name)
        createFolder(logs_folder_name)

        samples = yaml.load(f)

        # print("{folder}/skimming_submit.sh".format(folder=inputs_folder+yml_name))
        # print("input script: {folder}".format(folder=inputs_folder+yml_name))

        sub_script = open("{folder}/submit_all.sh".format(folder=process_folder_name),"w")
        sub_script.write("#!/bin/bash\n")

        for sample in samples:
            options = samples[sample]
            sample_path = options['path']
            mass        = options['mass']

            if len(masses) != 0 and not masses.count(mass):
                continue

            print("Processig: {sample}".format(sample=sample))

            list_of_files = glob.glob("{}/*".format(sample_path))
            if len(list_of_files) == 1:
                print(" - This sample contains to few Events (only 1 root file available)")
                print(" - Skipped!)")
                continue

            print(".... Number of jobs: {jobs}".format(jobs=len(list_of_files)))

            sample_scripts_folder = inputs_folder  + yml_name + "/" + tag + "/" + sample
            sample_logs_folder    = logs_folder  + yml_name + "/" + tag + "/" + sample
            output_scripts_folder = outputs_folder + "/" + yml_name + "/" + tag + "/" + sample + "/"
            createFolder(sample_scripts_folder)
            createFolder(sample_logs_folder)

            createFolder(output_scripts_folder)

            copy_cmd = "scp" if output_scripts_folder.find("pnfs") > -1 else "mv"
            copy_destination = "$USER@rw.iihe.ac.be:{}".format(output_scripts_folder) if output_scripts_folder.find("pnfs") > -1 else output_scripts_folder

            qsub_file = open("{folder}/Submit_{sample}.sh".format(folder=sample_scripts_folder, sample=sample),"w")
            qsub_file.write("#!/bin/bash\n")

            #print(list_of_files)
            for file in list_of_files:
                idx = re.findall(r"\d+",file)[-1]
                shell_filename = "bash_script_{idx}.sh".format(idx=str(idx))
                shell_script   = open("{folder}/{sh}".format(folder=sample_scripts_folder,sh=shell_filename), "w")
                command = "cmsRun {file_cfg} inputFile={path} outputFile={outFile} newIVF={flag}".format(file_cfg=cfgfile, path=file, outFile="output_"+str(idx)+".root", flag=("True" if newIVF else "False"))
                #print(command)
                script = bash_script_tmp.format(cmssw_src=cmssw_src_folder, command=command, copy_cmd=copy_cmd, outFile="output_"+str(idx)+".root", output_path=copy_destination)
                shell_script.write(script)
                shell_script.close()
                os.system("chmod 777 {folder}/{sh}".format(folder=sample_scripts_folder,sh=shell_filename))
                qsub_file.write("qsub -o {log} -e {log} {sh} \n".format(log=sample_logs_folder,sh=sample_scripts_folder+"/"+shell_filename))
            qsub_file.close()
            os.system("chmod 777 {folder}/Submit_{sample}.sh".format(folder=sample_scripts_folder, sample=sample))
            sub_script.write("bash {folder}/Submit_{sample}.sh\n".format(folder=sample_scripts_folder, sample=sample))
        sub_script.close()
