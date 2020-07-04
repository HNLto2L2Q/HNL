import os
import yaml
import glob
import re
from datetime import datetime


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--ymls", type=str, nargs='*', required=True, help="yml files" )
parser.add_argument("--sourcefile", type=str, required=True, help="path to source file to be used" )
parser.add_argument("--outputPath", type=str, required=True, help="output path" )
parser.add_argument("--onlyBackground", action='store_true', help='Submit skimming only for background samples')
parser.add_argument("--skipSignals", action='store_true', help='Submit skimming skipping signals samples')
parser.add_argument("--groups", type=str, choices=["dyjets","singletop","ttV","VVV","VV","data"], nargs='*',
                                help="Group of sample to process [dyjets, singletop, ttV, VVV, VV, data]" )

now = datetime.now()
time = now.strftime("%y%m%d")
tag = "skimming_{0}".format(time)


args = parser.parse_args()
ymls = args.ymls
groups = args.groups


if os.path.isfile(args.sourcefile):
    source_file = os.path.abspath(args.sourcefile)
else:
    print("Executeble file {source} does not exist".format(source=args.sourcefile))
    exit(11)

def createFolder(path):
    if os.path.isdir(path) == False:
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

        categories = yaml.load(f)

        # print("{folder}/skimming_submit.sh".format(folder=inputs_folder+yml_name))
        # print("input script: {folder}".format(folder=inputs_folder+yml_name))

        #sub_script = open("{folder}/submit_all.sh".format(folder=process_folder_name),"w")
        #sub_script.write("#!/bin/bash\n")

        for category in categories:
            isBkg  = True if category=='backgrounds' else False
            isData = True if category=='data' else False
            isSig  = True if category=='signals' else False
            isMC   = "false" if isData else "true"
            if isSig and args.skipSignals:
                continue
            if not isBkg and args.onlyBackground:
                continue

            processes = categories[category]
            for sample in processes.keys():

                if len(groups) > 0 :
                    if ("group" in processes[sample]):
                        if not groups.count(processes[sample]["group"]):
                            continue
                    else:
                        continue

                files_path = processes[sample]['path']
                print("Processig: {sample}".format(sample=sample))

                list_of_files = glob.glob("{}/*.root".format(files_path))

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
                    filename = file.split('/')[-1]
                    idx = re.findall(r"\d+",file)[-1]
                    shell_filename = "bash_script_{idx}.sh".format(idx=str(idx))
                    shell_script   = open("{folder}/{sh}".format(folder=sample_scripts_folder,sh=shell_filename), "w")
                    command = "{source} {path} {file} {mcflag}".format(source=source_file, path=files_path, file=filename, mcflag=isMC)
                    #command = "cmsRun {file_cfg} isMC=True isMCSignal=True hasLHE=True inputFile={path} outputFile={outFile} newIVF={flag}".format(file_cfg=cfgfile, path=file, outFile="output_"+str(idx)+".root", flag=("True" if newIVF else "False"))
                    #print(command)
                    script = bash_script_tmp.format(cmssw_src=cmssw_src_folder, command=command, copy_cmd=copy_cmd, outFile=filename, output_path=copy_destination)
                    shell_script.write(script)
                    shell_script.close()
                    os.system("chmod 777 {folder}/{sh}".format(folder=sample_scripts_folder,sh=shell_filename))
                    qsub_file.write("qsub -o {log} -e {log} {sh} \n".format(log=sample_logs_folder,sh=sample_scripts_folder+"/"+shell_filename))
                qsub_file.close()
                os.system("chmod 777 {folder}/Submit_{sample}.sh".format(folder=sample_scripts_folder, sample=sample))
                print("Please, "+bcolors.HEADER+"run"+bcolors.ENDC+":\n \t"+bcolors.OKGREEN+"bash "+sample_scripts_folder+"/Submit_"+sample+".sh "+bcolors.ENDC)

                #sub_script.write("bash {folder}/Submit_{sample}.sh\n".format(folder=sample_scripts_folder, sample=sample))
            #sub_script.close()
