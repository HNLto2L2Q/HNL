import os
import yaml

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--onlyBackground", action='store_true', help='Submit skimming only for background samples')
parser.add_argument("--skipSignals", action='store_true', help='Submit skimming skipping signals samples')

args = parser.parse_args()

def createFolder(path):
    if os.path.isdir(path) == False:
        os.system('mkdir -p ' + path)



inputs_folder  = "FARM/inputs"
logs_folder     = "FARM/logs"
outputs_folder = "FARM/outputs"


bash_script_tmp = """#!/bin/bash
set -e

FOLDER=\"{folder}\"

source $VO_CMS_SW_DIR/cmsset_default.sh
pushd $FOLDER
eval `scramv1 runtime -sh`
popd

{command}

mv *.root $FOLDER/{output}
"""

with open('samples.yml','r') as f:

    pwd = os.getcwd()

    createFolder(inputs_folder)
    createFolder(logs_folder)
    createFolder(outputs_folder)

    config = yaml.load(f) #  , Loader=yaml.FullLoader)
    path_bkg = config['path_bkg']
    path_sig = config['path_sig']
    categories = config['samples']

    print("{folder}/skimming_submit.sh".format(folder=inputs_folder))

    qsub_file = open("{folder}/skimming_submit.sh".format(folder=inputs_folder),"w")
    qsub_file.write("#!/bin/bash\n")

    for category in categories:
        isBkg  = True if category=='backgrounds' else False
        isData = True if category=='data' else False
        isSig  = True if category=='signals' else False
        isMC   = "false" if isData else "true"
        if isSig and args.skipSignals:
            continue
        if not isBkg and args.onlyBackground:
            continue
        path  = path_bkg if (isBkg or isData) else path_sig
        processes = categories[category]
        for process in processes.keys():
            shell_filename = "bash_script_{folder}.sh".format(folder=process)
            shell_script   = open("{folder}/{sh}".format(folder=inputs_folder,sh=shell_filename), "w")
            filename = processes[process]['filename']
            command = "$FOLDER/CloneTree.exe {path} {file} {mcflag}".format(path=path, file=filename, mcflag=isMC)
	    print(command)
            script = bash_script_tmp.format(folder=pwd, command=command, output=outputs_folder)
            shell_script.write(script)
            shell_script.close()
	    os.system("chmod 777 {folder}/{sh}".format(folder=inputs_folder,sh=shell_filename))
            qsub_file.write("qsub -o {log} -e {log} {sh} \n".format(log=logs_folder,sh=inputs_folder+"/"+shell_filename))
    qsub_file.close()
