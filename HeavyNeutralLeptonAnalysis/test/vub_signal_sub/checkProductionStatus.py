import glob
import re
from datetime import datetime

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--outputPath", type=str, required=True, help="output path" )

now = datetime.now()
time = now.strftime("%y%m%d")

args = parser.parse_args()
outputPath = args.outputPath

inputs_folder  = "FARM/inputs/"
logs_folder     = "FARM/logs/"

campaigns = glob.glob("{}/*".format(outputPath))



for campaign in campaigns:
    jobids   = glob.glob("{}/*".format(campaign))

    for jobid in jobids:
        samples  = glob.glob("{}/*/".format(jobid))
        input_jobid_path = jobid.replace(outputPath,inputs_folder)
        sub_script = open("{folder}/resubmit_failed.sh".format(folder=input_jobid_path),"w")
        sub_script.write("#!/bin/bash\n")

        for sample in samples:
            input_files_path = sample.replace(outputPath,inputs_folder)
            log_files_path   = sample.replace(outputPath,logs_folder)

            output_files = glob.glob("{}/*.root".format(sample))
            input_files  = glob.glob("{}/bash*.sh".format(input_files_path))
            n_outfiles   = len(output_files)
            n_infiles    = len(input_files)
            if n_outfiles == 0:
                print("Jobs not yet submitted")
                print(sample)
            elif  n_outfiles != n_infiles:
                outputs_id = set( [re.findall(r"\d+",ofile)[-1] for ofile in output_files] )
                inputs_id  = set( [re.findall(r"\d+",ifile)[-1] for ifile in input_files] )
                missing_ids = inputs_id - outputs_id
                for missing_id in missing_ids:
                    shell_filename = "bash_script_{id}.sh".format(id=missing_id)
                    sub_script.write("qsub -o {log} -e {log} {sh} \n".format(log=log_files_path,sh=input_files_path+shell_filename))

        sub_script.close()
                #print(n_outfiles,n_infiles)
                #print(campaign)
                #print(jobid)
                #print(sample)
                #print(inputs_id)
                #print(outputs_id)
                #print(missing_id)
