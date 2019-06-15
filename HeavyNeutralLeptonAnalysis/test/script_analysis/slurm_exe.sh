#!/bin/bash

files_list=$1
job_number=$2

set -e

./CloneTree_new.exe $files_list

hadd skimmed_Analysis_output_data_${job_number}.root skimmed_Analysis_output_prova_*.root
rm skimmed_Analysis_output_prova_*.root
rm ${files_list}
