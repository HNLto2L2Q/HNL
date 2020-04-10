#! /bin/bash
# The script generate a .yml file from a list a HNL displaced samples
# Example:
# bash creatYMLfile.sh  /pnfs/iihe/cms/store/user/tomc/heavyNeutrinoMiniAOD/Moriond17_aug2018_miniAODv3/displaced/ \
#                        HeavyNeutrino_lljj*_e_Dirac_massiveAndCKM_LO \
#                        Fall17_lljj_e

if [ $# -lt 3 ]
then
  echo "Missing arguments!!"
  echo "Usage: $0 [path] [regex] [output]"
  echo "[path]   : path containing the different signal folders"
  echo "[regex]  : regular expression used for matching the signal wanted"
  echo "[output] : output file name (without .yml)"
  exit
fi


signal_path=$1
regex=$2
output_file=$3


generateList (){

  full_path=$1$2
  list_of_dir=`ls -d $full_path`

  for sample in $list_of_dir
  do
    sample_name=`echo $sample | grep -E -o "Heavy.*" | tr . p`
    mass=`echo $sample | grep -E -o "M-[0-9]?[0-9]" | cut -d'-' -f2`
    coupling=`echo $sample | grep -E -o "V-[0-9].[0-9]*_" | cut -d'-' -f2 | cut -d'_' -f1`
    echo "$sample_name:"
    echo "  path: $sample"
    echo "  mass: $mass"
    echo "  coupling: $coupling"
  done
}


generateList $signal_path $regex &> $output_file.yml
