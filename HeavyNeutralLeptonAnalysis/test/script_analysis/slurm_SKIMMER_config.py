import os
import glob
import subprocess

from CP3SlurmUtils.Configuration import Configuration

config = Configuration()

num_files_per_job = 2

#--------------------------------------------------------------------------------
# 1. SLURM sbatch command options
#--------------------------------------------------------------------------------

config.sbatch_partition = 'cp3'
config.sbatch_qos = 'cp3'
config.sbatch_workdir = '.'
config.sbatch_time = '0-03:00'
config.sbatch_mem = '4096'
config.sbatch_output = '/dev/null'
config.sbatch_error = '/dev/null'
config.sbatch_additionalOptions = []

#--------------------------------------------------------------------------------
# 2. User batch script parameters that are same for all jobs
#--------------------------------------------------------------------------------

config.cmsswDir = os.environ.get('CMSSW_BASE')

config.inputSandboxDir = config.sbatch_workdir + '/slurm_input_sandboxes'

config.batchScriptsDir = config.sbatch_workdir + '/slurm_batch_scripts'

config.stageout = True
config.stageoutFiles = ['skimmed_*.root']
config.stageoutDir = '/nfs/user/ccaputo/HNL_Skim/2017/Data/job_array_${SLURM_ARRAY_JOB_ID}'
# We chose the filename of the outputs to be independent of the job array id number (but dependent on the job array task id number).
# So let's put the output files in a directory whose name contains the job array id number,
# so that each job array we may submit will write in a different directory.
#config.stageoutDir = '/nfs/user/ccaputo/HNL_Skim/Data_2017/job_array_${SLURM_ARRAY_JOB_ID}'
config.inputSandboxContent = ['CloneTree_new.exe','slurm_exe.sh']

config.writeLogsOnWN = True
config.separateStdoutStderrLogs = False
config.stageoutLogs = True
# The default filename of the slurm logs has already a job array id number and a job array task id number in it.
# So we can put all logs together (even from different job arrays we may submit) in a unique directory; they won't overwrite each other.
config.stageoutLogsDir = config.sbatch_workdir + '/slurm_logs'

config.useJobArray = True

# 2 jobs will be submitted, because the config parameter 'inputParams' has length 2.
config.numJobs = None

#--------------------------------------------------------------------------------
# 3 Job-specific input parameters and payload
#--------------------------------------------------------------------------------
config.inputParamsNames = ['inputFile','job_num']

# Get a list with all the input files.
inputFiles = glob.glob('/home/users/c/c/ccaputo/HNL/CMSSW_9_4_13_patch4/src/HNL/HeavyNeutralLeptonAnalysis/test/script_analysis/SingleMuon_2017*.txt')

#num_lines = sum(1 for line in open(inputFiles[0]))
#print num_lines
#tot_job = num_lines/num_files_per_job

# Now we will assume something that is almost always true: that the name of the input files
# contains an index and that we want to assign the same index to the output filename. 
# In this example the input filename has the index just before the '.root' file extension (my_input_file_<index>.root).
# We will loop over the list of input files and for each input file define the name of the output file,
# extracting from the input filename the desired index and adding it to the output filename.
config.inputParams = []
for inputFile in inputFiles:
    number = subprocess.check_output("cat "+inputFile+"| wc -l", shell = True)
    interval = 30
    number_of_jobs = (int(number.strip())/interval) + 1
    
    inputFileBasename = os.path.basename(inputFile) # here we removed the directory part and we are left with 'my_input_file_<index>.root'
    inputFileBasenameWithoutExt = os.path.splitext(inputFileBasename)[0] # here we removed the '.root' part and we are left with 'my_input_file_<index>'
    # print inputFileBasenameWithoutExt
    for x in range(1, int(number_of_jobs)+1):
        tmp_files_list = inputFileBasenameWithoutExt+"_"+str(x)+".txt"       
        os.system("sed '"+str(1+interval*(x-1))+","+str(interval*x)+"!d' "+inputFile+" > "+tmp_files_list)
        tmp_files_list_path = os.path.abspath(tmp_files_list)
        # print '\t '+tmp_files_list
        config.inputParams.append([tmp_files_list_path, '${SLURM_ARRAY_TASK_ID}'])

config.payload = \
"""
./slurm_exe.sh ${inputFile} ${job_num}
"""

