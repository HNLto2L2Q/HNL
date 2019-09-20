# coding: utf-8
import yaml
import subprocess
import re

with open('input_mc_2018.yml','r') as f:
    doc = yaml.load(f)
    # output_file = open('ext_list.txt','w')
    samples = doc['samples'].keys()
    for sample in samples:
        dataset = doc['samples'][sample]['dsname']
        dataset_re = re.sub(r'v15-v[0-9]', 'v*', dataset)
        ext = subprocess.check_output(['dasgoclient', '-query', 'dataset='+dataset_re ]).decode()
        
        if ext.find('ext') > -1:
            nevents = subprocess.check_output(['dasgoclient', '-query', 'file dataset='+dataset+'| sum(file.nevents)' ]).decode()
            print('------ {} - {}'.format(dataset, nevents))
            print('Extentions')
            print(ext)
	    #for ds in ext:
            #nevents_ext = subprocess.check_output(['dasgoclient', '-query', 'file dataset='+ext+'| sum(file.nevents)' ]).decode()
            print('   ------ {}'.format(ext))
      #  if site.find('T2_BE_UCL') > -1:
      #      print(site)
      #  else:
      #      output_file.write(dataset)
      #      print('--'*20)

