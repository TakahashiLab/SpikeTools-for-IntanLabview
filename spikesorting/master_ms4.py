def append_to_path(dir0): # A convenience function
    if dir0 not in sys.path:
        sys.path.append(dir0)

import argparse
import os, sys, json
import numpy as np
from matplotlib import pyplot as plt

# mountainlab imports
from mountainlab_pytools import mlproc as mlp
from mountainlab_pytools import mdaio


# imports from this repo
append_to_path('/nfs/home/tsusumu/github/mountainlab/packages/mountainsort_examples/python')
from mountainsort4_1_0 import sort_dataset as ms4_sort_dataset # MountainSort spike sorting

params_dir='/nfs/home/tsusumu/SpikeTools-for-IntanLabview/spikesorting/ms'

#######################################
# Initialize the pipeline object
#######################################
parser = argparse.ArgumentParser(
    description="Takahashilab mountainsort4 sorting ")


args = parser.parse_args()
work_path = 'ms'
#animal_day_output_path = args.output

# parse the epoch names from the input directory
epoch_names = [name for name in sorted(
    os.listdir(work_path)) if name.endswith('.mda')]



#dsdir=os.getcwd()+'/dataset'
output_base_dir=os.getcwd()+'/'+work_path
if not os.path.exists(output_base_dir):
    os.mkdir(output_base_dir)

#######################################
# RUN THE PIPELINE
#######################################
for epoch in epoch_names:
    dsdir=os.getcwd()+'/'+work_path+'/'
    output_dir=output_base_dir+'/'+epoch[0:-4]
    if not os.path.exists(output_base_dir):
        os.mkdir(output_base_dir)

    print(dsdir)
    print(epoch)
    print(params_dir)
    ms4_sort_dataset(dataset_dir=dsdir,dataset=epoch,params_dir=params_dir,output_dir=output_dir,adjacency_radius=-1,detect_threshold=5)

