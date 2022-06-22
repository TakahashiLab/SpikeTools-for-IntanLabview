def append_to_path(dir0): # A convenience function
    if dir0 not in sys.path:
        sys.path.append(dir0)

import argparse
import os, sys, json


# imports from this repo
append_to_path('/nfs/home/tsusumu/github/mountainlab/packages/mountainsort_examples/python')

params_dir='/nfs/home/tsusumu/SpikeTools-for-IntanLabview/spikesorting/ms'

#######################################
# Initialize the pipeline object
#######################################
parser = argparse.ArgumentParser(
    description="Takahashilab mountainview ")

args = parser.parse_args()
work_path = 'ms'


#animal_day_output_path = args.output

# parse the epoch names from the input directory
epoch_names = [name for name in sorted(
    os.listdir(work_path)) if name.endswith('.mda')]

#dsdir=os.getcwd()+'/dataset'
output_base_dir=os.getcwd()+'/'+work_path

#######################################
# RUN THE PIPELINE
#######################################
for epoch in epoch_names:
    dsdir=os.getcwd()+'/'+work_path+'/'
    output_dir=output_base_dir+'/'+epoch[0:-4]
    print(dsdir)
    print(epoch)
    cmd='qt-mountainview '+'--raw='+output_dir+'/'+epoch+' --samplerate=25000 --filt='+output_dir+'/filt.mda.prv'+' --pre='+output_dir+'/pre.mda.prv'+' --firings='+output_dir+'/firings_uncurated.mda --cluster_metrics='+output_dir+'/cluster_metrics.json --geom='+params_dir+'/geom.csv'
    print(cmd)
    os.system(cmd)


