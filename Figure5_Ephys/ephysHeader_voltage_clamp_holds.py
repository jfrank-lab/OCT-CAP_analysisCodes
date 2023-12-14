# -*- coding: utf-8 -*-
"""
Created on Thu Dec 14 08:57:33 2023

@author: howeca
"""

import os
import numpy as np
import pyabf

######## User inputs ########
# insert csX_cellX folder here. You need to run this for each experimental repeat
#e.g. folder=r'Z:\Labs\Frank Lab\Carmel\Ephys\HEK\230404_trans_noIncubation_OCT-CAP_x4\cs4_cell1_trans_noIncubation_OCT-CAP_5uM_UV-irradiation_CAP_1uM'
# !!!!!!the r before the file path is very important!!!!!!!
folder=r'INSERT CSX_CELLX FOLDER HERE'


# function to get the UV pulse indices so where the pulse starts and ends
# also have the UV pulse value if you need it but I never used it
def stepInformation(UV_led):
    # where the UV pulse starts
    # get the index where the UV pulse goes above 0.1 V i.e. a change in step
    step_start = next(x for x, val in enumerate(UV_led)
                                  if val > 0.1) 

    # postive current step amplitude 
    # in V
    step_value = round(UV_led[step_start] - UV_led[0])

    # where the UV pulse ends
    # get the index where the UV pulse goes below 0.1 V i.e. a change in step
    step_end = (next(x for x, val in enumerate(UV_led[step_start:len(UV_led)])
                                  if val < UV_led[step_start]-0.1)) + step_start 
    return step_start, step_end


"""
RUN this for voltage holds. UV holds and drug holds need the same code so just run twice 
for the recording you are analysing
Change the holdFolder below
"""
holdFolder =folder + '\\holds_UV'   # find the folder with the UV hold recordings
#holdFolder =folder + '\\holds_drug'  # find the folder with the drug hold recordings

filenames=np.array(os.listdir(holdFolder)) # list the files within that folder

abf = pyabf.ABF(holdFolder + '\\' + filenames[0]) # import the UV hold file

# assign variables from the abf file import
time = abf.sweepX    
current = abf.sweepY
voltage = abf.sweepC
UV_led= abf.data[2]
sample_rate=abf.sampleRate # samples per second


steadyStateCurrent=np.mean(current[0:30000]) # get the baseline current 

step_start,step_end = stepInformation(UV_led)

# end current is the average over the last 10 seconds of the UV pulse and 20 seconds after
# this is what I used for my negative controls
endCurrent = np.mean(current[step_end-(sample_rate*10):step_end+(sample_rate*20)]) 

# for results CAP/untethered/tethered I changed this to be the indices where the max response was
# save these indices
endCurrent = np.mean(current[2250000:2400000]) 


diff = endCurrent - steadyStateCurrent

##### here i manually saved the following ####
#file in final analysis "holding_current_values_CONDITION.csv"
# folder
# filename
# steadyStateCurrent
# endCurrent
# diff
# condition indices


  