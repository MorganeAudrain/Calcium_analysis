#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Morgane

This pipeline is use to chose settings for start the analyse
"""

#%% Importation and parameters

from Steps.run_steps_param import run_steps
import psutil
import caiman as cm

print("This pipeline is the one for choose the parameters that you will use for the rest of the mouse")
print('Choose the mouse, session, trial you want to analyse for decide global settings')
mouse_number = int(input("mouse number : "))
sessions = input(" session : ")
trial = input(" trial : ")

print('Choose which steps you want to run: 0 -> decoding, 1 -> cropping, 2 -> motion correction, 3 -> alignment, 4 -> equalization, 5 -> source extraction, 6 -> component evaluation, 7 -> registration ')
n_steps = input(' steps :')

# start a cluster

n_processes = psutil.cpu_count()
c, dview, n_processes = cm.cluster.setup_cluster(backend='local', n_processes=n_processes, single_thread=False)


process ='yes'
# Run steps that you want
while process == 'yes':

    run_steps(n_steps, mouse_number, sessions, trial, dview)
    print('This step is finish. Do you want to run another step with those settings ? (yes or no)')
    process=input("answer : ")
    if process == 'yes':
        print('Choose which steps you want to run: 0 -> decoding, 1 -> cropping, 2 -> motion correction, 3 -> alignment, 4 -> equalization, 5 -> source extraction, 6 -> component evaluation, 7 -> registration, all ->  every steps ')
        n_steps = input(' steps :')
        run_steps(n_steps, mouse_number, sessions, trial, dview)
    if process == 'no':
        print('Thanks for having used this pipeline. I hope everything went alright for you :)')
        dview.terminate()
    else:
        print('you did not choose a viable choice retry ')
        process = input("answer : ")


