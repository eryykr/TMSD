import os
import shutil
#import numpy as np

def make_input_file(directory,system_type,sim_ID):

   #making it slightly easier to write single lines
   def w(line):
      f.write(f'{line} \n') 


   f = open(f"input_{system_type}_{sim_ID}.in", "w")

   #Program Parameters
   w('backend = CPU')
   w(f'seed = 123{sim_ID}')
   w('interaction_type = DRH')
   w('use_average_seq = 0')
   w(f'seq_dep_file_DRH = {directory}/DRH_pars.txt')
   w(f'seq_dep_file_DNA = {directory}/DNA_pars.txt')
   w(f'seq_dep_file_RNA = {directory}/RNA_pars.txt')
   w('salt_concentration = 1')
   w('')

   #Simulation Options
   w('sim_type = VMMC')
   w('ensemble = NVT')
   w('steps = 1e9')
   w('check_energy_every = 100000')
   w('check_energy_threshold = 1.e-4')
   w('delta_translation = 0.1')
   w('delta_rotation = 0.2')
   w(f'T = 25C')
   w('verlet_skin = 1.00')
   w('umbrella_sampling = 1')
   w(f'op_file = {directory}/op.txt')
   w(f'weights_file = {directory}/wfile_{system_type}.txt') 
   w('maxclust = 12')
   w('small_system = 1')
   w('')
  
   #Input/Output
   w(f'last_hist_file = {directory}/Output/last_hist_{system_type}_{sim_ID}.dat')
   w(f'traj_hist_file = {directory}/Output/traj_hist_{system_type}_{sim_ID}.dat')
   w(f'topology = {directory}/{system_type}.top')
   w(f'conf_file = {directory}/tmsd_vmmc.dat')
   w(f'trajectory_file = {directory}/Output/trajectory_{system_type}_{sim_ID}.dat')
   w('no_stdout_energy = 0')
   w('restart_step_counter = 1')
   w(f'energy_file = {directory}/Output/energy_{system_type}_{sim_ID}.dat')
   w('print_conf_interval = 1e6')
   w('print_energy_every = 1e5')
   w('time_scale = linear')
   w('external_forces = 0')

   f.close()
  


system_types = ['fast', 'medium', 'slow', 'flat']

for i in range(10):
	for j in system_types:
		make_input_file('TMSD_kinetic_control',j,i)


#make_input_file('StrandDisplacement', 12)


