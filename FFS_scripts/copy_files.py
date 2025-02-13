import shutil



system_name = 'fast'

files_to_copy = ['input', f'{system_name}.top', f'tmsd_ffs.dat', 'DNA_pars.txt', 'RNA_pars.txt', 'DRH_pars.txt', 'op.txt']
new_names = ['input', 'unique.top', 'flux_initial.dat', 'DNA_pars.txt', 'RNA_pars.txt', 'DRH_pars.txt', 'op.txt']

for file, file_newname in zip(files_to_copy, new_names):
	
	shutil.copy(f'input_files/{file}', f'FLUX/{file_newname}')
	shutil.copy(f'input_files/{file}', f'IF1/{file_newname}')
	shutil.copy(f'input_files/{file}', f'IF2/{file_newname}')
	shutil.copy(f'input_files/{file}', f'IF3/{file_newname}')
	shutil.copy(f'input_files/{file}', f'IF4/{file_newname}')
	shutil.copy(f'input_files/{file}', f'IF5/{file_newname}')
	#shutil.copy(f'input_files/{file}', f'IF6/{file_newname}')
