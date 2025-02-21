#import numpy as np

#function for writing the contents of a list to a file
def list_to_file(file_name, list_name):
    file = open(file_name, 'w')
    for i in list_name:
        file.write(f'{str(i)} \n')
    file.close()

#for loading the contents of a trajectory file into a list
def file_to_list(filename):
	#reading trajectory file line by line
	file = open(filename, 'r')
	file_lines = file.readlines()
	file.close()

	return file_lines




'''
op1 = range(11)
op2 = range(15)
weights = [1,1,1] + [10**i for i in range(10)]
weights  = np.array(weights)

weights = weights/populations
weights = weights/min(weights)
'''

op1 = range(12)
op2 = range(16)

#weights = [250,50,5,1.0, 0.4, 0.3, 0.2, 0.1, 0.2, 0.3, 0.2, 0.3, 3.0, 3.1, 9.5] #fast
#weights = [250,50,5,1.0, 1.2, 2.1, 2.0, 3.1, 14.3, 61.6, 24, 10, 44.9, 17.3, 9.5] #medium
weights = [250,50,5,1.0, 3.6, 4, 10.2, 16.3, 13.8, 19.1, 24.1, 16.4, 6.5, 2.4, 1] #slow
#weights = [250,50,5,1, 1, 1, 1, 1, 1, 1, 3, 4, 1, 1, 1] #flat

#INVERTING BRANCH MIGRATION WEIGHTS FOR DNA-INVADING-HYBRID
weights = weights[0:4] + [1/weights[i] for i in range(len(weights)) if i > 3]

wfile_lines = []
for i in op2:
	for j in op1:
		if i>0 and i<16 and j>0:
			wfile_lines.append(f'{i} {j} {weights[i-1]}')
		else:
			wfile_lines.append(f'{i} {j} {0}')



list_to_file('wfile_slow.txt', wfile_lines)










