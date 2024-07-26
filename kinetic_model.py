
import numpy as np 







#DNA to RNA
def rev_comp_dna_to_rna(st):
    nn = {'A': 'U', 'C': 'G', 'G': 'C', 'T':'A'}
    return "".join(nn[n] for n in reversed(st))

#RNA to DNA
def rev_comp_rna_to_dna(st):
    nn = {'U': 'A', 'A': 'T', 'G': 'C', 'C':'G'}
    return "".join(nn[n] for n in reversed(st))

#SUGIMOTO MODEL (2021 ver.)
#model parameters (here we use the RNA sequence in the 3'->5' direction)
delta_H_pars = {
	'AA':-7800,
	'AC':-10100,
	'AG':-9400,
	'AU':-5800,
	'CA':-9800,
	'CC':-9500,
	'CG':-9000,
	'CU':-6100,
	'GA':-8600,
	'GC':-10600,
	'GG':-13300,
	'GU':-9300,
	'UA':-6600,
	'UC':-6500,
	'UG':-8900,
	'UU':-7400,
	#initiation parameters
	'C':0,
	'G':0,
	'A':0,
	'U':0
}

delta_S_pars = {
	'AA':-22.9,
	'AC':-27.7,
	'AG':-26.1,
	'AU':-17.4,
	'CA':-27.7,
	'CC':-25.1,
	'CG':-24.5,
	'CU':-18.4,
	'GA':-22.9,
	'GC':-27.7,
	'GG':-35.5,
	'GU':-25.5,
	'UA':-19.7,
	'UC':-16.4,
	'UG':-23.5,
	'UU':-24.5,
	#initiation parameters
	'C':-6.4,
	'G':-6.4,
	'A':-8.4,
	'U':-8.4
}

#in cal/mol
delta_G_pars = {
	'AA':-700,
	'AC':-1500,
	'AG':-1300,
	'AU':-400,
	'CA':-1200,
	'CC':-1700,
	'CG':-1400,
	'CU':-400,
	'GA':-1500,
	'GC':-2000,
	'GG':-2300,
	'GU':-1400,
	'UA':-500,
	'UC':-1400,
	'UG':-1600,
	'UU':200,
	#initiation parameters
	'C':2000,
	'G':2000,
	'A':2600,
	'U':2600
}
model_parameters = {
	'H':delta_H_pars,
	'S':delta_S_pars,
	'G':delta_G_pars
}

#SantaLucia & Hicks (5'->3' direction)
delta_G_pars_DNA = {
	'AA':-1.00,
	'AC':-1.44,
	'AG':-1.28,
	'AT':-0.88,
	'CA':-1.45,
	'CC':-1.84,
	'CG':-2.17,
	'CT':-1.28,
	'GA':-1.30,
	'GC':-2.24,
	'GG':-1.84,
	'GT':-1.44,
	'TA':-0.58,
	'TC':-1.30,
	'TG':-1.45,
	'TT':-1.00,

	#initiation parameter
	'init':1.96,

	#terminal AT correction
	'terminal_AT':0.05
}



#hybrids
def get_delta_(quantity, sequence, average=False, initiation_penalty=True):
	total = 0

	#initiation
	if initiation_penalty == True:
		if sequence[0] in 'AU' and sequence[-1] in 'AU':
			total = model_parameters[quantity]['A']
		else:
			total = model_parameters[quantity]['G']

	#pair contributions
	for i in range(len(sequence)-1):
		total += model_parameters[quantity][sequence[i:i+2]]

	#average case:
	if average == True:
		pair_avg = np.mean([model_parameters[quantity][i] for i in model_parameters[quantity]][0:16])
		init_avg = np.mean([model_parameters[quantity][i] for i in model_parameters[quantity]][16:20])
		total = init_avg+(len(sequence)-1)*pair_avg

	#converting to kcal/mol
	total = (total/1000)

	return total


#DNA
def get_delta_G_DNA(sequence, initiation_penalty=True):
	total = 0

	#initiation
	if initiation_penalty == True:
		total = delta_G_pars_DNA['init']

	#terminal AT penalty
	if 'AT' in sequence[0:2] or 'TA' in sequence[0:2] or 'AT' in sequence[-2:] or 'TA' in sequence[-2:]:
		total += delta_G_pars_DNA['terminal_AT']

	#pair contributions
	for i in range(len(sequence)-1):
		total += delta_G_pars_DNA[sequence[i:i+2]]

	return total






#NN estimation of free energy profile (enter RNA invader BMD sequence, 3' to 5')
def get_free_energy_profile(rna_seq):
	dG = np.zeros(len(rna_seq))
	dna_seq = rev_comp_rna_to_dna(rna_seq)[::-1]
	
	#sequence dependent changes
	for i in range(len(rna_seq)-1):
		dG_hybrid = get_delta_('G',rna_seq[i:i+2], initiation_penalty=False)
		dG_dna = get_delta_G_DNA(dna_seq[i:i+2], initiation_penalty=False)
		dG_dna = (0.63*dG_dna) - 1.667/len(rna_seq) #converting to 0.1M salt
		dG[i] = dG_hybrid - dG_dna

	#initiation penalty, so that net delta G is correct
	dG_hybrid = 2.3
	dG_dna = 1.96
	dG_dna = (0.63*dG_dna) - 1.667/len(rna_seq) #converting to 0.1M salt

	dG[-1] = dG_hybrid - dG_dna

	return 1 * dG/0.5922 #converting to kT





#KINETIC MODEL PARAMETERS
#free energy changes (units of kT)
dG_assoc = 2.5
dG_bp = 2.52
dG_bm = 9.3
dG_p = 3.5

#concentrations for a/di-ssociation
c = 10**(-6)
c_0 = 1

#rate constant defining absolute time
k_bp = 63200000


#RATES
def get_k_forward(bm_sequence, dG_rd, toehold_length):
	k_forward = []

	#association/dissociation
	k_forward.append(k_bp * np.exp( -(dG_assoc - np.log(c/c_0)) ))

	#toehold binding/unbinding
	for i in range(toehold_length-1):
		k_forward.append(k_bp)

	#initiating branch migration
	k_forward.append(k_bp * np.exp( -(dG_bm + dG_p + dG_rd[0]) )) 

	#continuing branch migration
	for i in range(1,len(bm_sequence)-1):
		k_forward.append(k_bp * np.exp( -(dG_bm + dG_rd[i]) ))

	#unbinding
	k_forward.append(k_bp * np.exp( -(-dG_assoc + np.log(c/c_0) - dG_bp + dG_rd[-1]) ))

	return k_forward


def get_k_back(bm_sequence,toehold_length):
	k_back = []

	#association/dissociation
	k_back.append(k_bp * np.exp(-dG_bp))

	#toehold binding/unbinding
	for i in range(toehold_length-1):
		k_back.append(k_bp * np.exp(-dG_bp))


	#initiating branch migration
	k_back.append(k_bp * np.exp( -(dG_bm) )) #double check that this is correct

	#continuing branch migration
	for i in bm_sequence[1:-1]:
		k_back.append(k_bp * np.exp( -(dG_bm) ))

	return k_back

def get_k_off(bm_sequence, toehold_length):

	#off rates during branch migration
	N = len(bm_sequence)
	k_off = []
	for n in range(len(bm_sequence)):
		k_off.append(k_bp * np.exp( -(N-n)*(dG_bp) ))
		
	#off rates during toehold binding
	k_off_toehold = [k_off[0] for i in range(toehold_length-1)]		
		
	return k_off_toehold + k_off




def get_j_term(k_forward, k_back, k_off, index):
	total_j = 1
	total_p = 1/k_forward[-1]
	
	if index > 1:
		for i in range(1,index):
			total_j = total_j + k_off[-i]*total_p 
			total_p = total_j/k_forward[-i-1] + total_p*(k_back[-i]/k_forward[-i-1])

	return total_j


def get_p_term(k_forward, k_back, k_off, index):
	total_j = 1
	total_p = 1/k_forward[-1]
	
	if index > 1:
		for i in range(1,index):
			total_j = total_j + k_off[-i]*total_p
			total_p = total_j/k_forward[-i-1] + total_p*(k_back[-i]/k_forward[-i-1])
			

	return total_p


def get_rate(seq, toehold_length):
	dG_rd = get_free_energy_profile(seq) 

	k_off = get_k_off(seq, toehold_length)
	k_forward = get_k_forward(seq, dG_rd, toehold_length)
	k_back = get_k_back(seq, toehold_length)

	T_pass = np.sum([get_p_term(k_forward, k_back, k_off, i) for i in range(1,len(k_forward)+1)]) /get_j_term(k_forward, k_back, k_off, len(k_forward))
	return 1/(T_pass*c)




#COMPUTING RATES
#Enter substrate strand displacement domain sequence
#The function 'get_rate()' takes sequence and toehold length as arguments
#Rates are given in units M^-1 s^-1
sequence = 'AAAATGTGTGTGTCCC'
sequence = rev_comp_dna_to_rna(sequence)[::-1]

print(get_rate(sequence, 4))

