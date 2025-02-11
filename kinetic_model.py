import numpy as np 


#DNA to RNA
def rev_comp_dna_to_rna(st):
    nn = {'A': 'U', 'C': 'G', 'G': 'C', 'T':'A'}
    return "".join(nn[n] for n in reversed(st))

#RNA to DNA
def rev_comp_rna_to_dna(st):
    nn = {'U': 'A', 'A': 'T', 'G': 'C', 'C':'G'}
    return "".join(nn[n] for n in reversed(st))


#PARAMETERS FOR DNA>DNA TMSD
#(coaxial stacking from Peyret N. 2000. Prediction of nucleic acid hybridization: parameters and algo- rithms.)

#in kcal/mol at 1M salt
#coaxial stacking interfaces (i.e. across the gap) in 3' to 5'
#this is from table 8.4 of the thesis; some parameters show up multiple times with different flanking sequence, in which case they are averaged

delta_G_pars_DNA_cxstck = {

	'AA':-2.55,
	'AC':-3.0,
	'AG':-2.4,
	'AT':-2.2,
	'CA':-1.7,
	'CC':-1.6,
	'CG':-0.8,
	'CT':-1.15,
	'GA':-2.9,
	'GC':-3.3,
	'GG':-2.6,
	'GT':-1.9,
	'TA':-1.6,
	'TC':-2.6,
	'TG':-1.7,
	'TT':-1.5,
}






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

	#applying 1M salt correction (NOTE: THIS SHOULD BE APPLIED TO A FULL DUPLEX )
	#total = (total + 1.667)/0.63

	return total


#DNA
def get_delta_G_DNA(sequence, average=False, initiation_penalty=True):
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

	#AVERAGE CASE; overriding previous total
	if average == True:
		pair_avg = np.mean([delta_G_pars_DNA[i] for i in delta_G_pars_DNA][0:16])
		init_avg = np.mean([delta_G_pars_DNA[i] for i in delta_G_pars_DNA][16])
		total = init_avg+(len(sequence)-1)*pair_avg

	return total




#OLD VERSION WHERE ONLY FORWARD STEP IS ALTERED

#NN estimation of free energy profile (enter RNA invader BMD sequence, 3' to 5')
def get_free_energy_profile_RDD(rna_seq):
	dG = np.zeros(len(rna_seq))

	#appending final base in toehold
	rna_seq = rna_seq
	dna_seq = rev_comp_rna_to_dna(rna_seq)[::-1]
	
	#sequence dependent changes
	for i in range(len(rna_seq)-1):

		#NEW KINETIC MODEL (RDD)
		dG_hybrid = get_delta_('G',rna_seq[i:i+2], initiation_penalty=False)
		#dG_hybrid = (dG_hybrid + (1.667/len(rna_seq)))/0.63 #converting to 1M salt
		dG_dna = get_delta_G_DNA(dna_seq[i:i+2], initiation_penalty=False)
		dG_dna = (0.63*dG_dna) - 1.667/len(rna_seq) #converting to 0.1M salt
		dG[i] = (dG_hybrid - dG_dna)


	#initiation penalty, so that net delta G is correct
	if (rna_seq[0] == 'C') or (rna_seq[0] == 'G'):
		dG_hybrid = 2
	else:
		dG_hybrid = 2.6
	#dG_hybrid = (dG_hybrid + (1.667/len(rna_seq)))/0.63 #converting to 1M salt
	dG_dna = 1.96
	dG_dna = (0.63*dG_dna) - 1.667/len(rna_seq) #converting to 0.1M salt

	dG[-1] = dG_hybrid - dG_dna   
		

	return 1 * dG/0.5922 #converting to kT





def get_free_energy_profile_DRD(rna_seq):
	dG = np.zeros(len(rna_seq))

	#appending final base in toehold
	rna_seq = rna_seq
	dna_seq = rev_comp_rna_to_dna(rna_seq)[::-1]
	
	#sequence dependent changes
	for i in range(len(rna_seq)-1):

		#NEW KINETIC MODEL (RDD)
		dG_hybrid = get_delta_('G',rna_seq[i:i+2], initiation_penalty=False)
		#dG_hybrid = (dG_hybrid + (1.667/len(rna_seq)))/0.63 #converting to 1M salt
		dG_dna = get_delta_G_DNA(dna_seq[i:i+2], initiation_penalty=False)
		dG_dna = (0.63*dG_dna) - 1.667/len(rna_seq) #converting to 0.1M salt
		dG[i] = (dG_hybrid - dG_dna) 



	#initiation penalty, so that net delta G is correct
	dG_hybrid = 2.3
	#dG_hybrid = (dG_hybrid + (1.667/len(rna_seq)))/0.63 #converting to 1M salt
	dG_dna = 1.96
	dG_dna = (0.63*dG_dna) - 1.667/len(rna_seq) #converting to 0.1M salt

	dG[-1] = dG_hybrid - dG_dna
	
		

	return 1 * -dG/0.5922 #converting to kT






#NN estimation of DNA>DNA free energy profiles based on stacking and coaxial stacking
def get_free_energy_profile_DDD(rna_seq, toehold_seq = 'G'):
	dG = np.zeros(len(rna_seq))

	dG_DD_1 = np.zeros(len(rna_seq))
	dG_DD_2 = np.zeros(len(rna_seq))

	#substrate strand, 5' to 3'
	sub_seq = rev_comp_rna_to_dna(rna_seq)[::-1]
	#adding a toehold base to compute stacking energies
	sub_seq = toehold_seq+sub_seq
	

	
	#converting invader seq to DNA (still 3' to 5')
	inv_seq = rna_seq.replace('U', 'T')
	
	#adding a toehold base to compute coaxial stacking energies
	inv_seq = rev_comp_dna_to_rna(toehold_seq).replace('U', 'T')+inv_seq
	

	#sequence dependent changes
	for i in range(len(rna_seq)-1):
		#SantaLucia NN pars
		dG[i] = -(get_delta_G_DNA(sub_seq[i+1:i+3]) - get_delta_G_DNA(sub_seq[i:i+2]))*0.63
		
		#coaxial stacking
		dG_cxstck = (delta_G_pars_DNA_cxstck[inv_seq[i+1:i+3]] - delta_G_pars_DNA_cxstck[inv_seq[i:i+2]])
		dG_cxstck = (0.63*dG_cxstck) #- 1.667/len(rna_seq)
		dG[i] += dG_cxstck 
		
   
	return 1 * dG/0.5922  #converting to kT


def get_dG_BM(rna_seq, average=False):
	dna_seq = rev_comp_rna_to_dna(rna_seq)[::-1]
	dG_hybrid = get_delta_('G',rna_seq, average=average)
	dG_dna = (0.63*get_delta_G_DNA(dna_seq, average=average)) - 1.667 #converting to 0.1M


	return (dG_hybrid - dG_dna)/0.5922



#KINETIC MODEL PARAMETERS
#free energy changes (units of kT)
dG_assoc = 2.5
dG_bp = 2.52

dG_bm_HYBRID = 9.48
dG_bm_DNA = 10.3

dG_p = 3.5


#concentrations for a/di-ssociation
c = 10**(-6)
c_0 = 1

#rate constant defining absolute time
k_bp = 64000000 



cas9_correction = np.array([6.1740569929793985, 2.292884965022965, -4.090985570212429, 2.2085630539898684, 0.17779482971361027, 0.7274571629898683, 3.167292265894292, -1.6118749919263764, -3.084192712010132, -0.12606526928638973, -1.0470872747771025, -0.6421837717771025, 3.104968555254306, -1.8225620009263763, -0.22861278411617703, 1.994250864979398, -0.11453906197703467, -2.6680194782124285, 1.0705685609898683, -8.729406661250254])


#RATES
#lines from old 4-par setup commented out
def get_k_forward(reaction_type, bm_sequence, dG_rd, toehold_length, second_toehold_length):
	#checking reaction type
	if (reaction_type == 'RDD') or (reaction_type == 'DRD'):
		dG_bm = dG_bm_HYBRID	
	elif reaction_type == 'DDD':
		dG_bm = dG_bm_DNA	


	k_forward = []

	#association/dissociation
	#k_forward.append(k_bp * np.exp( -(dG_assoc - np.log(c/c_0)) ))
	k_forward.append(k_bp * np.exp( -(dG_assoc - np.log(c/c_0)) ))

	#toehold binding/unbinding
	for i in range(toehold_length-1):
		k_forward.append(k_bp)

	#initiating branch migration
	k_forward.append(k_bp * np.exp( -(dG_bm + dG_p + dG_rd[0]/2) )) 

	#continuing branch migration
	for i in range(1,len(bm_sequence)-1):
		k_forward.append(k_bp * np.exp( -(dG_bm + dG_rd[i]/2) ))

	
	#dissociation/second toehold unbinding followed by dissociation
	if second_toehold_length == 0:
		#dissociation
		k_forward.append(k_bp * np.exp( -(-dG_assoc + np.log(c/c_0) + dG_rd[-1]/2) ))
	else:
		#final branch migration step
		k_forward.append(k_bp * np.exp( -(dG_bm - dG_p + dG_rd[-1]/2) )) 

		#breaking base pairs in second toehold
		for i in range(second_toehold_length-1):
			k_forward.append(k_bp * np.exp(-dG_bp))

		#dissociation
		k_forward.append(k_bp * np.exp( -(-dG_assoc + np.log(c/c_0) - dG_bp) ))
		


	return k_forward



def get_k_back(reaction_type, bm_sequence,dG_rd,toehold_length, second_toehold_length):
	#checking reaction type
	if (reaction_type == 'RDD') or (reaction_type == 'DRD'):
		dG_bm = dG_bm_HYBRID	
	elif reaction_type == 'DDD':
		dG_bm = dG_bm_DNA	



	k_back = []

	#association/dissociation
	#k_back.append(k_bp * np.exp(-dG_bp))
	k_back.append(k_bp * np.exp(-dG_bp))

	#toehold binding/unbinding
	for i in range(toehold_length-1):
		k_back.append(k_bp * np.exp(-dG_bp))


	#initiating branch migration
	k_back.append(k_bp * np.exp( -(dG_bm-dG_rd[0]/2) )) #double check that this is correct

	#continuing branch migration
	#for i in bm_sequence[1:-1]:
	for i in range(1,len(bm_sequence)-1):
		k_back.append(k_bp * np.exp( -(dG_bm-dG_rd[i]/2) ))

	if second_toehold_length != 0:
		#breaking base pairs in second toehold
		for i in range(second_toehold_length):
			k_back.append(k_bp)


	#dissociation (irreversible in this case)
	#k_back.append(0)

	return k_back


#check that everything here is working
def get_k_off(reaction_type, bm_sequence, toehold_length, second_toehold_length):

		
	#off rates during branch migration
	N = len(bm_sequence) + second_toehold_length
	k_off_DD = []
	k_off_RD = []
	for n in range(len(bm_sequence)):

	
		dummy_seq = 1000*'A'
		#SEQ-DEP DISSOCIATION FOR RNA INVADING dsDNA or DNA invading dsDNA
		remaining_bases = bm_sequence[n:].replace('U', 'T')
		dG_NN_avg = 0.63*(get_delta_G_DNA(dummy_seq, average=True)/(0.5922*len(dummy_seq)))
		
		
		dG_NN_seqdep = 0.63*(get_delta_G_DNA(remaining_bases)/(0.5922*len(remaining_bases)))
		#right here we are adding (as opposed to subtracting) dG_bp + dG_NN_avg because dG_bp is an absolute free energy
		adjusted_dG_bp = (dG_bp + dG_NN_avg) - dG_NN_seqdep

		
		#adjusted_dG_bp = abs(get_delta_G_DNA(remaining_bases, average=True))
		k_off_DD.append(k_bp * np.exp( -(N-n)*(adjusted_dG_bp) ))
	

		
		#SEQ-DEP DISSOCIATION FOR DNA invading hybrid
		remaining_bases = bm_sequence[n:]
		dG_NN_avg = (get_delta_('G', dummy_seq, average=True)/(0.5922*len(dummy_seq)))
	
		dG_NN_seqdep = (get_delta_('G', remaining_bases)/(0.5922*len(remaining_bases))) 
		#right here we are adding dG_bp + dG_NN_avg because dG_bp is an absolute free energy
		adjusted_dG_bp = (dG_bp + dG_NN_avg) - dG_NN_seqdep

		#adjusted_dG_bp = abs(get_delta_G_DNA(remaining_bases, average=True))
		k_off_RD.append(k_bp * np.exp( -(N-n)*(adjusted_dG_bp) ))
		
		
		

		
	#checking reaction type
	if (reaction_type == 'RDD') or (reaction_type == 'DDD'):
		k_off = k_off_DD	

	elif reaction_type == 'DRD':
		k_off = k_off_RD								

	#off rates during toehold binding
	k_off_toehold = [k_off[0] for i in range(toehold_length-1)]	

	#off rates during unbinding from second toehold (i.e., zero)
	k_off_second_toehold = list(np.zeros(second_toehold_length))


	#print(len(k_off_toehold+ k_off))
	return k_off_toehold + k_off + k_off_second_toehold






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







def get_rate(reaction_type, seq, toehold_length, second_toehold_length = 0, toehold_seq = 'G', mismatch_positions = False, cas9=False, cas9_dG = 0):
	
	#free energy profile depending on reaction type
	if reaction_type == 'DDD':
		dG_rd = get_free_energy_profile_DDD(seq, toehold_seq = toehold_seq)
	elif reaction_type == 'RDD':
		dG_rd = get_free_energy_profile_RDD(seq) 
	elif reaction_type == 'DRD':
		dG_rd = get_free_energy_profile_DRD(seq)


	#cas9 correction
	if cas9 == True:
			#cas9_correction = np.array([cas9_dG for i in range(20)])
			dG_rd = get_free_energy_profile_RDD(seq) + cas9_correction
		
	
	#optional mismatch
	if mismatch_positions != False:
		for i in mismatch_positions:
			dG_rd[i] = 5

	

	k_off = get_k_off(reaction_type, seq, toehold_length, second_toehold_length)

	k_forward = get_k_forward(reaction_type, seq, dG_rd, toehold_length, second_toehold_length)
	k_back = get_k_back(reaction_type, seq, dG_rd, toehold_length, second_toehold_length)

	
	T_pass = np.sum([get_p_term(k_forward, k_back, k_off, i) for i in range(1,len(k_forward)+1)]) /get_j_term(k_forward, k_back, k_off, len(k_forward))
	return 1/(T_pass*c)








#COMPUTING RATES
#Enter substrate strand displacement domain sequence
#The function 'get_rate()' takes reaction type ('DDD', 'RDD' or 'DRD', corresponding to DNA>DNA, RNA>DNA, DNA>RNA), sequence and toehold length as arguments
#Rates are given in units M^-1 s^-1
sequence = 'AAAATGTGTGTGTCCC'
sequence = rev_comp_dna_to_rna(sequence)[::-1]

print(get_rate('RDD', sequence, 4))

