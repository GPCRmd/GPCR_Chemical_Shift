import mdtraj as md
import os
import multiprocessing as mp
from json import loads, dump
import argparse as ap
import pandas as pd
import numpy as np
import math
import sys
import traceback
import resource
import time

# Incorporate Sparta and shiftx2 to the path, so mdtraj can work with them
filespath = ''
outpath = ''
os.environ['PATH'] = os.environ['PATH']+':.../software/shiftx2-linux/'
os.environ['PATH'] = os.environ['PATH']+':.../software/SPARTA+/'
os.environ["SPARTAP_DIR"] = ".../software/SPARTA+/"


### Functions
def json_dict(path):
	"""Converts json file to pyhton dict."""
	json_file=open(path)
	json_str = json_file.read()
	json_data = loads(json_str)
	return json_data

def cs_for_traj(dyn_id, overwrite, trajfiles):
	"""
	Calculate chemical shift for this trajectory file of this dynamic id
	"""

	# Residue one-letter to three-letter
	three_one = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
	 'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
	 'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
	 'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

	block_sizes=list(range(5,205,5))

	def parse_pdb_chains(pdb_file):
		"""
		Parses a PDB file and returns a dictionary with unique chain IDs as keys
		and their order of appearance as values.
		"""
		chain_dict = {}
		chain_order = 0

		with open(pdb_file, 'r') as file:
			for line in file:
				# Look for lines containing ATOM or HETATM entries
				if line.startswith(("ATOM", "HETATM")):
					chain_id = line[21]  # Chain ID is in column 22 (index 21 in 0-based indexing)
					if chain_id not in chain_dict:
						chain_dict[chain_id] = str(chain_order)
						chain_order += 1

		return chain_dict

	def prepare_traj(traj, compl_data):
		"""
		Filter out non-protein things, ligands and cappings from the system.
		We won't calculate the CS of any of these 
		"""

		# Identify the order of appearance of chain IDs in original pdb file (these numbers will be the new chainIDs in mdtraj objects)
		chain_dict = parse_pdb_chains(filespath+compl_data['dyn'+dyn_id]['struc_f'])

		# Take peptide ligand if exists. We will not include it in the calculation of chemical shift things
		peplig_sel = ''
		if compl_data['dyn'+dyn_id]['peplig']:
			peplig_chain_num = chain_dict[compl_data['dyn'+dyn_id]['peplig']]
			peplig_sel =  " and not chainid "+str(peplig_chain_num)

		# Filter out non-gpcr things
		print('Filtering trajectory of non-protein elements...')
		gpcr_chain = compl_data['dyn'+dyn_id]['gpcr_chain']
		gpcr_chain_num = chain_dict[gpcr_chain]
		atomprot = list(traj.top.select('protein'))
		trajprot = traj.atom_slice(atomprot)
		
		# Take capping residues atoms (first and last of every chain/segment)
		caps_atoms = ""
		prev_res = 0
		for res in trajprot.top.residues:
			stres = str(res.resSeq)
			if prev_res == 0:
				caps_atoms+=stres+" "
			elif res.resSeq != prev_res+1:
				caps_atoms+=str(prev_res)+" "
				caps_atoms+=stres+" "
			prev_res = res.resSeq
		caps_atoms+=stres

		# Remove cap atoms from trajectory object
		atomprot_nocap = list(trajprot.top.select('(not resSeq %s) and protein %s'%(caps_atoms, peplig_sel)))
		trajprot_nocap = trajprot.atom_slice(atomprot_nocap)

		return(trajprot_nocap)

	def zero_slope(data, chunksize = 7,max_slope = .0005):
		"""return the 'first' data point with zero slope
		data --> numpy ndarray - 2d [[x0,y0],[x1,y1],...]
		chunksize --> odd int
		returns numpy ndarray
		"""
		midindex = chunksize / 2
		for index in range(len(data) - chunksize):
			chunk = data[index : index + chunksize]
			# subtract the endpoints of the chunk
			# if not sufficient, maybe use a linear fit
			dy = abs(chunk[0] - chunk[-1])
			dx = 2
			if 0 <= dy / dx < max_slope:
				#print(dy/dx)
				return(chunk[int(midindex)])

	def ext_err(job,logfile):
		ext_err_list=list()
		output_log = open(logfile, "w")
		output_log.write("THESE RESIDUES FAILED WHEN CALCULATING THE CS BECAUSE THE MEAN IS THE SAME OVER THE MD\n")
		with open(job) as csv_f:
			lines = csv_f.readlines()[1:]
			for line in lines:
				line_s = line.strip()
				line_x = line_s.split(";")
				line_x = [x for x in line_x if x]# Remove empty strings
				cs_values=line_x[4:]
				est_err=list()
				overall_mean=(round(sum(list(map(float,cs_values)))/len(cs_values),4))
				for block_size in block_sizes:
					if block_size > len(cs_values)/2:#Avoid frame-block sizes bigger or equal than the entire simulation length 
						continue
					block_groups = [cs_values[n:n+block_size] for n in range(0, len(cs_values), block_size)]
					if len(block_groups[-1]) != block_size:
						block_groups.pop()
					block_average=list()
					for block in block_groups:
						block_average.append(round(sum(list(map(float, block)))/block_size,4))
					# Compute the square of all the average values in the colvar
					square_averages = list(map(lambda n: n ** 2, block_average))
					# Now compute the average over all the averages
					mean = sum(block_average)/len(block_average)
					# Compute the average of the squares of the individual averages
					mean2 = sum(square_averages)/len(square_averages)
					# Compute the population variance amongst the block averages
					population_variance = mean2 - mean*mean
					# Convert the population variance into a sample variance by multiplying by the bessel factor
					sample_variance = ( len(square_averages) / ( len(square_averages) - 1 ) )*population_variance
					# Print out the length of the blocks, the final average taken over all blocks and the square 
					# root of the sample variance divided by the number of data points that this estimate was
					# calcualted from.  This last term is a measure of the eror bar
					#print(line_x[:4])
					#print("block size=" + str(block_size) + "\nmean="+str(round(mean,4))+"\nsample variance="+str(abs(sample_variance))+"\nlen(square_averages)"+str(len(square_averages)))
					#print( block_size, round(mean,4), round(math.sqrt(abs(sample_variance)/len(square_averages)),4))
					est_err.append(round(math.sqrt(abs(sample_variance)/len(square_averages)),4))
					
					#print(npest_err)
					#print(line_x[:4])
				if any(est_err):
					if zero_slope(est_err) != None:
						ext_err_list.append([str(line_x[0]),str(line_x[1]),str(line_x[2]),str(line_x[3]),str(overall_mean),str(zero_slope(est_err))])
					else:
						ext_err_list.append([str(line_x[0]),str(line_x[1]),str(line_x[2]),str(line_x[3]),str(overall_mean),str(est_err[-1])])
				else:
					for failing_res in line_x[:4]:
						output_log.write(str(failing_res)+"\t")
					output_log.write('\n')
			output_log.close()
		return ext_err_list

	def sum_errors(ext_err):
		IntErr = {"C": 0.5330,"CA": 0.4412,"CB": 0.5163,"N": 1.1169,"HA": 0.1231,"H": 0.1711,"C_sc":0.9787,"H_sc":0.9482,}
		for atom in ext_err:
			if atom[1] in ["C","CA","CB","N","HA","H"]:
				atom.append(IntErr[atom[1]])
				atom.append(round(np.sqrt(IntErr[atom[1]]**2 + float(atom[-2])**2),4))
			else:
				if atom[1][0] == "H":
					atom.append(IntErr["H_sc"])
					atom.append(round(np.sqrt(IntErr["H_sc"]**2 + float(atom[-2])**2),4))
				if atom[1][0] == "C":
					atom.append(IntErr["C_sc"])
					atom.append(round(np.sqrt(IntErr["C_sc"]**2 + float(atom[-2])**2),4))
		return ext_err

	def format_output(total_err,outfile):
		output = open(outfile+".txt", "w") 
		output_users = open(outfile+"_usr.txt", "w")
		output.write("RES_ID\tRES_NAME\tATOM\tCS\tERROR\n")
		output_users.write("#Text file containing the mean chemical shift value for each atom predicted\n#Int_Error = Intrisical Error associated with the CS prediction made by SHIFTX+\n#Ext_Error = Extrinsical error associated to the fact that the MD is not infinite. The error is estimated by performing a block analysis for each atom.\n#Total_Error = Sum of squares of the Int_Error and the Ext_Error\n")
		output_users.write("ResidueID ResidueNAME Atom_Type Mean_CS(ppm) Int_Error(ppm) Ext_Error(ppm) Total_Error(ppm)\n")
		for item in total_err:
			output.write(str(item[0])+"\t"+str(item[2])+"\t"+str(item[1])+"\t"+str(item[4])+"\t"+str(item[7])+"\n")
			output_users.write(str(item[0])+" "+str(item[2])+" "+str(item[1])+" "+str(item[4])+" "+str(item[6])+" "+str(item[5])+" "+str(item[7])+"\n")
		output.close()
		output_users.close()

	def get_resname(t, row):
		"""
		Get residue names from atom names and residue ids in the sequence
		"""
		resid = row['resSeq']
		aindex = t.top.select('resSeq %d'%(resid))[0]
		resname = t.top.atom(aindex).residue.name
		return(resname)

	def run_CS(traj,outname,logfile,overwrite=False,cstype='shiftx2'):
		"""
		Run mdtraj CS analysis on trajectory
		"""
		outcsv = outname+'.csv'
		if os.path.exists(outcsv):
			print(cstype+' already computed for this system. Skipping...')
			return
		print('Running %s'%(cstype), outname)
		if (not os.path.exists(outcsv)) or overwrite:
			if cstype =='shiftx2':
				df = md.chemical_shifts_shiftx2(traj,pH=7.0,temperature=310)
			else:
				df = md.chemical_shifts_spartaplus(traj)
			df.dropna(inplace=True)
			# Add residue name column
			df_rei = df.reset_index() 
			df_rei.insert(2,'resname',df_rei.apply(lambda x: get_resname(traj, x), axis=1))
			df_rei.insert(3,'resname_s',df_rei.apply(lambda x: three_one[x['resname']], axis=1))
			df_rei.to_csv(outcsv, sep=';',index=False)

			# Calculate independent error of shift
			extrinsic_err=ext_err(outcsv,logfile)
			total_err=sum_errors(extrinsic_err)
			# Format independent error of shift
			format_output(total_err,outname)

	try: 

		# Skip if structure has a g protein (they are way too large, and crash Ori)
		if compl_data['dyn'+dyn_id]['gprot_chain_a']:	
			print('Dynamic %s has a G protein. SKipping...'%dyn_id)
			return

		# Skip if already computed
		for trajfile in trajfiles:
			filenum = trajfile.split('/')[-1].split('_')[0]
			shiftname = outfolder+"cs_shifty_dyn%s_%s"%(dyn_id, filenum)
			spartaname = outfolder+"cs_sparta_dyn%s_%s"%(dyn_id, filenum)
			shiftlog = logfolder+"cs_shifty_dyn%s_%s.log"%(dyn_id, filenum)
			spartalog = logfolder+"cs_sparta_dyn%s_%s.log"%(dyn_id, filenum)
			if os.path.exists(shiftname+'.csv') and os.path.exists(spartaname+'.csv') and not overwrite:
			# if os.path.exists(shiftname+'.csv'):
				print("Trajectory %s of dynid %s already has CS. Skipping..."%(trajfile, dyn_id))
				continue

			# Take number of file from trajectory file
			print('Computing CS for trajectory file '+trajfile)

			# Skip if no GPCR in system
			if not compl_data['dyn'+dyn_id]['gpcr_chain']:
				print("No GPCR found in system %s. Skipping"%dyn_id)
				continue

			# Append full path to trajectory and coordinates filenames
			trajfile = filespath+trajfile
			coordfile = filespath+compl_data['dyn'+dyn_id]['struc_f']

			# Try to do it using mdtraj shiftx+ implementation
			t = md.load(trajfile, top=coordfile)
			top = t.topology
			
			# Remove non-protein and other stuff from trajectory object
			# t = t[0]
			prep_traj = prepare_traj(t, compl_data)
			
			# Run shiftx2 and spartaplus on prepared trajectory (in paralel). Save its results in a CSV
			# pool = mp.Pool(processes=2)
			# for (name,log,cstype) in zip((shiftname,spartaname),(shiftlog,spartalog),('shiftx2',sparta)):
			# 	run_CS(prep_traj,name,log,overwrite,cstype=cstype)
			# 	x = pool.apply_async(run_CS,args=(prep_traj,name,log,overwrite,cstype))
			# 	# print(x.get()) # Print errors when activated, but also removes paralelization
			# pool.close()
			# pool.join() 

			run_CS(prep_traj,shiftname,shiftlog,overwrite,cstype='shiftx2')
			run_CS(prep_traj,spartaname,spartalog,overwrite,cstype='sparta')

	except Exception as e:
		print("%s chemical shift failed because of: "%dyn_id)
		print(traceback.format_exc())



############
## Main code
############

# Arguments
parser = ap.ArgumentParser(description="this calculates interaction frequencies for given simulation")
parser.add_argument(
	'--dynid',
	dest='dynid',
	action='store',
	nargs='+',
	default=False,
	help=' (int) dynamic id of the simulations to process, space-separated. If none specified, all dynamics avaliable will be processed'
)
parser.add_argument(
	'--overwrite',
	dest='overwrite',
	action='store_true',
	default=False,
	help='(bool) Repeat processing of this dynamics, even if results are already avaliable for it.'
)
parser.add_argument(
	'--cores',
	dest='cores',
	action='store',
	type=int,
	default=1,
	help='number of cores to use'
)
parser.add_argument(
	'--gpcrmd_only',
	dest='gpcrmd_only',
	action='store_true',
	default=False,
	help='Process GPCRmd community simulations only'
)
# parser.add_argument(
# 	'--memory_limit',
# 	dest='mem',
# 	action='store',
# 	default=20000,
# 	help='Max memory a given CS process is allowed to occupy (CS software is memory hungyr, and can easilty collaspe Ori)'
# )

# Arguments
args = parser.parse_args()
dyn_ids = args.dynid
overwrite = args.overwrite
cores = args.cores
gpcrmd_only = args.gpcrmd_only
# MEMORY_LIMIT_MB = int(args.mem)
# CHECK_INTERVAL = 10

if __name__ == '__main__':

	# Take original trajectory fileid from compl_data
	compl_data = json_dict(filespath+"Precomputed/compl_info.json")

	# Directory where formated CS files are to be stores
	outfolder = outpath+".../"
	os.makedirs(outfolder,exist_ok=True)
	logfolder = outpath+".../logs/"
	os.makedirs(logfolder,exist_ok=True)

	# Extract dynamics ID to process
	if dyn_ids:
		dyns = dyn_ids
	else:
		dyns = [ str(compl_data[a]['dyn_id']) for a in compl_data]

	#Calculate shifts
	pool = mp.Pool(processes=cores)
	# monitor = mp.Process(target=monitor_workers, args=(pool,))
	# monitor.start()
	for dyn_id in dyns:

		# Make directories for provisional result files
		trajfiles = compl_data['dyn'+dyn_id]['traj_f']
		# cs_for_traj(dyn_id, overwrite, trajfiles)
		x = pool.apply_async(cs_for_traj,args=(dyn_id, overwrite, trajfiles))
		# print(x.get()) # Print errors when activated, but also removes paralelization


	pool.close()
	pool.join() 
