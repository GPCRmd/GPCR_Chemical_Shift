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
import subprocess

import MDAnalysis as mda
from MDAnalysis.coordinates.PDB import PDBWriter
from MDAnalysis.topology.guessers import guess_types

from Bio.Data import IUPACData

# Incorporate Sparta and shiftx2 to the path, so mdtraj can work with them
filespath = './'
outpath = './example/'
example = './example/'
cspred_path = './CSpred/CSpred.py'  # Update with the actual path to CSpred.py

### Functions
# Function to convert 3-letter amino acid codes to 1-letter codes
def three_to_one(three_letter_code):
    return IUPACData.protein_letters_3to1.get(three_letter_code.capitalize(), 'X')

block_sizes=list(range(5,205,5))

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
			if cs_values == []:
				continue
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
	updated_ext_err = []
	for atom in ext_err:
		if atom[1] in ["C","CA","CB","N","HA","H"]:
			atom.append(IntErr[atom[1]])
			atom.append(round(np.sqrt(IntErr[atom[1]]**2 + float(atom[-2])**2),4))
			updated_ext_err.append(atom)
		else:
			if atom[1][0] == "H":
				atom.append(IntErr["H_sc"])
				atom.append(round(np.sqrt(IntErr["H_sc"]**2 + float(atom[-2])**2),4))
				updated_ext_err.append(atom)
			elif atom[1][0] == "C":
				atom.append(IntErr["C_sc"])
				atom.append(round(np.sqrt(IntErr["C_sc"]**2 + float(atom[-2])**2),4))
				updated_ext_err.append(atom)
	return updated_ext_err

def format_output(total_err,outfile):
	output = open(outfile+".txt", "w") 
	output_users = open(outfile+"_usr.txt", "w")
	output.write("RES_ID\tRES_NAME\tATOM\tCS\tERROR\n")
	output_users.write("#Text file containing the mean chemical shift value for each atom predicted\n#Int_Error = Intrisical Error associated with the CS prediction made by UCBShift\n#Ext_Error = Extrinsical error associated to the fact that the MD is not infinite. The error is estimated by performing a block analysis for each atom.\n#Total_Error = Sum of squares of the Int_Error and the Ext_Error\n")
	output_users.write("ResidueID ResidueNAME Atom_Type Mean_CS(ppm) Int_Error(ppm) Ext_Error(ppm) Total_Error(ppm)\n")
	for item in total_err:
		output.write(str(item[0])+"\t"+str(item[2])+"\t"+str(item[1])+"\t"+str(item[4])+"\t"+str(item[7])+"\n")
		output_users.write(str(item[0])+" "+str(item[2])+" "+str(item[1])+" "+str(item[4])+" "+str(item[6])+" "+str(item[5])+" "+str(item[7])+"\n")
	output.close()
	output_users.close()

def format_files(outcsv, outname, logfile):
	# Calculate independent error of shift
	extrinsic_err=ext_err(outcsv,logfile)
	total_err=sum_errors(extrinsic_err)
	# Format independent error of shift
	format_output(total_err,outname)

def get_atom_types():
	"""
	Create a dictionary where each amino acid (key) maps to a list of its unique atom types.

	Returns:
		dict: Dictionary with amino acids as keys and lists of unique atom types as values.
	"""
	partners = {
		"GLY": {"H": "N", "HA2": "CA", "HA3": "CA"}, 
		"ILE": {"H": "N", "HA": "CA", "HB": "CB", "HD1": "CD1", "HD2": "CD1", "HD3": "CD1", "HG12": "CG1", "HG13": "CG1", "HG21": "CG2", "HG22": "CG2", "HG23": "CG2"},
		"VAL": {"H": "N", "HA": "CA", "HG12": "CG1", "HG13": "CG1", "HG21": "CG2", "HG22": "CG2", "HB": "CB", "HG11": "CG1", "HG23": "CG2"}, 
		"ARG": {"H": "N", "HA": "CA", "HG2": "CG", "HB2": "CB", "HB3": "CB", "HD2": "CD", "HD3": "CD", "HE": "NE", "HG3": "CG","HH11": "NH1", "HH12": "NH1", "HH21": "NH2", "HH22": "NH2"},
		"TRP": {"H": "N", "HA": "CA", "HH2": "CH2", "HB2": "CB", "HB3": "CB", "HD1": "CD1", "HE1": "NE1", "HE3": "CE3", "HZ2": "CZ2",  "HZ3": "CZ3"},
		"GLU": {"H": "N", "HA": "CA", "HG2": "CG", "HG2": "CG", "HB2": "CB", "HB3": "CB", "HE2": "OE2", "HG3": "CG"},
		"GLN": {"H": "N", "HA": "CA", "HG2": "CG", "HB2": "CB", "HB3": "CB", "HE21": "NE2", "HE22": "NE2", "HG3": "CG"},
		"LYS": {"H": "N", "HA": "CA", "HG2": "CG", "HB2": "CB", "HB3": "CB", "HD2": "CD", "HD3": "CD", "HE2": "CE", "HE3": "CE", "HG3": "CG", "HZ1": "NZ", "HZ2": "NZ", "HZ3": "NZ"},
		"MET": {"H": "N", "HA": "CA", "HG2": "CG", "HB2": "CB", "HB3": "CB", "HE1": "CE", "HE2": "CE", "HE3": "CE", "HG3": "CG"}, 
		"PRO": {"H": "N", "HA": "CA", "HG2": "CG", "HB2": "CB", "HB3": "CB", "HD2": "CD", "HD3": "CD", "HG3": "CG"},
		"THR": {"H": "N", "HA": "CA", "HG21": "CG2", "HG22": "CG2", "HB": "CB", "HG1": "OG1", "HG23": "CG2"},
		"ALA": {"H": "N", "HA": "CA", "HB1": "CB", "HB2": "CB", "HB3": "CB"},
		"ASP": {"H": "N", "HA": "CA", "HB2": "CB", "HB3": "CB", "HD2": "OD2"},
		"ASN": {"H": "N", "HA": "CA", "HB2": "CB", "HB3": "CB", "HD21": "ND2", "HD22": "ND2"},
		"CYS": {"H": "N", "HA": "CA", "HB2": "CB", "HB3": "CB", "HG": "SG"}, 
		"HIS": {"H": "N", "HA": "CA", "HB2": "CB", "HB3": "CB", "HD1": "ND1", "HD2": "CD2", "HE1": "CE1", "HE2": "NE2"}, 
		"LEU": {"H": "N", "HA": "CA", "HB2": "CB", "HB3": "CB", "HD11": "CD1", "HD12": "CD1", "HD13": "CD1", "HD21": "CD2", "HD22": "CD2", "HD23": "CD2", "HG": "CG"}, 
		"PHE": {"H": "N", "HA": "CA", "HB2": "CB", "HB3": "CB", "HD1": "CD1", "HD2": "CD2", "HE1": "CE1", "HE2": "CE2", "HZ": "CZ"}, 
		"SER": {"H": "N", "HA": "CA", "HB2": "CB", "HB3": "CB", "HG": "OG"},  
		"TYR": {"H": "N", "HA": "CA", "HB2": "CB", "HB3": "CB", "HD1": "CD1", "HD2": "CD2", "HE1": "CE1", "HE2": "CE2", "HH": "OH"}
	}
	# Add general hbond acceptor partners for all residues
	hbond_acceptors_partners = {'N': 'CA', 'O': 'C', 'OE1': 'CD', 'OE2': 'CD', 'OD1': 'CG', 'OD2': 'CG', 'OH': 'CZ', 'ND1': 'CG', 'OG': 'CB', 'OG1': 'CB', 'NE2': 'CD'}
	for residue in partners:
		partners[residue].update(hbond_acceptors_partners)

	# Create a dictionary with amino acids as keys and lists of unique atom types as values
	atom_types_dict = {}
	for residue, atom_map in partners.items():
		unique_atoms = set(atom_map.keys()).union(atom_map.values())
		atom_types_dict[residue] = sorted(unique_atoms)
  
	# List of atoms that can be a part of the ring
	ring_atoms = {
		'PHE': ['CB', 'CD1', 'CD2', 'CE1', 'CE2', 'CG', 'CZ', 'HD1', 'HD2', 'HE1', 'HE2', 'HZ'],
		'TYR': ['CB', 'CD1', 'CD2', 'CE1', 'CE2', 'CG', 'CZ', 'HD1', 'HD2', 'HE1', 'HE2'],
		'TRP': ['CB', 'CD1', 'CD2', 'CE2', 'CG', 'HD1', 'HE1', 'HE3', 'CE3', 'CZ3', 'HZ3', 'CH2', 'HH2', 'CZ2', 'HZ2', 'NE1'],
		'HIS': ['CB', 'CD1', 'CD2', 'CE1', 'CE2', 'CG', 'HD1', 'HE1', 'HE2', 'CE3', 'CZ3', 'HZ3', 'CH2', 'HH2', 'CZ2', 'HZ2', 'ND1', 'NE2']
	}

	# Add ring atoms to partners if not already present
	for residue, atoms in ring_atoms.items():
		for atom in atoms:
			if atom not in atom_types_dict[residue]:
				atom_types_dict[residue].append(atom)
 
	return atom_types_dict


def xtctopdb(topology, trajectory, dyn_folder, filter_atoms=True):
    # Define the allowed atom types for each residue
	atom_types_dict = get_atom_types()

	u = mda.Universe(topology, trajectory)  # topology.pdb or .gro
	u.add_TopologyAttr('types', guess_types(u.atoms.names))
	u.add_TopologyAttr('elements', u.atoms.types)
	for ts in u.trajectory:
		selection = u.select_atoms("protein")
		frame = ts.frame  # Get the frame number from the timestep
		output = dyn_folder + "/" + tname[:-4] + "_" + str(frame) + ".pdb"
		if not os.path.exists(output):
			with PDBWriter(output, multiframe=True) as W:
				W.write(selection.atoms)
			if filter_atoms == True:
				# Filter the PDB file lines based on atom_types_dict
				filtered_lines = []
				with open(output, 'r') as pdb_file:
					for line in pdb_file:
						if line.startswith("ATOM"):
							atom_name = line[12:16].strip()
							residue_name = line[17:20].strip()
							if residue_name in atom_types_dict and atom_name in atom_types_dict[residue_name]:
								filtered_lines.append(line)
		
				# Sort the filtered lines by residue index
				filtered_lines.sort(key=lambda line: int(line[22:26].strip()))

				# Update atom numbering and write the filtered lines back to the PDB file
				updated_lines = []
				atom_index = 1
				for line in filtered_lines:
					if line.startswith("ATOM"):
						updated_line = line[:6] + f"{atom_index:5d}" + line[11:]
						updated_lines.append(updated_line)
						atom_index += 1
					else:
						updated_lines.append(line)
				with open(output, 'w') as pdb_file:
					pdb_file.writelines(updated_lines)
		else:
			print(f"Output file {output} already exists. Skipping...")

	return frame


def json_dict(path):
	"""Converts json file to pyhton dict."""
	json_file=open(path)
	json_str = json_file.read()
	json_data = loads(json_str)
	return json_data


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
	type=str,
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

parser.add_argument(
	'--get_outputs',
	dest='get_outputs',
	action='store_true',
	default=False,
	help='Get the outputs formatted with errors included'
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
cpu = args.cores
gpcrmd_only = args.gpcrmd_only
get_outputs = args.get_outputs
PH = "7"
# MEMORY_LIMIT_MB = int(args.mem)
# CHECK_INTERVAL = 10

if __name__ == '__main__':
    
	# Take original trajectory fileid from compl_data
	print("Reading GPCRmd information...")
	compl_data = json_dict(filespath+"/compl_info.json")

	# Directory where formated CS files are to be stores
	# outfolder = outpath+"Precomputed/chemical_shift/"
	# logfolder = outpath+"GPCRmd_precomputation/Chemical_shift/logs/"
 
	if get_outputs:
		print("Formatting output files with error calculations...")
		csv_files = [os.path.join(example, f) for f in os.listdir(example) if f.endswith('.csv')]
		for csv in csv_files:
			output_file = csv.replace('.csv', '').replace("cspred", "ucbshift") + ".txt"
			if not os.path.exists(output_file):
				print("Formatting file: ", csv)
				format_files(csv, csv.replace('.csv', '').replace("cspred", "ucbshift"), csv.replace('.csv', '_log.txt').replace("cspred", "ucbshift"))
			else:
				print(f"File {output_file} already exists. Skipping...")
		# format_files()
		sys.exit(0)

	# Extract dynamics ID to process
	if dyn_ids:
		dyns = dyn_ids
	else:
		dyns = [ str(compl_data[a]['dyn_id']) for a in compl_data]
 
	# Calculate shifts
	for dyn_id in dyns:

		print("Processing CSpred for dyn" + dyn_id)

		# Skip those not of GPCRmd community if so specified 
		is_community = compl_data['dyn' + dyn_id]['is_gpcrmd_community']
		if not is_community and gpcrmd_only:
			continue
		
		print("Checking directory example/dyn" + dyn_id)
		# Create a directory for the current dyn_id
		dyn_folder = os.path.join(example, f"dyn{dyn_id}")
		if not os.path.exists(dyn_folder):
			os.makedirs(dyn_folder)
		# Get trajectory and topology files
		trajectories = compl_data['dyn' + dyn_id]["traj_f"] 
		topology = filespath + compl_data['dyn' + dyn_id]["struc_f"] 
  
		# Create a directory for the output CSV files if it doesn't exist
		output_csv_dir = os.path.dirname(os.path.join(example, f"cspred_dyn{dyn_id}"))
		if not os.path.exists(output_csv_dir):
			os.makedirs(output_csv_dir)
   
		for trj in trajectories:
			tname = trj.split('/')[-1]
			output_csv = os.path.join(output_csv_dir, f"cs_cspred_dyn{dyn_id}_{tname.split('_')[0]}.csv")
			if not os.path.exists(output_csv):
				start_time = time.time()  # Start timing
				print("Processing CS for trajectory file... ", trj)
				subfolder_name = "trj" + tname.split('_')[0]
				subfolder_path = os.path.join(dyn_folder, subfolder_name)
				print("Checking subfolder... ", subfolder_path)
				if not os.path.exists(subfolder_path):
					os.makedirs(subfolder_path)
				# Transform xtc to pdb using MDAnalysis
				print("Generation of PDB files from trajectory...")
				last_frame = xtctopdb(topology, filespath + trj, subfolder_path, True)
	
				# Dictionary to hold aggregated data
				dict_data = {}
	
				# Iterate through all subdirectories and PDB files
				for root, dirs, files in os.walk(subfolder_path):
					# Collect all PDB files and their frame numbers
					pdb_files = []
					for file in files:
						if file.endswith('.pdb'):
							frame_number = int(file.split('_')[-1].split('.')[0])  # Extract frame number as integer
							pdb_files.append((frame_number, os.path.join(root, file)))

					# Sort PDB files by frame number
					pdb_files.sort(key=lambda x: x[0])

					# Limit to the first 3 frames for the example
					pdb_files = pdb_files

					# Iterate over sorted PDB files
					for frame_number, pdb_file in pdb_files:
						print("Running CSpred on file:", pdb_file)
						temp_output = os.path.join(subfolder_path, f'cspred_{frame_number}.csv')
						if not os.path.exists(temp_output):
							print(f"Analyzing the frame {frame_number} of {pdb_file}")
							subprocess.run(['python', cspred_path, pdb_file, '--output', temp_output, '-X', '--pH', PH, '--worker', cpu], check=True)

						# Read the generated CSV file
						print("Reading output CSV file:", temp_output)
						try:
							df = pd.read_csv(temp_output)
						except Exception as e:
							print(f"Error reading {temp_output}: {e}")
							continue

						print("Aggregate the information of each atom to the final csv file...")
						# Select columns containing 'UCBShift'
						ucbshift_columns = [col for col in df.columns if 'X' in col]
						filtered_df = df[['RESNUM', 'RESNAME'] + ucbshift_columns].copy()
						# Iterate over rows in the filtered DataFrame
						for i, row in filtered_df.iterrows():
							res_seq = row['RESNUM']
							res_name = row['RESNAME']
							res_name_s = three_to_one(res_name)
							ucbshift_values = row[ucbshift_columns].values
							for col, value in zip(ucbshift_columns, ucbshift_values):
								atom_name = col.split('_')[0]  # Extract atom name from column name
								key = f"{res_seq}_{atom_name}_{res_name_s}"
								if key not in dict_data:
									dict_data[key] = {
										'resSeq': res_seq,
										'name': atom_name,
										'resname': res_name,
										'resname_s': res_name_s
									}
								dict_data[key][frame_number] = round(value, 3)  # Round value to 3 decimals
						# # Remove the temporary output file after aggregation
						# if os.path.exists(temp_output):
						# 	os.remove(temp_output)

				end_time = time.time()  # End timing
				print(f"Processed {tname} in {end_time - start_time:.2f} seconds.")

				# Convert dict_data to a DataFrame directly
				print("Converting dict_data to DataFrame...")
				dict_df = pd.DataFrame.from_dict(dict_data, orient='index')

				# Drop rows with NaN values
				#print("Removing rows with NaN values...")
				#dict_df.dropna(inplace=True)

				# Save the DataFrame to a CSV file
				print(f"Saving data to {output_csv}...")
				dict_df.to_csv(output_csv, index=False, sep=';')
			else:
				print(f"File {output_csv} exist. Skip trajectory {trj} ")

