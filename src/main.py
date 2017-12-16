'''
This script is for generating Level 1 genbank files using an excel sheet with
maps showing which Level 0s belong in each Level 1, and in what order.

Template: RDJ Plasmids.xlsx

Post-processing:
 - A2 insulator annotation disappears
    (and probably others spanning an internal BsaI site too)
 - Annotations with spaces gain underscores
    (notable: gg overhang and gibson overhang)
 - Some primers span 1000s of bp for some reason, crossing parts
    (delete, then re-test for 1030F/R)

TODO:
 - Take Level 0 part names from DEFINITION field, and combine them with
   underscores to get pL1 definitions
    (eg: '[S2-3] A2_Inert_hEF1a_mKate2_Inert_Synth' comes from:
        S2-3 [NAME!!] | I A2 | P.2 hEF1a | 5.2 Inert | G mKate2 | 3 Inert | T Synth)
'''

import xlwings as xw
from pathlib import Path
from Bio.Restriction import BsaI	# Enzyme used for moclo
from Bio.Restriction import BsmBI	# For digesting FuRJ
import Bio.GenBank as gb
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

# Constants
LV_BB_NAME = 'pLV-RJ'				# Name of the LV cloning backbone
VECTOR_IDX = 2						# (1-based) Index of vector to clone into for quick check
SIZE_IDX = 4						# (1-based) Index of epected size in worksheet
L0_IDX = [2, 6, 7, 8, 9, 10, 11]	# (1-based) Index of pL0 parts in worksheet, including vector
L0_FRAG_SIZE = 2247					# L0 Backbone size (excluded during assembly)
GIBSON_VEC_FRAG_SIZE = 1130			# Gibson position vector backbone size (excluded during assembly)
FURJ_VEC_FRAG_SIZE = 1054			# FuRJ Lenti vector backbone size (excluded during assembly)
EXCLUDED_SIZES = [L0_FRAG_SIZE, GIBSON_VEC_FRAG_SIZE, FURJ_VEC_FRAG_SIZE]
G00101 = 'ATTACCGCCTTTGAGTGAGC'		# Forward L0 sequencing primer for additional help IDing L0 backbones
JH856R = 'CAGGGTTATTGTCTCATGAGC'	# Reverse L0 "                                      "

# Asus
L0_PATH = r'C:/Users/jones/Dropbox (MIT)/APE/Genbank pL0s/'		# Directory of L0 .gb files
L1_PATH = r'C:/Users/jones/Dropbox (MIT)/APE/Genbank pL1s/'		# Directory to deposit L1 .gb files
VV_PATH = r'C:/Users/jones/Dropbox (MIT)/APE/Genbank pVVs/'		# Directory to deposit VV .gb files
WORKBOOK = r'C:/Users/jones/Dropbox (MIT)/Team Herpes/Ross/RDJ Plasmids.xlsx' # Excel workbook with L0-->L1/VV mappings

# HP
# L0_PATH = r'C:/Users/Ross/Dropbox (MIT)/APE/Genbank pL0s/'
# L1_PATH = r'C:/Users/Ross/Dropbox (MIT)/APE/Genbank pL1s/'
# VV_PATH = r'C:/Users/Ross/Dropbox (MIT)/APE/Genbank pVVs/'
# WORKBOOK = r'C:/Users/Ross/Dropbox (MIT)/Team Herpes/Ross/RDJ Plasmids.xlsx'

# Debug status
DEBUG = False

def getFilename(file_ID):
	'''
	Finds a .gb, .ape, or .str file associated with the given ID.
	If none is found, the method returns <None>
	'''

	# Handle None-types
	if not file_ID.value: raise Exception('File ID is None')

	# Search for files in dir, accept the first match (prioritize gb > ape > str)
	for ext in ['.gb', '.gbk', '.ape', '.str']:
		if Path(L0_PATH + str(file_ID.value) + ext).exists():
			if DEBUG: print('\nFound file: ' + str(file_ID.value) + ext)

			return str(file_ID.value) + ext

	raise Exception('No file found for ID: ' + str(file_ID.value))


def digest(seq_record, enzyme):
	'''
	Performs a digestion on the given record using the given enzyme, assuming a 
	circular sequence. The fragments are then returned as sequence records.

	Returns None for any sequence inputs with fewer than 2 BsaI cut sites.
	'''

	# Problem w/ catalyze seems to be that it generates Seq
	#   instead of SeqRecord, which does not preserve annotations,
	#   so we do a manual search and "digestion" ourselves
	cut_pos = enzyme.search(seq_record.seq, linear = False)

	# There must be at least two cuts to generate fragments
	if len(cut_pos) < 2: return None

	# Auto-extract all fragments "internal" to the sequence
	#   There will be some sequence left over on the "edges",
	#   which must be manually dealt with, since python
	#   indexing doesn't work around corners.
	fragments = [seq_record[cut_pos[i] - 1 : cut_pos[i + 1] - 1] for i in range(len(cut_pos) - 1)]
	fragments.append(seq_record[cut_pos[-1] - 1 :] + seq_record[: cut_pos[0] - 1])

	if DEBUG: print(fragments)
	return fragments


def processMoclo(worksheet, path, start_idx = 1, end_idx = float('inf')):
	'''
	Takes in a worksheet of pL0-->pL1/pVV mappings and generates sequences 
	for the constructs, writing them to file. 
	'''

	for i, assembly_ID in enumerate(worksheet.range('A1:A10000')):
    	# | 1  |   2    |  3  |  4   | 5  | 6 | 7 | 8 | 9 | 10 | 11 | 12  | ....
    	# | ID | Vector | Glc | Size | ng | I | P | 5 | G | 3  | T  | Lig | ....
		
		if i < start_idx:
			continue
		if i > end_idx:
			break

		if DEBUG:
			print(assembly_ID)
			print(assembly_ID.value)

		# Skip empty and spacer rows
		if not assembly_ID.value or assembly_ID.value == 'ID': continue

		# Extract entire row - note indexing is in [ ] since enum i is 0-based idx
		assembly_row = worksheet[i, : L0_IDX[-1]]

		# Skip rows where an ID is assigned but not any parts
		if not assembly_row(VECTOR_IDX).value: continue

		# Prepare moclo product seq record
		moclo_product = SeqRecord(Seq('', IUPAC.ambiguous_dna))
		moclo_product_description = ""
		print('\nAssembling: ' + assembly_ID.value + \
			  '. Expected size: ' + str(assembly_row(SIZE_IDX).value))

		# Setup dictionary of outputs
		# moclo_products = dict()

		# Iterate over pL0 files for assembly
		missingPiece = False
		isLVRJ = False
		for j in L0_IDX:
			pL0_ID = assembly_row(j)
			# Find filename in the dir of pL0st
			try:
				pL0_filename = getFilename(pL0_ID)
			except Exception as e:
				missingPiece = True
				print(e)
				break

			# Read and parse file
			with open(L0_PATH + pL0_filename, 'rU') as f:
				pL0_record = SeqIO.read(f, format = 'gb')
				print('\tAdding: ' + pL0_ID.value)
				if (j > L0_IDX[0]): 	# Parts
					# In .gb files, the description comes in form "[Vector] Description [more description]."
					# 	We want to keep everything after [Vector] together and exclude the period.
					pL0_description = pL0_record.description[:-1].split(' ', 1)[1]
					if (j > L0_IDX[1]): # 2nd+ Part
						moclo_product_description += '_' +  pL0_description
					else: 				# 1st Part
						moclo_product_description += pL0_description
				else:					# Vector
					moclo_product_description += '[' + pL0_ID.value + '] ' # Use vector name from Excel Sheet

				# Cut w/ BsaI unless looking at LV backbone, which is cut w/ BsmBI
				if str(pL0_ID.value) == LV_BB_NAME:
					fragments = digest(pL0_record, BsmBI)
					isLVRJ = True
				else:
					fragments = digest(pL0_record, BsaI)

				for i, fragment in enumerate(fragments):
					if DEBUG: 
						print('Fragment {0}: {1} bp'.format(i, len(fragment)))
						print(type(fragment))
						print(type(fragment.seq))

					# Check if the fragment is the size of known vectors, or if it has
					# at least one of both G00101 and JH856R sites, implying backbone
					if not (len(fragment) in EXCLUDED_SIZES or \
							(fragment.seq.count(G00101) >= 1 and \
							 fragment.seq.count(JH856R) >= 1)):
						moclo_product += fragment
		
		# If missing a piece, break out of loop and don't assemble
		if missingPiece:
			continue

		# Save to dictionary
		# print(str(assembly_ID.value))
		# moclo_products[str(assembly_ID.value)] = moclo_product
		# print(moclo_products[str(assembly_ID.value)])
		
		# If a pLV-RJ, move starting point by 1000 bases so it looks 
		# nicer when linearized in Geneious
		if isLVRJ:
			moclo_product = moclo_product[1000:] + moclo_product[:1000]

		print("\t\tFinished assembling " + assembly_ID.value + \
		  '. Size = ' + str(len(moclo_product)))
		if DEBUG: 
			print("\nMoclo Product End:\n")
			print(moclo_product)

		writeAssembly(moclo_product, assembly_ID.value, moclo_product_description, path)

def writeAssembly(moclo_product, name, description, path):
	'''
    Writes the given moclo product to a genbank file with the given name 
    and description in the given path. The GenBank ID is set as the 
    filename as well.
    '''

	# Write file
	filename = path + name + '.gb'
	print('\t\tWriting file: ' + filename)
	with open(filename, 'w') as f:
		moclo_product.seq.alphabet = IUPAC.ambiguous_dna
		moclo_product.id = name
		moclo_product.name = name
		moclo_product.description = description
		SeqIO.write(moclo_product, f, 'gb')


##### Main Method #####

if __name__ == '__main__':
    
	# Defualt parameters
	startIdx = 1
	endIdx = 10**5

    # Handle argument inputs
	import sys
	argv = sys.argv[1:]
	if len(argv) > 0:
		startIdx = int(argv[0])
		if len(argv) > 1:
			endIdx = int(argv[1])
	
	# Open worksheet with mappings
	book = xw.Book(WORKBOOK)
    
	# Generate moclo products
	pL1_products = processMoclo(book.sheets['pL1s'], L1_PATH, start_idx = startIdx, end_idx = endIdx)
	# pVV_products = processMoclo(book.sheets['pVVs'], VV_PATH)
	

	# TODO Check size is equal to size in excel sheet / calculated from formula
