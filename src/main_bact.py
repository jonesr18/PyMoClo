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
from Bio.Restriction import BsaI
import Bio.GenBank as gb
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

# Constants
NUM_COLS = 6
VECTOR_IDX = 2
SIZE_IDX = 7
L0_IDX = [2, 3, 4, 5, 6] #includes vector
G00101 = 'GAGCGAGGAAGCGGAAGAGCG' # Level 0 forward seq primer
JH856R = 'GAAACCATTATTATCATGACATTAACC' # Level 0 reverse seq primer
SPEC_SEQ = 'CTGGAGAGAGCGAGATTCTCCGCGCTGTAGAAGTCACCATTGTTGTGCACGACGACATCATTCCGTGGCGTTATCCAGCTAAGCGCGAACTGCAATTTGGAGAATGGCAGCGCAATGACATTCTTGCAGGTATCTTCGAGCCAGC'

# Asus
L0_PATH = r'C:\\Users\\Weiss Lab\\Desktop\\Batch Export\\'
L1_PATH = r'C:\\Users\\Weiss Lab\\Desktop\\Batch Export\Level1s\\'
WORKBOOK = r'C:\\Users\\Weiss Lab\\Desktop\\PyClo_test.xlsx'

# HP
# L0_PATH = r'C:/Users/Ross/Dropbox (MIT)/APE/Genbank pL0s/'
# L1_PATH = r'C:/Users/Ross/Dropbox (MIT)/APE/Genbank pL1s/'
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


def writeL1(pL1_seq, file_ID):
    '''
    Writes the given pL1 record to a genbank file in the pL1 path with the given ID as a name.
    '''

    if DEBUG:
        print("\nMoclo Product End:\n")
        print(moclo_product)
        print('\tSize: ', len(moclo_product))

    with open(L1_PATH + file_ID.value + '.gb', 'w') as f:
        moclo_product.seq.alphabet = IUPAC.ambiguous_dna
        moclo_product.id = file_ID.value
        moclo_product.name = file_ID.value
        SeqIO.write(moclo_product, f, 'gb')


def digest(seq_record):
    '''
    Performs a digestion on the given record using BsaI, assuming a circular sequence. The fragments are then returned as sequence records.

    Returns None for any sequence inputs with fewer than 2 BsaI cut sites.
    '''
    # Problem w/ catalyze seems to be that it generates Seq
    #   instead of SeqRecord, which does not preserve annotations,
    #   so we do a manual search and "digestion" ourselves
    cut_pos = BsaI.search(seq_record.seq, linear = False)

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


if __name__ == '__main__':
    book = xw.Book(WORKBOOK)
    pL1_sheet = book.sheets['pL1s']
    for i, pL1_ID in enumerate(pL1_sheet.range('A1:A10000')):
        # | 1  |   2    |  3  |  4   | 5  | 6 | 7 | 8 | 9 | 10 | 11 | 12  | ....
        # | ID | Vector | Glc | Size | ng | I | P | 5 | G | 3  | T  | Lig | ....

        # print(ID)
        # print(ID.value)

        # Skip empty and spacer rows
        if not pL1_ID.value or pL1_ID.value == 'ID': continue

        # Extract entire row - note indexing is in [ ] since enum i is 0-based idx
        pL1_row = pL1_sheet[i, : NUM_COLS]

        # Skip rows where an ID is assigned but not any parts
        if not pL1_row(VECTOR_IDX).value: continue

        # Prepare pL1 output product
        moclo_product = SeqRecord(Seq('', IUPAC.ambiguous_dna))
        print('\nAssembling: ' + pL1_ID.value + \
              '. Expected size: ' + str(pL1_row(SIZE_IDX).value))

        # Iterate over pL0 files in pL1 part
        missingPiece = False
        for j in L0_IDX:
            pL0_ID = pL1_row(j)
            # Find filename in the dir of pL0st
            try:
               pL0_filename = getFilename(pL0_ID)
            except Exception as e:
               missingPiece = True
               print(e)
               break

            # Read and parse file
            # with open(r'../Test Files/' + filename, 'rU') as f:
            #     print(f.readline())
            with open(L0_PATH + pL0_filename, 'rU') as f:
                pL0_record = SeqIO.read(f, format = 'gb')
                print('\tAdding: ' + pL0_ID.value)

                fragments = digest(pL0_record)

                for i, fragment in enumerate(fragments):
                    # print('Fragment {0}: {1} bp'.format(i, len(fragment)))

                    # Check if the fragment is the size of known vectors, or if it has
                    # at least one of both G00101 and JH856R sites, implying backbone
                    if not (len(fragment) == 2247 or len(fragment) == 1130 or \
                            (fragment.seq.count(G00101) >= 1 and \
                            fragment.seq.count(JH856R) >= 1) or \
                            fragment.seq.count(SPEC_SEQ) >= 1):
                        moclo_product += fragment
        
        # If missing a piece, break out of loop and don't assemble
        if missingPiece:
           continue

        # Write to file
        writeL1(moclo_product, pL1_ID)
        print("\t\tFinished assembling " + pL1_ID.value + \
              '. Size = ' + str(len(moclo_product)))

        # TODO Check size is equal to size in excel sheet / calculated from formula
