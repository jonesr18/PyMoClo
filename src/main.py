import xlwings as xw
from pathlib import Path
from Bio.Restriction import BsaI
import Bio.GenBank as gb
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

NUM_COLS = 11
VECTOR_IDX = 2
G00101 = 'ATTACCGCCTTTGAGTGAGC'
JH856R = 'CAGGGTTATTGTCTCATGAGC'

def getFilename(fileID):
    '''
    Finds a .gb, .ape, or .str file associated with the given ID.
    If none is found, the method returns <None>
    '''

    # Handle None-types
    if not fileID.value: return None

    # Search for files in dir, accept the first match (prioritize gb > ape > str)
    for ext in ['.gb', '.gbk', '.ape', '.str']:
        if Path(r'../Test Files/' + str(fileID.value) + ext).exists():
            print('\n' + str(fileID.value) + ext)
            return str(fileID.value) + ext

    print('No file found for ID: ' + str(fileID.value))
    return None



if __name__ == '__main__':
    # print(dir(xw))
    book = xw.Book(r'C:/Users/Ross/Dropbox (MIT)/Team Herpes/Ross/RDJ Plasmids.xlsx')
    pL1s = book.sheets['pL1s']

    for i, ID in enumerate(pL1s.range('A1:A4')):
        # | 1  |   2    |  3  |  4   | 5  | 6 | 7 | 8 | 9 | 10 | 11 | 12  | ....
        # | ID | Vector | Glc | Size | ng | I | P | 5 | G | 3  | T  | Lig | ....

        # print(ID)
        # print(ID.value)

        # Skip empty and spacer rows
        if not ID.value or ID.value == 'ID': continue

        # Extract entire row - note indexing is in [ ] since enum i is 0-based idx
        row = pL1s[i, :NUM_COLS]

        # Skip rows where an ID is assigned but not any parts
        if not row(VECTOR_IDX).value:
            # print("ID: {0}\n".format(ID.value),
            #       "No vector assigned - skipping entry")
            continue

        # Prepare pL1 output product
        moclo_product = SeqRecord(Seq('', IUPAC.ambiguous_dna), id = ID.value, name = ID.value)
        print("\nMoclo Product Start:\n")
        print(moclo_product)

        # Iterate over pL0 files in pL1 part
        for pL0_ID in row[1:]:
            # Find filename in the dir of pL0s
            filename = getFilename(pL0_ID)
            if not filename: continue # Skip if no such filename

            # Read and parse file
            # with open(r'../Test Files/' + filename, 'rU') as f:
            #     print(f.readline())
            with open(r'../Test Files/' + filename, 'rU') as f:
                # rp = gb.RecordParser(2);
                # seq_record = rp.parse(f)
                seq_record = SeqIO.read(f, format = 'gb')
                print(seq_record.name)
                # print(seq_record.seq)
                # print(seq_record)

                # Problem w/ catalyze seems to be that it generates Seq
                #   instead of SeqRecord, which does not preserve annotations
                # fragments = BsaI.catalyze(seq_record.seq, linear = False)
                cut_pos = BsaI.search(seq_record.seq, linear = False)
                fragments = [seq_record[cut_pos[0] - 1 : cut_pos[1] - 1], \
                             seq_record[: cut_pos[0] - 1] + seq_record[cut_pos[1] - 1 :]]
                print(fragments)

                for i, fragment in enumerate(fragments):
                    # print('Fragment {0}: {1} bp'.format(i, len(fragment)))

                    # Check if the fragment is the size of known vectors, or if it has
                    # at least one of both G00101 and JH856R sites, implying backbone
                    if not (len(fragment) == 2247 or len(fragment) == 1130 or \
                            (fragment.seq.count(G00101) > 1 and \
                            fragment.seq.count(JH856R) > 1)):
                        moclo_product += fragment

        # Write to file
        with open(r'../Test Files/' + ID.value + '.gb', 'w') as f:
            print("\nMoclo Product End:\n")
            print(moclo_product)
            print('Size:\n', len(moclo_product))
            moclo_product.seq.alphabet = IUPAC.ambiguous_dna
            moclo_product.id = ID.value
            moclo_product.name = ID.value
            SeqIO.write(moclo_product, f, 'gb')


        # TODO Check size is equal to size in excel sheet / calculated from formula
