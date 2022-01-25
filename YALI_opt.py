# Yarrowia Codon optimization script

def Codon_optimize_YL(AA_seq, forbidden_seqs = ["GCGGCCGC", "GGTCTC"], add_stop = True, GC_limit_upper = 80, GC_limit_lower = 20, window_size_gc = 12, window_size_poly = 8, lenient_GC=False, quiet=True):
    """ Script to convert an AA sequence into a yarrowia lipolytica codon optimized DNA sequence.

    Args:
        AA_seq (str): Amino acid sequence
        forbidden_seqs (list, optional): List of sequences to avoid, such as restriction sites. Defaults to ["GCGGCCGC", "GGTCTC"].
        add_stop (bool, optional): Adds final stop codon if True. Defaults to True.
        GC_limit_upper (int, optional): Upper limit for local GC content. Defaults to 80.
        GC_limit_lower (int, optional): Lower limit for local GC content. Defaults to 20.
        window_size_gc (int, optional): Window size for local GC content. Defaults to 12.
        window_size_poly (int, optional): Shortest stretch of polynucleotides not allowed. Defaults to 8.
        lenient_GC (bool, optional): Tries to stay within GC constraints, but may proceed for with high local GC content is unavoidable. Defaults to False.
        quiet (bool, optional): Silences print comments during the script. Defaults to True.

    Returns:
        DNA_seq (str): Codon optimized DNA sequence
    """

    #Common restriction enzymes
    #NotI = "GCGGCCGC" 
    #BsmBI = "CGTCTC"
    #BsaI = "GGTCTC"

    preferred_codons = {
        'A': 'GCC', 'C': 'TGC', 'D': 'GAC', 
        'E': 'GAG', 'F': 'TTC', 'G': 'GGT', 
        'H': 'CAC', 'I': 'ATC', 'K': 'AAG', 
        'L': 'CTC', 'M': 'ATG', 'N': 'AAC', 
        'P': 'CCC', 'Q': 'CAG', 'R': 'CGA', 
        'S': 'TCC', 'T': 'ACC', 'V': 'GTC', 
        'W': 'TGG', 'Y': 'TAC', '-': 'TAA'}
    alternative_codons = {
        'A': ['GCT'], 'C': ['TGT'], 'D': ['GAT'], 
        'E': [], 'F': [], 'G': ['GGC'], 
        'H': [], 'I': ['ATT'], 'K': [], 
        'L': ['CTG', 'CTT'], 'M': [], 'N': [], 
        'P': ['CCT'], 'Q': [], 'R': [], 
        'S': ['TCT'], 'T': ['ACT'], 'V': ['GTT'], 
        'W': [], 'Y': [], '-': []}
    yl_codon_usage = {
        'AAA' : 2.7,	'CAA' : 4.9,	'GAA' : 6.4,	'TAA' : 83.0,
        'AAC' : 95.8,	'CAC' : 85.2,	'GAC' : 70.0,	'TAC' : 95.1,
        'AAG' : 97.3,	'CAG' : 95.1,	'GAG' : 93.6,	'TAG' : 10.0,
        'AAT' : 4.2,	'CAT' : 14.8,	'GAT' : 30.0,	'TAT' : 4.9,
        'ACA' : 1.9,	'CCA' : 1.9,	'GCA' : 2.2,	'TCA' : 2.1,
        'ACC' : 73.6,	'CCC' : 74.0,	'GCC' : 59.1,	'TCC' : 51.1,
        'ACG' : 0.9,	'CCG' : 0.9,	'GCG' : 0.8,	'TCG' : 2.9,
        'ACT' : 23.6,	'CCT' : 23.1,	'GCT' : 37.9,	'TCT' : 38.5,
        'AGA' : 4.2,	'CGA' : 86.1,	'GGA' : 16.2,	'TGA' : 7.0,
        'AGC' : 4.1,	'CGC' : 1.2,	'GGC' : 30.6,	'TGC' : 65.6,
        'AGG' : 1.5,	'CGG' : 2.6,	'GGG' : 1.4,	'TGG' : 100.0,
        'AGT' : 1.3,	'CGT' : 4.3,	'GGT' : 51.9,	'TGT' : 34.4,
        'ATA' : 1.2,	'CTA' : 0.8,	'GTA' : 1.2,	'TTA' : 0.8,
        'ATC' : 65.0,	'CTC' : 46.7,	'GTC' : 56.4,	'TTC' : 80.7,
        'ATG' : 100.0,	'CTG' : 25.4,	'GTG' : 8.5,	'TTG' : 1.8,
        'ATT' : 33.8,	'CTT' : 24.5,	'GTT' : 33.9,	'TTT' : 19.3}
    #if AA_seq.upper == Only AA letters
    AA_seq = AA_seq.upper()
    DNA_seq = ""
    
    #* Setting constraints
    # Adding reverse complement of forbidden sequences
    for seq in forbidden_seqs:
        if reverse_complement(seq) not in forbidden_seqs:
            forbidden_seqs.append(reverse_complement(seq))

    # Adding forbidden polynucleotide stretches. 
    if window_size_poly > 0:
        forbidden_seqs.insert(0,"A"*window_size_poly)
        forbidden_seqs.insert(0,"C"*window_size_poly)
        forbidden_seqs.insert(0,"G"*window_size_poly)
        forbidden_seqs.insert(0,"T"*window_size_poly)

    #* Building DNA sequence
    ignore_GC = False
    if not quiet: print("-- Initiating codon optimization --")
    for AA in AA_seq:
        codon_approval = False
        cycle_counter = 0
        while not codon_approval:

            # Evaluate different codon combinations
            #print("position:",int((len(DNA_seq)/3)+1),"start round", cycle_counter)
            if cycle_counter == 0: #Preferred new codon
                temp_DNA_seq = DNA_seq + preferred_codons[AA]

            if cycle_counter == 1: #Alternative new codon
                try: 
                    if isinstance(alternative_codons[AA][0], str): 
                        temp_DNA_seq = DNA_seq + alternative_codons[AA][0]
                        if not quiet: print("suggesting alternate codon at position", int(len(temp_DNA_seq)/3))
                except:
                    if not quiet: print("No alternative codons available at position", int(len(temp_DNA_seq)/3))
                    cycle_counter += 2
                    continue
            
            if cycle_counter == 2: #Alternative 2 new codon
                try: 
                    if isinstance(alternative_codons[AA][1], str):
                        temp_DNA_seq = DNA_seq + alternative_codons[AA][1]
                        if not quiet: print("suggesting alternate codon at position", int(len(temp_DNA_seq)/3))
                except:
                    if not quiet: print("No second alternative codon available at position", int(len(temp_DNA_seq)/3))
                    cycle_counter += 1
                    continue

            if cycle_counter == 3: #Alternative previous codon + Preferred new codon
                last_AA = Codon_lookup(DNA_seq[-3:])
                try: 
                    if isinstance(alternative_codons[last_AA][0], str): 
                        temp_DNA_seq = DNA_seq[:-3] + alternative_codons[last_AA][0] + preferred_codons[AA]
                        if not quiet: print("Suggesting alternate codon for previous codon at position", int(len(temp_DNA_seq)/3)-1)
                except:
                    if not quiet: print("No alternative codons available at position", int(len(temp_DNA_seq)/3)-1)
                    cycle_counter += 3
                    continue
            
            if cycle_counter == 4: #Alternative 2 previous codon + Preferred new codon
                #last_AA = Codon_lookup(DNA_seq[-3:])
                try: 
                    if isinstance(alternative_codons[last_AA][1], str): 
                        temp_DNA_seq = DNA_seq[:-3] + alternative_codons[last_AA][1] + preferred_codons[AA]
                        if not quiet: print("Suggesting alternate codon for previous codon at position", int(len(temp_DNA_seq)/3)-1)
                except:
                    if not quiet: print("No second alternative codon available at position", int(len(temp_DNA_seq)/3)-1)
                    cycle_counter += 1
                    continue

            if cycle_counter == 5: #Alternative previous codon + Alternative new codon
                #last_AA = Codon_lookup(DNA_seq[-3:])
                try: 
                    if isinstance(alternative_codons[last_AA][0], str) and isinstance(alternative_codons[AA][0], str): 
                        temp_DNA_seq = DNA_seq[:-3] + alternative_codons[last_AA][0] + alternative_codons[AA][0]
                        if not quiet: print("Suggesting alternate codons for previous and current codon at positions", int(len(temp_DNA_seq)/3)-1, "and", int(len(temp_DNA_seq)/3))
                except:
                    if not quiet: print("No alternative codons available at either position", int(len(temp_DNA_seq)/3), "or position", int(len(temp_DNA_seq)/3)-1)
                    cycle_counter += 1
                    continue

            if cycle_counter == 6: #Alternative second last codon + previously selected new codon + Preferred new codon
                second_last_AA = Codon_lookup(DNA_seq[-6:-3])
                try: 
                    if isinstance(alternative_codons[second_last_AA][0], str): 
                        temp_DNA_seq = DNA_seq[:-6] + alternative_codons[second_last_AA][0] + DNA_seq[-3:] + preferred_codons[AA]
                        if not quiet: print("Suggesting alternate codon for second last codon, at position", int(len(temp_DNA_seq)/3)-2)
                except:
                    if not quiet: print("No alternative codons available at position", int(len(temp_DNA_seq)/3)-2)
                    cycle_counter += 1
                    continue

            if cycle_counter == 7:
                if not ignore_GC:
                    if not quiet: print("DNA sequence could not be generated with the constraints given")
                    if lenient_GC:
                        if not quiet: print("proceeding with lenient GC restrictions")
                        cycle_counter = 0
                        ignore_GC = True
                        continue
                    else: 
                        print("Sequence generation aborted. Consider setting: lenient_GC=True")
                        print("PARTIAL DNA SEQ:", DNA_seq + preferred_codons[AA])
                        return "ERROR: DNA sequence could not be generated with the constraints given"
                if ignore_GC:
                    print("DNA sequence could not be generated with the constraints given")
                    print("Sequence generation aborted")
                    print("PARTIAL DNA SEQ:", DNA_seq + preferred_codons[AA])
                    return "ERROR: DNA sequence could not be generated with the constraints given"

            seq_length = len(temp_DNA_seq)

            #Check forbidden sequences
            forbidden_seq_status = True
            for seq in forbidden_seqs:
                if seq_length >= len(seq) +2:
                    if seq in temp_DNA_seq[-(len(seq) +2):]:
                        if not quiet: print("Position:", int(seq_length/3),"Forbidden seq detected:", seq)
                        if not quiet: print("...", temp_DNA_seq[-12:], sep="")
                        forbidden_seq_status = False
                elif seq_length >= len(seq):
                    if seq in temp_DNA_seq:
                        if not quiet: print("Position:", int(seq_length/3),"Forbidden seq detected:", seq)
                        if not quiet: print(temp_DNA_seq)
                        forbidden_seq_status = False

            #Check GC content
            GC_content_status = True
            if seq_length >= window_size_gc:
                GC_cont = GC_content(temp_DNA_seq[-window_size_gc:])
                if GC_cont > GC_limit_upper:
                    GC_content_status = False
                    if not quiet: print("Position:", int(seq_length/3), "GC content:", round(GC_cont, 1), "%")
                if GC_cont < GC_limit_lower:
                    GC_content_status = False
                    if not quiet: print("Position:", int(seq_length/3), "GC content:", round(GC_cont, 1), "%")
            if ignore_GC:
                GC_content_status = True #Consider adding stepwise GC limit increase instead

            #If OK, add codon to DNA seq
            if forbidden_seq_status and GC_content_status:
                DNA_seq = temp_DNA_seq
                codon_approval = True
            else:
                cycle_counter += 1

    #Final QC and if OK, return codon optimized DNA Sequence
    if not quiet: print("-- End of sequence --")
    if len(DNA_seq)/3 == len(AA_seq):
        new_AA_seq = translate(DNA_seq)
        if new_AA_seq == AA_seq:
            if not quiet: print("Sequence length and AA sequence verified.")
            #print(DNA_seq)
            if add_stop:
                if Codon_lookup(DNA_seq[-3:]) != "-":
                    DNA_seq += "TAA"
            return DNA_seq

        else: #Incorrect AA sequence
            print("ERROR, Incorrect AA sequence:", translate(DNA_seq))
            return "ERROR, Incorrect AA sequence"
    else: #Incorrect length
        print("ERROR, unexpected nucleotides inserted")
        print("Incorrect DNA:", DNA_seq)
        
        return "ERROR, unexpected nucleotides inserted"

#Check forbidden seqs, GC, and poly nucl, again.



def Codon_lookup(input, rev=False):
    """Returns the amino acid encoded by input codon, or if rev=True, returns the a list of all codons encoding the input amino acid.

    Args:
        input (str): codon or single letter amino acid
        rev (bool, optional): rev=True searches the reversed codon table. Defaults to False.

    Returns:
        if rev == False
        Str: Amino Acid
        if rev == True
        List: list of codons
    """
    codon_table = {
        "AAA" : "K", "AAC" : "N", "AAG" : "K", "AAT" : "N",
        "ACA" : "T", "ACC" : "T", "ACG" : "T", "ACT" : "T",
        "AGA" : "R", "AGC" : "S", "AGG" : "R", "AGT" : "S",
        "ATA" : "I", "ATC" : "I", "ATG" : "M", "ATT" : "I",
        "CAA" : "Q", "CAC" : "H", "CAG" : "Q", "CAT" : "H",
        "CCA" : "P", "CCC" : "P", "CCG" : "P", "CCT" : "P",
        "CGA" : "R", "CGC" : "R", "CGG" : "R", "CGT" : "R",
        "CTA" : "L", "CTC" : "L", "CTG" : "L", "CTT" : "L",
        "GAA" : "E", "GAC" : "D", "GAG" : "E", "GAT" : "D",
        "GCA" : "A", "GCC" : "A", "GCG" : "A", "GCT" : "A",
        "GGA" : "G", "GGC" : "G", "GGG" : "G", "GGT" : "G",
        "GTA" : "V", "GTC" : "V", "GTG" : "V", "GTT" : "V",
        "TAA" : "-", "TAC" : "Y", "TAG" : "-", "TAT" : "Y",
        "TCA" : "S", "TCC" : "S", "TCG" : "S", "TCT" : "S",
        "TGA" : "-", "TGC" : "C", "TGG" : "W", "TGT" : "C",
        "TTA" : "L", "TTC" : "F", "TTG" : "L", "TTT" : "F"}
    reversed_codon_table = {
        'A': ['GCA', 'GCC', 'GCG', 'GCT'], 
        'C': ['TGC', 'TGT'], 
        'D': ['GAC', 'GAT'], 
        'E': ['GAA', 'GAG'], 
        'F': ['TTC', 'TTT'], 
        'G': ['GGA', 'GGC', 'GGG', 'GGT'], 
        'H': ['CAC', 'CAT'], 
        'I': ['ATA', 'ATC', 'ATT'], 
        'K': ['AAA', 'AAG'], 
        'L': ['CTA', 'CTC', 'CTG', 'CTT', 'TTA', 'TTG'], 
        'M': ['ATG'], 
        'N': ['AAC', 'AAT'], 
        'P': ['CCA', 'CCC', 'CCG', 'CCT'], 
        'Q': ['CAA', 'CAG'], 
        'R': ['AGA', 'AGG', 'CGA', 'CGC', 'CGG', 'CGT'], 
        'S': ['AGC', 'AGT', 'TCA', 'TCC', 'TCG', 'TCT'], 
        'T': ['ACA', 'ACC', 'ACG', 'ACT'], 
        'V': ['GTA', 'GTC', 'GTG', 'GTT'], 
        'W': ['TGG'], 
        'Y': ['TAC', 'TAT'],
        '-': ['TAA', 'TAG', 'TGA']}
    if not rev:
        amino_acid = codon_table[input.upper()]
        return amino_acid
    if rev:
        list_of_codons = reversed_codon_table[input.upper()]
        return list_of_codons

    else:
        print("ERROR: rev expected to be True or False")


def GC_content(query):
    query = str(query).upper()
    c_count = query.count("C")
    g_count = query.count("G")
    GC_cont = 100*((c_count+g_count)/len(query))
    return GC_cont


def reverse_complement(DNA):
    """Returns the reverse complement to a given DNA string, allows for U and X as well as extended IUPAC

    Args:
        DNA (string): DNA sequence to be converted

    Returns:
        rev_comp_DNA (string): Reverse complement of DNA
    """
    Complement_DNA = { #Dictionary of complementary basepairs
    "A" : "T",
    "C" : "G",
    "G" : "C",
    "T" : "A",
    "U" : "X", #Unconventional - made up by me for USER cloning
    "X" : "U", #Unconventional - made up by me for USER cloning   
    "R" : "Y", #R = A or G; Purines
    "Y" : "R", #Y = T or C; Pyrimidines
    "S" : "S", #S = C or G; Strong bp
    "W" : "W", #W = A or T; Weak bp
    "K" : "M", #K = G or T; Keto
    "M" : "K", #M = A or C; Amino
    "B" : "V", #B = Not A
    "D" : "H", #D = Not C
    "H" : "D", #H = Not G
    "V" : "B", #V = Not T
    "N" : "N"  #N = any nucleotide
    }

    rev_comp_DNA = ""
    
    for nt in DNA.upper():
        rev_comp_DNA = Complement_DNA[nt] + rev_comp_DNA
    
    return rev_comp_DNA


def translate(ORF): 
    """Translate ORF into AA sequence, return AA sequence and molecular weight (MW)

    Args:
        ORF (str): String consisting of standard DNA nucleotides.

    Returns:
        AA_seq (str): Translated amino acid sequence
        MW (float): molecular weight for translated protein
    """
    Codon_table = {
        "AAA" : "K", "AAC" : "N", "AAG" : "K", "AAT" : "N",
        "ACA" : "T", "ACC" : "T", "ACG" : "T", "ACT" : "T",
        "AGA" : "R", "AGC" : "S", "AGG" : "R", "AGT" : "S",
        "ATA" : "I", "ATC" : "I", "ATG" : "M", "ATT" : "I",
        "CAA" : "Q", "CAC" : "H", "CAG" : "Q", "CAT" : "H",
        "CCA" : "P", "CCC" : "P", "CCG" : "P", "CCT" : "P",
        "CGA" : "R", "CGC" : "R", "CGG" : "R", "CGT" : "R",
        "CTA" : "L", "CTC" : "L", "CTG" : "L", "CTT" : "L",
        "GAA" : "E", "GAC" : "D", "GAG" : "E", "GAT" : "D",
        "GCA" : "A", "GCC" : "A", "GCG" : "A", "GCT" : "A",
        "GGA" : "G", "GGC" : "G", "GGG" : "G", "GGT" : "G",
        "GTA" : "V", "GTC" : "V", "GTG" : "V", "GTT" : "V",
        "TAA" : "-", "TAC" : "Y", "TAG" : "-", "TAT" : "Y",
        "TCA" : "S", "TCC" : "S", "TCG" : "S", "TCT" : "S",
        "TGA" : "-", "TGC" : "C", "TGG" : "W", "TGT" : "C",
        "TTA" : "L", "TTC" : "F", "TTG" : "L", "TTT" : "F"}

    i = 0
    ORF = ORF.upper()
    AA_seq = ""
    
    while i < len(ORF): #Translate into AA
        AA_seq = AA_seq + Codon_table[ORF[i:i+3]]
        i += 3    
    
    return AA_seq






#Usage example
my_AA_seq = "MVSRYVPDMGDLIWVDFDPTKGSEQAGHRPAVVLSPFMYNNKTGMCLCVPCTTQSKGYPFEVVLSGQERDGVALADQVKSIAWRARGATKKGTVAPEELQLIKAKINVLIG"

print(Codon_optimize_YL(my_AA_seq))

#! Arguments to consider:
# forbidden_seqs (list, optional): List of sequences to avoid, such as restriction sites. Defaults to ["GCGGCCGC", "GGTCTC"].
# GC_limit_upper (int, optional): Upper limit for local GC content. Defaults to 80.
# GC_limit_lower (int, optional): Lower limit for local GC content. Defaults to 20.
# window_size_gc (int, optional): Window size for local GC content. Defaults to 12.
# window_size_poly (int, optional): Shortest stretch of polynucleotides not allowed. Defaults to 8.
# lenient_GC (bool, optional): Tries to stay within GC constraints, but may proceed for with high local GC content is unavoidable. Defaults to False.
