from structures import NUCLEOTIDE_BASE, DNA_Codons, RNA_Codons
import random
from collections import Counter

# terminal cmd for GUI format change: pyuic5 "filename.ui" -o "filename.py"

class BioSeq:
    """
    DNA sequence class. Default value: ATCG, DNA, no label
    \nMany functions operate in tandem with the structures file containing dicts: NUCLEOTIDE_BASE, DNA_Codons, RNA_Codons
    \nThese are mainly for translation and validation
    """

    def __init__(self, seq= "ATCG", seq_type= "DNA", label='None'):
        """Sequence initialization, validation"""
        self.seq = seq.upper()
        self.label = label
        self.seq_type = seq_type

    # DNA Toolkit Functions

    def validate(self):
        return set(NUCLEOTIDE_BASE[self.seq_type]).issuperset(self.seq)


    def get_seq_biotype(self):
        return self.seq_type


    def get_seq_info(self):
        """returns 4 strings. Full sequence information"""
        return f"[Label]: {self.label} \n[Sequence]: {self.seq} \n[Type]: {self.seq_type} \n[Length]: {len(self.seq)}"


    def gen_random_seq(self, length= 10, seq_type="DNA"):
        """generate random DNA sequence: can use this for testing"""
        seq = ''.join((random.choice(NUCLEOTIDE_BASE[seq_type])
                       for _ in range(length)))
        self.__init__(seq, seq_type, "Random Gen. Seq.") # re-intializes / updates class with random seq.


    def nuc_freq(self):
        """counts the frequency of each nucleotide, returns dict"""
        return dict(Counter(self.seq))


    def transcription(self):
        """DNA -> RNA Transcription, replace T - U"""
        return self.seq.replace('T', 'U')


    def reverse_complement(self):
        """Complementary Base Pairing + Reversing String"""
        if self.seq_type == "DNA":
            mapping = str.maketrans("ATCG", "TAGC")
        else:
            mapping = str.maketrans("AUCG", "UAGC")
        return self.seq.translate(mapping)[::-1]


    def gc_content(self):
        """GC Content in DNA/RNA Sequence"""
        return round((self.seq.count('G') + self.seq.count('C')) / len(self.seq) * 100)


    def gc_subsec(self, k=20):
        """Calculates GC content in subsections of a specific length k"""
        sec_content = []
        subsec_list = []
        for i in range(0, len(self.seq) - k + 1, k): # move through sequence with step 'k'
            subsec = self.seq[i:i + k]
            subsec_list.append(subsec) # calculate GC as we go along, adding values to sec_content
            sec_content.append(round((subsec.count('G') + subsec.count('C')) / len(subsec) * 100)
)
        return sec_content


    def translate_seq(self , init_pos=0):
        """Translate from DNA into AA seq."""
        if self.seq_type == "DNA":
            return [DNA_Codons[self.seq[pos:pos + 3]] for pos in range(init_pos, len(self.seq) - 2, 3)]
        elif self.seq_type == "RNA":
            return[RNA_Codons[self.seq[pos:pos + 3]] for pos in range(init_pos, len(self.seq) - 2, 3)]


    def codon_frequency(self, aa):
        """Provides frequency of specific codons of a given amino acid"""
        tmpList = []  # stores codons for aas
        if self.seq_type == "DNA":
            for i in range(0, len(self.seq) - 2, 3):  # iterates through seq codon by codon
                if DNA_Codons[self.seq[i:i + 3]] == aa:  # translated codon slice is checked against amino acid search
                    tmpList.append(self.seq[i:i + 3])  # if matches, add it to temp list

        elif self.seq_type == "RNA":
            for i in range(0, len(self.seq) - 2, 3):
                if RNA_Codons[self.seq[i:i + 3]] == aa:
                    tmpList.append(self.seq[i:i + 3])

        freqDict = dict(Counter(tmpList))  # Counter makes a dictionary that counts key frequency
        total = sum(freqDict.values())  # Generates sum of all codon key freq. for that one aa search
        for target_aa in freqDict.keys():  # codon keys operations
            freqDict[target_aa] = round(freqDict[target_aa] / total, 2)
        return freqDict


    def gen_reading_frames(self):
        """Generate 6 reading frames of DNA seq, including reverse complement"""
        # when analyzing raw DNA, you might not know the correct orientation/frame
        # genes can occur on either strand of DNA, starting in any of the 3 possible frames
        # alternate ORFs, some genes use overlapping sequences that use alternate reading frames
        # viruses have space constraints and thus highly efficient coding
        # bioinformatics pipeline often analyzes seqs without knowing where translation begins
        frames = []
        frames.append(self.translate_seq(0))
        frames.append(self.translate_seq(1))
        frames.append(self.translate_seq(2))
        tmp_seq = BioSeq(self.reverse_complement(), self.seq_type)
        frames.append(tmp_seq.translate_seq(0))
        frames.append(tmp_seq.translate_seq(1))
        frames.append(tmp_seq.translate_seq(2))
        del tmp_seq
        return frames


    def reading_frame_proteins(self, seq): # can ignore this error
        """Compute all possible proteins in an aminoacid seq + return list of possible proteins"""
        current_prot = []
        proteins = []
        for aa in seq:  # iterates through string
            if aa == "_":  # STOP accumulating amino acids if _ is found
                if current_prot:  # if there's stuff in current prot
                    for p in current_prot:
                        proteins.append(p)  # add everything to the proteins[]
                    current_prot = []  # empty current_prot[]
            else:
                if aa == "M":  # START accumulating amino acids if M START was found
                    current_prot.append("")  # add empty space to current prot, making len = 1
                for i in range(len(current_prot)):  # checks if length of current protein is more than 0
                    current_prot[i] += aa  # if length > 0, then it adds aas but if not it skips
        return proteins


    def all_orf_proteins(self, startReadPos=0, endReadPos=0, ordered=False):
        """Compute all possible proteins for all open reading frames"""
        if endReadPos > startReadPos:
            tmp_seq = BioSeq(self.seq[startReadPos: endReadPos], self.seq_type)
            rfs = tmp_seq.gen_reading_frames()
        else:
            rfs = self.gen_reading_frames()

        res = []  # will accumulate all proteins we find
        for rf in rfs:
            proteins = self.reading_frame_proteins(rf)
            for p in proteins:
                res.append(p)


        if ordered:
            return sorted(res, key=len, reverse=True)  # orders chains by length
        return res
