"""
Manage the conversion of VCF variant calls to gdot format.

@author: letaw
"""
from typing import Tuple
from rinc.util.tx_eff_pysam import PysamTxEff
from rinc.util import chromosome_map

class HGVSNomenclature:

    def __init__(self, chrom: str, pos: int, ref: str, alt: str, pysam_tx: PysamTxEff) -> None:
        self.chrom = chrom
        self.pos = pos  # 1-based VCF position
        self.ref = ref
        self.alt = alt
        self.pysam = pysam_tx

        # Convenience stuff
        self.lref: int = len(self.ref)
        self.lalt: int = len(self.alt)
        self.padding_len: int = self._find_padding()
        self.start_coord: int = self.pos + self.padding_len
        self.raw_ref: str = ref[self.padding_len:]
        self.raw_alt: str = alt[self.padding_len:]
        self.raw_lref: int = len(self.raw_ref)
        self.raw_lalt: int = len(self.raw_alt)
        self.allele_size_gap: int = abs(self.raw_lref - self.raw_lalt)
        self.is_dup = self._check_dup_status()

    @staticmethod
    def _rev_comp(seq: str) -> str:
        """Reverse complements a DNA sequence """
        trans_table = str.maketrans('ATCG', 'TAGC')
        reversed_seq = seq[::-1].upper()
        return reversed_seq.translate(trans_table)

    def get_variant_type(self) -> str:
        """
        Classifies the core variant change (without padding) into a standard type.
        :return: A string representing the variant type (e.g., 'Substitution', 'Deletion', 'Insertion', 'Delins').
        """
        lref: int = self.raw_lref
        lalt: int = self.raw_lalt

        # Substitution - 1 base change (e.g., A > C)
        if lref == 1 and lalt == 1:
            return 'sub'

        # Inversion - at least two bases, reverese complemented (e.g., ATC > GAT)
        if self._rev_comp(self.raw_ref) == self.raw_alt:
            return 'inv'

        # Indel - complex deletions or insertions (e.g., ACC > CTG or C > TTT)
        if lref > 0 and lalt > 0:
            return 'delins'

        # Deletion - loss of sequence (e.g., TAA > T or AA > .)
        if lref > 0 and lalt == 0:
            return 'del'

        # Insertion - gain of sequence (e.g., T > TA or . > A)
        if lref == 0 and lalt > 0:
            if self.is_dup:
                return 'dup'
            return 'ins'

        raise Exception("Unknown variant type: " + self.chrom + ':' + str(self.pos) + self.ref + '>' + self.alt)

    def _check_dup_status(self):
        """
        Check to see if we should classify this as a duplication.  This happens when the sequence being inserted
        is the same as the sequence immediately following it.
        :return:
        """
        local_seq = self.pysam.faidx_query(self.chrom, self.pos)
        chk_start = self.pysam.size + self.padding_len
        chk_end = chk_start + self.raw_lalt
        return local_seq[chk_start:chk_end] == self.raw_alt

    def _find_padding(self) -> int:
        """
        If either the ref or alt allele is larger than length 1, and the first base of each matches,
        this is a padding base. Returns the length of the common prefix.

        :return: The length of the common prefix (padding sequence) as an integer.
        """
        padding: int = 0
        if self.lref > 1 or self.lalt > 1:
            try:
                while self.ref[padding] == self.alt[padding]:
                    padding += 1
            except IndexError:
                pass

        return padding

    def get_hgvs_coord(self) -> str:
        """Abstract method to be implemented by subclasses."""
        raise NotImplementedError("Subclass must implement abstract method 'get_hgvs_coord'")


class Indel(HGVSNomenclature):

    def __init__(self, chrom: str, pos: int, ref: str, alt: str, pysam_tx: PysamTxEff) -> None:
        super().__init__(chrom, pos, ref, alt, pysam_tx)
        self.delins_start = self.start_coord
        self.delins_end = self.delins_start + self.raw_lref - 1

    def get_hgvs_coord(self) -> str:
        # VCF: 10 TAG > AGC
        # HGVS: g.10delinsGCA
        if self.raw_lref == 1:
            return f"{self.delins_start}delins{self.raw_alt}"
        else:
            return f"{self.delins_start}_{self.delins_end}delins{self.raw_alt}"


class Insertion(HGVSNomenclature):

    def __init__(self, chrom: str, pos: int, ref: str, alt: str, pysam_tx: PysamTxEff) -> None:
        super().__init__(chrom, pos, ref, alt, pysam_tx)
        self.local_seq = self.pysam.faidx_query(self.chrom, self.pos)
        self.offset, self.new_ins_seq = self._find_3prime_offset()
        self.ins_start = self.start_coord + self.offset - 1
        self.ins_end = self.ins_start + 1

    def _find_3prime_offset(self) -> Tuple[int, str]:
        """
        Use a sliding window to move 3' down the sequence (including inserted allele) and compare resulting sequences
        after removing a sequence of length equal to the deleted sequence.  If the resulting sequences are identical,
        continue until they are not.  Collect final offset value to be used in adjusting g. position for insertions.
        Collect the new inserted sequence as well in case of realignment.
        :return: offset, or how far can we move the deletion and result in identical reference
        """
        size: int = self.pysam.size
        start_pos: int = size + self.padding_len
        seq_with_ins: str = self.local_seq[:start_pos] + self.raw_alt + self.local_seq[start_pos:]
        comp_seq: str = seq_with_ins[:size + self.padding_len] + seq_with_ins[size + self.raw_lalt + self.padding_len:]
        new_seq: str = seq_with_ins[:size + self.padding_len] + seq_with_ins[size + self.raw_lalt + self.padding_len:]
        offset: int = 0
        while comp_seq == new_seq and size + self.raw_lalt + self.padding_len + offset < (size * 2) + 1:
            offset += 1
            new_seq = (seq_with_ins[:size + self.padding_len + offset] +
                       seq_with_ins[size + self.raw_lalt + self.padding_len + offset:])
        else:
            offset -= 1
            new_ins_start = size + self.padding_len + offset
            new_ins_seq: str = seq_with_ins[new_ins_start:new_ins_start + self.raw_lalt]
        return offset, new_ins_seq

    def get_hgvs_coord(self) -> str:
        # VCF: 10 A > AGCA
        # HGVS: g.11_13insGCA
        return f"{self.ins_start}_{self.ins_end}ins{self.new_ins_seq}"


class Duplication(Insertion):

    def __init__(self, chrom: str, pos: int, ref: str, alt: str, pysam_tx: PysamTxEff) -> None:
        super().__init__(chrom, pos, ref, alt, pysam_tx)
        self.dup_start = self.ins_start - self.raw_lalt + 1
        self.dup_end = self.ins_start
        # self.dup_start = self.ins_start
        # self.dup_end = self.dup_start + self.raw_lalt - 1

    def get_hgvs_coord(self) -> str:
        # VCF: 10 TAG > AGC
        # HGVS: g.10delinsGCA
        if self.allele_size_gap == 1 and self.raw_lalt == 1:
            return f"{self.dup_start}dup"
        else:
            return f"{self.dup_start}_{self.dup_end}dup"


class Inversion(HGVSNomenclature):

    def __init__(self, chrom: str, pos: int, ref: str, alt: str, pysam_tx: PysamTxEff) -> None:
        super().__init__(chrom, pos, ref, alt, pysam_tx)
        self.inv_start: int = self.start_coord
        self.inv_end: int = self.start_coord + self.raw_lref - 1

    def get_hgvs_coord(self) -> str:
        # VCF: 10 ATCG > CGAT
        # HGVS: g.10_13inv
        return f"{self.inv_start}_{self.inv_end}inv"


class Substitution(HGVSNomenclature):

    def __init__(self, chrom: str, pos: int, ref: str, alt: str, pysam_tx: PysamTxEff) -> None:
        super().__init__(chrom, pos, ref, alt, pysam_tx)

    def get_hgvs_coord(self) -> str:
        # VCF: 10 A > C
        # HGVS: g.10A>C
        return f"{self.pos}{self.ref}>{self.alt}"


class Deletion(HGVSNomenclature):

    def __init__(self, chrom: str, pos: int, ref: str, alt: str, pysam_tx: PysamTxEff) -> None:
        super().__init__(chrom, pos, ref, alt, pysam_tx)
        self.local_seq: str = self.pysam.faidx_query(self.chrom, self.pos)
        self.offset_3prime: int = self._find_3prime_offset()
        self.del_start: int = self.start_coord + self.offset_3prime
        self.del_end: int = self.del_start + self.raw_lref - 1

    def _find_3prime_offset(self) -> int:
        """
        Use a sliding window to move 3' down the sequence and compare resulting sequences after removing
        a sequence of length equal to the deleted sequence.  If the resulting sequences are identical, continue
        until they are not.  Collect final offset value to be used in adjusting g. position for deletions.
        :return: offset, or how far can we move the deletion and result in identical reference
        """
        size: int = self.pysam.size
        comp_seq: str = self.local_seq[:size + self.padding_len] + self.local_seq[size + self.raw_lref
                                                                                  +self.padding_len:]
        new_seq: str = self.local_seq[:size + self.padding_len] + self.local_seq[size + self.raw_lref
                                                                                 +self.padding_len:]
        offset: int = 0
        while comp_seq == new_seq and size + self.raw_lref + self.padding_len + offset < (size * 2) + 1:
            offset += 1
            new_seq = (self.local_seq[:size + self.padding_len + offset] +
                       self.local_seq[size + self.raw_lref + self.padding_len + offset:])
        else:
            offset -= 1
        return offset

    def get_hgvs_coord(self) -> str:
        if self.allele_size_gap == 1:
            # Single base deletion (e.g., pos+1 del T)
            return f"{self.del_start}del"
        else:
            # Multi-base deletion (e.g., pos+1_pos+3 del ATG)
            return f"{self.del_start}_{self.del_end}del"

def get_gdot_plus(chrom: str, pos: int, ref: str, alt: str, pysam_tx: PysamTxEff) -> tuple[str,str]:
    """
    Returns g_dot and the variant type
    1) Return the full g.dot include chromsome andg. prefix
    2) Returns the variant type (eg del, delins, etc)
    """
    vrnt_type_map = {'del': Deletion,
                     'delins': Indel,
                     'dup': Duplication,
                     'ins': Insertion,
                     'inv': Inversion,
                     'sub': Substitution
                     }

    vrnt_type = HGVSNomenclature(chrom, pos, ref, alt, pysam_tx).get_variant_type()
    gdot_change = vrnt_type_map[vrnt_type](chrom, pos, ref, alt, pysam_tx).get_hgvs_coord()
    refseq_chromosome = chromosome_map.get_refseq(chrom)
    
    full_g_dot = refseq_chromosome + ':g.' + gdot_change
    
    return full_g_dot, vrnt_type
    
    
    