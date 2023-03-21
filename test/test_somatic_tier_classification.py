import argparse
import os
import shutil
import sys
import unittest

sys.path.append(os.getcwd())
from somatic_tiers import CatSalutSomaticTiers, InvalidAlteration, MissingGene


# https://stackoverflow.com/questions/1305532/convert-nested-python-dict-to-object
class Struct:
    def __init__(self, **entries):
        self.__dict__.update(entries)


class SomaticTiersTester(unittest.TestCase):
    """ """

    @classmethod
    def setUpClass(cls) -> None:
        """ """
        variants_tsv = os.path.join(
            os.getcwd(), "test", "tiers", "Tiers_cancer.v1.tsv"
        )
        print(variants_tsv)
        cls._somatic_tier = CatSalutSomaticTiers(variants_tsv)

    def test_somatic_kras_g13_1(self) -> None:
        """
        KRAS codon 13 mutations
        """
        tier = self._somatic_tier.classify(
            gene="KRAS",
            variant_type="SNV",
            hgvsp="G13V",
            exon="2",
            csq="missense_variant",
        )
        self.assertTrue(tier)

    def test_somatic_kras_g13_1(self) -> None:
        """
        KRAS codon 13 mutations
        """
        tier = self._somatic_tier.classify(
            gene="KRAS",
            variant_type="SNV",
            hgvsp="p.G13C",
            exon="2",
            csq="missense_variant",
        )
        self.assertTrue(tier)

    def test_somatic_tier_egfr_1(self) -> None:
        """
        Exon 21 L858R
        """
        tier = self._somatic_tier.classify(
            gene="EGFR",
            variant_type="SNV",
            hgvsp="L858R",
            exon="21",
            csq="Missense_variant",
        )
        self.assertTrue(tier)

    def test_somatic_tier_egfr_2(self) -> None:
        """
        Exon 19 deletions
        """
        tier = self._somatic_tier.classify(
            gene="EGFR", variant_type="Deletion", exon="19", csq="Inframe_deletion"
        )
        self.assertTrue(tier)

    def test_somatic_tier_egfr_3(self) -> None:
        """
        Exon 20 insertions
        """
        tier = self._somatic_tier.classify(
            gene="EGFR",
            variant_type="Insertion",
            hgvsp="p.D770_N771insG",
            exon="20",
            csq="inframe_insertion",
        )
        self.assertTrue(tier)

    def test_somatic_tier_met_amplification(self) -> None:
        """
        Exon 20 insertions
        """
        tier = self._somatic_tier.classify(
            gene="ERBB2", variant_type="Amplification", hgvsp=".", exon=".", csq="."
        )
        self.assertTrue(tier)

    def test_somatic_braf_v600e(self) -> None:
        """
        BRAF V600E
        """
        tier = self._somatic_tier.classify(
            gene="BRAF",
            variant_type="SNV",
            hgvsp="V600E",
            exon="15",
            csq="missense_variant",
        )
        self.assertTrue(tier)

    def test_somatic_met_ex14_skipping1(self) -> None:
        """
        MET exon 14 skipping mutations (exon 14)
        """
        tier = self._somatic_tier.classify(
            gene="MET", variant_type="Mutation", exon="14", csq="splice_region_variant"
        )
        self.assertTrue(tier)

    def test_somatic_met_ex14_skipping2(self) -> None:
        """
        MET exon 14 skipping mutations (intron 13)
        """
        tier = self._somatic_tier.classify(
            gene="MET", variant_type="Mutation", intron="13", csq="intron_variant"
        )
        self.assertTrue(tier)

    def test_somatic_met_ex14_skipping3(self) -> None:
        """
        MET exon 14 skipping mutations (intron 15)
        """
        tier = self._somatic_tier.classify(
            gene="MET", variant_type="Mutation", intron="15", csq="intron_variant"
        )
        self.assertTrue(tier)

    def test_somatic_met_ex14_skipping4(self) -> None:
        """
        MET exon 14 skipping mutations (intron 13)
        """
        tier = self._somatic_tier.classify(
            gene="MET",
            variant_type="Deletion",
            intron="13",
            csq="splice_acceptor_variant&coding_sequence_variant&splice_polypyrimidine_tract_variant",
        )
        self.assertTrue(tier)

    def test_somatic_met_ex14_skipping5(self) -> None:
        """
        MET exon 14 skipping mutations (intron 15)
        """
        tier = self._somatic_tier.classify(
            gene="MET",
            variant_type="Deletion",
            exon="14",
            csq="splice_acceptor_variant,coding_sequence_variant,splice_polypyrimidine_tract_variant",
        )
        print(tier)
        self.assertTrue(tier)

    def test_somatic_kras_g12c(self) -> None:
        """
        KRAS G12C
        """
        tier = self._somatic_tier.classify(
            gene="KRAS",
            variant_type="SNV",
            hgvsp="G12C",
            exon="2",
            csq="missense_variant",
        )
        self.assertTrue(tier)

    def test_somatic_kras_other_g12c(self) -> None:
        """
        KRAS Other codon 12 mutations
        """
        hgvsp = "G12V"
        tier = self._somatic_tier.classify(
            gene="KRAS",
            variant_type="SNV",
            hgvsp=hgvsp,
            exon="2",
            csq="missense_variant",
        )
        self.assertTrue(tier)

    def test_somatic_erbb2_exon20_mutations(self) -> None:
        """
        ERBB2 exon 20 mutations
        """
        tier = self._somatic_tier.classify(
            gene="ERBB2",
            variant_type="Mutation",
            exon="20",
            csq="splice_region_variant",
        )
        self.assertTrue(tier)

    def test_somatic_tier_fusion_ret(self) -> None:
        """
        Fusion RET
        """
        tier = self._somatic_tier.classify(gene="RET", variant_type="Fusion")
        self.assertTrue(tier)

    def test_somatic_tier_fusion_alk(self) -> None:
        """
        Fusion ALK
        """
        tier = self._somatic_tier.classify(gene="ALK", variant_type="Fusion")
        self.assertTrue(tier)

    # Negative tests
    def test_non_actionable_gene(self) -> None:
        """
        Test that a non-actionable gene raises a custom exception.
        Force error with force_gene=True
        """
        with self.assertRaises(MissingGene):
            self._somatic_tier.classify(
                gene="SCN5A", variant_type="Mutation", force_gene=True
            )

    def test_invalid_alteration_definition(self) -> None:
        """
        Test that an non-supported alteration definition
        """
        with self.assertRaises(InvalidAlteration):
            self._somatic_tier.classify(gene="EGFR", variant_type="INDEL")

    def test_invalid_met_ex14_1(self) -> None:
        """
        Test a METex14 non actionable mutation
        """
        tier = self._somatic_tier.classify(gene="MET", variant_type="Mutation", exon=12)
        self.assertFalse(tier)


if __name__ == "__main__":
    unittest.main()
