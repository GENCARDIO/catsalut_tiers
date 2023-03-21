import csv
import os
import argparse
from collections import defaultdict
from typing import Optional
import os

class InvalidTsvSchema(Exception):
    pass


class MissingGene(Exception):
    pass


class InvalidAlteration(Exception):
    pass


class CatSalutSomaticTiers:
    """ """

    """
    Class to read and classify clinically actionable variants supported by CatSalut
    """

    def __init__(self, variants_tsv: str):
        self._variants_tsv = variants_tsv

        if not os.path.isfile(self._variants_tsv):
            msg = f"Missing input tsv file {self._variants_tsv}"
            raise FileNotFoundError(msg)

        self._tier_variants = self.load_known_tier_variants()

    def load_known_tier_variants(self) -> list():
        """
        Read clinically accionable variants supported by CatSalut
        """
        header = [
            "Gene",
            "Exon",
            "Intron",
            "Alteration",
            "Alteration Comments",
            "Consequence",
            "HGVSp",
            "Exclusion Criteria",
            "ESCAT",
            "Treatments",
            "Treatment Lines",
            "Comments",
            "Clinical Trials",
            "Tier",
            "Automatized",
            "Skip",
            "Version",
            "Date",
        ]

        tier_variants = defaultdict(dict)
        with open(self._variants_tsv, mode="r") as tsv_file:
            tsv_reader = csv.DictReader(tsv_file, delimiter="\t")
            for idx, row in enumerate(tsv_reader):
                if idx == 0:
                    for key in row:
                        if key not in header:
                            msg = f"{key} not found between:{','.join(header)}"
                            raise InvalidTsvSchema(msg)
                if not row["Gene"] in tier_variants:
                    tier_variants[row["Gene"]] = []
                tier_variants[row["Gene"]].append(row)
        tsv_file.close()
        return tier_variants

    def _normalize_input(self, value: Optional[str]) -> Optional[str]:
        """
        Normalizes the input value by removing '.' and 'p.'
        """
        if value:
            if value == ".":
                return ""
            else:
                if type(value) == str:
                    if value.startswith("p."):
                        value = value.replace("p.", "")
        return value

    def classify(
        self,
        gene: str,
        variant_type: str,
        exon: Optional[str] = None,
        intron: Optional[str] = None,
        csq: Optional[str] = None,
        hgvsp: Optional[str] = None,
        force_gene: bool = False,
    ) -> str:
        """
        Note that this method does not check that input genes
        and variant codes are correct since it is assumed that they
        come from the annotation process
        """

        exon = self._normalize_input(exon)
        intron = self._normalize_input(intron)
        csq = self._normalize_input(csq)
        hgvsp = self._normalize_input(hgvsp)

        valid_alterations = [
            "SNV",
            "MNV",
            "Insertion",
            "Deletion",
            "Fusion",
            "Amplification",
            "Loss",
            "Mutation",
            "SV",
        ]
        if variant_type:
            if variant_type not in valid_alterations:
                msg = f"Invalid {variant_type}. \
                    Accepted variant types are {','.join(valid_alterations)}"
                raise InvalidAlteration(msg)

        tier = None
        if gene not in self._tier_variants:
            msg = f"Gene {gene} not found on {self._variants_tsv}"
            if force_gene:
                raise MissingGene(msg)
            else:
                return tier

        for item in self._tier_variants[gene]:

            if item["Automatized"] == "NO":
                continue

            if item["Skip"] == "YES":
                continue

            is_ok = True
            if item["Alteration"] != "Mutation":
                if variant_type.lower() not in item["Alteration"].lower():
                    is_ok = False

            if intron:
                valid_introns = item["Intron"].split(",")
                for elem in valid_introns:
                    if elem != item["Intron"]:
                        is_ok = False
                    else:
                        is_ok = True

            # Now check exon (if applies)
            if exon:
                if exon != item["Exon"]:
                    is_ok = False
            else:
                if item["Exon"]:
                    is_ok = False
            if csq:
                if item["Consequence"]:
                    csq_list = csq.split(",")
                    if not csq_list:
                        csq_list = csq.split("&")
                    for csq_item in csq_list:
                        if csq_item.lower() not in item["Consequence"].lower():
                            is_ok = False

            # check protein change 
            if hgvsp:
                # Assuming that comes formatted in VEP style
                if hgvsp.startswith("p."):
                    hgvsp = hgvsp.replace("p.", "")
                if item["Exclusion Criteria"]:
                    if hgvsp == item["Exclusion Criteria"]:
                        is_ok = False
                    if item["HGVSp"] not in hgvsp:
                        is_ok = False
                else:
                    if is_ok:
                        if item["HGVSp"] == "":
                            pass
                        else:
                            if item["HGVSp"] in hgvsp:
                                is_ok = True
                            else:
                                is_ok = False
                    else:
                        if hgvsp != item["HGVSp"]:
                            is_ok = False
            if is_ok is True:
                tier = item["Tier"]
                break
        return tier



def main():
    """ """
    parser = argparse.ArgumentParser(description='CatSalut tier assignment')

    parser.add_argument( '--gene', type=str, help='Gene in HUGO nomenclature', required=True)
    parser.add_argument( '--variant_type', type=str, help='Variant type:', required=True)
    parser.add_argument( '--exon', type=str, help='Exon number', required=False)
    parser.add_argument( '--intron', type=str, help='Intron Number', required=False)
    parser.add_argument( '--csq', type=str, help='VEP consequence field', required=False)
    parser.add_argument( '--hgvsp', type=str, help='VEP HGVSp field', required=False)
    parser.add_argument( '--force_gene', type=bool, help='Force input gene', required=False)
    args = parser.parse_args()
    
    print(args.gene, args.variant_type, args.exon, args.intron)
    curr_dir = os.path.dirname(os.path.realpath(__file__))

    variants_tsv = os.path.join(curr_dir, "test/tiers/Tiers_cancer.v1.tsv")

    cs = CatSalutSomaticTiers(variants_tsv)
    tier = cs.classify(gene=args.gene, variant_type=args.variant_type, exon=args.exon)
    print(tier)


if __name__ == "__main__":
    main()
