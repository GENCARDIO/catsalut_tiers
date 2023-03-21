# CatSalutSomaticTiers
This Python script provides a class `CatSalutSomaticTiers` to classify clinically actionable variants supported by CatSalut 
The official instruction is available at: [Perfil genetic tumors sòlids 2022](https://scientiasalut.gencat.cat/bitstream/handle/11351/8424.2/determinacions_perfil_genetic_tumors_solids_adult_2022.pdf?sequence=5&isAllowed=y)

## Installation
No installation is required. Simply download the `cat_salut_somatic_tiers.py` file and import it into your Python project.

## Usage
To use the `CatSalutSomaticTiers` class, you need to provide a TSV (tab-separated values) file containing the list of known tier variants.

There's a tsv file with variant-specific rules file at `test/tiers/Tiers_cancer.v1.tsv`

Here's an example of how to use the class:

```python
from cat_salut_somatic_tiers import CatSalutSomaticTiers

variants_tsv = "path/to/your/variants.tsv"
cat_salut = CatSalutSomaticTiers(variants_tsv)

gene = "BRAF"
variant_type = "SNV"
exon = "15"
csq = "missense_variant"
hgvsp = "p.V600E"

tier = cat_salut.classify(gene, variant_type, exon=exon, csq=csq, hgvsp=hgvsp)
print(f"Tier classification: {tier}")
Each row should represent a known tier variant with the corresponding information.
```
## Class Methods
The main methods provided by the CatSalutSomaticTiers class are:
    `__init__(self, variants_tsv: str)`: Initializes the class with the given TSV file containing known tier variants.
    `load_known_tier_variants(self)`: Reads the TSV file and stores the tier variants in a dictionary.
    `classify(self, gene: str, variant_type: str, exon: Optional[str] = None, intron: Optional[str] = None, csq: Optional[str] = None, hgvsp: Optional[str] = None, force_gene: bool = False)`: Classifies a given variant based on the gene, variant type, and other optional parameters. Returns the tier classification.

## Exceptions

The module defines the following custom exceptions:
    `InvalidTsvSchema`: Raised when the provided TSV file has an incorrect schema.
    `MissingGene`: Raised when the provided gene is not found in the known tier variants.
    `InvalidAlteration`: Raised when the provided variant type is not valid.

## TSV File Format
The input TSV file should have the following header:
`Gene, Exon, Intron, Alteration, Alteration Comments, Consequence, HGVSp, Exclusion Criteria, ESCAT, Treatments, Treatment Lines, Comments, Clinical Trials, Tier, Automatized, Skip, Version, Date`

Each row should represent a known tier variant with the corresponding information.

## Tests

```
python -m pytest test/ 
```

## Limitations
Currently it is only supported the analysis of LUAD.
In addition, the current version does not support the automatic annotation of all available mutation classes:
- EGFR uncommon mutations.
- ALK mutations of resistance.
- RET mutations of resistance.
- NTRK mutations of resistance.

