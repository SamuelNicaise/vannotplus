from cyvcf2 import cyvcf2
import numpy as np

from vannotplus.commons import get_variant_id


def get_gmc_header(gene_field: str) -> dict[str, str | int]:
    """
    To be used with cyvcf2.VCF.set_format()
    """
    return {
        "ID": "GMC",
        "Number": 1,
        "Type": "Integer",
        "Description": f"Gene Mutations Count, i.e. how many variants were called in the current gene (where gene is defined by the {gene_field} field)",
    }


def genotypes_to_counts(genotypes: np.ndarray) -> list[int]:
    """
    cyvcf2.Variant.gt_types returns a numpy array indicating the type of each sample.
    HOM_REF=0, HET=1. For gts012=True HOM_ALT=2, UNKNOWN=3
    """
    # return [1 if 1 <= v <= 2 else 0 for v in genotypes]
    return np.where(1 <= genotypes <= 2, 1, 0)


def get_gmc_by_variant(vcf_path: str, gene_field: str) -> dict[str, np.ndarray]:
    # gts012=True is extremely important for genotypes_to_counts
    vcf = cyvcf2.VCF(vcf_path, gts012=True)
    variant_gene_dict = {}
    gene_gmc_dict = {}

    for variant in vcf:
        key = get_variant_id(variant)
        try:
            gene = variant.INFO[gene_field]
        except KeyError:
            # variant is not in a gene
            continue
        variant_gene_dict[key] = gene

        if gene in gene_gmc_dict:
            gene_gmc_dict[gene] = np.add(
                gene_gmc_dict[gene], genotypes_to_counts(variant.gt_types)
            )
        else:
            gene_gmc_dict[gene] = genotypes_to_counts(variant.gt_types)

    # overwrite genes with GMC in variant_gene_dict to save RAM
    for variant, gene in variant_gene_dict.items():
        if "/" in gene:
            raise NotImplementedError(f"gene field '{gene_field}' needs to contain only one gene. Got '{gene}' for variant {variant}")
        variant_gene_dict[variant] = gene_gmc_dict[gene]

    return variant_gene_dict
