from cyvcf2 import cyvcf2
import logging as log
import numpy as np


from vannotplus.commons import get_variant_id


def get_gmc_header(gene_field: str, do_filtered_gmc) -> list[dict[str, str | int]]:
    """
    To be used with cyvcf2.VCF.set_format()
    """
    headers = [
        {
            "ID": "GMC",
            "Number": 1,
            "Type": "Integer",
            "Description": f"Gene Mutations Count, i.e. how many variants were called in the current gene (where gene is defined by the {gene_field} field)",
        }
    ]
    if do_filtered_gmc:
        headers.append(
            {
                "ID": "GMC_FILTERED",
                "Number": 1,
                "Type": "Integer",
                "Description": f"Filtered Gene Mutations Count, i.e. number of variants in the current gene (defined by the {gene_field} field) that pass Cutevariant's AR htz filter (from DIAG_v3.yml). The intended use of this column is to add GMC_FILTERED >=2 to said AR htz filter, so that only genes with at least 2 variants passing the filter are visible.",
            }
        )
    return headers


def genotypes_to_counts(genotypes: np.ndarray) -> np.ndarray:
    """
    cyvcf2.Variant.gt_types returns a numpy array indicating the type of each sample.
    HOM_REF=0, HET=1. For gts012=True HOM_ALT=2, UNKNOWN=3
    """
    # return [1 if 1 <= v <= 2 else 0 for v in genotypes]
    # return np.where(1 <= genotypes <= 2, 1, 0)
    mask = (genotypes >= 1) & (genotypes <= 2)
    return mask.astype(np.int32)


def variant_to_filtered_counts(variant: cyvcf2.Variant, default_empty_array: np.ndarray, gmc_config: dict, eps: float = 1e-8) -> np.ndarray:
    """
    default_empty_array is an array of zeros with length equal to the number of samples in the VCF
    It is computed once and passed as an arg to avoid recomputing it for each variant
    
    If filters don't pass, return default_empty_array
    Otherwise, return an array of 1 or 0 for each sample depending on if filter passed for each sample

    For maximum computing speed, checks that fails the most should be done first to avoid further checks.

    There is one last filtering step that is not done in this func for performance optimization: unfiltered GMC should be >= 1, but it is not computed yet when this func is executed.
    This last filter is implemented in filter_gmc_by_gmc(), called by get_gmc_by_variant() after this func.

    Note: this relies on numpy float comparisons. Due to floating point precision issues,  a variant with allele frequency exactly equal to the threshold might be considered as above or below the threshold (e.g. 0.01 might be stored as 0.009999999 or 0.010000001).
    To avoid this, a small epsilon (default: 1e-8) is added to float thresholds conservatively.
    """

    for pop_freq in gmc_config["pop_freq_fields"]:
        try:
            if variant.INFO[pop_freq] > gmc_config["pop_freq_threshold"] + eps:
                return default_empty_array
        except KeyError:
            # field not present = variant passes filter
            # print("KeyError for pop_freq field:", variant.CHROM, variant.POS, variant.REF, variant.ALT)
            pass

    for pop_homcount_field in gmc_config["pop_homcount_fields"]:
        try:
            if variant.INFO[pop_homcount_field] > gmc_config["pop_homcount_threshold"]:
                return default_empty_array
        except KeyError:
            pass

    # OMIM ID can't be null
    try:
        variant.INFO[gmc_config["omim_id_field"]]
    except KeyError:
        # /!\ opposite of pop filters: absence of omim_id means fail
        return default_empty_array

    try:
        if "AR" not in variant.INFO[gmc_config["omim_inheritance_field"]]:
            return default_empty_array
    except KeyError:
        pass

    for inner_freq in gmc_config["allelefreq_fields"]:
        try:
            if variant.INFO[inner_freq] > gmc_config["allelefreq_threshold"] + eps:
                return default_empty_array
        except KeyError:
            pass

    for inner_count in gmc_config["homcount_fields"]:
        try:
            if variant.INFO[inner_count] > gmc_config["homcount_threshold"]:
                return default_empty_array
        except KeyError:
            pass

    #last ones are more technical : they are sample based
    result_array = np.ones(len(variant.gt_types), dtype=bool)

    # genotype filter
    if gmc_config["gt"] == 1:
        # keep only heterozygous variants
        result_array &= (variant.gt_types == 1)
    else:
        # default method: keep both HET and HOM_ALT
        result_array = (variant.gt_types >= 1) & (variant.gt_types <= 2)

    # VAF filter
    try:
        vaf_array = variant.format("VAF")
    except KeyError:
        vaf_array = None
    vaf_threshold = gmc_config.get("vaf_threshold", None)

    if vaf_array is not None and vaf_threshold is not None:
        # this manipulation is required because cyvcf2 returns a 2D array even for single ALT alleles
        if vaf_array.ndim == 2 and vaf_array.shape[1] == 1:
            vaf_array = vaf_array[:, 0]
        else:
            raise ValueError(f"Unexpected VAF array shape: {vaf_array.shape} for variant {variant.CHROM} {variant.POS} {variant.REF} {variant.ALT}\\nThis code assumes no multiallelic records. You should be able to fix this with bcftools norm -m- <input.vcf>")
        result_array &= ~np.isnan(vaf_array) & (vaf_array >= vaf_threshold)

    log.debug(f"Variant passed all filters: {variant.CHROM} {variant.POS} {variant.REF} {variant.ALT} {vaf_array}")
    return result_array.astype(np.int32)

def filter_gmc_by_gmc(gmc: np.ndarray, filtered_gmc: np.ndarray) -> np.ndarray:
    """
    Filter GMC on itself: if GMC < 2, then filtered GMC should be 0. This is because genes containing only one variant should not appear in the AR hom filter.

    Example:
    >>> filter_gmc_by_gmc(np.array([0, 1, 4, 5]), np.array([0, 1, 2, 3]))
    array([0, 0, 2, 3])
    """
    return np.where(gmc < 2, 0, filtered_gmc)

def get_gmc_by_variant(
    vcf_path: str, gmc_config: dict, do_filtered_gmc: bool = False
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    """
    Returns two dictionaries:
    - variant_gene_dict: keys are variant IDs, values are numpy arrays of GMC for each sample
    - filtered_variant_gene_dict: keys are variant IDs, values are numpy arrays of filtered GMC for each sample

    Steps:
    1) build a gene to GMC mapping
        If do_filtered_gmc is True, apply the same logic but only for variants passing the filter
    2) build a variant to GMC mapping using the gene to GMC mapping

    If do_filtered_gmc is True, a second dictionary is built along the way, counting only variants passing the filter


    If a variant is not in a gene, it is not included in the output dictionaries
    If a gene field contains multiple genes (e.g. "GENE1/GENE2"), raise NotImplementedError
    GMC is computed as the sum of genotypes (HET=1, HOM_ALT=2) for each sample
    The filtered GMC is computed similarly but only for variants passing the filter
    """
    # gts012=True is extremely important for genotypes_to_counts
    vcf = cyvcf2.VCF(vcf_path, gts012=True)
    variant_gene_dict = {}
    filtered_variant_gene_dict = {}
    gene_gmc_dict = {}
    gene_filtered_gmc_dict = {}
    # We keep np.int32 type to be able to set its value to "." (minimal value of np.int32 in cyvcf2)
    default_filtered_gmc = np.zeros(len(vcf.samples), dtype=np.int32)

    log.debug(f"do_filtered_gmc: {do_filtered_gmc}")
    for variant in vcf:
        key = get_variant_id(variant)
        try:
            gene = variant.INFO[gmc_config["gene_field"]]
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

        # if gene == "DVL1":
        #     print(variant)
        #     print(variant.CHROM, variant.POS, variant.REF, variant.ALT, variant.gt_types)
        #     print(gene_gmc_dict[gene])

        # same for filtered GMC
        if do_filtered_gmc:
                if gene in gene_filtered_gmc_dict:
                    gene_filtered_gmc_dict[gene] = np.add(
                        gene_filtered_gmc_dict[gene],
                        variant_to_filtered_counts(variant, default_filtered_gmc, gmc_config),
                    )
                else:
                    gene_filtered_gmc_dict[gene] = variant_to_filtered_counts(variant, default_filtered_gmc, gmc_config)

    # overwrite genes with GMC in variant_gene_dict to save RAM
    for variant, gene in variant_gene_dict.items():
        if "/" in gene:
            raise NotImplementedError(
                f"gene field '{gmc_config['gene_field']}' needs to contain only one gene. Got '{gene}' for variant {variant}"
            )
        variant_gene_dict[variant] = gene_gmc_dict[gene]

        if do_filtered_gmc:
            # Unlike gmc that is always null (if variant isn't in a gene) or >= 1 (if variant is in a gene),
            # filtered_gmc can be 0 if the variant is in a gene but no variant in that gene passed the filter

            #there is one ultimate step that can't be done before gmc is computed: filter filtered_gmc on gmc itself
            final_filtered_gmc = filter_gmc_by_gmc(gene_gmc_dict[gene], gene_filtered_gmc_dict[gene])

            filtered_variant_gene_dict[variant] = final_filtered_gmc

    return variant_gene_dict, filtered_variant_gene_dict
