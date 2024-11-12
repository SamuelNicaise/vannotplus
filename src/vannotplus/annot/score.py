from cyvcf2 import cyvcf2

from vannotplus.annot.gmc import get_gmc_by_variant, get_gmc_header
from vannotplus.annot.splicing import get_splicing_score
from vannotplus.commons import get_variant_id, get_variant_info


def main_annot(
    input_vcf_path: str,
    output_vcf_path: str,
    config: dict,
    do_vannotscore: bool = False,
) -> None:
    """
    VANNOT score has been replaced by PZTScore_transcript computed by howard
    if do_vannotscore == True: compute it and add it to vcf
    By default it is False and VANNOT score is not computed

    This code is still required in VANNOT to add GMC to the final VCF
    """
    input_vcf = cyvcf2.VCF(input_vcf_path)

    if do_vannotscore:
        input_vcf.add_info_to_header(
            {
                "ID": "vannotscore",
                "Number": 1,
                "Type": "Integer",
                "Description": "VANNOT score, an emulation of VaRank score",
            }
        )

    input_vcf.add_format_to_header(get_gmc_header(config["gene_field"]))
    variant_gmc_dic = get_gmc_by_variant(input_vcf_path, config["gene_field"])

    output_vcf = cyvcf2.Writer(output_vcf_path, input_vcf)

    for variant in input_vcf:
        # vannotscore
        if do_vannotscore:
            variant.INFO["vannotscore"] = get_score(variant, config)

        # GMC
        variant_id = get_variant_id(variant)
        try:
            variant.set_format("GMC", variant_gmc_dic[variant_id])
        except KeyError:
            # variant is not in a gene
            pass

        output_vcf.write_record(variant)

    output_vcf.close()


def get_score(variant: cyvcf2.Variant, config: dict) -> int:
    """
    missing:
    S_EssentialSplice : significant effect on splicing -> 90
    S_CloseSplice: Mutation outside of the canonical splice sites (donor site is -3 to +6', acceptor site -12 to +2) -> 70
    S_LSEstrong : Strong local splice effect -> 40
    S_LSEweak: Weak local splice activation -> 35
    S_DeepSplice : Intronic mutation resulting in a significant effect on splicing -> 25

    # Add SIFT & PPH2 bonus for missense only.
    incr score [addBonus $siftPred $siftMed $thePPH2 $phastcons]
    """
    score_config = config["score_config"]

    clinvar = get_variant_info(variant, "CLINVAR_clnsig").lower()
    snpeff_annotation = get_variant_info(variant, "snpeff_annotation").lower()
    outcome = get_variant_info(variant, "outcome").lower()

    if "pathogenic" in clinvar:  # covers both Pathogenic and Probably pathogenic
        score = score_config["S_Known"]
    elif any([v in outcome for v in ["frameshift", "stop_gained", "stop gained"]]):
        score = score_config["S_StopGain"]
        score += get_phastcons_bonus(variant, score_config)
    elif any(
        [
            v in outcome
            for v in [
                "startloss",
                "start_lost",
                "start loss",
                "stoploss",
                "stop_lost",
                "stop loss",
            ]
        ]
    ):
        score = score_config["S_StartStopLoss"]
        score += get_phastcons_bonus(variant, score_config)
    elif "missense" in snpeff_annotation:
        score = score_config["S_Missense"]
        score += get_bonus_score(variant, score_config)
    elif any([v in outcome for v in ["in-frame", "inframe"]]):
        score = score_config["S_Inframe"]
    elif any([v in outcome for v in ["synonymous", "start_retained", "stop_retained"]]):
        score = score_config["S_Synonymous"]
    # S_ExonIntron / varlocation
    elif any([v in snpeff_annotation for v in ["intron", "exon"]]):
        score = score_config["S_ExonIntron"]
        score += get_phastcons_bonus(variant, score_config)
    elif "utr" in snpeff_annotation:
        score = score_config["S_UTR"]
    else:
        score = 0

    score = get_splicing_score(variant, score, score_config)

    return score


def get_bonus_score(variant: cyvcf2.Variant, score_config) -> int:
    """
    Polyphen2 prediction based on HumDiv :
    - D (probably damaging HDIV score in [0.9571] or rankscore in [0.558590.91137])
    - P (possibly damaging  HDIV score in [0.4540.956] or rankscore in [0.370430.55681])
    - B (benign  HDIV score in [00.452] or rankscore in [0.030610.36974]).

    SIFT : same system

    Phastcons: based on a score threshold

    Add: alphamissense
    Add: provean
    Add: FATHMM # ignore
    Add: CADD
    Add: mistic
    """
    bonus = 0

    # What's the difference between SIFT4G_pred and SIFT_pred?
    # SIFT4G is more recent. We keep only SIFT4G by Jean's decision
    sift = get_variant_info(variant, "SIFT_pred")
    if sift == "D":
        bonus += 5

    pph2 = get_variant_info(variant, "Polyphen2_HDIV_pred")
    if pph2 == "D":
        bonus += 5

    # Which phastcons? Only one is left, use that one
    bonus += get_phastcons_bonus(variant, score_config)

    return bonus


def get_phastcons_bonus(variant: cyvcf2.Variant, score_config: dict):
    phastcons = get_variant_info(variant, "phastCons100way")
    if phastcons != "" and float(phastcons) > score_config["Threshold_Phastcons"]:
        return score_config["B_phastCons"]
    return 0


if __name__ == "__main__":
    # input_vcf = "/home1/L_PROD/NGS/BAS/HOWARD/data/nicaises/cut.vcf"
    input_vcf = "/home1/L_PROD/NGS/BAS/HOWARD/data/nicaises/KLA2403985.final.vcf"
    output_vcf = "/home1/L_PROD/NGS/BAS/HOWARD/data/nicaises/score/prio.vcf"
    config = {}
    main_annot(input_vcf, output_vcf, config)
