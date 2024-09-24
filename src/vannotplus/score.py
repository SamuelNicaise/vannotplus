from cyvcf2 import cyvcf2


def main_vannotscore(input_vcf: str, output_vcf: str) -> None:
    input_vcf = cyvcf2.VCF(input_vcf)
    input_vcf.add_info_to_header(
        {
            "ID": "vannotscore",
            "Number": 1,
            "Type": "Integer",
            "Description": "VANNOT score, an emulation of VaRank score",
        }
    )

    output_vcf = cyvcf2.Writer(output_vcf, input_vcf)

    for variant in input_vcf:
        variant.INFO["vannotscore"] = get_score(variant)
        output_vcf.write_record(variant)

    output_vcf.close()


def get_variant_info(variant: cyvcf2.Variant, field: str) -> str:
    """
    Do not deal with types until needed for maximum performance
    """
    try:
        return variant.INFO[field]
    except KeyError:
        return ""


def get_score(variant: cyvcf2.Variant) -> int:
    """
    missing:
    S_EssentialSplice : significant effect on splicing -> 90
    S_CloseSplice: Mutation outside of the canonical splice sites (donor site is -3 to +6', acceptor site -12 to +2) -> 70
    S_LSEstrong : Strong local splice effect -> 40
    S_LSEweak: Weak local splice activation -> 35
    S_DeepSplice : Intronic mutation resulting in a significant effect on splicing -> 25
    S_ExonIntron

    #verify that missense is correct

    # Add SIFT & PPH2 bonus for missense only.
    incr score [addBonus $siftPred $siftMed $thePPH2 $phastcons]
    """

    clinvar = get_variant_info(variant, "CLINVAR_clnsig").lower()
    snpeff_annotation = get_variant_info(variant, "snpeff_annotation").lower()
    outcome = get_variant_info(variant, "outcome").lower()

    if "pathogenic" in clinvar:
        score = 110
    elif any([v in outcome for v in ["frameshift", "stop_gained", "stop gained"]]):
        score = 100
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
        score = 80
    elif "missense" in snpeff_annotation:
        score = 50
        score += get_bonus_score(variant)
    elif any([v in outcome for v in ["in-frame", "inframe"]]):
        score = 30
    elif any([v in outcome for v in ["synonymous", "start_retained", "stop_retained"]]):
        score = 10
    # S_ExonIntron / varlocation
    elif any([v in snpeff_annotation for v in ["intron", "exon"]]):
        score = 2
    elif "utr" in snpeff_annotation:
        score = 1
    else:
        score = 0
    return score


def get_bonus_score(variant: cyvcf2.Variant) -> int:
    """
    Polyphen2 prediction based on HumDiv :
    - D (probably damaging HDIV score in [0.9571] or rankscore in [0.558590.91137])
    - P (possibly damaging  HDIV score in [0.4540.956] or rankscore in [0.370430.55681])
    - B (benign  HDIV score in [00.452] or rankscore in [0.030610.36974]).
    """
    bonus = 0

    # TODO: what's the difference between SIFT4G_pred and SIFT_pred?
    # TODO: what is siftMed ?
    sift = get_variant_info(variant, "SIFT_pred")
    if sift == "D":
        bonus += 5

    pph2 = get_variant_info(variant, "Polyphen2_HDIV_pred")
    if pph2 == "D":
        bonus += 5

    # TODO: which phastcons? 4 possibilities
    phastcons = get_variant_info(variant, "phastCons30way_mammalian")
    if phastcons != "" and float(phastcons) > 0.95:
        bonus += 5

    return bonus


if __name__ == "__main__":
    # input_vcf = "/home1/L_PROD/NGS/BAS/HOWARD/data/nicaises/cut.vcf"
    input_vcf = "/home1/L_PROD/NGS/BAS/HOWARD/data/nicaises/KLA2403985.final.vcf"
    output_vcf = "/home1/L_PROD/NGS/BAS/HOWARD/data/nicaises/score/prio.vcf"
    app = ""
    config = {}
    main_vannotscore(input_vcf, output_vcf, app, config)