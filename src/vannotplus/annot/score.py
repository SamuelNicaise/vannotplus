from cyvcf2 import cyvcf2

import numpy as np

from vannotplus.annot.gmc import get_gmc_by_variant, get_gmc_header
from vannotplus.commons import get_variant_id


def main_annot(input_vcf_path: str, output_vcf_path: str, config: dict) -> None:
    input_vcf = cyvcf2.VCF(input_vcf_path)
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
        variant.INFO["vannotscore"] = get_score(variant)
        variant_id = get_variant_id(variant)
        try:
            variant.set_format("GMC", variant_gmc_dic[variant_id])
        except KeyError:
            # variant is not in a gene
            pass
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

    # Add SIFT & PPH2 bonus for missense only.
    incr score [addBonus $siftPred $siftMed $thePPH2 $phastcons]
    """

    clinvar = get_variant_info(variant, "CLINVAR_clnsig").lower()
    snpeff_annotation = get_variant_info(variant, "snpeff_annotation").lower()
    outcome = get_variant_info(variant, "outcome").lower()

    if "pathogenic" in clinvar:  # covers both Pathogenic and Probably pathogenic
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

    score = get_splicing_score(variant, score)

    return score


def get_splicing_score(variant: cyvcf2.Variant, score: int) -> int:
    score = get_spliceai_score(variant, score)
    # spip: do the same as spliceai?
    return score


def get_spliceai_score(variant: cyvcf2.Variant, score: int) -> int:
    SPLICEAI_THRESHOLD = 0.5
    S_ESSENTIALSPLICE = 90
    S_CLOSESPLICE = 70
    S_DEEPSPLICE = 25

    spliceai_symbols = get_variant_info(variant, "SpliceAI_SYMBOL").split(",")
    index = spliceai_symbols[0]
    # TODO: change when HOWARD's transcript prioritization is done
    # tnomen = get_variant_info(variant, "TNOMEN")
    # index =  spliceai_symbols.index(tnomen)
    spliceai_scores_fields = [
        "SpliceAI_DS_AG",
        "SpliceAI_DS_AL",
        "SpliceAI_DS_DG",
        "SpliceAI_DS_DL",
    ]
    spliceai_pos_fields = [
        "SpliceAI_DP_AG",
        "SpliceAI_DP_AL",
        "SpliceAI_DP_DG",
        "SpliceAI_DP_DL",
    ]
    splice_ai_scores = [
        float(get_variant_info(variant, f)) for f in spliceai_scores_fields
    ]
    splice_ai_pos = [int(get_variant_info(variant, f)) for f in spliceai_pos_fields]
    max_score = max(splice_ai_scores)
    index = splice_ai_scores.index(max_score)
    dist = splice_ai_pos[index]

    if max_score > SPLICEAI_THRESHOLD:
        # do the pos thing
        if score < S_ESSENTIALSPLICE:
            if dist in (1, 2):
                score = S_ESSENTIALSPLICE  # + phastcons
        elif score < S_CLOSESPLICE:
            if (index in (2, 3) and -3 < dist < 6) or (
                index in (0, 1) and -12 < dist < 2
            ):
                # index in (0,1) = acceptor = 3'
                # index in (2,3) = donor = 5'
                score = S_CLOSESPLICE  # + phastcons
        elif score < S_DEEPSPLICE:
            snpeff_annotation = get_variant_info(variant, "snpeff_annotation").lower()
            if "intron" in snpeff_annotation:
                score = S_DEEPSPLICE  # +phastcons
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
    config = {}
    main_annot(input_vcf, output_vcf, config)
