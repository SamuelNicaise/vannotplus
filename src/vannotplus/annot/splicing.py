from cyvcf2 import cyvcf2

from vannotplus.commons import get_variant_info


def get_splicing_score(variant: cyvcf2.Variant, score: int, score_config: dict) -> int:
    """
    Returns the maximum possible score from either spip or spliceAI
    If score is already higher than the maximum possible with splicing, return it directly.
    """

    # do not bother computing score if already > max
    if score >= score_config["S_EssentialSplice"]:
        return score
    else:
        score = get_spip_score(variant, score, score_config)

    # same if SPIP already got max score
    if score >= score_config["S_EssentialSplice"]:
        return score
    else:
        score = get_spliceai_score(variant, score, score_config)

    return score


def get_info_from_tuple(
    variant: cyvcf2.Variant, info_field: str, tnomen_index: int
) -> str:
    """
    Many fields have one value per transcript
    Only return the one that corresponds to the prioritized transcript
    """
    info = get_variant_info(variant, info_field)
    if isinstance(info, tuple):
        return info[tnomen_index]
    elif isinstance(info, str):
        # str lists are not loaded as tuples
        return info.split(",")[tnomen_index]
    return info


def get_spip_score(variant: cyvcf2.Variant, score: int, score_config: dict) -> int:
    tnomen_index = 0
    # TODO: change when HOWARD's transcript prioritization is done
    # spip_transcripts = get_variant_info(variant, "SPiP_transcript").split(",")
    # tnomen = get_variant_info(variant, "TNOMEN")
    # tnomen_index =  spip_transcripts.index(tnomen)

    interpretation = get_info_from_tuple(variant, "SPiP_Interpretation", tnomen_index)

    # TODO: specify a list of effect? A threshold?
    if interpretation not in ("", ".", "NTR"):
        dist = int(get_info_from_tuple(variant, "SPiP_DistSS", tnomen_index))
        nearest_ss = get_info_from_tuple(variant, "SPiP_NearestSS", tnomen_index)
        score = get_generalized_splicing_score(
            variant, score, dist, nearest_ss, score_config
        )

    return score


def get_spliceai_score(variant: cyvcf2.Variant, score: int, score_config: dict) -> int:
    SPLICEAI_THRESHOLD = score_config["Threshold_SpliceAI"]

    tnomen_index = 0
    # TODO: change when HOWARD's transcript prioritization is done
    # spliceai_symbols = get_variant_info(variant, "SpliceAI_SYMBOL").split(",")
    # tnomen = get_variant_info(variant, "TNOMEN")
    # tnomen_index =  spliceai_symbols.index(tnomen)

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
        get_info_from_tuple(variant, f, tnomen_index) for f in spliceai_scores_fields
    ]
    splice_ai_pos = [
        get_info_from_tuple(variant, f, tnomen_index) for f in spliceai_pos_fields
    ]

    if splice_ai_scores == [""] * 4:
        # no spliceai score
        return score
    max_delta_score = max([float(v) for v in splice_ai_scores if v != ""])
    if max_delta_score > SPLICEAI_THRESHOLD:

        index = splice_ai_scores.index(max_delta_score)
        dist = int(splice_ai_pos[index])
        if index in (2, 3):
            nearest_ss = "donor"  # 5'
        elif index in (0, 1):
            nearest_ss = "acceptor"  # 3'
        else:
            raise ValueError(
                f"Unexpected index when computing nearest_ss; index should  be in [0,3] ; got instead: index:{index} ; splice_ai_scores:{splice_ai_scores}"
            )

        score = get_generalized_splicing_score(
            variant, score, dist, nearest_ss, score_config
        )

    return score


def get_generalized_splicing_score(
    variant: cyvcf2.Variant, score: int, dist: int, nearest_ss: str, score_config: dict
) -> int:
    S_ESSENTIALSPLICE = score_config["S_EssentialSplice"]
    S_CLOSESPLICE = score_config["S_CloseSplice"]
    S_DEEPSPLICE = score_config["S_DeepSplice"]

    if score < S_ESSENTIALSPLICE:
        if dist in (1, 2):
            score = S_ESSENTIALSPLICE + get_bonus(variant, score_config)
    elif score < S_CLOSESPLICE:
        if (nearest_ss == "donor" and -3 < dist < 6) or (
            nearest_ss == "acceptor" and -12 < dist < 2
        ):
            score = S_CLOSESPLICE + get_bonus(variant, score_config)
    elif score < S_DEEPSPLICE:
        snpeff_annotation = get_variant_info(variant, "snpeff_annotation").lower()
        if "intron" in snpeff_annotation:
            score = S_DEEPSPLICE + get_bonus(variant, score_config)

    # TODO: adjust condition with bonuses
    return score


def get_bonus(variant: cyvcf2.Variant, score_config: dict):
    """for splicing: phastcons

    Planned new bonuses:  + phylop + cadd
    """
    bonus = 0
    phastcons = get_variant_info(variant, "phastCons100way")
    if phastcons != "" and float(phastcons) > score_config["Threshold_Phastcons"]:
        bonus += score_config["B_phastCons"]

    return bonus
