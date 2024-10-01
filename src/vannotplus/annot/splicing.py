from cyvcf2 import cyvcf2

from vannotplus.commons import get_variant_info


def get_splicing_score(variant: cyvcf2.Variant, score: int, config: dict) -> int:
    """
    Returns the maximum possible score from either spip or spliceAI
    If score is already higher than the maximum possible with splicing, return it directly.
    """

    score_matrix = config["score_matrix"]

    # do not bother computing score if already > max
    if score >= score_matrix["S_ESSENTIALSPLICE"]:
        return score
    else:
        score = get_spip_score(variant, score, score_matrix)

    # same if SPIP already got max score
    if score >= score_matrix["S_ESSENTIALSPLICE"]:
        return score
    else:
        score = get_spliceai_score(variant, score, score_matrix)

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
    return info


def get_spip_score(variant: cyvcf2.Variant, score: int, score_matrix: dict) -> int:
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
            variant, score, dist, nearest_ss, score_matrix
        )

    return score


def get_spliceai_score(variant: cyvcf2.Variant, score: int, score_matrix: dict) -> int:
    SPLICEAI_THRESHOLD = score_matrix["SPLICEAI_THRESHOLD"]

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
            variant, score, dist, nearest_ss, score_matrix
        )

    return score


def get_generalized_splicing_score(
    variant: cyvcf2.Variant, score: int, dist: int, nearest_ss: str, score_matrix: dict
) -> int:
    S_ESSENTIALSPLICE = score_matrix["S_ESSENTIALSPLICE"]
    S_CLOSESPLICE = score_matrix["S_CLOSESPLICE"]
    S_DEEPSPLICE = score_matrix["S_DEEPSPLICE"]

    if score < S_ESSENTIALSPLICE:
        if dist in (1, 2):
            score = S_ESSENTIALSPLICE  # + phastcons + phylop + cadd
    elif score < S_CLOSESPLICE:
        if (nearest_ss == "donor" and -3 < dist < 6) or (
            nearest_ss == "acceptor" and -12 < dist < 2
        ):
            score = S_CLOSESPLICE  # + phastcons  + phylop + cadd
    elif score < S_DEEPSPLICE:
        snpeff_annotation = get_variant_info(variant, "snpeff_annotation").lower()
        if "intron" in snpeff_annotation:
            score = S_DEEPSPLICE  # + phastcons  + phylop + cadd

    # TODO: adjust condition with bonuses
    return score
