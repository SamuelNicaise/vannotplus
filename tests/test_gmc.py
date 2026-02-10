import filecmp
import os
from os.path import join as osj
import tempfile

from vannotplus.commons import load_config, set_log_level
from vannotplus.annot.score import main_annot


def test_gmc(debug_mode=False):
    current_dir = os.path.dirname(os.path.abspath(__file__))
    config = load_config(osj(current_dir, "data", "config.yml"))
    input_vcf = osj(current_dir, "data", "gmc_mini_input.vcf")
    control_vcf = osj(current_dir, "controls", "gmc_mini_control.vcf")

    tmp_dir = tempfile.TemporaryDirectory()
    output_vcf = osj("gmc_mini_out.vcf")
    config["gmc"]["do_filtered_gmc"] = False
    main_annot(input_vcf, output_vcf, config, do_vannotscore=False, do_filtered_gmc=False)
    print("output_vcf", output_vcf)

    assert os.path.exists(output_vcf), "Output VCF file was not created."
    assert filecmp.cmp(
        output_vcf, control_vcf
    ), f"Output {output_vcf} does not match control {control_vcf}."
    print("exists", os.path.exists(output_vcf))

    if not debug_mode:
        tmp_dir.cleanup()

def test_filtered_gmc(debug_mode=False):
    current_dir = os.path.dirname(os.path.abspath(__file__))
    config = load_config(osj(current_dir, "data", "config.yml"))
    input_vcf = osj(current_dir, "data", "filtered_gmc_input.vcf")
    control_vcf = osj(current_dir, "controls", "filtered_gmc_control.vcf")

    tmp_dir = tempfile.TemporaryDirectory()
    output_vcf = osj(tmp_dir.name, "filtered_gmc_control.vcf")
    main_annot(input_vcf, output_vcf, config, do_vannotscore=False, do_filtered_gmc=False)
    print("output_vcf", output_vcf)

    assert os.path.exists(output_vcf), "Output VCF file was not created."
    assert filecmp.cmp(
        output_vcf, control_vcf
    ), f"Output {output_vcf} does not match control {control_vcf}."
    print("exists", os.path.exists(output_vcf))

    if not debug_mode:
        tmp_dir.cleanup()

if __name__ == "__main__":
    set_log_level("DEBUG")
    test_gmc(debug_mode=True)
