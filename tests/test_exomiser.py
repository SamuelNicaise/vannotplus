import filecmp
import os
from os.path import join as osj
import tempfile

from vannotplus.commons import load_config, set_log_level
from vannotplus.exomiser.exomiser import main_exomiser

def test_exomiser(debug_mode=False):
    current_dir = os.path.dirname(os.path.abspath(__file__))
    config = load_config(osj(current_dir, "data", "config.yml"))
    input_vcf = osj(current_dir, "data", "exomiser.vcf")
    app = "FAKE_APP"
    control_vcf = osj(current_dir, "controls", "exomiser.vcf")

    tmp_dir = tempfile.TemporaryDirectory()
    output_vcf = osj(tmp_dir.name, "exomiser_out.vcf")
    main_exomiser(input_vcf, output_vcf, app, config)
    print("output_vcf", output_vcf)

    assert os.path.exists(output_vcf), "Output VCF file was not created."
    assert filecmp.cmp(
        output_vcf, control_vcf
    ), f"Output {output_vcf} does not match control {control_vcf}."
    print("exists", os.path.exists(output_vcf))

    if not debug_mode:
        tmp_dir.cleanup()


if __name__ == "__main__":
    set_log_level("DEBUG") # will also keep exomiser's container's temporary files
    test_exomiser(debug_mode=True)

    # current_dir = os.path.dirname(os.path.abspath(__file__))
    # config = load_config(osj(current_dir, "data", "config.yml"))
    # # config = load_config("/home1/HUB/bin/vannotplus/vannotplus/src/vannotplus/config.yml")
    # input_vcf = osj(current_dir, "data", "exomiser.vcf")
    # app = "FAKE_APP"
    # output_vcf = osj(current_dir, "exomiser_out.vcf")
    # control_vcf = osj(current_dir, "controls", "exomiser.vcf")

    # main_exomiser(input_vcf, output_vcf, app, config)
    # assert filecmp.cmp(
    #     output_vcf, control_vcf
    # ), f"Output {output_vcf} does not match control {control_vcf}."