import filecmp
import logging as log
import os
from os.path import join as osj

from vannotplus.commons import load_config
from vannotplus.family.barcode import main_barcode_fast
import tempfile


def test_main_barcode_fast_family(debug_mode=False):
    current_dir = os.path.dirname(os.path.abspath(__file__))
    config = load_config(osj(current_dir, "data", "config.yml"))
    input_vcf = osj(current_dir, "data", "test_family.vcf")
    app = "FAKE_APP"
    control_vcf = osj(current_dir, "controls", "test_family_barcoded.vcf")

    if debug_mode:
        output_vcf = "test_family_barcoded.vcf"
        main_barcode_fast(input_vcf, output_vcf, app, config)
        return

    tmp_dir = tempfile.TemporaryDirectory()
    output_vcf = osj(tmp_dir.name, "test_family_barcoded.vcf")
    main_barcode_fast(input_vcf, output_vcf, app, config)

    assert os.path.exists(output_vcf), "Output VCF file was not created."
    assert filecmp.cmp(
        output_vcf, control_vcf
    ), "Output VCF does not match control VCF."
    tmp_dir.cleanup()


def test_main_barcode_fast_pools(debug_mode=False):
    current_dir = os.path.dirname(os.path.abspath(__file__))
    config = load_config(osj(current_dir, "data", "config.yml"))
    input_vcf = osj(current_dir, "data", "test_pools.vcf")
    app = "FAKE_APP"
    control_vcf = osj(current_dir, "controls", "test_pools_barcoded.vcf")

    if debug_mode:
        output_vcf = "test_pools_barcoded.vcf"
        main_barcode_fast(input_vcf, output_vcf, app, config)
        return

    tmp_dir = tempfile.TemporaryDirectory()
    output_vcf = osj(tmp_dir.name, "test_pools_barcoded.vcf")
    main_barcode_fast(input_vcf, output_vcf, app, config)

    assert os.path.exists(output_vcf), "Output VCF file was not created."
    assert filecmp.cmp(
        output_vcf, control_vcf
    ), "Output VCF does not match control VCF."
    tmp_dir.cleanup()


if __name__ == "__main__":
    log.basicConfig(level=log.DEBUG)
    test_main_barcode_fast_family(debug_mode=True)
    test_main_barcode_fast_pools(debug_mode=True)

    # main_barcode_fast(
    #     "/home1/HUB/bin/vannotplus/vannotplus/tinyexome.vcf",
    #     "/home1/HUB/bin/vannotplus/vannotplus/exome.vcf",
    #     "WES_AGILENT",
    #     load_config("/home1/HUB/bin/vannotplus/vannotplus/src/vannotplus/config.yml"),
    # )
