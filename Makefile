install:
	pip install -e .

run:
	python -m vannotplus barcode -i /home1/L_PROD/NGS/BAS/HOWARD/data/nicaises/ped/trio_base.vcf -o /home1/HUB/bin/vannotplus/vannotplus/output.vcf -p /home1/data/WORK_DIR_SAM/Ped_raw/Data/ -a WES_AGILENT