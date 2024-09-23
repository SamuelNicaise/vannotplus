install:
	pip install -e .

barcode:
	python -m vannotplus barcode -i /home1/L_PROD/NGS/BAS/HOWARD/data/nicaises/ped/trio_base.vcf -o /home1/HUB/bin/vannotplus/vannotplus/output.vcf -c src/vannotplus/config.yml -a WES_AGILENT

exomiser:
	python -m vannotplus exomiser -i /home1/L_PROD/NGS/BAS/HOWARD/data/nicaises/ped/jb_trio_output.vcf -o /home1/L_PROD/NGS/BAS/HOWARD/data/nicaises/test/merged.vcf -c src/vannotplus/config.yml -a WES_AGILENT -v debug

send:
	rm -rf /home1/L_PROD/NGS/BAS/HOWARD/data/nicaises/vannotplus
	cp -r /home1/HUB/bin/vannotplus/vannotplus /home1/L_PROD/NGS/BAS/HOWARD/data/nicaises