install:
	pip install -e .

barcode:
	python -m vannotplus barcode -i /home1/L_PROD/NGS/BAS/HOWARD/data/nicaises/ped/trio_base.vcf -o /home1/HUB/bin/vannotplus/vannotplus/output.vcf -c src/vannotplus/config.yml -a WES_AGILENT

exomiser:
	python -m vannotplus exomiser -i /home1/L_PROD/NGS/BAS/HOWARD/data/nicaises/ped/jb_trio_output.vcf -o /home1/L_PROD/NGS/BAS/HOWARD/data/nicaises/test/merged.vcf -c src/vannotplus/config.yml -a WES_AGILENT -v debug

annot:
	python -m vannotplus annot -i /home1/L_PROD/NGS/BAS/HOWARD/data/nicaises/SGT151797.final.clean.full-annotation.v4.noprio.vcf.gz -o /home1/L_PROD/NGS/BAS/HOWARD/data/nicaises/score/annot.vcf -c src/vannotplus/config.yml -v debug

annot_wes:
	python -m vannotplus annot -i /home1/L_PROD/NGS/BAS/HOWARD/data/nicaises/KLA2403985.final.vcf -o /home1/L_PROD/NGS/BAS/HOWARD/data/nicaises/score/annotKLA.vcf -c src/vannotplus/config.yml

debug:
	python -m vannotplus exomiser -i /home1/L_PROD/NGS/BAS/DOCKER_STARK_MAIN_FOLDER/data/diebrant/SGT2301773.vcf.gz -o /home1/L_PROD/NGS/BAS/HOWARD/data/nicaises/test/SGT2301773.vcf -c src/vannotplus/config.yml -a WES_AGILENT -v debug

send:
	rm -rf /home1/L_PROD/NGS/BAS/HOWARD/data/nicaises/vannotplus
	cp -r /home1/HUB/bin/vannotplus/vannotplus /home1/L_PROD/NGS/BAS/HOWARD/data/nicaises