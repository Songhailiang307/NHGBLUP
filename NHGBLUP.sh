#!/bin/bash
module load python/3.8.3
module load plink/1.9
module load conda/3-2020.07
module load R/4.3.2
source activate clean_env


PHE=$1
SNP=$2
REL=$3
VAL=$4


cp ${REL} relphe.txt
cp ${VAL} valphe.txt
./bin/plink  --bfile ${SNP}  --recodeA --out simu --allow-extra-chr

awk '{print $1,$3}' valphe.txt >val_tbv
awk '{print $1,$3}' relphe.txt >pheno.txt

awk 'NR>1 {print $2}' simu.raw  > geno_id

./bin/hiblup --make-xrm --code-method 1 --bfile ${SNP} --add --dom --step 10000 --thread 50 --out demo
./bin/hiblup --trans-xrm --xrm demo.GA --out demoA
./bin/hiblup --trans-xrm --xrm demo.GD --out demoD	

python ./bin/GINV.py demoD.txt DMA.txt_GINV
python ./bin/Gma_to_3lineID.py DMA.txt_GINV geno_id Ginv
python ./bin/GINV.py demoA.txt GMA.txt_GINV
python ./bin/Gma_to_3lineID.py GMA.txt_GINV geno_id Ginv

#GBLUP and GDBLUP
awk '{print $1,"1",$2}' pheno.txt >relphe.txt
cp ./bin/gblup.DIR ../../bin/gdblup.DIR  ../../bin/g2dblup.DIR ./
./bin/r_dmuai gblup
./bin/r_dmuai gdblup


python ./bin/ggebv.py gdblup.SOL val_tbv val_gdebv
python ./bin/gebv.py gblup.SOL val_tbv val_gebv

python ./bin/COR_REG_used2.py val_gdebv GDBLUP.txt 1
python ./bin/COR_REG_used2.py val_gebv GBLUP.txt 1	


####################################
# 2. NHGBLUP and NHGDBLUP
########################################

for i in $(awk 'BEGIN{for(x=0.1;x<=1.0001;x+=0.1) printf "%.1f\n",x}')
do
	echo "Running AD dBLUP with parameter = $i"

	################ tanh ################
	python ./bin/dBLUP.py tanh $i

	awk 'BEGIN{OFS="	"} NR==FNR{a[$1]=$2;next} FNR>1 && ($1 in a){print $1,a[$1],$2}' val_tbv Validation_GEBV.txt > gblup_gebv
	awk 'BEGIN{OFS="	"} NR==FNR{a[$1]=$2;next} FNR>1 && ($1 in a){print $1,a[$1],$3}' val_tbv Validation_GEBV.txt > gvblup_gebv

	python ./bin/COR_REG_used2.py gblup_gebv GBLUP_${i}.txt 1
	python ./bin/COR_REG_used2.py gvblup_gebv D2GBLUP_${i}.txt 1

	python ./bin/GINV.py Ghybrid.txt Ghybrid.txt_GINV
	python ./bin/Gma_to_3lineID.py Ghybrid.txt_GINV geno_id Ginv
	./bin/r_dmuai g2dblup
	python ./bin/ggebv.py g2dblup.SOL val_tbv val_g2debv
	python ./bin/COR_REG_used2.py val_g2debv D2GDBLUP_${i}.txt 1

done
