DIR="$(dirname "${BASH_SOURCE[0]}")"
DIR="$(realpath "${DIR}")"


# Download the SILVA 138 pre-trained classifier for 16s sequences (99% OTUs from V4 region of sequences) 
# from https://docs.qiime2.org/2023.5/data-resources/ (accessed on 14 July 2023)
ls $DIR/16s/ | ( grep silva-138-99-515-806-nb-classifier.qza > /dev/null && echo "SILVA 138 pre-trained classifier already downloaded" ) \
|| ( echo "Download the SILVA 138 pre-trained classifier" && wget -P $DIR/16s/ https://data.qiime2.org/2023.5/common/silva-138-99-515-806-nb-classifier.qza ) \
|| ( echo "Can't download the SILVA 138 pre-trained classifier." && exit )

# Download the UNITE pre-trained classifier for ITS sequences (99% OTUs from ITS1f/ITS2 region of sequences)
# from https://github.com/colinbrislawn/unite-train/releases (accessed on 14 July 2023)
ls $DIR/ITS/ | ( grep unite_ver9_99_all_29.11.2022-Q2-2023.5.qza > /dev/null && echo "UNITE pre-trained classifier already downloaded" ) \
|| ( echo "Download the UNITE pre-trained classifier" && wget -P $DIR/ITS/ https://github.com/colinbrislawn/unite-train/releases/download/9.0-qiime2-2023.5-demo/unite_ver9_99_all_29.11.2022-Q2-2023.5.qza ) \
|| ( echo "Can't download the UNITE pre-trained classifier." && exit )


# Download SRA files used in this article (16s data)
mkdir $DIR/16s/fastq

PRJNA1127768_16s=(SRR29528551 SRR29528550 SRR29528459 SRR29528536 SRR29528477 SRR29528450 \
SRR29528527 SRR29528468 SRR29528441 SRR29528494 SRR29528549 SRR29528514 SRR29528483 SRR29528482 \
SRR29528481 SRR29528480 SRR29528479 SRR29528462 SRR29528461 SRR29528460 SRR29528458 SRR29528457 \
SRR29528456 SRR29528455 SRR29528542 SRR29528541 SRR29528540 SRR29528539 SRR29528538 SRR29528537 \
SRR29528535 SRR29528510 SRR29528509 SRR29528508 SRR29528507 SRR29528506 SRR29528505 SRR29528504 \
SRR29528503 SRR29528478 SRR29528476 SRR29528475 SRR29528474 SRR29528473 SRR29528472 SRR29528471 \
SRR29528454 SRR29528453 SRR29528452 SRR29528451 SRR29528449 SRR29528448 SRR29528447 SRR29528534 \
SRR29528533 SRR29528532 SRR29528531 SRR29528530)

all_SRA=(${PRJNA1127768_16s[*]})

for sra in ${all_SRA[*]}
do
ls $DIR/16s/fastq | grep $sra > /dev/null || ( prefetch $sra -p -O $DIR/16s/fastq/ && fasterq-dump $DIR/16s/fastq/$sra/$sra.sra -O $DIR/16s/fastq \
&& gzip $DIR/16s/fastq/$sra"_1.fastq" \
&& gzip $DIR/16s/fastq/$sra"_2.fastq" \
&& rm -R $DIR/16s/fastq/$sra/ )
done

for sra in ${all_SRA[*]}
do
ls $DIR/16s/fastq | grep $sra > /dev/null 
    if [[ $? != 0 ]]
        then
        echo "Can't download $sra. Try run script again." && exit 111
        break
    fi
done

if [[ $? == 111 ]]
    then
    exit
fi

# Download SRA files used in this article (ITS data)

mkdir $DIR/ITS/fastq

PRJNA1127768_ITS=(SRR29528529 SRR29528528 SRR29528502 SRR29528501 SRR29528500 SRR29528499 \
SRR29528498 SRR29528497 SRR29528496 SRR29528495 SRR29528470 SRR29528469 SRR29528467 SRR29528466 \
SRR29528465 SRR29528464 SRR29528463 SRR29528446 SRR29528445 SRR29528444 SRR29528443 SRR29528442 \
SRR29528440 SRR29528439 SRR29528526 SRR29528525 SRR29528524 SRR29528523 SRR29528522 SRR29528521 \
SRR29528520 SRR29528519 SRR29528493 SRR29528492 SRR29528491 SRR29528490 SRR29528489 SRR29528488 \
SRR29528487 SRR29528438 SRR29528437 SRR29528436 SRR29528548 SRR29528547 SRR29528546 SRR29528545 \
SRR29528544 SRR29528543 SRR29528518 SRR29528517 SRR29528516 SRR29528515 SRR29528513 SRR29528512 \
SRR29528511 SRR29528486 SRR29528485 SRR29528484)

all_SRA=(${PRJNA1127768_ITS[*]})

for sra in ${all_SRA[*]}
do
ls $DIR/ITS/fastq | grep $sra > /dev/null || ( prefetch $sra -p -O $DIR/ITS/fastq/ && fasterq-dump $DIR/ITS/fastq/$sra/$sra.sra -O $DIR/ITS/fastq \
&& gzip $DIR/ITS/fastq/$sra"_1.fastq" \
&& gzip $DIR/ITS/fastq/$sra"_2.fastq" \
&& rm -R $DIR/ITS/fastq/$sra/ )
done

for sra in ${all_SRA[*]}
do
ls $DIR/ITS/fastq | grep $sra > /dev/null 
    if [[ $? != 0 ]]
        then
        echo "Can't download $sra. Try run script again." && exit 111
        break
    fi
done

if [[ $? == 111 ]]
    then
    exit
fi

# Rename SRA for input in qiime2 (16s data)
data_16s_1=($(ls $DIR/16s/fastq | grep _1.fastq.gz))
for oldname in ${data_16s_1[*]}
do
newname=$(echo $oldname | sed "s/_1.fastq.gz/_S1_L001_R1_001.fastq.gz/g")
mv $DIR/16s/fastq/$oldname $DIR/16s/fastq/$newname
done

data_16s_2=($(ls $DIR/16s/fastq | grep _2.fastq.gz))
for oldname in ${data_16s_2[*]}
do
newname=$(echo $oldname | sed "s/_2.fastq.gz/_S1_L001_R2_001.fastq.gz/g")
mv $DIR/16s/fastq/$oldname $DIR/16s/fastq/$newname
done

# Rename SRA for input in qiime2 (ITS data)
data_ITS_1=($(ls $DIR/ITS/fastq | grep _1.fastq.gz))
for oldname in ${data_ITS_1[*]}
do
newname=$(echo $oldname | sed "s/_1.fastq.gz/_S1_L001_R1_001.fastq.gz/g")
mv $DIR/ITS/fastq/$oldname $DIR/ITS/fastq/$newname
done

data_ITS_2=($(ls $DIR/ITS/fastq | grep _2.fastq.gz))
for oldname in ${data_ITS_2[*]}
do
newname=$(echo $oldname | sed "s/_2.fastq.gz/_S1_L001_R2_001.fastq.gz/g")
mv $DIR/ITS/fastq/$oldname $DIR/ITS/fastq/$newname
done