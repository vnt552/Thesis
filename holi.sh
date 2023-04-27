for infile in $(pwd)/*sga4.fq
do
bname=$(basename $infile)
echo $bname
#bname2=$(echo $bname | sed 's/.fq*/_holi/')
basepath=$(pwd)/
#basefolder=$basepath
echo $basepath
#echo $bname2
#mkdir $basepath$bname2
#cd $basepath$bname2
#pwd

## Mapping remaining reads against downloaded databases
for DB in /willerslev/datasets/ncbi_nt_Nov2020/nt.fa.?
do
echo Mapping $bname against $DB
nice -n 5 bowtie2 --threads 60 -k 5000 -x $DB -U $bname --no-unal | samtools view -bS - > $bname.$(basename $DB).bam
done

for DB in /willerslev/datasets/refseq_23dec2020/vertebrate_other/vertebrate_other.?
do
echo Mapping $bname against $DB
nice -n 5 bowtie2 --threads 60 -k 5000 -x $DB -U $bname --no-unal | samtools view -bS - > $bname.$(basename $DB).bam
done

for DB in /willerslev/datasets/refseq_23dec2020/vertebrate_mammalian/vertebrate_mammalian.?
do
echo Mapping $bname against $DB
nice -n 5 bowtie2 --threads 60 -k 5000 -x $DB -U $bname --no-unal | samtools view -bS - > $bname.$(basename $DB).bam
done

for DB in /willerslev/edna/database/Homo_sapiens/homo_sapiens.fa?
do
echo Mapping $bname against $DB
nice -n 5 bowtie2 --threads 60 -k 5000 -x $DB -U $bname --no-unal | samtools view -bS - > $bname.$(basename $DB).bam
done

for DB in /willerslev/datasets/refseq_23dec2020/invertebrate/invertebrate.?
do
echo Mapping $bname against $DB
nice -n 5 bowtie2 --threads 60 -k 5000 -x $DB -U $bname --no-unal | samtools view -bS - > $bname.$(basename $DB).bam
done

for DB in /willerslev/datasets/refseq_23dec2020/archaea_fungi_virus/archaea_fungi_virus.fa
do
echo Mapping $bname against $DB
nice -n 5 bowtie2 --threads 60 -k 5000 -x $DB -U $bname --no-unal | samtools view -bS - > $bname.$(basename $DB).bam
done

for DB in /willerslev/datasets/ycw/polar_animals/arctic_animals.fa
do
echo Mapping $bname against $DB
nice -n 5 bowtie2 --threads 60 -k 5000 -x $DB -U $bname --no-unal | samtools view -bS - > $bname.$(basename $DB).bam
done

for DB in /willerslev/datasets/ycw/polar_animals/arctic_animals_other.fa
do
echo Mapping $bname against $DB
nice -n 5 bowtie2 --threads 60 -k 5000 -x $DB -U $bname --no-unal | samtools view -bS - > $bname.$(basename $DB).bam
done

for DB in /willerslev/datasets/ycw/polar_animals/arctic_animal_b3.fa
do
echo Mapping $bname against $DB
nice -n 5 bowtie2 --threads 60 -k 5000 -x $DB -U $bname --no-unal | samtools view -bS - > $bname.$(basename $DB).bam
done

for DB in /willerslev/datasets/refseq_23dec2020/plant/plant.?
do
echo Mapping $bname against $DB
nice -n 5 bowtie2 --threads 60 -k 5000 -x $DB -U $bname --no-unal | samtools view -bS - > $bname.$(basename $DB).bam
done

for DB in /willerslev/datasets/refseq_23dec2020/protozoa/protozoa.fa
do
echo Mapping $bname against $DB
nice -n 5 bowtie2 --threads 60 -k 5000 -x $DB -U $bname --no-unal | samtools view -bS - > $bname.$(basename $DB).bam
done

for DB in /willerslev/datasets/ycw/norway_plant_ctgenoms/complete_genomes/final_db/norPlantCom.?
do
echo Mapping $bname against $DB
nice -n 5 bowtie2 --threads 60 -k 5000 -x $DB -U $bname --no-unal | samtools view -bS - > $bname.$(basename $DB).bam
done

for DB in /willerslev/edna/antonio/geogenetics/DBs/gtdb/r202/bowtie2/gtdb-r202.?
do
echo Mapping $bname against $DB
nice -n 5 bowtie2 --threads 60 -k 5000 -x $DB -U $bname --no-unal | samtools view -bS - > $bname.$(basename $DB).bam
done


## Merging all alignment files
samtools merge --verbosity 5  $bname.merged.sam.gz $bname*.bam -@ 60

## Sorting the merged sam.gz file
echo Printing header
nice -n 5 time samtools view --threads 60  -H $bname.merged.sam.gz | gzip > $bname.merged.Header.sam.gz
echo Printing alignment
nice -n 5 time samtools view --threads 60 $bname.merged.sam.gz | gzip > $bname.merged.alignment.sam.gz
echo Sorting alignment file
time /willerslev/software/gz-sort/gz-sort -S 30G -P 10 $bname.merged.alignment.sam.gz $bname.merged.alignment.sort.sam.gz
echo Merging Header and sorted alignment
time zcat $bname.merged.Header.sam.gz $bname.merged.alignment.sort.sam.gz | samtools view -h -o $bname.merged.sort.sam.gz
rm $bname.merged.Header.sam.gz $bname.merged.alignment.sam.gz $bname.merged.alignment.sort.sam.gz

ls -lh *bam

rm *nt.fa*
rm *vert_other.?*
rm *vert_mam.?*
rm *highres-db-micro_20201125*
rm *invert.?*
rm *norPlantCom*
rm *viral_fungi_archaea*
rm *arctic_animal*
rm *plant*
rm $bname.merged.sam.gz

## Running ngsLCA and metaDamage
for file in *.merged.sort.sam.gz
do
nam=/willerslev/edna/ncbi_taxonomy3dec2020/names.dmp
nod=/willerslev/edna/ncbi_taxonomy3dec2020/nodes.dmp
ac2tax=/willerslev/edna/ncbi_taxonomy3dec2020/combined_taxid_accssionNO_20201120.gz
/willerslev/edna/metadamage/metadamage lca -simscorelow 0.95 -simscorehigh 1.0 -names $nam -nodes $nod -acc2tax $ac2tax -bam $file -outnames $file.allrank -howmany 15
done

cd $basepath
done
