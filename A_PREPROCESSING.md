# PRE-PROCESSING

## Reads QC
```{bash}
module add fastqc/0.11.8

WORKDIR="path_to_working_directory"
FQDIR="path_to_fastq_files"

echo "Running FASTQC"

cd ${FQDIR}

ls *.fastq.gz | xargs -I % -n 1 -P 48 sh -c 'echo %; fastqc -t 48 -q %'

cd ${WORKDIR}

echo "Finished FASTQC"
```

&nbsp;

## Reads Trimming
```{bash}
WORKDIR="path_to_working_directory"
FQDIR="path_to_fastq_files"

echo "Running TRIMMOMATIC"

cd ${FQDIR}

cat samples.txt | while read line
    do
        echo "Processing" ${line}
        fq1=`echo ${line}"_R1_001.fastq.gz"`
        fq2=`echo ${line}"_R2_001.fastq.gz"`
        echo ${fq1} "|" ${fq2}
        
        cat ${line}"_R1_001_fastqc/fastqc_data.txt" | grep Over -A 100 | grep 'Illumina\|TruSeq' | grep -P '^[A-Z]' | nl | awk '{print ">" $1 "_adapter\n" $2}' > ${line}"_R1_adapters.fa"
        cat ${line}"_R2_001_fastqc/fastqc_data.txt" | grep Over -A 100 | grep 'Illumina\|TruSeq' | grep -P '^[A-Z]' | nl | awk '{print ">" $1 "_adapter\n" $2}' > ${line}"_R2_adapters.fa"
        cat ${line}"_R1_adapters.fa" ${line}"_R2_adapters.fa" > ${line}"_adapters.fa"
        rm ${line}"_R1_adapters.fa"
        rm ${line}"_R2_adapters.fa"
        
        paired1=`basename ${fq1} | sed -e "s/.fastq.gz/.NoAdapt.Trim.fastq.gz/"`
        paired2=`basename ${fq2} | sed -e "s/.fastq.gz/.NoAdapt.Trim.fastq.gz/"`
        unpaired1=`basename ${fq1} | sed -e "s/.fastq.gz/.Unp.NoAdapt.Trim.fastq.gz/"`
        unpaired2=`basename ${fq2} | sed -e "s/.fastq.gz/.Unp.NoAdapt.Trim.fastq.gz/"`
        echo $paired1 "|" $paired2
        echo $unpaired1 "|" $unpaired2
        
        java -jar /work/RESOURCES/TOOLS/Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 $fq1 $fq2 "$paired1" "$unpaired1" "$paired2" "$unpaired2" ILLUMINACLIP:${line}"_adapters.fa":2:30:10 -threads 48 SLIDINGWINDOW:4:18 LEADING:7 TRAILING:7 MINLEN:35
    done

cd ${WORKDIR}

echo "Finished TRIMMOMATIC"
```

&nbsp;

## Reads Alignment
```{bash}
module add star/2.5.2b

WORKDIR="path_to_working_directory"
FQDIR="path_to_fastq_files"
TRIMDIR="path_to_trimmed_fastq_files"

echo "Running STAR"

cd ${TRIMDIR}

cat samples.txt | while read line
    do
        echo "Processing" ${line}
        fq1=`echo ${line}"_R1_001.NoAdapt.Trim.fastq.gz"`
        fq2=`echo ${line}"_R2_001.NoAdapt.Trim.fastq.gz"`
        outputname=`basename ${fq1} | sed -e "s/_R1_001.NoAdapt.Trim.fastq.gz/_OUTPUT/"`
        echo "Processing" ${line} "-->" ${outputname}
        
        STAR --runThreadN 48 --genomeDir /work/RESOURCES/DATABASES/MM10_GRCm38p6_GENCODEvM17_STAR --readFilesIn ${fq1} ${fq2} --readFilesCommand zcat --sjdbGTFfile /work/RESOURCES/DATABASES/MM10_GRCm38p6_GENCODEvM17_STAR/gencode.vM17.protein_coding.gtf --outFilterType BySJout --outFilterMultimapNmax 10 --alignSJoverhangMin 10 --alignSJDBoverhangMin 1 --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outFilterMismatchNmax 3 --twopassMode Basic --outFileNamePrefix ${outputname}
    done

cd ${WORKDIR}

echo "Running STAR"
```

&nbsp;

## BAM Filtering
```{bash}
module add samtools
module add RSeQC/2.6.4

WORKDIR="path_to_working_directory"
FQDIR="path_to_fastq_files"
TRIMDIR="path_to_trimmed_fastq_files"
STARDIR="path_to_star_alignment_output"

echo "Filtering BAM"

cd ${STARDIR}

## Multimapped reads
echo "=====> Listing multimapped reads"
mkdir MULTIMAPPEDREADS/
for file in `ls *Aligned.sortedByCoord.out.bam`
 do
  echo ${file}
  outputname=`basename ${file} | sed -e "s/Aligned.sortedByCoord.out.bam/.MultiMappedReads.txt/"`
  echo ${outputname}
  samtools view ${file} | grep -v NH:i:1 | perl -pe 's/AS.+(NH:i:\d+)/\$1/' | cut -f1,10,12 | perl -pe 's/NH:i://' | sort -u -k3,3nr > MULTIMAPPEDREADS/${outputname}
done
echo "=====> Finished listing multimapped reads"


## Fetch only primary alignment (remove unmapped,chimeric etc etc)
echo "=====> Fetching primary alignment reads"
for file in `ls *Aligned.sortedByCoord.out.bam`
 do
  echo ${file}
  outputname=`basename ${file} | sed -e "s/Aligned.sortedByCoord.out.bam/.Primary.bam/"`
  echo ${outputname}
  samtools view -F 256 -b ${file} > ${outputname}
 done
echo "=====> Finished fetching primary alignment"


## Split Primary alignment in rRNA (in.bam) and non rRNA (ex.bam)
echo "=====> Splitting rRNA reads"
for file in `ls *.Primary.bam`
 do
  samname=`basename ${file} | sed -e "s/.Primary//"`
  split_bam.py -i ${file} -r /work/RESOURCES/DATABASES/MM10_GRCm38p6_GENCODEvM17_STAR/mm10_rRNA.bed -o "${samname}"
  echo $samname
 done
echo "=====> Finished splitting rRNA reads"

rm *.junk.bam


## Get uniquely mapped
echo "=====> Fetching uniquely mapped reads"
 for file in `ls *.ex.bam`
  do
   echo ${file}
   outputname=`basename $file | sed -e "s/.ex.bam/.MAPQ.bam/"`
   echo ${outputname}
   (samtools view -H $file; samtools view -F 2308 $file | grep -w 'NH:i:1') | samtools view -bS - > "$outputname"
  done
echo "=====> Finished fetching uniquely mapped reads"


## Counting the number of reads per sample with and without (i) MAPQ filtering, (ii) rRNA removal, (iii) mapped and unmapped
for file in `ls *.MAPQ.bam`
    do
        basename=`echo ${file} | sed "s/.bam.MAPQ.bam//g"`
        starout=`echo ${basename}"Aligned.sortedByCoord.out.bam"`
        primary=`echo ${basename}".Primary.bam"`
        worrna=`echo ${basename}".bam.ex.bam"`
        mapqout=`echo ${basename}".bam.MAPQ.bam"`

        count1=`samtools view -c ${starout}`
        count2=`samtools view -c ${primary}`
        count3=`samtools view -c ${worrna}`
        count4=`samtools view -c ${mapqout}`

        echo ${basename} ${count1} ${count2} ${count3} ${count4} | sed "s/ /\t/g"
    done

cd ${WORKDIR}

```

&nbsp;


## Counting reads per gene
```{bash}
module add HTSeq/0.6.1

WORKDIR="path_to_working_directory"
FQDIR="path_to_fastq_files"
TRIMDIR="path_to_trimmed_fastq_files"
STARDIR="path_to_star_alignment_output"

echo "Counting reads per gene"
cd ${STARDIR}

ls *.MAPQ.bam | xargs -I % -n 1 -P 36 sh -c 'echo %; htseq-count -f bam -r pos -s reverse -t gene -i gene_name -m intersection-strict % /work/RESOURCES/DATABASES/MM10_GRCm38p6_GENCODEvM17_STAR/gencode.vM17.protein_coding.gtf > %.count;'

echo "Finished COUNTS"
cd ${WORKDIR}




## Remove list lines from HTseq count

CNTDIR="path_to_counts_directory"

cd ${CNTDIR}

for file in `ls *.count`
    do
        newname=`basename $file | sed -e "s/.bam.MAPQ.bam.count/.LastLinesRem.txt/"`
        head -n -5 $file > "$newname"
        echo $file
        echo $newname
    done

```

