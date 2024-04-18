**PRE-PROCESSING**

**Reads QC**
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

**Reads Trimming**
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

**Reads Alignment**
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

**BAM Filtering**
```{bash}

```