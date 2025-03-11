#Thomas Tekle
#March 11, 2025
#Bioinformatics for Genomics - Assignment 2

#This analysis was conducted in the cedar under Dr. Lewis Luken's groups

#loading the sra-toolkit module
module load sra-toolkit 


#retrieving the SRA files regarding the study conducted by Mardanov et al (2020).
SRA_files=('SRR10551658' 'SRR10551660' 'SRR10551662' 'SRR10551664' 'SRR10551657' 'SRR10551659' 'SRR10551661' 'SRR10551663' 'SRR10551665')
for sra in "${SRA_files[@]}"; do
    prefetch "$sra"
done

#Retrieving the fastq files with respect to SRA files, using fasterq-dump tool
for i in *.sra; do base = $(basename "$i" .sra); fasterq-dump "$i" --outfile ${base}.fq --outdir ~/scratch/Assignment_2/files/FQ_files 

#Conducting a quality check for the reads used in this analysis
module load fastqc

cd ~/scratch/Assignment_2/files/QC_reports/pretrimmed

fastqc ../FQ_files/*.fq

unzip *fastqc.zip

multiqc ./ -o ./


#The mulitqc report shows poor quality, the reads were trimmed using the trimmomatic tool

module load trimmomatic

cd  ~/scratch/Assignment_2/files/FQ_files

for fq in *.fq; do base=$(basename "$fq" .fq);  java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar SE -threads 6 "$fq" ${base}_trimmed.fq ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 HEADCROP:20

#Conducting an additional quality check on the trimmed reads

cd ~/scratch/Assignment_2/files/QC_reports/trimmed

fastqc ../FQ_files/*._trimmed.fq

unzip *fastqc.zip

multiqc ./ -o ./

#Aligning reads to reference genome using STAR splice-aware aligner

cd ~/scratch/Assignment_2/files/Alignment

module load star

for fq in ../FQ_files/*_trimmed.fq; do base=$(basename "$fq" _trimmed.fq); STAR --genomeDir ../../Assignment_2_Genome/ --runThreadN 6 --readFileIn "$fq" --outFileNamePrefix ${base}; done

#Converting Sam files to Bam files then sorting those files

module load samtools

cd ~/scratch/Assignment_2/files/Bam_sorted

for i in ../Alignment/*.sam; do base=$(basename "$i" .sam); samtools view -b $i >  ${base}.bam; done

for i in *.bam; do base = $(basename "$i" .bam); do samtools sort "$i" > ${base}_sorted.bam; done

#Counting the number of reads that match to reference genome features using the featureCounts tool from the subread module

Bam_list=('SRR10551657_sorted.bam' 'SRR10551658_sorted.bam' 'SRR10551659_sorted.bam' 'SRR10551660_sorted.bam' 'SRR10551661_sorted.bam' 'SRR10551662_sorted.bam' 'SRR10551663_sorted.bam' 'SRR10551664_sorted.bam' 'SRR10551665_sorted.bam')

featureCounts -T 6 -a ../../Assignment_2_Genome/*.gtf -o count_output.txt ../Bam_sorted/"${Bam_list[@]}"

#The count_out.txt file was used for the remainder of the RNA-seq analysis in R.
