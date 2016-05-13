Assembly Tutorial
=================

### Dependencies 
- Kraken (with local download of minikrakenDB)
- Samtools 1.2
- BWA
- SPAdes
- Prokka
- Mauve
- Tablet 
- Fastqc
- Quast
- Bandage
 
# Summary

For this tutorial we will be rigourously QCing and assembling a bacterial genome. 

The example genomes we will be using are *Renibacterium salmoninarum* isolates taken from this publication - [Microevolution of *Renibacterium salmoninarum*: evidence for intercontinental dissemination associated with fish movements](http://www.nature.com/ismej/journal/v8/n4/full/ismej2013186a.html)

We will be using the finished reference genome ATCC_33209 from - [Genome sequence of the fish pathogen *Renibacterium salmoninarum* suggests reductive evolution away from an environmental Arthrobacter ancestor.](http://www.ncbi.nlm.nih.gov/pubmed/18723615). This has been reannottated for you with Prokka. 

The workflow will be as follows:
- Taxonomic assignment of reads.
- QC of raw reads including adapter removal and trimming and filtering.
- Removal of known contaminants.
- Assembly using SPAdes.
- Contig QC including length filtering, read coverage filtering and taxonomic identification. 
- Contig annotation, visualisation and QC.


# Assembly and Annotation

For this exercise we will assemble a bacterial genome and rigorously QC it. 

Each of you will be provided with a selection of *Renibacterium* genomes to assemble and annotate. 

You will have copied the files for the project into your work directory in the previous session. Each group will be assigned a sample to assemble from the Wednesday18thMay/Renibacterium/Fastq/ folder.


| Group Number | Student A | Student B | Student C |
| ------------ | ----------| --------- | --------- |
| 1 | Jenni Lehtimäki | Katariina Pärnänen | Tero Tuomivirta | 
| 2 | Jing Cheng | Lijuan Yan | Kirsi Hyytiäinen |
|3|  Alejandra Culebro  |  Rouger Amelie | Mohammad Jaber Alipour|
|4|      Julija Svirskaite     |                    Hanna Castro              |                Anniina Jaakkonen|
|5|      Seyed Abdollah Mousavi  |      Sajan Raju                             |       Outi-Maaria Sietiö|
|6|      Noora Ottman             |            Ann-Katrin Llarena               |      Minna Santalahti|
|7|      Egle   Kudirkiene         |             Jani Halkilahti |


*Make a directory for your sample of interest and download/move the fastq.gz files to this directory.* 

## Read QC

The first stage in assembling a genome is to assess the quality of your fastq file (raw reads). 

I will be using the example of a file called ERR327970. **You will have to adapt the code for your problems**. 

It is important to remember - *'rubbish in, rubbish out'*. 

There are a number of factors to considered at this stage. Reads will contain many confounding factors that can cause missassembly. These include:
* Sequencing adapters/primers from the library preparation stage.
* phiX calibration spike DNA.
* Poor quality sequence. 
* DNA from other samples from the same run. 
* Overabundant k-mers.

We will keep things tidy, so make a directory for the Read QC. 

```
mkdir Read_QC
cd Read_QC
```
Copy your reads to this directory 

#### FastQC

We need to be able to assess the quality of our fastq file. We will use a program called Fastqc.  

Fastqc can be run with a user interface or programatically. 

If you want to use the UI then type fastqc and load the files you are interested in, you must have X11 forwarding active if you are working remotely. 

If you programmatically want to generate a report:

```
fastqc -t 2 ERR327970.fastq.gz ERR327970.fastq.gz 
```

This will generate a html report for each file.

Open these files and take note of the overrepresented sequences, Illumina adapter content and k-mer content fields. Also, what do you notice about the quality of the forward and reverse reads?

We can see that there has been sequence identified as Illumina Single End PCR Primer 1 present in the reads. 

#### Trimomatic 
Our reads have adapter sequence. Lets tackle that and poor quality sequence first.

This can be removed using a number of programs. Popular choices include sickle, cutadapt and trimmomatic. For this exercise I have chosen trimmomatic as it is fast, flexible and has easily understandable syntax. 

*Note that trimmomatic will run your commands sequentially so if you filter your reads by minimum length before trimming adapter sequence you will get a different result than the other way around.*

```
trimmomatic PE ERR327970_1.fastq.gz ERR327970_2.fastq.gz ERR327970_1.paired.fastq.gz ERR327970_1.unpaired.fastq.gz ERR327970_2.paired.fastq.gz ERR327970_2.unpaired.fastq.gz ILLUMINACLIP:../PE_All.fasta:2:30:10 LEADING:20 TRAILING:3 SLIDINGWINDOW:4:20 CROP:94 MINLEN:30
```
You will be left with paired and unpaired reads for both forward and reverse files. 

[Trimmomatic manual](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf)

So what I have I done here:
- Our reads are paired reads so I used PE - pair end mode.
- LEADING and TRAILING : commands are commonly used to remove poor quality seqience from the beggining and end of reads (<phred score 2)
- ILLUMINACLIP: Removes adapter sequences specified by the user. 
- SLIDINGWINDOW : Removes bases below this quality in a moving window from the end --> beggining of the read.
- CROP: I noticed overabundant K-mers at the end of reads which may confound assembly, this crops reads to this length.
- MINLEN: Reads below this side are useless for assembly as we will not be using k-mers below this size. 

*Now check your paired reads again using fastqc.*

#### K-mer correction
We will not be performing any k-mer correction here. However, be aware that this is a possible step during read filtering and may improve assembly quality. Possible software includes Quake and Musket.

#### Taxonomic assignment using Kraken
We can assign the raw reads to a database of genomes in order to see if they are what we think they are. 

Kraken will try to assign an each reads to a sequence from a database (in this case the MinKrakenDB provided with the software, a reduced version of the NCBI reference genome collection)

```
# If it is your first time using/installig kraken then you need to follow the next few steps
# Navigate to a directory you want to install the kraken scripts.
mkdir kraken && cd kraken
wget ccb.jhu.edu/software/kraken/dl/kraken-0.10.5-beta.tgz && tar -vxzf kraken-0.10.5-beta.tgz 
cd kraken-0.10.5-beta
./install_kraken.sh ../
cd ../ && wget ccb.jhu.edu/software/kraken/dl/minikraken.tgz && tar -vxzf minikraken.tgz # MinkrakenDB

export PATH=$PATH:./ # May not work

# When I call kraken in the commands below you will need to enter the path to the correct executable in the ../kraken/bin directory. The same goes for the MiniKrakenDB (directory not file)

```

```
# Preload database and assign reads 
kraken --preload --db /mnt/data/bioinformatics/Databases/MiniKraken/ ERR327970_1.paired.fastq.gz ERR327970_2.paired.fastq.gz --threads 6 --paired --output kraken_result 

# Convert IDs to usable taxonomic labels.
kraken-translate --db /mnt/data/bioinformatics/Databases/MiniKraken/ kraken_result > sequences.labels

# Produce  a summary of the kraken output 
awk '{print $4}' sequences.labels | sort | uniq -c | sort -nr
```

Next to Renibacterium we have a large number of reads classified as phage in our sample. Open the file and see what these are. 

# Remove Contaminants
From the kraken output it is clear we have some phiX carryover from the sequencing. This is used as a calibration spike in sequencing runs. We will want to remove these reads before assembly. 
```
# Index reference fasta file
bwa index ../phiX.fasta

# Map your reads to the phiX genome 
bwa mem -t 6 ../phiX.fasta ERR327970_1.paired.fastq.gz ERR327970_2.paired.fastq.gz > ERR327970.phiX_mapping.sam

# Generate some mapping stats
samtools flagstat ERR327970.phiX_mapping.sam > contam.stats
grep "^SN" contam.stats

# Extract only reads that DID NOT map to the phiX (-f 12) 
samtools view -@ 6 -f 12 -T ../phiX.fasta -bSu ERR327970.phiX_mapping.sam | samtools sort -@ 6 -T temp.bwa -O bam - | samtools bam2fq -s temp.fq - > ERR327970.phiX_filtered.fastq

# Filter your reads to seperate files. 
grep -A3 "^@.*/1$" ERR327970.phiX_filtered.fastq > ERR327970_1.filtered.fastq
grep -A3 "^@.*/2$" ERR327970.phiX_filtered.fastq > ERR327970_2.filtered.fastq
gzip ERR327970_1.filtered.fastq ERR327970_2.filtered.fastq # Compress them to save space.
```
Now check you have an identical number of reads in each file.

Re-run kraken and analyse the output. 

We still have a reasonable amount of Vibrio DNA. We will deal with that after assembly.  

#### Assembly 

We will assemble our genomes with SPAdes, arguably the best assembler for bacterial genomes available. 

```
cd ../
mkdir Assembly

cp Read_QC/ERR327970_1.filtered.fastq.gz ./
cp Read_QC/ERR327970_2.filtered.fastq.gz ./

# Assemble using spades
spades.py -t 6 --pe1-1 ERR327970_1.filtered.fastq.gz --pe1-2 ERR327970_2.filtered.fastq.gz --careful -o Assembly/
```

SPAdes was run using 6 cores (-t 6) on paired end data (--pe1-1/--pe-2) with error correction on (--careful) which correction the contigs by remapping the reads to output folder Assembly (-o)

In the results folder we have a number of important files. 

Fistly we will look at the contig.fasta and scaffolds.fasta file.
* contigs.fasta - this contains the contigs (contiguous sequences) as generated by SPAdes  
* scaffolds.fasta - the contigs generated by SPAdes after repeat resolution and scaffolding using the pair end information of the PE reads. 
* assembly_graph.fastg - the assembly graph generate by SPAdes. This can be used to visualise the assembly. 

We can compare the contigs and scaffolds quantitative way using a program called Quast. 

```
mkdir Contig_QC

quast.py -t 6 -R ATCC_33209.fasta Assembly/contigs.fasta Assembly/scaffolds.fasta -o Contig_QC/RawContigs
```

Quast generates a number of assembly statistic quickly and compares the assemblies to the reference genome (-R). It can be run without a reference genome to generate just the assembly stats. 

#### Assembly QC

*** Explain N50 and where files are 
In some situations (often if there is residual adapter or phiX sequence) the scaffold.fasta file will contain a greater number of misassemblies than contig.fasta. 
However, the scaffold.fasta will often have a larger N50 or fewer contigs and therefore is 'better' in that respect. 
In my view fewer misasseblies are preferred so we will use the contigs.fasta. 

The header of each contig will have the format NODE_X_length_Y_cov_Z 
-The node number is arbitary but they will often be assorted in decending size order.
- length in bp
- k-mer coverage (this is the coverage of the last K-mer used by spades and does NOT represent read coverage but can often be used as a proxy). 


#### Taxonomic assignment of raw contigs 

Again, we will use kraken:
```
kraken --preload --db /mnt/data/bioinformatics/Databases/MiniKraken/ Assembly/contigs.fasta --threads 6 --output kraken_result_contigs
kraken-translate --db /mnt/data/bioinformatics/Databases/MiniKraken/ kraken_result_contigs > sequence_contigs.labels
awk '{print $4}' sequence_contigs.labels | sort | uniq -c | sort -nr # Summary (Save with > file.txt appended to end)
```

Contamination of sample with even a small amount of exogenous DNA is an ever present problem in WGS. [Example]
Very sensitive assemblers such as SPAdes, which assembles over a range of k-mer and read coverages it can be especially problematic. 
However, with a little post processing of the contigs we can usually deal with small to moderate amounts of contamination.

#### Visualise assembly graph with bandage 

Bandage allows us to visualise the assembly graph. 


This can be used to visualise connected sequences that nevertheless get labelled as separate contigs because their position cannot be resolved by the assembler. 

Bandage can also be used to resolve assembly error, insertion of large tracks of DNA and DNA repeats. It is a useful tool. 

```
Bandage load assembly_graph.fastg
click 'draw graph'
```

We can see that there is a large amount of 'connected' contigs, which most likely represents our genome. There is a lot of fragmented sequence, which is most likely contamination. We will need to remove this in the following steps.

## Filtering Contigs 



#### Filtering on Read Coverage

It is often useful to map reads to contigs to get a good idea of the coverage of the original reads to the resulting contigs.
We can generate some useful statistics including how many reads map back the the resulting contig, how may are paired correctly etc.  

```
mkdir Mapping

# Copy fasta file to working directory and make them into one line multifastas (this is useful for some downstream applications )
awk '/^>/{print (NR==1)?$0:"\n"$0;next}{printf "%s", $0}END{print ""}' < Assembly/contigs.fasta > Mapping/contigs.fasta # Make working fasta one line fasta 

cd Mapping

# Index the contigs 
bwa index contigs.fasta # index contig sequences

# Map and sort the reads 
bwa mem -t 6 contigs.fasta ../ERR327970_1.filtered.fastq.gz ../ERR327970_2.filtered.fastq.gz | samtools view -@ 6 -T contigs.fasta -bS - | samtools sort -@ 6 -T contigs.bwa -o contigs.bwa.bam -

# Generate mapping and coverage stats. 
samtools stats contigs.bwa.bam > contigs.stats # Mapping stats
samtools depth contigs.bwa.bam > contigs.coverage # Coverage per base 
	
# Summarise 
grep "^SN" contigs.stats
```
How many of our filtered reads mapped to our assembly? How many of them were correctly paired? 

Now lets calculate the average coverage of reads across our whole assembly. 

```
awk '{cnt+=$3; n+=1} END{ if (n > 0){ print cnt/n } else { print "0" } }' < contigs.coverage
```
We need to identify an acceptable lower bound for coverage. Lets set 10% of our estimated average coverage. 

All contigs below this threshold we can consider are likely to be contaminants. **Beware this is threshold is entirely arbitrary, it will depend on the data/sample**
```
THR_COV=9

# Calculate average coverage per contig 
awk '{a[$1]+=$3;++c[$1]}END{for(i in a)printf "%s\t%.1f\n", i, a[i]/c[i]}' < contigs.coverage > contigs.contig_coverage 
		
# Filter per contigs on coverage threshold.
cat contigs.coverage | awk '{print $1}' | uniq > contigs.contig.list # get all contig IDs
echo -n "" > contigs.c_filtered.fasta
echo -n "" > temp.count
cat contigs.contig_coverage | while read cID av_cov; do pf=$(echo "$av_cov>$THR_COV" | bc); if [ $pf == 1 ]; then seq=$(grep -m1 -A1 $cID contigs.fasta); printf "%s\n" $seq >> contigs.c_filtered.fasta ; else N=$(cat temp.count);N=$((N+1));echo $N > temp.count ; fi ; done
c_removed=$(cat temp.count); echo "Number of removed contigs = $c_removed"
```
How many contigs did we remove? How many remain? 

*Run kraken on our filtered set of contigs*

#### Length Filter
It looks like we still have some contamination. Lets remove the very small contigs. They are not contributing in a meaningful way to our analysis and may just be k-mer carry-over from SPAdes. We will set a threshold of 500 bp. Quast ignores contigs below this size anyway. 

```
# Count Contigs below threshold
awk -v L=500 'BEGIN {OFS = "\n"} {header = $0 ; getline seq ; if (length(seq) < L){ count++ }} END{ print count }'  < contigs.c_filtered.fasta 

# Remove contigs below length threshold 
awk -v L=500 'BEGIN {OFS = "\n"} {header = $0 ; getline seq ; if (length(seq) > L) {print header "\n" seq }}' < contigs.c_filtered.fasta > contigs.c_l_filtered.fasta 
```
How many contigs remain?
```
grep "^>" contigs.c_l_filtered.fasta | wc -l
```
What species do these contigs represent? *Run Kraken again*
```
kraken --preload --db /mnt/data/bioinformatics/Databases/MiniKraken/ contigs.c_l_filtered.fasta --threads 6 --output kraken_result_contigs_postfilter
kraken-translate --db /mnt/data/bioinformatics/Databases/MiniKraken/ kraken_result_contigs_postfilter > sequence_contigs_postfilter.labels
awk '{print $4}' sequence_contigs_postfilter.labels | sort | uniq -c | sort -nr # Summary (Save with > file.txt appended to end)
```

Our assembly should be looking much better at this point. If it isn't then we may need to use stricted thresholds or filter our reads to a better quality. 

#### Adding Annotation
Having a 'final' assembly is good, but having it annotated is even better for comparison. This will let us know what parts of the genome are missing or different between our strains and reference strains. 

For this purpose we will use Prokka, which annotates our strain in a few minutes from reference databases.  
```
cd ../

awk '/^>/{print ">" ++i; next}{print}' < Mapping/contigs.fasta > ERR327970.spades_contigs.fasta  # avoids a prokka error with long file names. 
prokka --outdir Annotation --force --cpus 8 --addgenes --prefix ERR327970 --genus Renibacterium --species salmoninarum ERR327970.spades_contigs.fasta
```
This should take a few minutes and will produce a folder in the output directory containing annotation information. We willbe using the .gff and.gbk files for further QC. 

#### Comparing the contigs to a reference. 
Lets run quast again and see what our 'final' assembly compares to the unfiltered one. 
```
quast.py -t 6 -o Contig_QC/Final -R ATCC_33209.fasta -G Annotation/ATCC-33209.gff ERR327970.spades_contigs.fasta
```
The use of -G *.ggf  also allows for the comparison of our annotated genome to the annotation of the reference. Which CDSs were not detected in our sample? Is there a trend in the missing genes?

We can also visualise the contigs relative to a reference to aid comparison. 

This can allow for the visualisation of large genomic structural variants such as translocations, insertions and deletions. It also lets us see if the genome we have assembled looks sane, i.e. do we largely see that contigs are syntenic (in the same order) as in a closely related reference genome?  

This can be performed using a number of tools (Artemis/IGV/Mauve/ACT/Tablet). We will use an alignment toolbox called Mauve:

```
# If Mauve is not already installed 
wget darlinglab.org/mauve/snapshots/2015/2015-02-13/linux-x64/mauve_linux_snapshot_2015-02-13.tar.gz
tar -vxzf mauve_linux_snapshot_2015-02-13.tar.gz

# Use th path of the directory in which you just downloaded the jar file in the command below. 

java -Djava.awt.headless=true -Xmx8000m -cp /opt/mauve/mauve_snapshot_2015-02-13/Mauve.jar org.gel.mauve.contigs.ContigOrderer -output "Ordered_contigs" -ref ../ATCC_33209.gb -draft contigs.c_l_filtered.fasta
```
Find the highest iterated in the Ordered_contigs folder. 
```
progressiveMauve --output=alignment ../ATCC_33209.gb Ordered_contigs/alignment3/contigs.c_l_filtered.fasta

# Visualise
mauve alignment

OR

java -jar ~/Path2/mauve_snapshot_2015-02-13/Mauve.jar alignment
```

The window shows the alignment between our contigs and the reference  genome. Zoom in on the alignment sections where there is no alignments between our contigs and the reference. Is there any commonality between these regions? Do any contigs not align? Blast them and see what they are. 

### Map the reads to get final statistics 

Now we have our 'final' assembly we want to remap our reads to it to generate some final statistics. 
```
mkdir ../Mapping_Final
cp Ordered_contigs/alignment3/contigs.c_l_filtered.fasta ../Mapping_Final/contigs.fasta

cd ../Mapping_Final/
bwa index contigs.fasta
bwa mem -t 6 contigs.fasta ../ERR327970_1.paired.fastq.gz ../ERR327970_2.paired.fastq.gz | samtools view -@ 6 -T contigs.fasta -bS - | samtools sort -@ 6 -T contigs.bwa -o contigs.bwa.bam -
	
samtools stats contigs.bwa.bam > contigs.stats # save stats
samtools depth contigs.bwa.bam > contigs.coverage # get per base coverage 
	
# Get some mapping stats
grep "^SN" contigs.stats

samtools index contigs.bwa.bam
tablet contigs.bwa.bam contigs.fasta
```

You might want to download and run tablet locally, otherwise is it likely to be prohibitively slow. 

This and other similar software allows us to visualise the read pile-up on top of our contigs. Are there any notable features in your Mauve alignment that you can find an explanation for in the Tablet window?
