# Assemblies for newbies
## Steps

1. Assemble
2. Evaluate
3. Improve
4. Compare
5. 1.




## Context

Imagine there is a cholera outbreak in Krumlov. “Luckily” you have access to nearby sequencers and we produced datasets from bacterial isolates to sequence. So you have both a PacBio run and an Illumina run. Now the task is to assemble the data in order to do various analyses later: check the phylogeny of that strain, what makes it different from other strains in terms of gene content, SNPs, structural variants, etc.. And remember: don’t drink the water, it’s a known vector of contamination. Czech beer is perfectly safe though.

You are given several sequencing datasets from the same organism:  V. cholerae, which has a genome size of ~4 Mbp. The goal is to perform an assembly of the V. cholerae genome. It is known that the genome has 2 chromosomes and has around 3,800 annotated genes.

The following datasets are provided in the Workshop AMI (in ~/workshop_materials/assembly/) but can also be downloaded by using the ERRxxxxxx identifiers on the Sequence Read Archive. See how at the end of this tutorial.

PacBio sequencing (these are raw noisy “CLR” reads, not “CCS” or high fidelity  reads):
ERR1716491.fastq
Illumina sequencing (paired-end interleaved, insert size 150 bp):
SRR531199.fastq

The Part 1 consist to play a bit with the two datasets and have a rough idea of the characteristic of the available datasets.

The goal of Part 2 is to obtain a quick, initial assembly of the V. cholerae genome using any of these datasets.
You could use any of the assemblers which have been pre-installed on your Instance. “But, which one”, you ask? Well, the point of this part is to let you to take initiatives and pick one!
So, for this first attempt, if the assembler asks for any parameter, try to guess a reasonable value but do not over-think it. At this point, it does not matter if the assembly is of poor quality.
In the step you can try to work with a sub sampling of your datasets to accelerate the process.

In the Part 3 , you will measure the quality of your initial assembly and recognize that it could possibly be improved. Once you have generated your first assembly, move to Part 2. If you are stuck, don’t panic, either ask a TA or follow the detailed steps below.

## PART 1: Play with datasets

*How many bases is there in the datasets?*

*What is you mean read length?*

*Find your longest read!*

*What is the estimated coverage?*

*Make datasets containing 10% of the initial coverage to make "quick" assembly test.*

HINT: Use seqkit

My commands :
```
  seqkit stats ERR1716491.fasta
  seqkit sort -l ERR1716491.fasta -r --line-width 0 > LRS.fasta
  seqkit sample -p 0.1 ERR1716491.fasta --line-width 0 > LR_0.1.fa
```
## PART 2: Your first assembly

### Available assemblers
##### Long reads assemblers
1. Miniasm
2. Raven
3. Flye


1. Miniasm

Miniasm is a rather special long read assembler as it does not include a consensus step.
The resulting contigs are just merged erroneous long reads and still contains many sequencing errors.
Produced contigs are structurally correct, but at the nucleotide level, there are many many mismatches and indels.
For most application the contigs need to be polished. E.g. using the Racon software
or minipolish than use it aproprietly to correct the contigs errors.

  Miniasm Work in two steps:
  1. Find overlaps (minimap2)
  2. Second generate contigs (miniasm)
  3. Miniasm do not include a polishing step, you can try minipolish(https://github.com/rrwick/Minipolish)

  Website: https://github.com/lh3/miniasm
```
minimap2 -x ava-pb -t8 ERR1716491.fastq.gz ERR1716491.fastq.gz | gzip -1 > reads.paf.gz
 miniasm -f ERR1716491.fastq.gz reads.paf.gz > assembly.gfa
 minipolish -t 8 ERR1716491.fastq.gz assembly.gfa > polished.gfa
```

2. Raven

Raven is an improved version of Ra that is based on existing components: minimap, racon (polishing), but with layout stage.
```
raven LR_0.1.fa  --graphical-fragment-assembly raven_graph.gfa -t 8 > raven_contigs.fa
```
Website: https://github.com/lbcb-sci/raven


3. Flye

Flye is a long read assembler based on a original paradigm. It build a repeat graph that is conceptually similar to  De Bruijn graph but using approximate matches.
Flye is also able to assemble metagenomes.

```
 flye --pacbio-raw LR_0.1.fa --threads 8 --out-dir flye01  --genome-size 4000000
 ```
 Website: https://github.com/fenderglass/Flye


4. Canu

Canu is one of the reference assemblers for long reads data. It gives very good results, however in the context of this workshop it might take a while to run and
will not be used for this workshop.


 Website: https://github.com/marbl/canu

###### Short reads assemblers
1. Spades
2. MEGAHIT
3. Minia



1. Spades

  Spades is an Illumina assembler designed for prokaryotic and small eukaryotic genomes. It does an excellent job at assembling bacteria with short read sequencing, either multi-cell or single-cell data, and also small metagenomes. It uses multiple size of k to construct the best possible contigs. It takes generally longer time and memory than other assemblers. It can also improve its contigs using long reads.

  Website: https://github.com/ablab/spades
  ```
  spades.py --12 SRR531199.fastq -o work_dir
  spades.py --12 SRR531199.fastq -o work_dir --pacbio ERR1716491.fastq
  ```
  Spades can be quite long, you can use a single size of k to get results faster
  ```
  spades.py --12 SRR531199.fastq -o work_dir -k 33
  ```
   Tips: If you run the multiple K assembly, you can take a look at the successive contigs files produced from different value of K in the work_dir/K21/final_contigs.fasta,  work_dir/K33/final_contigs.fasta before the whole process is finished.

2. MEGAHIT

  MEGAHIT is an Illumina assembler designed to assemble large metagenomic experiment. It is very conservative and is able to assemble even low coverage regions. It is very fast and memory efficient despite the fact that it use several kmer size.

  Website: https://github.com/voutcn/megahit
  ```
  megahit --12 SRR531199.fastq -o out -t 8
  ```
  Like Spades, you can specify your own kmer sizes to use with
  --k-list  parameter
  ```
  megahit --12 SRR531199.fastq -o out -t 8 --k-list 21,41,61,81,101,121,141
```

3. Minia

  Minia is an Illumina assembler designed to be ressources efficient and able to assemble very large genomes.
  You can run the GATB pipeline that correct the reads (Bloocoo), generate contigs (Minia) and scafold them (BESST)
  ```
  ./gatb --12 interleaved_reads.fastq
  ```

  Website: https://github.com/GATB/gatb-minia-pipeline


### Other assemblers (long reads)

+ Canu (https://github.com/marbl/canu/commits/master)
+ Redbean (https://github.com/ruanjue/wtdbg2)
+ FALCON (https://github.com/PacificBiosciences/FALCON)
+ SMARTdenovo (https://github.com/ruanjue/smartdenovo)

### Other assemblers (short reads)

+ AByss (https://github.com/bcgsc/abyss)
+ Unicycler (https://github.com/rrwick/Unicycler)
+ Discovardenovo (https://software.broadinstitute.org/software/discovar/blog/)

## PART 3: Assembly Evaluation

### Evaluate your assembly with Quast
To Evaluate your assembly you will run  Quast on your contigs file.

 A (very nice) manual can be found here http://quast.bioinf.spbau.ru/manual.html.
```
./gatb quast.py -o output_directory assembly.fa
```

 Move into your QUAST output directory and examine the report.txt file. You may also take a look at the HTML file.

 Note only "large enough" contigs will be considered, contigs smaller than 500bp are ignored by Quast.

Now we will compare your contigs with the reference genome.
With the -r option, Quast will align your contigs on your reference genome and estimate their accuracy

```
./gatb quast.py -o output_directory assembly.fa -r vcholerae_h1.fasta
```


In real life you will likely not have a reference genome, but if you have a closely related genome you can use it as reference and get some nice stats from QUAST such as ‘genome fraction’ (=genome coverage).

*How many large / small misassemblies were made by the assembler?*

*What part of your genome is covered by your contigs?*

If your contigs contain too many errors quast will not align them on the reference. In such cases (ie miniasm unpolished assembly) you can use the --min-identity quast parameter.
```
./gatb quast.py -o output_directory assembly.fa -r vcholerae_h1.fasta  --min-identity 0.8
```
### Visualize your assembly graph (optionnal)

Use the Bandage software to visualize the assembly graph.
It is generally the file that ends with “.gfa”.

If you do not have a gfa files you may skip this step.
```
Bandage assembly.gfa
```

*Is your graph well connected?*

*How many disjoint component is there?*

This provides some indications of whether you had sufficient coverage to perform the assembly. It could also mean that the sequencer produced some sequences that did not overlap others.

*Does the graph has many branching nodes?*

This is indicative of variants or repetitions that could not have been resolved by the assembler.




### Annotate your assembly
You may use the tool Prokka to annotate the assembly.
```
prokka assembly_file
```

*How many annotated proteins were predicted by Prokka?*
```
wc -l output.faa
```



Across all reported assemblies in the Google Docs, what is the variability of the number of annotated genes?
### Call SNP (Optionnal)
Actually the reference genome, the Illumina and  Pacbio datasets aren’t from the same exact strain so they may differ one from another.
To visualize thoses differences, you may align your Illumina reads to the reference genome and/or a Pacbio assembly and call SNPs.
To do so you may use  samtools mpileup: http://samtools.sourceforge.net/mpileup.shtml.

*How many SNPs were found between the two strains?*

As a check, do the same procedure by aligning reads from the Illumina dataset to an Illumina-only assembly.

*Does SNP calling of a dataset on its own assembly report “false positive” SNPs?*

## Part 4 Assembly comparison

*Report your assembly on google doc and compare with the results from the other participants.*


At this point, you may be tempted to re-run your assembly with better parameters, with a higher coverage or to test another assembler. This is the time to do so!
