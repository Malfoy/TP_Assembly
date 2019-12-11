# Assembly for newbies
## Steps

1. Assemble
2. Evaluate
3. Improve
4. Compare
5. 1.


## Your first Assembly!

#### Context

#### Play with datasets

How many bases is there in the datasets?

What is you mean read length?

Find your longest read!

What is the estimated coverage?

Make datasets containing 10% of the initial coverage to make "quick" assembly test.

HINT: Use seqkit

My commands :
```
  seqkit stats ERR1716491.fasta
  seqkit sort -l ERR1716491.fasta -r --line-width 0 > LRS.fasta
  seqkit sample -p 0.1 ERR1716491.fasta --line-width 0 > LR_0.1.fa
```

#### Installed assemblers
###### Long reads assembler
1. Miniasm (https://github.com/lh3/miniasm)
2. Raven (https://github.com/lbcb-sci/raven)
3. Flye (https://github.com/fenderglass/Flye)


1. Miniasm
Work in two steps:
  1. Find overlaps (minimap2)
  2. Second generate contigs (miniasm)
  3. Miniasm do not include a polishing step, you can try minipolish(https://github.com/rrwick/Minipolish)
```
minimap2 -x ava-pb -t8 ERR1716491.fastq.gz ERR1716491.fastq.gz | gzip -1 > reads.paf.gz
 miniasm -f ERR1716491.fastq.gz reads.paf.gz > assembly.gfa
 minipolish -t 8 ERR1716491.fastq.gz assembly.gfa > polished.gfa
```

2. Raven
```
raven LR_0.1.fa  --graphical-fragment-assembly raven_graph.gfa -t 8 > raven_contigs.fa
```

3. Flye
```
 flye --pacbio-raw LR_0.1.fa --threads 8 --out-dir flye01  --genome-size 4000000
 ```

###### Short reads assembler
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

+ Canu
+ Redbean
+ FALCON
+ SMARTdenovo

### Other assemblers (short reads)

+ AByss
+ Unicycler
+ Discovardenovo

## Assembly Evaluation

Use Quast. A manual can be found here http://quast.bioinf.spbau.ru/manual.html.
```
./gatb quast.py -o output_directory assembly.fa
```

Now use Quast to align your contigs on your reference genome and estimate their accuracy

```
./gatb quast.py -o output_directory assembly.fa -r vcholerae_h1.fasta
```
How many large / small misassemblies were made by the assembler?

What part of your genome is covered by your contigs?

If your contigs contain too many errors quast will not align them on the reference. In such cases (ie miniasm unpolished assembly) you can use the --min-identity quast parameter.
```
./gatb quast.py -o output_directory assembly.fa -r vcholerae_h1.fasta  --min-identity 0.8
```

Use the Bandage software to visualize the assembly graph.
It is generally the file that ends with “.gfa”
```
Bandage assembly.gfa
```

Is your graph well connected?

How many disjoint component is there?

Does the graph has many branching nodes?


Annotate the assembly
You may use the tool Prokka to annotate the assembly.
```
prokka assembly_file
```

How many annotated proteins were predicted by Prokka?
```
wc -l output.faa
```
Across all reported assemblies in the Google Docs, what is the variability of the number of annotated genes?


## Assembly comparison

Report your assembly on google doc and compare with the results from the other participants.


At this point, you may be tempted to re-run your assembly with better parameters, with a higher coverage or to test another assembler. This is the time to do so!

## Improve your assembly
