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

How many nucleotides is there in the datasets ?

What is you mean length?

Find your longest read!

What is the estimated coverage ?

HINT: Use seqkit

My commands :
```
  seqkit sort -l ERR1716491.fasta -r --line-width 0 > LRS.fasta
  ```

#### Installed assemblers
###### Long reads assembler
1. Miniasm
2. Raven
3. Flye


###### Short reads assembler
1. Spades
2. MEGAHIT
3. Unicycler


### Other assemblers (long reads)

+ Canu
+ Redbean
+ FALCON
+ SMARTdenovo

### Other assemblers (short reads)

+ AByss
+ Minia
+ Discovardenovo

## Assembly Evaluation

Use quast

## Improve your assembly

## Assembly comparison
