---
title: "Explore Alteromonas assembly strategy 3"

author: Nina Dombrowski
date-modified: last-modified
date: 05/03/2024

engine: knitr 

execute:
  eval: false
format:
  html:
    toc: true
    embed-resources: true

bibliography: references.bib
csl: american-society-for-microbiology.csl
---

## Description

Workflow to assemble the Alteromonas J4 genome.

## Software installation

Software versions:

-   Porechop v0.2.4, [link to github](https://github.com/rrwick/Porechop), no paper to cite
-   Filtlong v0.2.1, [link to github](https://github.com/rrwick/Filtlong), no paper to cite
-   NanoPlot v1.42.0, [@decoster2023]
-   SeqKit 2.7.0, [@shen]
-   Trycycler v0.5.5, [@wick2021]
    -   Flye v2.9.3 [@kolmogorov2019]
    -   Miniasm v0.3-r179 [@li2016]
        -   Minipolish 0.1.3 [@wick2019]
        -   Racon 1.5.0 [@vaser2017]
        -   Minimap2 2.28-r1209 [@li2018]
        -   any2fasta 0.4.2, [link to github](https://github.com/tseemann/any2fasta), no paper to cite
    -   Raven 1.8.3 [@vaser2021]
    -   Unicycer 0.5.0 [@wick2017]
-   Qualimap 2.3.0 [@garcía-alcalde2012]
-   Quast v5.2.0 [@quast2013]
-   CheckM2 v1.0.2 [@chklovskil]
-   Prokka 1.14.6 [@seemann2014]
-   Pseudofinder 1.1.0 [@syberg-olsen2022]
-   Homopolish 0.4.1 [@huang2021]
-   Gtdbtk v2.3.2 [@chaumeil2019] together with the GTDB r214 database [@parks2020]
- Circlator 1.5.5


## Set working directory

```{bash}
wdir="/zfs/omics/projects/sargo/j4_assembly_analysis"
cd $wdir
```

## Filter reads

-   To remove adapters, porechop was run with `--discard_middle` on the file provided by GM called Nanopore_J4_SuperBasecalling.fastq.gz (as described in the notes for the first trial, in the notebook: 01_chopper_flye.qmd)
-   After removing the adapters, Filtlong was used on the output of the porechop analysis to remove the 5% worst reads and reads shorter than 1000 bp (settings: `--min_length 1000 --keep_percent 95`)

```{bash}
mkdir -p results/v3_trycycler/filtlong/p95

conda activate filtlong_0.2.1

#seq and quality have different sequence lengths with and without `--mean_q_weight 10`
filtlong --min_length 1000 --keep_percent 95 results/v1_chopper_flye/filtered_reads/01_porechop/Nanopore_J4_SuperBasecalling_ad.fastq.gz | gzip > results/v3_trycycler/filtlong/p95/Nanopore_J4_filtlong.fastq.gz

seqkit stats -e -a -To results/v3_trycycler/filtlong/p95/seqkit_stats.tsv results/v3_trycycler/filtlong/p95/Nanopore_J4_filtlong.fastq.gz

conda deactivate
```

## Assemble with Trycycler

### Subsample reads

Notice:

-   During the first assembly approach I used jellyfish + genomescope to estimate the genome size (notes on how this was done are in 01_chopper_flye.qmd), which is a bit smaller than the assembled flye genome. But according to the trycycler documentation, \~10% deviation is fine.

```{bash}
mkdir -p results/v3_trycycler/subsampling

mamba activate trycycler_0.5.5

trycycler subsample \
    --reads results/v3_trycycler/filtlong/p95/Nanopore_J4_filtlong.fastq.gz \
    --out_dir results/v3_trycycler/subsampling \
    --count 12 \
    --genome_size 4.2m

```

**Comments**

-   subset_depth = 54.7x
-   reads per subset: 72,629

### Run Assemblies

Trycycler was run 3x each using the following assemblers:

-   Raven
-   Unicycler
-   Flye
-   Minipolish

Other assemblers to consider for future runs are: Canu, NECAT and NextDenovo/NextPolish (not tested)

```{bash}
mkdir assemblies

#submit all 12 assemblies: 48169
sbatch scripts/trycycler.sh

#check number of contigs: 
grep -c ">" assemblies/*fasta

#cleanup 
mv assemblies results/v3_trycycler/
```

```{bash}
#| code-fold: true
#| code-summary: "Code for `scripts/trycycler.sh`, click to see"

#!/bin/bash
#SBATCH --job-name=Trycycler
#SBATCH --output=logs/Trycycler_%j.out
#SBATCH --error=logs/Trycycler_%j.err
#SBATCH --cpus-per-task=30
#SBATCH --mem=200G
#SBATCH --time=48:00:00

#activate env
source ~/.bashrc
mamba activate trycycler_0.5.5

threads=30

#first set
flye --nano-hq results/v3_trycycler/subsampling/sample_01.fastq --threads "$threads" --out-dir assembly_01 && cp assembly_01/assembly.fasta assemblies/assembly_01.fasta && cp assembly_01/assembly_graph.gfa assemblies/assembly_01.gfa && rm -r assembly_01

miniasm_and_minipolish.sh results/v3_trycycler/subsampling/sample_02.fastq "$threads" > assemblies/assembly_02.gfa && any2fasta assemblies/assembly_02.gfa > assemblies/assembly_02.fasta

raven --threads "$threads" --disable-checkpoints --graphical-fragment-assembly assemblies/assembly_03.gfa results/v3_trycycler/subsampling/sample_03.fastq > assemblies/assembly_03.fasta

unicycler --long results/v3_trycycler/subsampling//sample_04.fastq --threads "$threads" --out assembly_04 && cp assembly_04/assembly.fasta assemblies/assembly_04.fasta && rm -r assembly_04


#second set
flye --nano-hq results/v3_trycycler/subsampling/sample_05.fastq --threads "$threads" --out-dir assembly_05 && cp assembly_05/assembly.fasta assemblies/assembly_05.fasta && cp assembly_05/assembly_graph.gfa assemblies/assembly_05.gfa && rm -r assembly_05

miniasm_and_minipolish.sh results/v3_trycycler/subsampling/sample_06.fastq "$threads" > assemblies/assembly_06.gfa && any2fasta assemblies/assembly_06.gfa > assemblies/assembly_06.fasta

raven --threads "$threads" --disable-checkpoints --graphical-fragment-assembly assemblies/assembly_07.gfa results/v3_trycycler/subsampling/sample_07.fastq > assemblies/assembly_07.fasta

unicycler --long results/v3_trycycler/subsampling//sample_08.fastq --threads "$threads" --out assembly_08 && cp assembly_08/assembly.fasta assemblies/assembly_08.fasta && rm -r assembly_08


#third set
flye --nano-hq results/v3_trycycler/subsampling/sample_09.fastq --threads "$threads" --out-dir assembly_09 && cp assembly_09/assembly.fasta assemblies/assembly_09.fasta && cp assembly_09/assembly_graph.gfa assemblies/assembly_09.gfa && rm -r assembly_09

miniasm_and_minipolish.sh results/v3_trycycler/subsampling/sample_10.fastq "$threads" > assemblies/assembly_10.gfa && any2fasta assemblies/assembly_10.gfa > assemblies/assembly_10.fasta

raven --threads "$threads" --disable-checkpoints --graphical-fragment-assembly assemblies/assembly_11.gfa results/v3_trycycler/subsampling/sample_11.fastq > assemblies/assembly_11.fasta

unicycler --long results/v3_trycycler/subsampling//sample_12.fastq --threads "$threads" --out assembly_12 && cp assembly_12/assembly.fasta assemblies/assembly_12.fasta && rm -r assembly_12
```

**Comments**

Contigs per assembly:

-   assemblies/assembly_01.fasta:1
-   assemblies/assembly_02.fasta:1
-   assemblies/assembly_03.fasta:1
-   assemblies/assembly_04.fasta:1
-   assemblies/assembly_05.fasta:1
-   assemblies/assembly_06.fasta:1
-   assemblies/assembly_07.fasta:1
-   assemblies/assembly_08.fasta:1
-   assemblies/assembly_09.fasta:1
-   assemblies/assembly_10.fasta:1
-   assemblies/assembly_11.fasta:2
-   assemblies/assembly_12.fasta:1

### Cluster contigs

```{bash}
mkdir results/v3_trycycler/trycycler

trycycler cluster \
    --assemblies results/v3_trycycler/assemblies/*.fasta \
    --reads results/v3_trycycler/filtlong/p95/Nanopore_J4_filtlong.fastq.gz \
    --out_dir results/v3_trycycler/trycycler

#check the assignment of assemblies per cluster:
ls results/v3_trycycler/trycycler/cluster_00*/1*/*fasta
```

**Comments**

Distribution of assemblies per cluster:

Cluster1:

-   results/v3_trycycler/trycycler/cluster_001/1_contigs/A_contig_1.fasta
-   results/v3_trycycler/trycycler/cluster_001/1_contigs/B_utg000001c.fasta
-   results/v3_trycycler/trycycler/cluster_001/1_contigs/D_1.fasta
-   results/v3_trycycler/trycycler/cluster_001/1_contigs/E_contig_1.fasta
-   results/v3_trycycler/trycycler/cluster_001/1_contigs/F_utg000001c.fasta
-   results/v3_trycycler/trycycler/cluster_001/1_contigs/G_Utg1404.fasta
-   results/v3_trycycler/trycycler/cluster_001/1_contigs/H_1.fasta
-   results/v3_trycycler/trycycler/cluster_001/1_contigs/I_contig_1.fasta
-   results/v3_trycycler/trycycler/cluster_001/1_contigs/J_utg000001c.fasta
-   results/v3_trycycler/trycycler/cluster_001/1_contigs/K_Utg1354.fasta
-   results/v3_trycycler/trycycler/cluster_001/1_contigs/L_1.fasta

Cluster2:

-   results/v3_trycycler/trycycler/cluster_002/1_contigs/C_Utg1376.fasta

Cluster3:

-   results/v3_trycycler/trycycler/cluster_003/1_contigs/K_Utg1356.fasta

After running trycycler cluster the newick file was downloaded and explored in FigTree.

-   Assembly cluster_3_K_Utg1356_144930_bp_118.8x was a clear outlier (even though the node was not extremely long) and belonged to the only assembly with 2 contigs
-   Assembly cluster_2_C_Utg1376_4497983_bp_119.0x showed a bit a distance between the remaining 10 assemblies and could also be considered for removal
-   We will proceed with the assemblies in cluster1 (and thus not further work with the 2 assemblies described above).

### Reconcile contigs

```{bash}
trycycler reconcile \
    --reads results/v3_trycycler/filtlong/p95/Nanopore_J4_filtlong.fastq.gz \
     --cluster_dir results/v3_trycycler/trycycler/cluster_001

trycycler dotplot --cluster_dir results/v3_trycycler/trycycler/cluster_001
```

After first run, we got an error "Error: the following sequences failed to circularise: K_Utg1354."

I decided to remove this assembly to proceed with analysis (after checking the Trycycler documentation for recommendations):

```{bash}
mv results/v3_trycycler/trycycler/cluster_001/1_contigs/K_Utg1354.fasta results/v3_trycycler/trycycler/cluster_001/1_contigs/K_Utg1354.fasta.bad 

trycycler reconcile \
    --reads results/v3_trycycler/filtlong/p95/Nanopore_J4_filtlong.fastq.gz \
     --cluster_dir results/v3_trycycler/trycycler/cluster_001

```

**Comments**

-   Trycycler was unable to find a suitable known starting sequence (it looks for dnaA)
-   Trycycler uses the edlib aligner to get global alignments between all pairs of sequences. The manual gives the following advice: This can help you to spot any problematic sequences that should be excluded before continuing. If you see any sequences with notably worse identities, you can remove them.
-   Based on this advice , I chose to remove the following assemblies due to a low identity around 50. Notice, that the rest of assemblies have identities ranging from 70-95, so more subsampling + additional assemblers could be used to push the number of final assemblies that can be reconciled
    -   D_1
    -   G_Utg1404

```{bash}
#cleanup of assemblies before proceeding
mv results/v3_trycycler/trycycler/cluster_001/1_contigs/D_1.fasta results/v3_trycycler/trycycler/cluster_001/1_contigs/D_1.fasta.bad 
mv results/v3_trycycler/trycycler/cluster_001/1_contigs/G_Utg1404.fasta results/v3_trycycler/trycycler/cluster_001/1_contigs/G_Utg1404.fasta.bad 

trycycler reconcile \
    --reads results/v3_trycycler/filtlong/p95/Nanopore_J4_filtlong.fastq.gz \
     --cluster_dir results/v3_trycycler/trycycler/cluster_001
```

### Run MSA

```{bash}
srun --cpus-per-task 20 --mem=10G trycycler msa \
    --cluster_dir results/v3_trycycler/trycycler/cluster_001 \
    --threads 20
```

### Partition reads

```{bash}
srun --cpus-per-task 20 --mem=10G  trycycler partition \
    --reads results/v3_trycycler/filtlong/p95/Nanopore_J4_filtlong.fastq.gz \
    --cluster_dirs results/v3_trycycler/trycycler/cluster_001 \
    --threads 20
```

### Generate a consensus

```{bash}
trycycler consensus --cluster_dir results/v3_trycycler/trycycler/cluster_001

#check nr of contigs
grep -c ">" results/v3_trycycler/trycycler/cluster_001/7_final_consensus.fasta

conda deactivate
```

### Check raw assembly

#### Run prokka to identify CDS

```{bash}
conda activate prokka

srun --cpus-per-task 10 --mem=10G prokka \
    --outdir results/v3_trycycler/trycycler/cluster_001/prokka \
    --prefix trycycler_raw \
    results/v3_trycycler/trycycler/cluster_001/7_final_consensus.fasta \
    --cpus 10 --kingdom Bacteria --compliant

#calculate the average protein length
seqtk seq results/v3_trycycler/trycycler/cluster_001/prokka/trycycler_raw.faa | awk 'NR % 2 == 0' | awk '{sum += length($0); n++} END {print sum/n;}'

conda deactivate
```

#### Run pseudofinder to identify pseudogenes

```{bash}
mkdir results/v3_trycycler/trycycler/cluster_001/pseudofinder/

#run pseduofinder without refrence
conda activate pseudofinder

pseudo_db="/home/ndombro/personal/projects/alteromonas_j4/db/uniprot_sprot.fasta"

srun --cpus-per-task 20 --mem=10G pseudofinder.py annotate \
    --genome  results/v3_trycycler/trycycler/cluster_001/prokka/trycycler_raw.gbk \
    --outprefix results/v3_trycycler/trycycler/cluster_001/pseudofinder/pseudo_ \
    --database $pseudo_db --threads 20

conda deactivate
```

**Comments**

-   Prokka
    -   Found 4806 CDS
    -   Found 334 potential pseudo-genes
    -   Mean length: 272.944 \<-- this is better as compared to the other approaches, which are around 250-260
-   Pseudofinder:
    -   Pseudogenes (total): 984
    -   Intact genes: 3599

## Polishing the assembly with medaka

```{bash}
mkdir -p results/v3_trycycler/medaka/r1

conda activate medaka_1.11.3 

#run medaka : 46079
srun --cpus-per-task 20 --mem=50G medaka_consensus \
    -i results/v3_trycycler/filtlong/p95/Nanopore_J4_filtlong.fastq.gz \
    -d results/v3_trycycler/trycycler/cluster_001/7_final_consensus.fasta \
    -o results/v3_trycycler/medaka/r1 \
    -t 20 \
    -m r941_min_fast_g507

conda deactivate
```

### Check medaka assembly

#### Run prokka

```{bash}
conda activate prokka

srun --cpus-per-task 10 --mem=10G prokka \
    --outdir results/v3_trycycler/medaka/r1/prokka \
    --prefix trycycler_med \
    results/v3_trycycler/medaka/r1/consensus.fasta \
    --cpus 10 --kingdom Bacteria --compliant

#calculate average protein length
seqtk seq results/v3_trycycler/medaka/r1/prokka/trycycler_med.faa | awk 'NR % 2 == 0' | awk '{sum += length($0); n++} END {print sum/n;}'

conda deactivate
```

#### Run pseudofinder

```{bash}
mkdir results/v3_trycycler/medaka/r1/pseudofinder

#run pseduofinder without refrence
conda activate pseudofinder

pseudo_db="/home/ndombro/personal/projects/alteromonas_j4/db/uniprot_sprot.fasta"

srun --cpus-per-task 20 --mem=10G pseudofinder.py annotate \
    --genome results/v3_trycycler/medaka/r1/prokka/trycycler_med.gbk \
    --outprefix results/v3_trycycler/medaka/r1/pseudofinder/pseudo_ \
    --database $pseudo_db --threads 20

conda deactivate
```

**Comments**

-   Prokka
    -   Found 4937 CDS
    -   Found 385 potential pseudo-genes
    -   Mean length: 265.041 \<-- this is shorter than the initial assembly
-   Pseudofinder:
    -   Pseudogenes (total): 1093
    -   Intact genes: 3610

**In this case, medaka seem to make the assembly worse, so I wont proceed with this polished assembly but I will use the trycycler consensus directly!**

## Polishing with homopolish

```{bash}
mkdir results/v3_trycycler/homopolish 

conda activate homopolish

#no genomes found with default mash_threshold, tried lowering to:
#0.9,0.85, 0.8
#only 0.8 gave results
srun --cpus-per-task 20 --mem=20G python3 /zfs/omics/projects/bioinformatics/software/homopolish/homopolish.py polish \
  -a results/v3_trycycler/trycycler/cluster_001/7_final_consensus.fasta \
  -m R9.4.pkl \
  --mash_threshold 0.80 \
  -s /zfs/omics/projects/bioinformatics/databases/homopolish/290424/bacteria.msh \
  -o results/v3_trycycler/homopolish \
  -t 20

conda deactivate
```

**Comment:**

-   To find a set of genomes to compare the assembly to mash_threshold needed to be lower to 0.8 This is consistent with gtdbtk suggesting that J4 is a new genus with only very distant genomes in its database based on the ANI. Its unclear if this is due to J4 being indeed a new genus or the sequence quality increasing the distance to reference genomes. Regardless, it is something to keep in mind as the more similar the reference genomes are, the better the polishing will go

### Check homopolish assembly

#### Run prokka

```{bash}
conda activate prokka

srun --cpus-per-task 10 --mem=10G prokka \
    --outdir results/v3_trycycler/homopolish/prokka \
    --prefix trycycler_homopol \
     results/v3_trycycler/homopolish/7_final_consensus_homopolished.fasta \
    --cpus 10 --kingdom Bacteria --compliant

#estimate average protein length
seqtk seq results/v3_trycycler/homopolish/prokka/trycycler_homopol.faa | awk 'NR % 2 == 0' | awk '{sum += length($0); n++} END {print sum/n;}'

conda deactivate
```

#### Run pseudofinder

```{bash}
mkdir results/v3_trycycler/homopolish/pseudofinder

#run pseduofinder without refrence
conda activate pseudofinder

pseudo_db="/home/ndombro/personal/projects/alteromonas_j4/db/uniprot_sprot.fasta"

srun --cpus-per-task 20 --mem=10G pseudofinder.py annotate \
    --genome results/v3_trycycler/homopolish/prokka/trycycler_homopol.gbk \
    --outprefix results/v3_trycycler/homopolish/pseudofinder/pseudo_ \
    --database $pseudo_db --threads 20

conda deactivate
```

#### Run quast

```{bash}
#run quast
conda activate quast_5.2.0

srun --cpus-per-task 20 --mem=10G quast.py \
  results/v3_trycycler/homopolish/7_final_consensus_homopolished.fasta \
  -r  db/references/GCF_016756315.1_ASM1675631v1_genomic.fna \
  -g results/v3_trycycler/homopolish/prokka/trycycler_homopol.gff \
  --nanopore data/Nanopore_J4_SuperBasecalling.fastq.gz \
  -o results/v3_trycycler/homopolish//quast_report \
  -t 20

conda deactivate 
```

#### Run qualimap

##### Do read mapping

```{bash}
#read mapping
mkdir results/v3_trycycler/homopolish/minimap2/

mamba activate minimap2

srun --cpus-per-task 20 --mem=20G minimap2 -ax map-ont -t 20 \
    results/v3_trycycler/homopolish/7_final_consensus_homopolished.fasta \
    data/Nanopore_J4_SuperBasecalling.fastq.gz | \
    samtools sort -o results/v3_trycycler/homopolish/minimap2/homopolished_mapped.bam

#filter out mapped only reads
samtools view -F 260  results/v3_trycycler/homopolish/minimap2/homopolished_mapped.bam \
    -o  results/v3_trycycler/homopolish/minimap2/homopolished_mapped_only.bam

#cleanup
rm results/v3_trycycler/homopolish/minimap2/homopolished_mapped.bam

conda deactivate
```

##### Run qualimap

```{bash}
#check quality
mamba activate qualimap_2.3.-0

srun --cpus-per-task 20 --mem=50G  qualimap bamqc \
  -nt 20 \
  -nw 50000 \
  -bam  results/v3_trycycler/homopolish/minimap2/homopolished_mapped_only.bam \
  -outdir  results/v3_trycycler/homopolish/qualimap

conda deactivate
```

#### Check taxonomy with Gtdbtk

```{bash}
mkdir results/v3_trycycler/homopolish/gtdb 

conda activate gtdbtk_2.3.2

srun --cpus-per-task 20 --mem=100G gtdbtk classify_wf --genome_dir  results/v3_trycycler/homopolish/ \
    --extension fasta \
    --out_dir results/v3_trycycler/homopolish//gtdb \
    --mash_db /zfs/omics/projects/metatools/DB/GTDB_tk/GTDB_tk_r214/mash/ \
    --cpus 20

conda deactivate
```

#### Check completeness with Checkm2

```{bash}
mkdir results/v3_trycycler/homopolish/checkm2  

#run checkm2
conda activate checkm2_1.0.2

srun --cpus-per-task 20 --mem=10G checkm2 predict --threads 30 \
  --input  results/v3_trycycler/homopolish/  \
  -x fasta --force \
  --output-directory results/v3_trycycler/homopolish/checkm2 

conda deactivate
```


**Comments on the quality control steps**

-   Prokka
    -   4720 CDS
    -   278.014 average gene length
-   Pseudofinder
    -   Pseudogenes (total): 905
    -   Intact genes: 3629
-   GTDBtk:
    -   assigned to: d\_\_Bacteria;p\_\_Pseudomonadota;c\_\_Gammaproteobacteria;o\_\_Enterobacterales_A;f\_\_Alteromonadaceae;g\_\_Alteromonas
    -   closest genome: d\_\_Bacteria;p\_\_Pseudomonadota;c\_\_Gammaproteobacteria;o\_\_Enterobacterales_A;f\_\_Alteromonadaceae;g\_\_Alteromonas;s\_\_Alteromonas sp009811495 (85.27 ANI)
-   Checkm2
    -   98.1 compl \<-- that is \~10% better as compared to the NCBI genome
    -   5.11 cont
    -   0.877 coding density
    -   4,492,914 genome size
    -   0.44 GC


## Run Circlator to rotate

```{bash}
mkdir results/v3_trycycler/circlator/

conda activate circlator

circlator get_dnaa results/v3_trycycler/circlator/circlator

circlator fixstart results/v3_trycycler/homopolish/7_final_consensus_homopolished.fasta results/v3_trycycler/circlator/

conda deactivate
```

**Comment**

- Notice: circlator get_dnaa  does not work since the link to uniprot is broken
- Since we don't have a dnaA database, dnaA was not found, so  a gene predicted by prodigal was used - the match was on the forward strand. 
- this is not ideal, so I checked for more tools that we can use 



## Run dnaapler v0.7.0 to rotate

Link to [github](https://github.com/gbouras13/dnaapler)

```{bash}
mamba create -n dnaapler_env
conda activate dnaapler_env
mamba install -c bioconda dnaapler

mkdir results/v3_trycycler/dnaapler/

dnaapler chromosome -i results/v3_trycycler/homopolish/7_final_consensus_homopolished.fasta \
    -o results/v3_trycycler/dnaapler/ -p J4 -t 8
```

Notice:

- The top blastx hit for dnaA is sp|A8G7Q2|DNAA_SERP5, which has length of 463 AAs.  The alignment did not begin with a valid start codon.
- 422 AAs were covered in the best hit, with an overall coverage of 91.14%
- 258 AAs were identical, with an overall identity of 61.14%.
- The CDS most overlapping the tophit has a start coordinate of 955704


### CheckM on dnaapler output

```{bash}
mkdir results/v3_trycycler/dnaapler/checkm2  

#run checkm2
conda activate checkm2_1.0.2_old

srun --cpus-per-task 20 --mem=10G checkm2 predict --threads 30 \
  --input  results/v3_trycycler/dnaapler/  \
  -x fasta --force \
  --output-directory results/v3_trycycler/dnaapler/checkm2

conda deactivate
```

#### Run prokka

```{bash}
conda activate prokka

srun --cpus-per-task 10 --mem=10G prokka \
    --outdir results/v3_trycycler/dnaapler/prokka \
    --prefix dnaapler \
     results/v3_trycycler/dnaapler/J4_reoriented.fasta \
    --cpus 10 --kingdom Bacteria --compliant

grep -c ">" results/v3_trycycler/dnaapler/prokka/dnaapler.faa

#estimate average protein length: 278.112
seqtk seq results/v3_trycycler/dnaapler/prokka/dnaapler.faa | awk 'NR % 2 == 0' | awk '{sum += length($0); n++} END {print sum/n;}'

conda deactivate
```

#### Run pseudofinder

```{bash}
mkdir results/v3_trycycler/dnaapler/pseudofinder

#run pseduofinder without refrence
conda activate pseudofinder

pseudo_db="/home/ndombro/personal/projects/alteromonas_j4/db/uniprot_sprot.fasta"

srun --cpus-per-task 20 --mem=10G pseudofinder.py annotate \
    --genome results/v3_trycycler/dnaapler/prokka/dnaapler.gbk \
    --outprefix results/v3_trycycler/dnaapler/pseudofinder/pseudo_ \
    --database $pseudo_db --threads 20

conda deactivate
```

**Comments**

-   Prokka
    -   4719 CDS
    -   278.112 average gene length
-   Pseudofinder
    -   Pseudogenes (total): 903
    -   Intact genes: 3633
-   Checkm2
    -   97.98 compl
    -   5.07 cont
    -   0.877 coding density
    -   4,492,914 genome size
    -   0.44 GC


## Summary

**Genome stats when looking at the genome downloaded from [NCBI, GCF_036870955:](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_036870955.1/)**

-   Prokka
    -   5277 CDS
    -   247.029 average CDS length
-   Pseudofinder
    -   Pseudogenes (total): 1346
    -   Intact genes: 3589

**Genome stats after trycycler + homopolish (approach 3):**

-   Prokka
    -   4720 CDS
    -   278.014 average CDS length \<-- this is a bit longer compared to only the trycycler assembly, thus I will keep the homopolish results
-   Pseudofinder
    -   Pseudogenes (total): 905
    -   Intact genes: 3629


## References

::: {#refs}
:::