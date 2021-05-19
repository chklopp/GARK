# GARK
Check genome assembly repeat content using kmers

It is now common to compare genome assembly and read kmer content to verify assembly completness with kat comp (fig 1. of https://www.biorxiv.org/content/10.1101/064733v2.full.pdf). This gives a quick overview of missing kmers in the homozygote or heterozygote sections of the genome. 

Hereunder I present an extension of this analysis for repeated parts of the genome. The idea emerged while comparing ONT and PacBio HiFi assemblies. HiFi assemblies where larger than ONT assemblies with close to the same busco scores. HiFi assembly lengths were also closer to the expected genome size infered from kmer content using genomescope2 (http://qb.cshl.edu/genomescope/genomescope2.0/). 

kat comp compares kmer counts in the genome versus in the reads. Zero counts correspond to kmers missing in the genome but present in the reads. In a haploid genome assembly, counts of 1 correspond to homozygous or heterozogous sections. In a diploid genome assembly, kmer counts of 1 correspond to heterozygous sections and counts of 2 to homozygous sections. Counts above 2 correspond to section which are several times in the genome and can arise due to polyploidy or non resolved ancient duplications or recent duplications or local duplications. 

For repeats, the idea is first to threshold kmers from the genome. We could for example decide that we call a repeated kmer a kmer present at least 50 times in the genome. Then in the reads this kmer should be 50 * coverage (homozygous coverage peak of genomescope2) time in the reads. If this is the case the ration (read kmer count)/(genome kmer count) should be equal to the coverage for all repeated kmers. Kmers having a ratio far from the expected coverage show repeat collapse or expension events in the genome compared to the reads. 

The procedure chains seven steps :
- setting the genome kmer count threshold
- calculating the read kmer count threshold
- creating a kmer database from the genome with kmc using the thresold
- creating a kmer database from the reads with kmc using the thresold
- comparing both databases and creating an inserction database
- extraction the counts from the genome and the reads for the intersecting kmer and calculating the ratio
- ploting the ratio histogram 

We are going to compare two genome assemblies performed with ONT and HiFi reads with the HiFi reads.
- ont_genome.fasta is the ONT genome fasta file
- hifi_genome.fasta is the HiFi genome fasta file
- @hifi.list contains the urls of the HiFi read fasta files

Here are the corresponding command lines

```
mkdir tmp
# -ci2000 corresponds to 20 time de repeat threshold valuer, here you should chose a value close to coverage * repead threshold
kmc -v -k21 -m60 -t4 -ci2000 -cs1000001 @hifi.list  hifi_read.k21 tmp
# -ci2000 is the repeat threshold 
kmc -v -k21 -m60 -t4 -ci100 -cs1000001 -fm ont_genome.fasta ont_genome.k21 tmp
kmc -v -k21 -m60 -t4 -ci100 -cs1000001 -fm hifi_genome.fasta hifi_genome.k21 tmp
kmc_tools simple ont_genome.k21 hifi_read.k21 intersect common_hifi_read_ont_genome.k21 -ocleft
kmc_tools simple hifi_genome.k21 hifi_read.k21 intersect common_hifi_read_hifi_genome.k21 -ocleft
kmc_tools transform common_hifi_read_ont_genome.k21 dump common_hifi_read_ont_genome.k21.txt
kmc_tools transform common_hifi_read_hifi_genome.k21 dump common_hifi_read_hifi_genome.k21.txt
kmc_tools transform hifi_read.k21 dump hifi_read.k21.txt

sort -k1,1 common_hifi_read_ont_genome.k21.txt > common_hifi_read_ont_genome.k21.txt.sort
sort -k1,1 common_hifi_read_hifi_genome.k21.txt > common_hifi_read_hifi_genome.k21.txt.sort
sort -k1,1 hifi_read.k21.txt > hifi_read.k21.txt.sort

join -1 1 -2 1 common_hifi_read_ont_genome.k21.txt.sort hifi_read.k21.txt.sort | awk '{print $1"\t"$2"\t"$3"\t"($3/$2)"\t"(int(log($3/$2)/log(10)));}' > common_hifi_read_ont_genome.k21.txt.sort.ratio
join -1 1 -2 1 common_hifi_read_hifi_genome.k21.txt.sort hifi_read.k21.txt.sort | awk '{print $1"\t"$2"\t"$3"\t"($3/$2)"\t"(int(log($3/$2)/log(10)));}' > common_hifi_read_hifi_genome.k21.txt.sort.ratio
```

Once you have the ratio files you can load them into R and produce the plots.
```
library(ggplot2)
hifi <- read.table("common_hifi_read_hifi_genome.k21.txt.sort.ratio")
ont <- read.table("common_hifi_read_ont_genome.k21.txt.sort.ratio")
names(hifi) <- c("kmer","genome","reads","ratio","class")
names(ont) <- c("kmer","genome","reads","ratio","class")
hifi$origin <- "hifi"
ont$origin <- "ont"
all <- rbind(hifi,ont)
ggplot(all[all$origin == "ont" & all$reads != 1000001,] , aes(ratio, fill = origin)) + geom_histogram(binwidth =1) 
ggplot(all[all$origin == "hifi" & all$reads != 1000001,] , aes(ratio, fill = origin)) + geom_histogram(binwidth =1)
```

If the genome contains all the repeated kmers of the reads the histogram should look like : [!(https://user-images.githubusercontent.com/4004727/118785312-0b0b7100-b891-11eb-95ce-d16c30d61c64.png)]


