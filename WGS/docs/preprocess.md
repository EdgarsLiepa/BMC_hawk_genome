
# Preprocessing

### IN: 

forward and reverese reads from sequencer

- <SAMPLE_ID>\_1.fq.gz 
- <SAMPLE_ID>\_2.fq.gz

### OUT:

**DIR=** /home_beegfs/edgars01/Ineta/WGS/starpfaili/trimmed_sequences

Files for creted sample:

**CutAdapt:**
- *<SAMPLE_ID>*\_1_cutadapt.fq.gz
- *<SAMPLE_ID>*\_2_cutadapt.fq.gz

**Trimmomatic:**
- *<SAMPLE_ID>*\_F_paired.fq.gz
- *<SAMPLE_ID>*\_R_paired.fq.gz
- *<SAMPLE_ID>*\_F_unpaired.fq.gz
- *<SAMPLE_ID>*\_R_unpaired.fq.gz


*Prints processed sample names*

``` bash

ls /home_beegfs/edgars01/Ineta/WGS/starpfaili/trimmed_sequences/*.fq.gz | sed -e 's/_1_cutadapt.fq.gz//' -e 's/_2_cutadapt.fq.gz//' | sed -e 's/_F_paired.fq.gz//' -e 's/_R_paired.fq.gz//' | sed -e 's/_F_unpaired.fq.gz//' -e 's/_R_unpaired.fq.gz//' |sed -e 's/\/home_beegfs\/edgars01\/Ineta\/WGS\/starpfaili\/trimmed_sequences\///'  | sort -u

```

 