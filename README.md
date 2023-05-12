# comSVis
Visualizing complex SVs 

# Installation
This will clone the repository. You can run the comSVis within this directory.
```
git clone https://github.com/vanessa-porter/comSVis.git
```

# Input Files
Edit these paths in events.yaml
- Sniffles VCF File
- TXT file of read names frmo a SV of interest
- Whole genome BAM file (used to create the CNV files)
- Ploidy file from PloiDetect (models.txt)
- CN file from PloiDetect (cna_condensed.txt)

### **Run snakemake**
This is the command to run it with conda. The `-c` parameter can be used to specify maximum number of threads. 

```
snakemake -c 30 --use-conda 
```
