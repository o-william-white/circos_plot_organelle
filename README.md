
## Circos plot organelle 

Python wrapper script to create an annotated circos plot. 

### Install dependencies using conda
```
# create environment
conda create -n circos_plot_organelle -c bioconda samtools trf circos python

# activate environment
conda activate circos_plot_organelle
```

### Clone repo

```
git clone https://github.com/o-william-white/circos_plot_organelle
cd circos_plot_organelle
```

### Run the python script using the example data
```
python circos_plot_organelle.py \
   --fasta example_data/010210013_StH_H3.fasta \
   --bam example_data/010210013_StH_H3.bam \
   --bed example_data/010210013_StH_H3.bed \
   --window 50
```

### Output

See the image circos.png, which should look like the image below. The plot shows the following attributes from outside to inside:
- Sequence position
- Annotation names
- Annotations on the + strand
- Annotations on the - strand
- Coverage
- GC content
- Repeat content

![alt text](example_data/circos.png?)

