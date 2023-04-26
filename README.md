
## Circos plot organelle 

Python wrapper script to create an annotated circos plot. 

Install dependencies using conda
```
# create environment
conda create -n circos_plot_organelle -c bioconda samtools trf circos python

# activate environment
 conda activate circos_plot_organelle
```

Run the python script using the example data
```
python circos_plot_organelle.py \
   --fasta example_data/010210013_StH_H3.fasta \
   --bam example_data/010210013_StH_H3.bam \
   --bed example_data/010210013_StH_H3.bed \
   --window 50
```


