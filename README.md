



## fig1 : circos plot

simulation performed on the full human genome. Download it with:
```bash
wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz
gunzip hg38.fa.gz
```

the configfile is provided in the relative figure folder

```bash
# generate the synthetic VCF
pixi run inSVert simulate fig1-circos/config.yaml data/Homo_sapiens.GRCh38.dna.primary_assembly.fa --seed 123 -o fig1-circos/simulated.vcf

# plot the synthetic variants
pixi run Rscript scripts/plot_circlize.R fig1-circos/simulated.vcf data/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai fig1-circos/fig1 --format png
```

## fig2 : distributions

```bash
# Generate the normal and Pareto distributions.
pixi run inSVert simulate fig2-distributions/config_normal.yaml data/Homo_sapiens.GRCh38.dna.primary_assembly.fa --seed 123 -o fig2-distributions/simulated_normal.vcf

pixi run inSVert simulate fig2-distributions/config_pareto.yaml data/Homo_sapiens.GRCh38.dna.primary_assembly.fa --seed 123 -o fig2-distributions/simulated_pareto.vcf

# Plot the distributions.
pixi run Rscript scripts/plot_distributions.R fig2-distributions/simulated_normal.vcf fig2-distributions/simulated_pareto.vcf -o fig2-distributions
```

## fig3 : dotplot

```bash
# Generate synthetic VCF
pixi run inSVert simulate fig3-dotplot/config.yaml data/cerevisiae_test.fa --seed 123 -o fig3-dotplot/simulated.vcf

# Insert variant in the genome
pixi run inSVert insert data/cerevisiae_test.fa fig3-dotplot/simulated.vcf --ploidy 1 -o fig3-dotplot/simulated.fa --truth-vcf fig3-dotplot/inserted.vcf

# Align the edited genome and the unchanged reference to the reference.
# -c records gaps within alignments, which the plot needs to show insertions and deletions.
pixi run minimap2 -cx asm5 -t 4 data/cerevisiae_test.fa fig3-dotplot/simulated.fa > fig3-dotplot/simulated-vs-reference.paf

pixi run minimap2 -cx asm5 -t 4 data/cerevisiae_test.fa data/cerevisiae_test.fa > fig3-dotplot/reference-vs-reference.paf

# Dotplot of the edited genome against reference
pixi run Rscript scripts/dotplot.R fig3-dotplot/simulated-vs-reference.paf -o fig3-dotplot/simulated-dotplot --format png --display-contigs

# Dotplot of reference against reference (control)
pixi run Rscript scripts/dotplot.R fig3-dotplot/reference-vs-reference.paf -o fig3-dotplot/control-dotplot --format png --display-contigs

# Control and edited genome side by side
pixi run Rscript scripts/compare_dotplots.R fig3-dotplot/reference-vs-reference.paf fig3-dotplot/simulated-vs-reference.paf -o fig3-dotplot/control-vs-simulated-dotplot.png
```
