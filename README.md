



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

## fig 2 : dot plot


```bash
# generate the synthetic VCF
pixi run inSVert simulate fig2-dotplot/config.yaml data/cerevisiae_test.fa --seed 123 -o fig2-dotplot/simulated.vcf
# insert variants into a synthetic genome
pixi run inSVert insert data/cerevisiae_test.fa fig2-dotplot/simulated.vcf --ploidy 1 -o fig2-dotplot/simulated.fa
```