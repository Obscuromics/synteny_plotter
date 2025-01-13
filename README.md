# synteny_plotter

A tool to visualise synteny relationships between genomes.

Installation

```
conda create synteny_plotter
conda activate synteny_plotter
conda install conda-forge::r-dplyr
conda install -c conda-forge r-argparse
```

## Testing runs

Required parameters are:
- `--busco_list`: list of the BUSCO `full_table.tsv` files (top genome is reference for now), see example `test_data/busco_list.tsv`
- `--chrom_list`: list of teh chromosome tsv files, see example `test_data/chrom_list.tsv` (columns: chrom, length, order, invert, annot)
- `-o`: prefix of the output PDF file containing the plot

Default run is with colouring by chromosomes based on the reference genome (top file in the busco list):
```
Rscript scripts/generate_synteny_plot.R --busco_list test_data/busco_list.tsv --chrom_list test_data/chrom_list.tsv -o example_output/example_three_fly_species_plot
```

You can also provide a file with ancestral linkage groups (ALGs), see example `test_data/ALGs_n5_72_all_8`. The plot will be coloured based on which ALG a busco gene belongs to (use ALG colours ecoded in the script):
```
Rscript scripts/generate_synteny_plot.R --busco_list test_data/busco_list.tsv --chrom_list test_data/chrom_list.tsv -o example_output/example_three_fly_species_plot --colour_by algs --alg_file test_data/ALGs_n5_72_all_8
```

Colour by ALGs and provide ALG colours manually as a table, see example `test_data/ALGs_colours.tsv`:
```
Rscript scripts/generate_synteny_plot.R --busco_list test_data/busco_list.tsv --chrom_list test_data/chrom_list.tsv -o example_output/example_three_fly_species_plot --colour_by algs --alg_file test_data/ALGs_n5_72_all_8 --alg_colours test_data/ALGs_colours.tsv
```

See the example of the plot here: `example_output/example_three_fly_species_plot.pdf`
