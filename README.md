# PlotDiversiVision

This purpose of this repo is to play with ideas to automate measures of plant diversity in images.  The images may be found here: https://huggingface.co/datasets/imageomics/NEON-plant-subplot-pilot.  These are images taken from two sites and several plots within the sites.  Species identified in these images are provided in this repo. 

This repo includes code to 1) call a foundational model BioCLIP 2 to determine probabilities of species appearing in eachimage, 2) process the  BioCLIP 2 results to create a structured dataset of specific species predictions and 3) assess quality of predictions.

## Run BioCLIP 2

This repository includes scripts for preparing plant species labels and running
BioCLIP 2 grid predictions. The workflow is:

1. Create TaxonoPy-passed species lists for NEON plots to get ground-truth
   labels and for species lists from different sources.
2. Map downstream labels by lookup from that resolved species list.
3. Run BioCLIP 2 over image grid crops using the resolved species list.

### Directory Layout

- `assets/NEON_plotData.csv`: source NEON plot data used to extract plot-level
  species labels.
- `assets/neonSiteSpeciesList.csv`: source plant list from NEON
  used for site-level species lists.
- `assets/<name>_BONAPlist.csv`: source BONAP plant list used for region-level
  species lists.
- `assets/conus_plant_lists_accepted.csv`: source CONUS accepted plant list
  used for state-level species lists.
- `assets/species_list/`: TaxonoPy-passed species lists for BioCLIP 2 label sets.
- `assets/test_labels/`: subplot-level label files used for image
  benchmarking.

### TaxonoPy/GNVerifier Setup

TaxonoPy uses GNVerifier for name resolution. Install `taxonopy` from
`requirements.txt`, then make sure a `gnverifier` executable is available.

You can find different versions of gnverifier from the [GitHub release](https://github.com/gnames/gnverifier/releases).
On macOS ARM64, download and extract the GNVerifier release:

```bash
mkdir -p outputs/tools/gnverifier/

curl -L \
  -o outputs/tools/gnverifier/gnverifier-v1.3.7-mac-arm64.tar.gz \
  https://github.com/gnames/gnverifier/releases/download/v1.3.7/gnverifier-v1.3.7-mac-arm64.tar.gz

tar -xzf outputs/tools/gnverifier/gnverifier-v1.3.7-mac-arm64.tar.gz \
  -C outputs/tools/gnverifier \
  --strip-components 1

chmod +x outputs/tools/gnverifier/gnverifier
```

The species-list scripts prepend that directory to `PATH` automatically and
redirect `HOME` into the run-specific `outputs/species_list/<name>/taxonopy/`
directory. That keeps GNVerifier config files out of the user home directory
and avoids permission issues in sandboxed runs.

### Create Resolved Species Lists

Use `scripts/create_taxonopy_neon_species_list.py` when starting from
`assets/NEON_plotData.csv`. It extracts unique species labels for a plot,
runs TaxonoPy, and writes only the final resolved list to
`assets/species_list/`.

Example for one SCBI plot:

```bash
python scripts/create_taxonopy_neon_plot_species_list.py \
  --source-csv assets/NEON_plotData.csv \
  --plot-id SCBI_008 \
  --name SCBI_008
```

To create a merged species list for multiple plots:

```bash
python scripts/create_taxonopy_neon_plot_species_list.py \
  --source-csv assets/NEON_plotData.csv \
  --plot-id SCBI_005 --plot-id SCBI_008 --plot-id SCBI_015 --plot-id SCBI_021 \
  --name SCBI_plot
```

This writes intermediate files to `outputs/species_list/SCBI_plot/`
and the final label file to `assets/species_list/SCBI_plot_labels.csv`.

For a single state from the CONUS accepted plant list, use
`scripts/create_taxonopy_conus_species_list.py`:

```bash
python scripts/create_taxonopy_conus_species_list.py \
  --source-csv assets/conus_plant_lists_accepted.csv \
  --state Colorado \
  --name CPER_state
```

This writes the intermediate state species CSV to
`outputs/species_list/CPER_state/` and the final resolved list to
`assets/species_list/CPER_state_labels.csv`.

For other plant lists without the need to filter states, use
`scripts/create_taxonopy_other_species_list.py`:

```bash
python scripts/create_taxonopy_other_species_list.py \
  --source-csv assets/CPER_BONAPlist.csv \
  --column scientificName \
  --name CPER_BONAPlist
```

This writes the intermediate state species CSV to
`outputs/species_list/CPER_BONAPlist/` and the final resolved list to
`assets/species_list/CPER_BONAPlist_labels.csv`.

### Create Test Label Files

Once a TaxonoPy-passed species list exists, downstream label files should be
mapped by lookup instead of resolving the same labels again. Use
`scripts/map_labels_from_resolved_species_list.py` to create one subplot-level
test label file for each benchmark plot:

```bash
python scripts/map_labels_from_resolved_species_list.py \
  --plot-id SCBI_005
```

By default this reads `assets/NEON_plotData.csv`, looks up labels in
`assets/species_list/<site_name>_plot_labels.csv`, and writes
`assets/test_labels/<plot_id>_subplot_labels.csv`. The output has one row per
subplot, with original NEON labels, resolved BioCLIP labels, resolved
scientific names, TaxonoPy taxonomy strings, TaxonoPy resolution statuses, and
any labels that could not be mapped.

### Grid-Based BioCLIP 2 Predictions

After preparing a species list, use `scripts/predict_grid_species.py` to split
an image into a `3x3` or `4x4` grid and run BioCLIP 2 on each crop. The script
uses the Python API and saves the full per-grid species probability list;
integration across crops is left for a later step.

Example:

```bash
python scripts/predict_grid_species.py \
  --data-root data \
  --plot-id SCBI_008 \
  --grid-size 3 \
  --output-csv outputs/grid_predictions/SCBI_008_predictions.csv
```

When `--species-list` is omitted, the script infers the site-level label set
from the matched images, for example `assets/species_list/SCBI_NEON_labels.csv` for
SCBI plots and `assets/species_list/CPER_labels.csv` for CPER plots. Pass
`--species-list` explicitly only when you want to override that behavior with a
plot-specific or custom label set. The script can also be pointed at explicit
image files, directories, or globs with `--images`.

By default, inference labels come from the `resolved_taxonomic_labels` column.
Pass `--species-column` only when testing a different label form.

The first run for a species list computes BioCLIP text embeddings and caches
them under `outputs/text_embeddings/`. Later runs with the same model, species
list, species column, and label order reuse that cache. Use
`--no-text-embedding-cache` to force recomputation.

The prediction CSV records the image path, parsed plot metadata, subplot ID,
image date, grid position, crop bounds, and one probability column per species.
Probability column names come from the species list's `resolved_labels` column;
if that column is unavailable, the selected inference label column is used. The
year remains available as parsed metadata, but it is not a required filtering
layer.

### Project Predictions Onto Original Images

Use `scripts/project_grid_predictions.py` to create a PNG for each original
image represented by a prediction CSV. The figure keeps the source image's
native pixel dimensions, draws the grid boundaries, and annotates each grid
cell with every species and its probability when the probability is greater
than `0.1`.

```bash
python scripts/project_grid_predictions.py \
  --predictions temp_results \
  --output-dir temp_results/projected_figures \
  --threshold 0.1
```

The output directory mirrors the prediction CSV directory structure, preventing
files generated with different species lists or grid sizes from overwriting one
another. Source images are resolved first from `image_path`, then from
`relative_image_path` under `--data-root` (which defaults to `data`).

## Process BioCLIP 2 Results
Outcomes from Bioclip 2 include one probability per grid crop per candidate species (4x4 = 16 corps). The number of probabilities depends on the chose grid size to analyze each image and the number of species depends on the list of possible labels to be considered.  Scripts in R were written to transform matrices of probabilities to predictions of specific species names in each image.    Also, ground truth data are imported to a) calculate true positive, true negative, false positive, false negative, and false discover rates (aka success rates).  The R scripts are stored in the folder playWithBioClipResults.

### Transform probabilities to species names and calculate success/failure rates

The  script, accessBioClipOutcomes.R,  transforms probabilities to a predicted species list and calculates success rates for images at one site,  subplot, with a grid specification, and labels file.  The  script, accessBioClipOutcomesALLsite.R,  does the same for all sites and all conditions on which BioCLIP 2 was run.  In turn, the following files are output from the code: tprFpr.csv and allPredAndTrue_except3.csv.  The file allPredAndTrue_except3.csv lists both predicted and ground truth species across all images considered, except for 3 bioclip runs.  There is a small bug in the code that causes it to fail.  

To run the scripts, simply update the workspace and data pathnames for the required files. Requited files are found in assets/ and temp_results/, but listed explicitly in the file, filesForBioclipPredictions.csv.  For example, the first row of  filesForBioclipPredictions.csv provides the prediction file name, labels file name and file for ground truth to assess bioclip results from CPER_001_bonap_grid2.csv. 

## Assess BioCLIP 2 Results: 

### Graph BioCLIP species predictions. 
The script in R, initialLook_dtb.R, creates graph to explore how well the current pipeline does relative to ground truth. 

### Phylo Script Description

Phylo script is to be run after 'PlaywithBioClip' script where data variables are put into format which feeds species labels into a df ready to be run through 'UPhyloMaker2' which is a rooted tree wrapper that gives a rough approximate of species phylogeny based on global vascular species models.

### Funding

This work was supported by both the [Imageomics Institute](https://imageomics.org) as part of [FloraPalooza](https://github.com/Imageomics/FloraPalooza-2026). The Imageomics Institute is funded by the US National Science Foundation's Harnessing the Data Revolution (HDR) program under [Award #2118240](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2118240) (Imageomics: A New Frontier of Biological Information Powered by Knowledge-Guided Machine Learning). Any opinions, findings and conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the National Science Foundation.
