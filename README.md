# AIRR 2026 Dowser Tutorial

## Overview

This is a tutorial for a basic walk through of making and analyzing B cell clonal families using [10x Genomics](https://www.10xgenomics.com/products/single-cell-immune-profiling) BCR (B cell receptor) sequencing data. The tutorial in it's entirey can be found in `airr_2026_tutorial.rmd`. Installing Dowser and its various dependencies locally is encouraged, but a [Docker container for this tutorial is also available](https://hub.docker.com/repository/docker/colejensen/tyche-dowser/general).

Knowledge of using R is assumed. Additionally, various preprocessing steps have been assumed to be done. Further instruction for how to do said sets can be explained in other software tutorials such AIRR-flow or in be found in the [Immcantation tutorials](https://immcantation.readthedocs.io/en/stable/getting_started/10x_tutorial.html#assign-v-d-and-j-genes-using-igblast).

## Getting started

We will be using a simulated dataset generated using [simble](https://github.com/hoehnlab/simble/tree/main#). You can find our example dataset in `data/all_samples_airr.tsv`. Additional files will be needed in the tutorial. The UCA inference steps require various model files and they can be found in [OLGA](https://github.com/statbiophys/OLGA/tree/master/olga/default_models). With TyCHE, template files are needed to run BEAST2 through Dowser. They can be found in `data/xml-templates` for this tutorial, but the most update to date templates can be found [here](https://github.com/hoehnlab/xml-templates).

#### Install and load libraries and programs
In this tutorial we will use various programs/languages. Below are instructions on how to install the various needed tools.

### Option 1: Docker
Pull and run the container:
```bash
docker pull colejensen/tyche-dowser:latest
docker run -d --rm -p 8788:8787 -e USER=rstudio -e PASSWORD=rstudio colejensen/tyche-dowser:latest
```
Then go to http://localhost:8788 in your browser and follow along! The username and password are `rstudio`. The tutorial can be found in `/data`.

#### Troubleshooting
Depending on your setup, you may run into some problems launching the container. Below are troubleshooting tips based on your system.

##### Linux
```bash
# Check what's on 8788
sudo ss -tlnp | grep 8788

# If it's a running container
sudo docker ps
sudo docker stop <container_id>

# If it's native RStudio Server
sudo systemctl stop rstudio-server

# If you get a permission error, add your user to the docker group
# (then log out and back in)
sudo usermod -aG docker $USER

# Then start your container
sudo docker run -d --rm -p 8788:8787 -e USER=rstudio -e PASSWORD=rstudio colejensen/tyche-dowser:latest
```

##### Mac
```bash
# Check what's on 8788
lsof -i :8788

# If it's a running container
docker ps
docker stop <container_id>

# Then start your container
docker run -d --rm -p 8788:8787 -e USER=rstudio -e PASSWORD=rstudio colejensen/tyche-dowser:latest
```

##### Windows
```powershell
# Check what's on 8788
netstat -ano | findstr :8788

# If it's a running container
docker ps
docker stop <container_id>

# Then start your container
docker run -d --rm -p 8788:8787 -e USER=rstudio -e PASSWORD=rstudio colejensen/tyche-dowser:latest
```
### Option 2: Local installation

#### Dowser
Dowser is an R package with various dependencies. Below are instructions on how get Dowser to install locally.

Dowser has various Bioconductor dependencies that don't always install without specifying them.
```r
# This will check if the package is installed and install it if not
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

biopackages <- c('ggtree', 'treeio', 'pwalign', 'Biostrings', 'GenomicAlignments',
                  'IRanges', 'S4Arrays', 'Rsamtools', 'SparseArray', 'DelayedArray',
                  'SummarizedExperiment')
new_biopackages <- biopackages[!(biopackages %in% installed.packages())]
if (length(new_biopackages)) {BiocManager::install(new_biopackages)}
```

Now the latest version of Fowser can be installed from Github using the devtools package.
```r
devtools::install_github("immcantation/dowser")
```

#### IgPhyML
We will be making trees using IgPhyML, a B cell specific lineage tree model.

##### Compiling from source on Linux
```bash
apt-get install automake autoconf libblas-dev liblapack-dev libatlas-base-dev
git clone https://github.com/immcantation/igphyml
cd igphyml
./make_phyml_omp
```
##### Compiling from source on Mac
This is a bit more tricky. [Please see instructions here](https://igphyml.readthedocs.io/en/latest/install.html)

##### Compiling from source on Windows
This is not available.

#### IMGT References
IMGT references will be needed. Run the following:
```bash
curl -O https://raw.githubusercontent.com/immcantation/immcantation/master/scripts/fetch_imgtdb.sh
chmod +x fetch_imgtdb.sh
./fetch_imgtdb.sh
```

#### Python
First, make sure that Python is installed and accessible from the command line:
```bash
python3 --version
```

On some systems, the command may be `python` instead of `python3`:
```bash
python --version
```

Next, install the required Python packages.
```bash
pip install biopython logomaker matplotlib numpy olga pandas scipy
```

Depending on your setup, you may need to use:
```bash
python3 -m pip install biopython logomaker matplotlib numpy olga pandas scipy
```

or
```bash
python -m pip install biopython logomaker matplotlib numpy olga pandas scipy
```

It is important to verify that the Python executable used by R is the same Python installation where these packages were installed. This is a common source of errors, especially when using virtual environments or multiple Python installations.

To check which Python executable R is calling, in R run:
```r
system2("python3", "--version")
system2("python3", "-c", "import sys; print(sys.executable)")
```

If your system uses `python` instead of `python3`, replace `python3` with `python` in the commands above.

You can also test whether the required packages are available from within R:
```r
system2("python3", args = c(
  "-c", shQuote('import Bio, logomaker, matplotlib, numpy, olga, pandas, scipy; print("All required Python packages are installed and accessible.")')))
```

If this command runs without errors, then R is using a Python installation with all required packages available.

If you do find you are missing a package, you can install it through R using this command:
```r
system2("python3", args = c("-m", "pip", "install", "missing package name"))
```

#### BEAST2 and TyCHE

##### On Mac and Windows machines, we recommend:
1. Download the appropriate version for your device from [the v2.7.7 release](https://github.com/CompEvol/beast2/releases/tag/v2.7.7)
2. Open BEAUti, click on the “File” menu, and select “Manage Packages…”.
3. In the package manager, find and install the “BEAST Classic” package.
4. Follow this [tutorial to add the “extra packages” package repository](www.beast2.org/managing-packages). Use https://github.com/CompEvol/CBAN/blob/master/packages-extra-2.7.xml as the package repository URL.
5. In the package manager, find and install the “TyCHE” package. This tutorial relies on version v0.0.10 or later.
6. In the package manager, find and install the “rootfreqs” package.

If using BEAUti's graphical interface for the package manager does not work, these 
alternative steps work on Mac:
```bash
# alterative for step 3, install BEAST Classic
/Applications/BEAST\ 2.7.7/bin/packagemanager -add BEAST_CLASSIC

# alternative for step 4, add "extra packages" package repo
# first confirm that this is the correct path for your beauti.properties file
ls ~/Library/Application\ Support/BEAST/2.7/beauti.properties
# add the extra packages repo 
echo "packages.url=https\://raw.githubusercontent.com/CompEvol/CBAN/master/packages-extra-2.7.xml" >> ~/Library/Application\ Support/BEAST/2.7/beauti.properties

# alterative for steps 5 and 6
# install TyCHE package
/Applications/BEAST\ 2.7.7/bin/packagemanager -add TyCHE

# install rootfreqs package
/Applications/BEAST\ 2.7.7/bin/packagemanager -add rootfreqs
```

##### For Linux machines, we recommend running:
```bash
# Choose appropriate version for your architecture (x86 or aarch64)
BEAST=BEAST.v2.7.7.Linux.x86.tgz # or BEAST=BEAST.v2.7.7.Linux.aarch64.tgz

# download file and uncompress
curl -O https://github.com/CompEvol/beast2/releases/download/v2.7.7/$BEAST
tar -xvzf $BEAST

# optionally remove the compressed file
rm $BEAST

# run BEAST, at least with help, to allow it to set up its directories
~/beast/bin/beast -help

# install BEAST Classic package
~/beast/bin/packagemanager -add BEAST_CLASSIC

# add "extra packages" package repo
echo "packages.url=https\://raw.githubusercontent.com/CompEvol/CBAN/master/packages-extra-2.7.xml" >> ~/.beast/2.7/beauti.properties

# install TyCHE package
~/beast/bin/packagemanager -add TyCHE

# install rootfreqs package
~/beast/bin/packagemanager -add rootfreqs
```
