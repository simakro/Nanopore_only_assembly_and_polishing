# Nanopore_only_assembly_and_polishing
Genome assembly and polishing from (long read) Nanopore-only data; originally for HCMV (should be good for bacteria as well)

# Installation

## GTDBTK
For classification with GTDBTK first a database has to be downloaded.
For this:
    1. create a new conda environment on the same machine you intend to run the workflow:
        conda create -n gtdbtk
    2. change into this env
        conda activate gtdbtk
    3. install the tool:
        conda install -c bioconda gtdbtk=2.3.2

When installation completes the following message can be seen:
"""
    GTDB-Tk v2.3.2 requires ~78G of external data which needs to be downloaded
    and extracted. This can be done automatically, or manually.

    Automatic:

        1. Run the command "download-db.sh" to automatically download and extract to:
            /home/simon/mambaforge/envs/gtdbtk/share/gtdbtk-2.3.2/db/

    Manual:

        1. Manually download the latest reference data:
            wget https://data.gtdb.ecogenomic.org/releases/release214/214.0/auxillary_files/gtdbtk_r214_data.tar.gz

        2. Extract the archive to a target directory:
            tar -xvzf gtdbtk_r214_data.tar.gz -C "/path/to/target/db" --strip 1 > /dev/null
            rm gtdbtk_r214_data.tar.gz

        3. Set the GTDBTK_DATA_PATH environment variable by running:
            conda env config vars set GTDBTK_DATA_PATH="/path/to/target/db
"""

Use the automatic versionwith the shell script.
Within the snakefile in rule classify_gtdbtk, in the shell directive, the first
line has to be edited to fit to the path given in (1.) e.g.:

"mamba env config vars set GTDBTK_DATA_PATH=/homes/simon/.conda/envs/gtdbtk/share/gtdbtk-2.3.2/db && "
 would be wrong here and had to be changed to 
"mamba env config vars set GTDBTK_DATA_PATH=/home/simon/mambaforge/envs/gtdbtk/share/gtdbtk-2.3.2/db && "

## PlasClass & PlasFlow
Both packages can not be installed via conda without trouble.
Therefore I decided to place them as packages into the pkgs folder.
THe PlasClass package is small enough so I could include it in my github repository.
However, PlasFlow is so big (65-85Mb) that I would not want to include it.
It can be downloaded by cloning the PlasFlow repository into the pkgs folder.
To avoid stacking repo within repo, the .git folder within PlasFlow package should be removed.
!!!It is extremely important to be careful to only delete the .git in PlasFLow and not the main workflow!!!
Thus change dir into PlasFLow and run rm .git only there.


## Create DAG graph
snakemake --dag | dot > DAG.dot
dot -Tpng DAG.dot -o DAG2.png