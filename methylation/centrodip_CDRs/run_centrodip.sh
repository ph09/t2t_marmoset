#!/bin/bash

#SBATCH --job-name=marmoset_centrodip
#SBATCH --partition=long
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=7-00:00:00
#SBATCH --output=logs/marmoset_centrodip.%j.log

set -euo pipefail
# set -x

THREADS=16

# smoothing parameters
WINDOW_SIZE=10000

# dip detection parameters
PROMINENCE=0.25
HEIGHT=0.1
BROADNESS=0.9

# filtering parameters
MIN_SIZE=1000
MIN_SCORE=500
CLUSTER_SIZE=-1

#########################################
# === run centrodip on Baguette_pri === #
#########################################
centrodip \
    --plot \
    --threads "$THREADS" \
    --mod-code "m" \
    --window-size "$WINDOW_SIZE" \
    --prominence "$PROMINENCE" \
    --height "$HEIGHT" \
    --broadness "$BROADNESS" \
    --min-size "$MIN_SIZE" \
    --min-score "$MIN_SCORE" \
    --cluster-distance "$CLUSTER_SIZE" \
    /private/groups/migalab/jmmenend/marmoset/data/modkit/Baguette_pri.CpG.pileup.bed \
    /private/groups/migalab/jmmenend/marmoset/data/cenSat/GCA_049354715.1.cenSat.v1.0.active_alpha.bed \
    GCA_049354715.1_calJac240_Baguette_pri.centrodip.bed

#########################################
# === run centrodip on Baguette_alt === #
#########################################
centrodip \
    --plot \
    --threads "$THREADS" \
    --mod-code "m" \
    --window-size "$WINDOW_SIZE" \
    --prominence "$PROMINENCE" \
    --height "$HEIGHT" \
    --broadness "$BROADNESS" \
    --min-size "$MIN_SIZE" \
    --min-score "$MIN_SCORE" \
    --cluster-distance "$CLUSTER_SIZE" \
    /private/groups/migalab/jmmenend/marmoset/data/modkit/Baguette_alt.CpG.pileup.bed \
    /private/groups/migalab/jmmenend/marmoset/data/cenSat/GCA_049354655.1.cenSat.v1.0.active_alpha.bed \
    GCA_049354655.1_calJac240_Baguette_alt.centrodip.bed

#######################################
# === run centrodip on Maisie_pri === #
#######################################
centrodip \
    --plot \
    --threads "$THREADS" \
    --mod-code "m" \
    --window-size "$WINDOW_SIZE" \
    --prominence "$PROMINENCE" \
    --height "$HEIGHT" \
    --broadness "$BROADNESS" \
    --min-size "$MIN_SIZE" \
    --min-score "$MIN_SCORE" \
    --cluster-distance "$CLUSTER_SIZE" \
    /private/groups/migalab/jmmenend/marmoset/data/modkit/Maisie_pri.CpG.pileup.bed \
    /private/groups/migalab/jmmenend/marmoset/data/cenSat/GCA_049354665.1.cenSat.v1.0.active_alpha.bed \
    GCA_049354665.1_calJac220_Maisie_pri.centrodip.bed

#######################################
# === run centrodip on Maisie_alt === #
#######################################
centrodip \
    --plot \
    --threads "$THREADS" \
    --mod-code "m" \
    --window-size "$WINDOW_SIZE" \
    --prominence "$PROMINENCE" \
    --height "$HEIGHT" \
    --broadness "$BROADNESS" \
    --min-size "$MIN_SIZE" \
    --min-score "$MIN_SCORE" \
    --cluster-distance "$CLUSTER_SIZE" \
    /private/groups/migalab/jmmenend/marmoset/data/modkit/Maisie_alt.CpG.pileup.bed \
    /private/groups/migalab/jmmenend/marmoset/data/cenSat/GCA_049354675.1.cenSat.v1.0.active_alpha.bed \
    GCA_049354675.1_calJac220_Maisie_alt.centrodip.bed
