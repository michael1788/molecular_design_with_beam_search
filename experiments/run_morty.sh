#!/bin/bash

CONFIGFILE=$1

if [ $# -eq 0 ] ; then
    echo "Configfile path not supplied."
else

sh run_processing.sh $CONFIGFILE &&
sh run_training.sh $CONFIGFILE &&
sh run_beam_generation.sh $CONFIGFILE &&
sh run_final_ranking.sh $CONFIGFILE &&
sh run_boxplot.sh $CONFIGFILE &&
sh run_tsne_plot.sh $CONFIGFILE

fi
