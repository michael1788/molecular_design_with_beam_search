#!/bin/bash

CONFIGFILE=$1

if [ $# -eq 0 ] ; then
    echo "Configfile path not supplied."
else

python do_tsne_plot.py --configfile $CONFIGFILE --verbose True

fi