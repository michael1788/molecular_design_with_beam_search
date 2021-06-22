#!/bin/bash

conda env create -f environment.yml
git submodule init
git submodule update
