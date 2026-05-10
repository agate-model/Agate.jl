#!/bin/bash

#make a clean temp folder
mkdir /tmp/Agate.jl
#copy folder but exclude unneeded folders
rsync -av --progress /home/joost/Agate.jl/ /tmp/Agate.jl/ --exclude .git \
--exclude data --exclude pre-processing --exclude analysis
#copy temp folder to bp
scp -r /tmp/Agate.jl/ ba18321@bp1-login.acrc.bris.ac.uk:/user/work/ba18321/
#if copy ok, remove temp folder
rm -rf /tmp/Agate.jl
