#!/bin/bash
#

label=`ls -d *_??`
for dirn in $label ; do
    cd $dirn/
    mkdir tmp
    cd tmp/
    if [ -f ../$dirn.xml ]; then mv ../$dirn.xml . ; fi
    if [ -f ../INFO.OUT ]; then mv ../INFO.OUT . ; fi
    if [ -f ../info.xml ]; then mv ../info.xml . ; fi
    if [ -f ../input.xml ]; then mv ../input.xml . ; fi
    if [ -f ../RMSDVEFF.OUT ]; then mv ../RMSDVEFF.OUT . ; fi
    if [ -f ../GEOMETRY.OUT ]; then mv ../GEOMETRY.OUT . ; fi
    if [ -f ../TOTENERGY.OUT ]; then mv ../TOTENERGY.OUT . ; fi
    if [ -f ../DTOTENERGY.OUT ]; then mv ../DTOTENERGY.OUT . ; fi
    cd ../
    ls -l | grep ^- | awk '{print $9}' | xargs rm
    mv tmp/* .
    rm -r tmp
    cd ../
done
