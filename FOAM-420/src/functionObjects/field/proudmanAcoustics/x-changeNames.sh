#!/bin/bash

oldName=wallShearStress
newName=proudmanAcoustics

mv $oldName.H $newName.H
mv $oldName.C $newName.C

sed -i s/$oldName/$newName/g $newName.H
sed -i s/$oldName/$newName/g $newName.C
