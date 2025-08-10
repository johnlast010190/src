#!/bin/bash

oldName=patchObserver
newName=vibroacousticObserver

mv "$oldName".H     "$newName".H
mv "$oldName".C     "$newName".C

sed -i s/$oldName/$newName/g    "$newName".H
sed -i s/$oldName/$newName/g    "$newName".C
