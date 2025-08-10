#!/bin/bash

oldName=uniformFixedValue
newName=fixedValueZeroGradient

mv "$oldName"FaPatchFields.H    "$newName"FaPatchFields.H
mv "$oldName"FaPatchFields.C    "$newName"FaPatchFields.C
mv "$oldName"FaPatchFieldsFwd.H "$newName"FaPatchFieldsFwd.H
mv "$oldName"FaPatchField.H     "$newName"FaPatchField.H
mv "$oldName"FaPatchField.C     "$newName"FaPatchField.C

sed -i s/$oldName/$newName/g    "$newName"FaPatchFields.H
sed -i s/$oldName/$newName/g    "$newName"FaPatchFields.C
sed -i s/$oldName/$newName/g    "$newName"FaPatchFieldsFwd.H
sed -i s/$oldName/$newName/g    "$newName"FaPatchField.H
sed -i s/$oldName/$newName/g    "$newName"FaPatchField.C
