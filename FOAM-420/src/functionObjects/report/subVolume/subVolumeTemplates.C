/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  FOAM (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : 4.2.0
|    o     o     |  ESI Ltd. <http://esi.com/>
|       o        |
\*---------------------------------------------------------------------------
License
    This file is part of FOAMcore.
    FOAMcore is based on OpenFOAM (R) <http://www.openfoam.org/>.

    FOAMcore is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    FOAMcore is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with FOAMcore.  If not, see <http://www.gnu.org/licenses/>.

Copyright
    (c) 2010-2012 Esi Ltd.
    (c) 1991-2007 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "subVolume/subVolume.H"
#include "fields/volFields/volFields.H"
#include "db/IOstreams/IOstreams/IOmanip.H"
#include "fields/fvPatchFields/basic/zeroGradient/zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class T>
void Foam::functionObjects::subVolume::findFields(wordList& typeFieldNames)
{
    typeFieldNames.setSize(fieldNames_.size());
    label typeFieldI = 0;

    forAll(fieldNames_, fieldI)
    {
        const word& fldName = fieldNames_[fieldI];

        if (foundObject<T>(fldName))
        {
            typeFieldNames[typeFieldI++] = fldName;
        }
    }

    typeFieldNames.setSize(typeFieldI);
}


template <class T>
void Foam::functionObjects::subVolume::writeSubset(const word& fieldName)
{
    const GeometricField<T, fvPatchField, volMesh>& fld =
        lookupObject<GeometricField<T, fvPatchField, volMesh>>
        (
            fieldName
        );

    GeometricField<T, fvPatchField, volMesh> subFields( subsetMesh_().interpolate(fld,true) );

    subFields.rename(fieldName);

    subFields.write();
}


template <class T>
void Foam::functionObjects::subVolume::writeSubset(const wordList& typeFields)
{
    forAll(typeFields, i)
    {
        writeSubset<T>(typeFields[i]);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
