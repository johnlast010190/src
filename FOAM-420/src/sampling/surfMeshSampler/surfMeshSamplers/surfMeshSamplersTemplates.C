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
    (c) 2016 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "surfMeshSampler/surfMeshSamplers/surfMeshSamplers.H"
#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::wordList
Foam::surfMeshSamplers::acceptType() const
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;

    return mesh_.names<VolFieldType>(fieldSelection_);
}

#if 0
template<class Type>
Foam::wordList
Foam::surfMeshSamplers::acceptType
(
    const IOobjectList& objects,
    bool fromFiles
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;

    if (fromFiles_)
    {
        // This should actually be in the caller:
        // IOobjectList objects1 = objects.lookup(fieldSelection_);

        return objects.names(VolFieldType::typeName, fieldSelection_);
    }
    else
    {
        return mesh_.names<VolFieldType>(fieldSelection_);
    }
}
#endif


// ************************************************************************* //
