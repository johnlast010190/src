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

#include "ensight/part/ensightPart.H"
#include "ensight/part/ensightParts.H"
#include "ensight/type/ensightPTraits.H"
#include "primitives/direction/direction.H"
#include "db/typeInfo/typeInfo.H"

// * * * * * * * * * * Static Private Member Functions * * * * * * * * * * * //

template<template<typename> class FieldContainer, class Type>
void Foam::ensightSerialOutput::writeFieldContent
(
    const word& key,
    const FieldContainer<Type>& fld,
    ensightFile& os
)
{
    if (fld.size())
    {
        os.writeKeyword(key);

        for (direction d=0; d < pTraits<Type>::nComponents; ++d)
        {
            const label cmpt = ensightPTraits<Type>::componentOrder[d];

            os.writeList(fld.component(cmpt));
        }
    }
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class Type>
bool Foam::ensightSerialOutput::writeField
(
    const Field<Type>& fld,
    const ensightPartFaces& part,
    ensightFile& os,
    const bool nodeValues
)
{
    if (part.size() && fld.size())
    {
        if (nodeValues)
        {
            os.beginPart(part.index());

            os.writeKeyword("coordinates");
            for (direction d=0; d < pTraits<Type>::nComponents; ++d)
            {
                const label cmpt = ensightPTraits<Type>::componentOrder[d];

                os.writeList(fld.component(cmpt));
            }
        }
        else
        {
            os.beginPart(part.index());

            for (label typei=0; typei < ensightFaces::nTypes; ++typei)
            {
                const ensightFaces::elemType what =
                    ensightFaces::elemType(typei);

                writeFieldContent
                (
                    ensightFaces::key(what),
                    Field<Type>(fld, part.faceIds(what)),
                    os
                );
            }
        }
    }

    return true;
}


template<class Type>
bool Foam::ensightSerialOutput::writeField
(
    const GeometricField<Type, fvPatchField, volMesh>& fld,
    const ensightPartCells& part,
    ensightFile& os
)
{
    if (part.size() && fld.size())
    {
        os.beginPart(part.index());

        for (label typei=0; typei < ensightCells::nTypes; ++typei)
        {
            const ensightCells::elemType what = ensightCells::elemType(typei);

            writeFieldContent
            (
                ensightCells::key(what),
                Field<Type>(fld, part.cellIds(what)),
                os
            );
        }
    }

    return true;
}


template<class Type>
bool Foam::ensightSerialOutput::writeField
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const ensightParts& list,
    ensightFile& os
)
{
    forAllConstIter(ensightParts::StorageType, list, iter)
    {
        if (isA<ensightPartFaces>(*iter))
        {
            const ensightPartFaces& part =
                dynamicCast<const ensightPartFaces&>(*iter);

            const label patchi = part.patchIndex();
            if (patchi >= 0 && patchi < vf.boundaryField().size())
            {
                writeField
                (
                    vf.boundaryField()[patchi],
                    part,
                    os,
                    false
                );
            }
        }
        else
        {
            const ensightPartCells& part =
                dynamicCast<const ensightPartCells&>(*iter);

            writeField(vf, part, os);
        }
    }

    return true;
}


// ************************************************************************* //
