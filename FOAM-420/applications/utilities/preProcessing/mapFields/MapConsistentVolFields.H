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
    (c) 2011-2017 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#ifndef MapConsistentVolFields_H
#define MapConsistentVolFields_H

#include "fields/GeometricFields/GeometricField/GeometricField.H"
#include "meshToMesh0/meshToMesh0.H"
#include "db/IOobjectList/IOobjectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type, class CombineOp>
void MapConsistentVolFields
(
    const IOobjectList& objects,
    const meshToMesh0& meshToMesh0Interp,
    const meshToMesh0::order& mapOrder,
    const CombineOp& cop
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    const fvMesh& meshSource = meshToMesh0Interp.fromMesh();
    const fvMesh& meshTarget = meshToMesh0Interp.toMesh();

    IOobjectList fields = objects.lookupClass(fieldType::typeName);

    forAllIter(IOobjectList, fields, fieldIter)
    {
        Info<< "    interpolating " << fieldIter()->name()
            << endl;

        // Read field
        fieldType fieldSource(*fieldIter(), meshSource);

        IOobject fieldTargetIOobject
        (
            fieldIter()->name(),
            meshTarget.time().timeName(),
            meshTarget,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        );

        if (fieldTargetIOobject.typeHeaderOk<fieldType>(true))
        {
            // Read fieldTarget
            fieldType fieldTarget(fieldTargetIOobject, meshTarget);

            // Interpolate field
            meshToMesh0Interp.interpolate
            (
                fieldTarget,
                fieldSource,
                mapOrder,
                cop
            );

            // Write field
            fieldTarget.write();
        }
        else
        {
            fieldTargetIOobject.readOpt() = IOobject::NO_READ;

            // Interpolate field
            fieldType fieldTarget
            (
                fieldTargetIOobject,
                meshToMesh0Interp.interpolate
                (
                    fieldSource,
                    mapOrder,
                    cop
                )
            );

            // Write field
            fieldTarget.write();
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
