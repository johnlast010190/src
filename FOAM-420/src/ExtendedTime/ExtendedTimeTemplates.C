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
    (c) 2011-2012 OpenFOAM Foundation
    (c) 2017 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "ExtendedTime.H"

namespace Foam
{

// * * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
void ExtendedTime::addCurrentFieldsToCP(label cpi)
{
    typedef GeometricField<Type, PatchField, GeoMesh> GeometricFieldType;

    const fvMesh& mesh = lookupObject<fvMesh>(Foam::fvMesh::defaultRegion);

    forAll(cpTable_, tbi)
    {
        bool fieldIsActive =
            mesh.foundObject<GeometricFieldType>(cpTable_[tbi]);

        if (fieldIsActive)
        {
            bool fieldIsInCP =
                checkpoint_[cpi].foundObject<GeometricFieldType>(cpTable_[tbi]);

            if (fieldIsInCP)
            {
                // Copy field to cp
                GeometricFieldType& field =
                    const_cast<GeometricFieldType&>
                    (
                        checkpoint_[cpi].lookupObject<GeometricFieldType>
                        (
                            cpTable_[tbi]
                        )
                    );
                field.forceAssign(mesh.lookupObject<GeometricFieldType>(cpTable_[tbi]));
            }
            else
            {
                // Create field in cp
                autoPtr<GeometricFieldType> field
                (
                    new GeometricFieldType
                    (
                        IOobject
                        (
                            cpTable_[tbi],
                            checkpoint_[cpi].instance(),
                            checkpoint_[cpi],
                            IOobject::NO_READ,
                            IOobject::NO_WRITE
                        ),
                        mesh.lookupObject<GeometricFieldType>(cpTable_[tbi])
                    )
                );
                field.ptr()->store();
            }
        }
    }
}

template<class Type, template<class> class PatchField, class GeoMesh>
void ExtendedTime::restoreFieldsFromCP(label cpi)
{
    typedef GeometricField<Type, PatchField, GeoMesh> GeometricFieldType;

    const fvMesh& mesh = lookupObject<fvMesh>(Foam::fvMesh::defaultRegion);

    const wordList& nameList =
        checkpoint_[cpi].names<GeometricFieldType>();

    forAll(nameList, ni)
    {
        if (debug)
        {
            Info<< "Restoring field " << nameList[ni]
                << " from checkpoint " << cpi << endl;
        }

        GeometricFieldType& field =
            const_cast<GeometricFieldType&>
            (
                mesh.lookupObject<GeometricFieldType>
                (
                    nameList[ni]
                )
            );
        const GeometricFieldType& cpfield =
            checkpoint_[cpi].lookupObject<GeometricFieldType>
            (
                nameList[ni]
            );

        // Rename field so as to avoid the triggering of storeOldTimes function
        field.rename(field.name() + "renamed_0");
        field.forceAssign(cpfield);
        field.rename(nameList[ni]);
    }
}

} // End namespace FOAM

// ************************************************************************* //
