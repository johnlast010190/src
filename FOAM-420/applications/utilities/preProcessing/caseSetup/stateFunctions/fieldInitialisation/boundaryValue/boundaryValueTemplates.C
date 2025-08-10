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
    (c) 2016 Esi Ltd.

Class
    Foam::fieldInits::boundaryValue

Description
    boundaryValue initialization method.

SourceFiles
    boundaryValue.C

\*--------------------------------------------------------------------*/

#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * //

template<class Type>
bool Foam::fieldInitializations::boundaryValue::setBoundaryValue
(
    label patchID
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> GeoField;

    if (localDb().foundObject<GeoField>(name()))
    {
        GeoField& f
            = const_cast<GeoField&>
            (localDb().lookupObject<GeoField>(name()));

        Switch initialiseBoundaries
        = initDict().lookupOrDefault<Switch>("initialiseBoundaries", false);

        if (initialiseBoundaries)
        {
            Info<< "   " << "Initialising " << name()
                 << " boundary conditions" <<endl;
            f.correctBoundaryConditions();
        }


        f.primitiveFieldRef() = gAverage(f.boundaryField()[patchID]);

        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



// ************************************************************************* //
