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

\*---------------------------------------------------------------------------*/

#include "relativeHumidity.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fieldInitializations
{
    defineTypeNameAndDebug(relativeHumidity, 0);
    addToRunTimeSelectionTable(fieldInit, relativeHumidity, initMethod);
}
}

// * * * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fieldInitializations::relativeHumidity::relativeHumidity
(
    const fvMesh& mesh,
    const fvSolutionRegistry& localDb,
    const dictionary& fieldDict,
    const word& fN
)
:
    fieldInit(mesh, localDb, fieldDict, fN)
{
    //requires temperature field
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::fieldInitializations::relativeHumidity::correct()
{
    //if already initialized stop
    if (initialised())
    {
        return;
    }

    fieldInit::initMsg();

    word TName = initDict().lookupOrDefault<word>("T", "T");
    volScalarField& T
        = const_cast<volScalarField&>
        (localDb().lookupObject<volScalarField>(TName));

    scalar pAbs = initDict().lookupOrDefault<scalar>("Pabs", 1e5);
    scalar Mvap = initDict().lookupOrDefault<scalar>("Mvap", 18.02);
    scalar Mmix = initDict().lookupOrDefault<scalar>("Mmix", 28.96);
    scalar relativeHumidity
        = readScalar(initDict().lookup("relativeHumidity"));

    //hardcoded saturation pressure for water vapour
    scalarField Ppartial
    (
        (2337 * Foam::exp(6879*(-1/T.primitiveField() + 1/293.15)
        - 5.031*Foam::log(T.primitiveField()/293.15)))
        * relativeHumidity
    );

    volScalarField& f
        = const_cast<volScalarField&>
        (localDb().lookupObject<volScalarField>(name()));

    f.primitiveFieldRef()
        = Ppartial * Mvap / (Ppartial * Mvap + (pAbs - Ppartial)*Mmix);

    // the field has been initialized
    initialised() = true;

    initCoupledBoundaries();
}



// ************************************************************************* //
