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
    (c) 2013 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "submodels/MPPIC/CorrectionLimitingMethods/absolute/absolute.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

namespace Foam
{
namespace CorrectionLimitingMethods
{
    defineTypeNameAndDebug(absolute, 0);

    addToRunTimeSelectionTable
    (
        CorrectionLimitingMethod,
        absolute,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CorrectionLimitingMethods::absolute::absolute(const dictionary& dict)
:
    CorrectionLimitingMethod(dict),
    e_(readScalar(dict.lookup("e")))
{}


Foam::CorrectionLimitingMethods::absolute::absolute(const absolute& cl)
:
    CorrectionLimitingMethod(cl),
    e_(cl.e_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::CorrectionLimitingMethods::absolute::~absolute()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::vector Foam::CorrectionLimitingMethods::absolute::limitedVelocity
(
    const vector uP,
    const vector dU,
    const vector uMean
) const
{
    const vector uRelative = uP - uMean;

    return minMod
    (
        dU,
      - (1.0 + this->e_)*uRelative
       *mag(uP)/max(mag(uRelative), SMALL)
    );
}


// ************************************************************************* //
