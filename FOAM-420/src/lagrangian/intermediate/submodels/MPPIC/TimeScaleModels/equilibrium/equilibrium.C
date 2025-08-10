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
    (c) 2013-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "submodels/MPPIC/TimeScaleModels/equilibrium/equilibrium.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace TimeScaleModels
{
    defineTypeNameAndDebug(equilibrium, 0);

    addToRunTimeSelectionTable
    (
        TimeScaleModel,
        equilibrium,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::TimeScaleModels::equilibrium::equilibrium
(
    const dictionary& dict
)
:
    TimeScaleModel(dict)
{}


Foam::TimeScaleModels::equilibrium::equilibrium
(
    const equilibrium& hc
)
:
    TimeScaleModel(hc)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::TimeScaleModels::equilibrium::~equilibrium()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::FieldField<Foam::Field, Foam::scalar>>
Foam::TimeScaleModels::equilibrium::oneByTau
(
    const FieldField<Field, scalar>& alpha,
    const FieldField<Field, scalar>& r32,
    const FieldField<Field, scalar>& uSqr,
    const FieldField<Field, scalar>& f
) const
{
    static const scalar a =
        16.0/sqrt(3.0*constant::mathematical::pi)
       *0.25*(1.0 - e_*e_);

    return
        a
       *alpha*sqrt(max(uSqr, scalar(0)))/max(r32, SMALL)
       *alphaPacked_/max(alphaPacked_ - alpha, SMALL);
}


// ************************************************************************* //
