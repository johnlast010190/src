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

#include "multiphaseEulerFoamRASState.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
namespace stateFunctions
{

defineTypeNameAndDebug(multiphaseEulerFoamRASState, 0);
addToRunTimeSelectionTable(stateFunction, multiphaseEulerFoamRASState, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

multiphaseEulerFoamRASState::multiphaseEulerFoamRASState
(
    word region,
    const dictionary& input,
    const dictionary& defaults,
    const stateFunction& master,
    const stateIndex& index,
    word meshName
)
:
    multiphaseEulerState
    (
        region,
        input,
        defaults,
        master,
        index,
        meshName
    )
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void multiphaseEulerFoamRASState::finalise()
{
    //reverse order for finalise compared to other dictionary manipulators
    multiphaseEulerState::finalise();

    // clean-up
    // remove SIMPLE subdict
    if (system().found("fvSolution"))
    {
        if (system().subDict("fvSolution").found("SIMPLE"))
        {
            system().subDict("fvSolution").remove("SIMPLE");
        }
    }
}

} // End namespace stateFunction
} // End namespace Foam

// ************************************************************************* //
