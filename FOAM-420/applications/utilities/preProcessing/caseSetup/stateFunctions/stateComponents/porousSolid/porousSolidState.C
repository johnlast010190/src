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
    (c) 2016-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "porousSolidState.H"
#include "stateIndex/stateIndex.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
namespace stateFunctions
{

defineTypeNameAndDebug(porousSolidState, 0);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

porousSolidState::porousSolidState
(
    word region,
    const dictionary& input,
    const dictionary& defaults,
    const stateFunction& master,
    const stateIndex& index,
    word meshName
)
:
    solidState
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

void porousSolidState::finalise()
{
    solidState::finalise();

    label otherRegioni = -1;
    const PtrList<stateFunction>& regions(index_.regions());
    forAll(regions, regioni)
    {
        if
        (
            &regions[regioni] != this
         && regions[regioni].meshName() == meshName()
        )
        {
            if (otherRegioni != -1)
            {
                FatalErrorInFunction
                    << "Only two regions sharing the same mesh are supported "
                    << "for 'porous' state."
                    << nl << exit(FatalError);
            }
            else
            {
                otherRegioni = regioni;
            }
        }
    }
    if (otherRegioni == -1)
    {
        FatalErrorInFunction
            << "Two regions sharing the same mesh are required for 'porous' "
            << "state."
            << nl << exit(FatalError);
    }

    const dictionary& matDict = system().subDict("materialProperties");
    const word matType(matDict.lookup("materialType"));
    if (matType != "solid")
    {
        FatalErrorInFunction
            << "The 'porous' state should only be applied to a solid "
            << "material" << nl << exit(FatalError);
    }
    const dictionary& otherMatDict =
        index_.regions()[otherRegioni].system().subDict("materialProperties");
    const word otherMatType(otherMatDict.lookup("materialType"));
    if (otherMatType != "fluid" && otherMatType != "reactingFluid")
    {
        FatalErrorInFunction
            << "The 'porous' state should be applied to a solid material "
            << "which shares the same mesh as a fluid material"
            << nl << exit(FatalError);
    }

    // fvOptions
    if
    (
        system().found("fvOptions")
     && system().subDict("fvOptions").found("solidHeatConductionSolver")
    )
    {
        system().subDict("fvOptions").subDict("solidHeatConductionSolver").add
        (
            word("fluidRegion"),
            index_.regions()[otherRegioni].regionName()
        );
    }
}

} // End namespace stateFunction
} // End namespace Foam

// ************************************************************************* //
