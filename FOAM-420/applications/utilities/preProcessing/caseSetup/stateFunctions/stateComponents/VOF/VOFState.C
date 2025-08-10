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

#include "VOFState.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "multiphaseThermo/multiphaseThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
namespace stateFunctions
{

defineTypeNameAndDebug(VOFState, 0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

VOFState::VOFState
(
    word region,
    const dictionary& input,
    const dictionary& defaults,
    const stateFunction& master,
    const stateIndex& index,
    word meshName
)
:
    multiphaseMultiComponentState
    (
        region,
        input,
        defaults,
        master,
        index,
        meshName
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


void VOFState::initialise()
{
    word specifiedPassive =
        input().subOrEmptyDict("materialProperties").lookupOrDefault
        (
            "passivePhase", word::null
        );

    forAll(phaseNames(), pi)
    {
        // Don't create a fieldMap for the passive phase
        if
        (
            phaseNames()[pi] == specifiedPassive
         || (
                specifiedPassive == word::null
             && phaseNames().size() == 2 && pi == 1
            )
        )
        {
            continue;
        }
        // fieldMaps
        stateDict_.subDict("fieldMaps").add
        (
            word("alpha." + phaseNames()[pi]),
            input().subDict("fieldMaps").lookup("alpha")
        );
    }
    // Remove old one from input to prevent it being merged later
    input().subDict("fieldMaps").remove("alpha");

    multiphaseMultiComponentState::initialise();

    // Check that we have a multiphase material model
    dictionary& matDict = system().subDict("materialProperties");
    const word matType(matDict.lookup("materialType"));
    if (matType != "multiphase")
    {
        FatalErrorInFunction
            << "The 'VOF' state requires a multiphase material"
            << nl << exit(FatalError);
    }

    if (specifiedPassive != word::null)
    {
        matDict.add("passivePhase", specifiedPassive);
    }
    else if (phaseNames().size() == 2)
    {
        matDict.add("passivePhase", phaseNames()[1]);
    }

    forAll(phaseNames(), pi)
    {
        // find grad(alpha) in defaults and write to state using proper phase name
        addAlphaTypeGrad(phaseNames()[pi]);
    }

    multiphaseThermo::sigmaTable sigmas;
    if (matDict.found("sigmas"))
    {
        multiphaseThermo::sigmaTable newSigmas
        (
            matDict.lookup("sigmas"), multiphaseThermo::iNewSigma()
        );
        sigmas.transfer(newSigmas);
    }

    multiphaseThermo::sigmaTable defSigmas
    (
        defaults().subDict("materialProperties").subDict
        (
            "binaryPhaseData"
        ).lookup("sigmas"),
        multiphaseThermo::iNewSigma()
    );

    forAll(phaseNames(), pi)
    {
        for (label pj = pi+1; pj < phaseNames().size(); pj++)
        {
            Foam::interfacePair pij(phaseNames()[pi], phaseNames()[pj]);
            multiphaseThermo::sigmaTable::iterator iter = defSigmas.find(pij);
            Function1<scalar>* ptr;
            if (iter.found())
            {
                ptr = defSigmas.remove(iter);
            }
            else
            {
                ptr =
                    new Function1Types::Constant<scalar>
                    (
                        "sigma",
                        0.0
                    );
            }
            if (!sigmas.insert(pij, ptr))
            {
                delete ptr;
            }
        }
    }
    matDict.add("sigmas", sigmas, true);
}


void VOFState::finalise()
{
    multiphaseMultiComponentState::finalise();

    // clean-up
    if (system().found("fvSchemes"))
    {
        if (system().subDict("fvSchemes").found("gradSchemes"))
        {
            system().subDict("fvSchemes")
                .subDict("gradSchemes").remove("grad(alpha)");
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace stateFunction
} // End namespace Foam

// ************************************************************************* //
