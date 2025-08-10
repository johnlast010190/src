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

#include "multiphaseFluidState.H"
#include "stateComponents/singlePhaseFluid/singlePhaseFluidState.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
namespace stateFunctions
{

defineTypeNameAndDebug(multiphaseFluidState, 0);

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void multiphaseFluidState::addAlphaTypeGrad(const word& phaseName)
{
    if (input().found("system"))
    {

        if (input().subDict("system").found("fvSchemes"))
        {
            if
            (
                input().subDict("system")
                .subDict("fvSchemes").found("gradSchemes")
            )
            {
                if
                (
                    input().subDict("system")
                    .subDict("fvSchemes").subDict("gradSchemes")
                    .found("grad(alpha)")
                )
                {
                    system().add(word("fvSchemes"), dictionary(), false);
                    system().subDict("fvSchemes")
                        .add(word("gradSchemes"), dictionary(), false);

                    system().subDict("fvSchemes").subDict("gradSchemes").add
                    (
                        word("grad(alpha."+phaseName +")"),
                        input().subDict("system").subDict("fvSchemes")
                        .subDict("gradSchemes").lookup("grad(alpha)"),
                        false
                    );
                }
            }
        }
    }
}


Xfer<HashTable<dictionary, stateFunction::interfacePair, stateFunction::interfacePair::hash>>
stateFunctions::multiphaseFluidState::assembleBinaryDictionaries
(
    word property,
    const dictionary& defaultValue,
    const List<dictionary>& matprops
)
{
    label nPhases = phases_.size();
    label nBinaryPairs = (nPhases*nPhases - nPhases)/2;

    HashTable<dictionary, stateFunction::interfacePair, stateFunction::interfacePair::hash>
        binaryData(nBinaryPairs);


    for (label n = 0; n < nPhases-1; n++)
    {
        for (label m = n+1; m < nPhases; m++)
        {
            //check for reciprocal data definitions, primary phase wins
            //if phases define incompatible reciprocals
            dictionary cDat(defaultValue);

            //do not insert entry at all if not available
            bool found = false;

            if (matprops[n].found("binaryPhaseData"))
            {
                if (matprops[n].subDict("binaryPhaseData").found(phases_[m]))
                {
                    const dictionary& bpdDict
                    (
                        matprops[n].subDict("binaryPhaseData")
                        .subDict(phases_[m])
                    );

                    if (bpdDict.found(property))
                    {
                        cDat = bpdDict.subDict(property);

                        found = true;
                    }
                }
            }
            if (matprops[m].found("binaryPhaseData") && !found)
            {
                if (matprops[m].subDict("binaryPhaseData").found(phases_[n]))
                {

                    const dictionary& bpdDict
                    (
                        matprops[m].subDict("binaryPhaseData")
                        .subDict(phases_[n])
                    );

                    if (bpdDict.found(property))
                    {
                        cDat = bpdDict.subDict(property);
                    }
                }
            }

            binaryData.insert(interfacePair(phases_[n], phases_[m]), cDat);
        }
    }

    return binaryData.xfer();
}



Xfer<List<dictionary>> multiphaseFluidState::assembleMaterialProperties()
{

    List<dictionary> matprops(phases_.size(), dictionary());

    forAll(phases_, pi)
    {
        word phase(phases_[pi]);

        //construct merged input dictionary
        dictionary& materialProperties(matprops[pi]);

        materialProperties.add("materialName", phase);

        materialProperties.merge
        (
            defaults().subDict("materialProperties").subOrEmptyDict(phase)
        );

        if (input().found("materialProperties"))
        {
            materialProperties.merge
            (
                input().subDict("materialProperties").subOrEmptyDict(phase)
            );
        }
    }


    return matprops.xfer();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

multiphaseFluidState::multiphaseFluidState
(
    word region,
    const dictionary& input,
    const dictionary& defaults,
    const stateFunction& master,
    const stateIndex& index,
    word meshName
)
:
    turbulenceModelState
    (
        region,
        input,
        defaults,
        master,
        index,
        meshName
    ),
    phases_(0)
{
    //check that more than one material is specified

    phases_ = wordList(stateFunction::input().lookup("materials"));

    if (phases_.size() < 2)
    {
        FatalErrorInFunction
            << "Single material specified for a multi-phase state: "
            << phases_ << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void multiphaseFluidState::correct()
{
    // clean-up, remove generic alpha entry, since it has been specialised
    dictionary& fields(stateDict_.subDict("fields"));

    if (fields.found("alpha"))
    {
        fields.remove("alpha");
    }

    turbulenceModelState::correct();

}


void multiphaseFluidState::finalise()
{
    //reverse order for finalise compared to other dictionary manipulators
    turbulenceModelState::finalise();

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
