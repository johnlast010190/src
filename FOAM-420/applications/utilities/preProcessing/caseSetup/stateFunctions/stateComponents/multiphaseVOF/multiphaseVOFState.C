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

#include "multiphaseVOFState.H"
#include "stateComponents/singlePhaseFluid/singlePhaseFluidState.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
namespace stateFunctions
{

defineTypeNameAndDebug(multiphaseVOFState, 0);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void multiphaseVOFState::initialiseTwoPhase()
{
    // fieldMaps
    stateDict_.subDict("fieldMaps").add
    (
        word("alpha."+phases_[0]),
        input().subDict("fieldMaps").lookup("alpha")
    );

    stateDict_.subDict("fieldMaps").remove("alpha");

    // find grad(alpha) in defaults and write to state using proper phase name
    addAlphaTypeGrad(phases_[0]);

    // multiphaseVOF specific transportProperties entries
    if (compressibility() == ctIncomp)
    {
        dictionary& tpp(constant().subDict("transportProperties"));
        tpp.add("phases", phases_, false);
    }


    //assemble material properties
    List<dictionary> materialProperties(assembleMaterialProperties());

    // inject material properties
    forAll(phases_, pi)
    {
        word phase(phases_[pi]);

        if (compressibility() == ctIncomp)
        {
            //add phase sub-dicts
            dictionary& tpp(constant().subDict("transportProperties"));
            tpp.add
            (
                phase,
                singlePhaseFluidState::incompressibleProperties
                (
                    phase,
                    materialProperties[pi]
                )
            );
        }
        else if (compressibility() == ctComp)
        {
            //add thermophysicalProperties if not present
            constant().add
            (
                word("thermophysicalProperties."+phase),
                dictionary(),
                false
            );
            dictionary& tpp
            (
                constant().subDict("thermophysicalProperties."+phase)
            );


            tpp.merge
            (
                singlePhaseFluidState::compressibleProperties
                (
                    phase,
                    materialProperties[pi]
                )
            );
        }
        else
        {
            FatalErrorInFunction
                << "Invalid compressibility type: "
                << compressibility() << exit(FatalError);
        }
    }

    //add sigmas
    HashTable<scalar, interfacePair, interfacePair::hash> sigmas
    (
        assembleBinaryData<scalar>("sigma", 0.0, materialProperties)
    );

    //post entries
    if (compressibility() == ctIncomp)
    {
        dictionary& tpp(constant().subDict("transportProperties"));

        tpp.add
        (
            "sigma",
            sigmas[interfacePair(phases_[0], phases_[1])],
            false
        );
    }
    else if (compressibility() == ctComp)
    {
        dictionary& tpp(constant().subDict("thermophysicalProperties"));

        tpp.add
        (
            "sigma",
            sigmas[interfacePair
            (phases_[0],
            phases_[1])],
            false
        );
    }

}


void multiphaseVOFState::initialiseMultiPhase()
{

    //assemble material properties
    List<dictionary> materialProperties(assembleMaterialProperties());

    //add sigmas
    HashTable<scalar, interfacePair, interfacePair::hash> sigmas
    (
        assembleBinaryData<scalar>("sigma", 0.0, materialProperties)
    );

    forAll(phases_, pi)
    {
        // fieldMaps
        stateDict_.subDict("fieldMaps").add
        (
            word("alpha." + phases_[pi]),
            input().subDict("fieldMaps").lookup("alpha")
        );

        // find grad(alpha) in defaults and write to state using proper phase name
        addAlphaTypeGrad(phases_[pi]);
    }


    // multiphaseVOF specific transportProperties entries
    if (compressibility() == ctIncomp)
    {
        dictionary& tpp(constant().subDict("transportProperties"));

        HashTable<dictionary,word> phaseDicts(phases_.size());

        forAll(phases_, pi)
        {
            phaseDicts.insert
            (
                phases_[pi],
                singlePhaseFluidState::incompressibleProperties
                (
                    phases_[pi],
                    materialProperties[pi]
                )
            );
        }

        tpp.add("phases", phaseDicts, false);

        tpp.add("sigmas", sigmas);
    }
    else if (compressibility() == ctComp)
    {
        // inject material properties
        forAll(phases_, pi)
        {
            word phase(phases_[pi]);

            //add thermophysicalProperties if not present
            constant().add
            (
                word("thermophysicalProperties."+phase),
                dictionary(),
                false
            );

            dictionary& tpp
            (
                constant().subDict("thermophysicalProperties."+phase)
            );

            tpp.merge
            (
                singlePhaseFluidState::compressibleProperties
                (
                    phases_[pi],
                    materialProperties[pi]
                )
            );
        }

        dictionary& tpp(constant().subDict("thermophysicalProperties"));
        tpp.add("sigmas", sigmas);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

multiphaseVOFState::multiphaseVOFState
(
    word region,
    const dictionary& input,
    const dictionary& defaults,
    const stateFunction& master,
    const stateIndex& index,
    word meshName
)
:
    multiphaseFluidState
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


void multiphaseVOFState::initialise()
{
    // combined entries
    if (compressibility() == ctIncomp)
    {
        //add transportProperties if not present
        constant().add("transportProperties", dictionary(), false);
    }
    else if (compressibility() == ctComp)
    {
        //add thermophysicalProperties if not present
        constant().add("thermophysicalProperties", dictionary(), false);
        dictionary& tpp(constant().subDict("thermophysicalProperties"));

        //add phases if not already added
        tpp.add("phases", phases_, false);

        //add default pMin if not already added
        tpp.add("pMin", 10000, false);
    }
    else
    {
        FatalErrorInFunction
            << "Invalid compressibility type: "
            << compressibility() << exit(FatalError);
    }


    if (phases_.size() == 2)
    {
        initialiseTwoPhase();
    }
    else
    {
        initialiseMultiPhase();
    }

    multiphaseFluidState::initialise();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace stateFunction
} // End namespace Foam

// ************************************************************************* //
