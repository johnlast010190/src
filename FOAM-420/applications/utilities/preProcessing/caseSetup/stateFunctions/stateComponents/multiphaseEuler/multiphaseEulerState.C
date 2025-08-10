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

#include "multiphaseEulerState.H"
#include "stateComponents/singlePhaseFluid/singlePhaseFluidState.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
namespace stateFunctions
{

defineTypeNameAndDebug(multiphaseEulerState, 0);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


Xfer<dictionary> multiphaseEulerState::eulerPhaseProperties
(
    const word& materialName,
    const dictionary& mat
) const
{
    dictionary matDict
    (
        singlePhaseFluidState::incompressibleProperties(materialName, mat)
    );

    word diameterModel
    (
        mat.lookupOrDefault<word>
        (
            "diameterModel",
            "constant"
        )
    );
    matDict.add
    (
        "diameterModel",
        diameterModel
    );

    word diameterModelCoeffs(diameterModel + "Coeffs");
    dictionary constantCoeffsDef;
    constantCoeffsDef.add("d", 1e-3);

    if (mat.found(diameterModelCoeffs))
    {
        matDict.add
        (
            diameterModelCoeffs,
            mat.subDict(diameterModelCoeffs)
        );

    }
    else
    {
        matDict.add
        (
            diameterModelCoeffs,
            constantCoeffsDef
        );
    }

    return matDict.xfer();
}

void multiphaseEulerState::initialiseMultiPhase()
{

    //assemble material properties
    List<dictionary> materialProperties(assembleMaterialProperties());

    //add sigmas
    HashTable<scalar, interfacePair, interfacePair::hash> sigmas
    (
        assembleBinaryData<scalar>("sigma", 0.0, materialProperties)
    );

    //interfaceCompression
    HashTable<scalar, interfacePair, interfacePair::hash> intfcomp
    (
        assembleBinaryData<scalar>("interfaceCompression", 1.0, materialProperties)
    );

    //virtualMass
    HashTable<scalar, interfacePair, interfacePair::hash> virtmass
    (
        assembleBinaryData<scalar>("virtualMass", 0.0, materialProperties)
    );

    //drag
    HashTable<dictionary, interfacePair, interfacePair::hash> drag
    (
        assembleBinaryDictionaries("drag", dictionary(), materialProperties)
    );

    Info<< "U input: " << input().subDict("fieldMaps") << endl;

    forAll(phases_, pi)
    {
        // fieldMaps
        stateDict_.subDict("fieldMaps").add
        (
            word("alpha."+phases_[pi]),
            input().subDict("fieldMaps").lookup("alpha")
        );

        stateDict_.subDict("fieldMaps").add
        (
            word("U."+phases_[pi]),
            input().subDict("fieldMaps").lookup("U")
        );

        // find grad(alpha) in defaults and write to state using proper phase name
        addAlphaTypeGrad(phases_[pi]);
    }

    //  specific transportProperties entries
    if (compressibility() == ctIncomp)
    {
        dictionary& tpp(constant().subDict("transportProperties"));

        HashTable<dictionary,word> phaseDicts(phases_.size());

        forAll(phases_, pi)
        {
            phaseDicts.insert
            (
                phases_[pi],
                eulerPhaseProperties(phases_[pi], materialProperties[pi])
            );
        }

        tpp.add("phases", phaseDicts, false);

        tpp.add("sigmas", sigmas);

        tpp.add("interfaceCompression", intfcomp);

        tpp.add("virtualMass", virtmass);

        tpp.add("drag", drag);

        // This is a dummy to support the Smagorinsky model
        tpp.add("transportModel", word("Newtonian"));
        tpp.add("nu", dimensionedScalar("nu", dimViscosity/dimDensity, 0.0));
    }
    else if (compressibility() == ctComp)
    {
        FatalErrorInFunction
            << "Invalid compressibility type: "
            << compressibility() << exit(FatalError);
    }

    //add dummy Newtonian for LES models
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

multiphaseEulerState::multiphaseEulerState
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


void multiphaseEulerState::initialise()
{
    // combined entries
    if (compressibility() == ctIncomp)
    {
        //add transportProperties if not present
        constant().add("transportProperties", dictionary(), false);
    }
    else if (compressibility() == ctComp)
    {
        FatalErrorInFunction
            << "Invalid compressibility type: "
            << compressibility() << exit(FatalError);
    }
    else
    {
        FatalErrorInFunction
            << "Invalid compressibility type: "
            << compressibility() << exit(FatalError);
    }


    if (phases_.size() == 2)
    {
        initialiseMultiPhase();
    }
    else
    {
        initialiseMultiPhase();
    }

    multiphaseFluidState::initialise();
}

void multiphaseEulerState::correct()
{
    //remove generic U entry from fields as it has been specialised
    dictionary& fields(stateDict_.subDict("fields"));

    if (fields.found("U"))
    {
        fields.remove("U");
    }


    multiphaseFluidState::correct();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace stateFunction
} // End namespace Foam

// ************************************************************************* //
