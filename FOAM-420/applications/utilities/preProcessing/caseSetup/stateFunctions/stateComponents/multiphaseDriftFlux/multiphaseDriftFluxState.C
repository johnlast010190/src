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
    (c) 2021 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "multiphaseDriftFluxState.H"
#include "stateComponents/singlePhaseFluid/singlePhaseFluidState.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
namespace stateFunctions
{

defineTypeNameAndDebug(multiphaseDriftFluxState, 0);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


Xfer<dictionary> multiphaseDriftFluxState::eulerPhaseProperties
(
    const word& materialName,
    const dictionary& mat
) const
{
    dictionary matDict
    (
        singlePhaseFluidState::incompressibleProperties(materialName, mat)
    );

    word relativeVelocityModel
    (
        mat.lookupOrDefault<word>
        (
            "relativeVelocityModel",
            "general"
        )
    );
    matDict.add
    (
        "relativeVelocityModel",
        relativeVelocityModel
    );

    Switch useTotalSolids
    (
        mat.lookupOrDefault<Switch>
        (
            "useTotalSolids",
            false
        )
    );
    matDict.add
    (
        "useTotalSolids",
        useTotalSolids
    );

    word relativeVelocityModelCoeffs(relativeVelocityModel + "Coeffs");
    dictionary generalCoeffsDef;
    generalCoeffsDef.add("V0", vector(0,-0.00135,0));
    generalCoeffsDef.add("a", 607.6);
    generalCoeffsDef.add("a1", 6250);
    generalCoeffsDef.add("residualAlpha", 6.9e-06);

    if (mat.found(relativeVelocityModelCoeffs))
    {
        matDict.add
        (
            relativeVelocityModelCoeffs,
            mat.subDict(relativeVelocityModelCoeffs)
        );

    }
    else
    {
        matDict.add
        (
            relativeVelocityModelCoeffs,
            generalCoeffsDef
        );
    }

    return matDict.xfer();
}

void multiphaseDriftFluxState::initialiseMultiPhase()
{

    //assemble material properties
    List<dictionary> materialProperties(assembleMaterialProperties());

    Info<< "field input: " << input().subDict("fieldMaps") << endl;

    forAll(phases_, pi)
    {
        // fieldMaps
        stateDict_.subDict("fieldMaps").add
        (
            word("alpha."+phases_[pi]),
            input().subDict("fieldMaps").lookup("alpha")
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

        scalar alphaMax = 1.0;
        word contiPhase = "water";
        word gasPhase = "air";
        word transMod = "BokilBewtra";

        if (input().subDict("constant").isDict("transportProperties"))
        {
            const dictionary& intpp = input().subDict("constant").subDict("transportProperties");

            if (intpp.found("alphaMax"))
            {
                alphaMax = readScalar(intpp.lookup("alphaMax"));
            }
            if (intpp.found("continuousPhase"))
            {
                contiPhase = word(intpp.lookup("continuousPhase"));
            }
            if (intpp.found("gasPhase"))
            {
                gasPhase = word(intpp.lookup("gasPhase"));
                tpp.add("gasPhase", gasPhase);
            }

            if (intpp.found("mixture"))
            {
                if (intpp.subDict("mixture").found("transportModel"))
                {
                    transMod = word(intpp.subDict("mixture").lookup("transportModel"));
                }
            }
        }

        word transModCoeffs(transMod + "Coeffs");
        dictionary transModCoeffsDef;
        transModCoeffsDef.add("coeff", 0.00327);
        transModCoeffsDef.add("exponent", 191.4);
        transModCoeffsDef.add("muMax", 1);
        transModCoeffsDef.add("alphac", 0.00048276);

        dictionary viscModel;
        viscModel.add("transportModel", transMod);

        if (input().subDict("constant").isDict("transportProperties"))
        {
            const dictionary& intpp = input().subDict("constant").subDict("transportProperties");

            if (intpp.found("mixture"))
            {
                if (intpp.subDict("mixture").isDict(transModCoeffs))
                {
                    viscModel.add
                    (
                        transModCoeffs,
                        intpp.subDict("mixture").subDict(transModCoeffs)
                    );
                }
            }
        }

        if (!viscModel.found(transModCoeffs))
        {
            viscModel.add(transModCoeffs, transModCoeffsDef);
        }

        tpp.add("continuousPhase", contiPhase);
        tpp.add("alphaMax", alphaMax);
        tpp.add("mixture", viscModel, false);
    }
    else if (compressibility() == ctComp)
    {
        FatalErrorInFunction
            << "Invalid compressibility type: "
            << compressibility() << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

multiphaseDriftFluxState::multiphaseDriftFluxState
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


void multiphaseDriftFluxState::initialise()
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

void multiphaseDriftFluxState::correct()
{
//    //remove generic U entry from fields as it has been specialised
//    dictionary& fields(stateDict_.subDict("fields"));

    multiphaseFluidState::correct();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace stateFunction
} // End namespace Foam

// ************************************************************************* //
