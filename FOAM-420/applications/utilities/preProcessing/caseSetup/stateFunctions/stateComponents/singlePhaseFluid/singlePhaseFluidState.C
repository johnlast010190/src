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

#include "singlePhaseFluidState.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "stateComponents/multiphaseMultiComponent/multiphaseMultiComponentState.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
namespace stateFunctions
{

defineTypeNameAndDebug(singlePhaseFluidState, 0);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Xfer<dictionary> singlePhaseFluidState::incompressibleProperties
(
    const word& matName,
    const dictionary& materialProperties
)
{
    dictionary matDict;
    matDict.add("materialName", matName);

    word tpmdl = materialProperties.lookup("transportModel");

    matDict.add("transportModel", tpmdl, false);

    if (materialProperties.found(word(tpmdl+"Coeffs")))
    {
        const dictionary& coeffDict(materialProperties.subDict(word(tpmdl+"Coeffs")));
        matDict.add(word(tpmdl+"Coeffs"), coeffDict, false);
    }
    else
    {
        matDict.add(word(tpmdl+"Coeffs"), dictionary(), false);
    }

    dimensionedScalar rho(word::null, materialProperties.lookup("rho"));
    matDict.add("rho", rho, false);

    const wordList properties = {"nu", "TRef", "beta", "Cp", "Prt", "kappa", "lambda"};
    forAll(properties, i)
    {
        const word& propName = properties[i];
        if (materialProperties.found(propName))
        {
            dimensionedScalar propI(word::null, materialProperties.lookup(propName));
            matDict.add(propName, propI, false);
        }
    }

    if (!materialProperties.found("nu"))
    {
        dimensionedScalar mu(word::null, materialProperties.lookup("mu"));
        dimensionedScalar nu
        (
            word::null,
            mu.dimensions()/rho.dimensions(),
            mu.value()/rho.value()
        );
        matDict.add("nu", nu, false);
    }

    return matDict.xfer();
}


Xfer<dictionary> singlePhaseFluidState::compressibleProperties
(
    const word& matName,
    const dictionary& materialProperties
)
{
    dictionary matDict;
    matDict.add("materialName", matName);
    matDict.add("thermoType", materialProperties.subDict("thermoType"));

    //add mixture properties (pure mixture)
    matDict.add("mixture", dictionary());
    dictionary& mixture(matDict.subDict("mixture"));

    mixture.add("specie", materialProperties.subDict("specie"));
    mixture.add("thermodynamics", materialProperties.subDict("thermodynamics"));
    mixture.add("transport", materialProperties.subDict("transport"));
    mixture.add("equationOfState", materialProperties.subOrEmptyDict("equationOfState"));

    return matDict.xfer();
}


Xfer<dictionary> singlePhaseFluidState::setupMaterialProperties
(
    const word& matName,
    const dictionary& materialProperties
)
{
    dictionary matDict;
    matDict.add("materialName", matName);
    const wordList properties
    ({
        "materialType",
        "energy",
    });
    forAll(properties, i)
    {
        matDict.add(properties[i], materialProperties.lookup<word>(properties[i]));
    }

    matDict.merge
    (
        multiphaseMultiComponentState::speciesProperties
        (
            materialProperties.lookup<word>("materialType"),
            materialProperties
        )
    );

    // Buoyancy model is optional entry
    if (input().lookup<wordList>("state").found("buoyant"))
    {
        if (input().isDict("buoyantSteady"))
        {
            matDict.add
            (
                "buoyancyModel",
                input().subDict("buoyantSteady").subDict("materialProperties").lookupOrDefault<word>("buoyancyModel", "rhoModel")
            );
        }
        else
        {
            matDict.add("buoyancyModel", "rhoModel");
        }
    }

    // Dictionaries to add
    const wordList dictsToAdd
    ({
        "referenceFields"
    });
    forAll(dictsToAdd, i)
    {
        matDict.add(dictsToAdd[i], materialProperties.subOrEmptyDict(dictsToAdd[i]));
    }

    if (matDict.found("thermoType"))
    {
        matDict.remove("thermoType");
    }

    return matDict.xfer();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

singlePhaseFluidState::singlePhaseFluidState
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
    materialName_(word::null),
    componentNames_(word::null)
{
    wordList materials(stateFunction::input().lookup("materials"));

    if (materials.size() > 1)
    {
        // compressible multi component
        componentNames_ = materials;

    }
    else if (materials.size() == 1)
    {
        // single component incompressible and compressible flows
        materialName_ = materials.first();
    }
    else
    {
        FatalErrorInFunction
            << "At least one material must be defined"
            <<" for this state." << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void singlePhaseFluidState::initialiseMaterials
(
    dictionary& materialProperties,
    const dictionary& matDefaults,
    const dictionary& matInputs,
    const bool isMatInInput
)
{
    // This name is use to make difference between thermo and material library
    const word matName(materialName_ + "Material");

    // Pick up referenceFields
    materialProperties.merge(matDefaults.subDict("referenceValues"));
    materialProperties.merge(matDefaults.subOrEmptyDict(matName));
    if (isMatInInput)
    {
        materialProperties.merge(matInputs.subOrEmptyDict(materialName_));
    }

    const dictionary& stateDefaluts =
        defaults().subDict("solvers").subDict(stateType());
    if (stateDefaluts.found("materialProperties"))
    {
        materialProperties.merge(stateDefaluts.subDict("materialProperties"));
    }

    if (compressibility() == ctComp || compressibility() == ctIncomp)
    {
        system().add("materialProperties", dictionary(), false);
        system().subDict("materialProperties").merge
        (
            setupMaterialProperties(materialName_, materialProperties)
        );
        if (isMatInInput)
        {
            system().subDict("materialProperties").merge
            (
                matInputs.subOrEmptyDict(materialName_)
            );
        }
    }
    else
    {
        FatalErrorInFunction
            << "Invalid compressibility type: "
            << compressibility() << exit(FatalError);
    }
}


void singlePhaseFluidState::initialiseThermo
(
    dictionary& materialProperties,
    const dictionary& matDefaults,
    const dictionary& matInputs,
    const bool isMatInInput
)
{
    materialProperties.merge(matDefaults.subOrEmptyDict(materialName_));
    if (isMatInInput)
    {
        materialProperties.merge(matInputs.subOrEmptyDict(materialName_));
    }
    if (compressibility() == ctIncomp)
    {
        constant().add("transportProperties", dictionary(), false);
        constant().subDict("transportProperties").merge
        (
            incompressibleProperties(materialName_, materialProperties)
        );
        if (isMatInInput)
        {
            if (constant().found("transportProperties"))
            {
                constant().subDict("transportProperties").merge
                (
                    matInputs.subOrEmptyDict(materialName_)
                );
            }
        }
    }
    else if (compressibility() == ctComp)
    {
        constant().add("thermophysicalProperties", dictionary(), false);

        constant().subDict("thermophysicalProperties").merge
        (
            compressibleProperties(materialName_, materialProperties)
        );
        // Merge other user-defined specifications inside the
        // thermophysical properties. Some of the dictionaries must be
        // removed from the top level because they have been already added to
        // the mixture
        if (isMatInInput)
        {
            if (constant().found("thermophysicalProperties"))
            {
                dictionary croppedInputDict = matInputs.subOrEmptyDict(materialName_);
                const wordList checkRemove = {"specie", "thermodynamics", "transport", "equationOfState"};
                forAll(checkRemove, i)
                {
                    if (croppedInputDict.found(checkRemove[i]))
                    {
                        croppedInputDict.remove(checkRemove[i]);
                    }
                }
                constant().subDict("thermophysicalProperties").merge(croppedInputDict);
            }
        }
    }
    else
    {
        FatalErrorInFunction
            << "Invalid compressibility type: "
            << compressibility() << exit(FatalError);
    }
}


void singlePhaseFluidState::initialise()
{
    if (materialName() == word::null)
    {
        FatalErrorInFunction
            << "Only one material can be defined"
            <<" for this state. Current materials defined are: "
            << componentNames() << exit(FatalError);
    }

    const bool isMaterials =
        defaults().lookup<wordList>("materialLibStates").found(stateType());

    // merge materials
    const bool isMatInInput(input().found("materialProperties"));
    const dictionary& matDefaults = defaults().subDict("materialProperties");
    // This is already input joined with the state defaults
    const dictionary matInputs = input().subOrEmptyDict("materialProperties");

    dictionary materialProperties;
    if (isMaterials)
    {
        initialiseMaterials
        (
            materialProperties,
            matDefaults,
            matInputs,
            isMatInInput
        );
    }
    else
    {
        initialiseThermo
        (
            materialProperties,
            matDefaults,
            matInputs,
            isMatInInput
        );
    }

    turbulenceModelState::initialise();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace stateFunction
} // End namespace Foam

// ************************************************************************* //
