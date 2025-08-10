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

#include "solidState.H"
#include "db/dictionary/functionEntries/includeEtcEntry/includeEtcEntry.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
namespace stateFunctions
{

defineTypeNameAndDebug(solidState, 0);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

solidState::solidState
(
    word region,
    const dictionary& inputDict,
    const dictionary& defaults,
    const stateFunction& master,
    const stateIndex& index,
    word meshName
)
:
    regionState
    (
        region,
        inputDict,
        defaults,
        master,
        index,
        meshName
    ),
    materialName_(word::null)
{
    //check that only one material is specified

    wordList materials(input().lookup("materials"));

    if (materials.size() != 1)
    {
        FatalErrorInFunction
            << "Multiple materials specified for a single phase state: "
            << materials << exit(FatalError);
    }

    materialName_ = materials[0];
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void solidState::initialise()
{
    const wordList matLibStates =
        stateFunction::defaults().lookup<wordList>("materialLibStates");
    const bool isMaterials = matLibStates.found(stateFunction::stateType());

    //merge materials
    dictionary materialProperties;

    materialProperties.merge
    (
        defaults().subDict("materialProperties").subOrEmptyDict
        (
            isMaterials ? materialName_ + "Material" : materialName_
        )
    );

    if (input().found("materialProperties"))
    {
        materialProperties.merge
        (
            input().subDict("materialProperties").subOrEmptyDict(materialName_)
        );
    }
    if (isMaterials)
    {
        system().add("materialProperties", dictionary(), false);
        dictionary& matDict(system().subDict("materialProperties"));
        matDict.add("materialName", materialName_);

        matDict.merge(materialProperties);

        // Check that required entries are present
        word matType =
            materialProperties.lookupOrDefault<word>("materialType", "solid");
        if (matType != "solid")
        {
            FatalIOErrorInFunction(input().subDict("materialProperties"))
                << "Expected material type 'solid' but got '"
                << matType << "' instead." << nl << exit(FatalIOError);
        }
        matDict.add
        (
            "materialType",
            matType
        );

        matDict.add
        (
            "equationOfState",
            materialProperties.lookup<word>("equationOfState")
        );
        matDict.add
        (
            "thermodynamics",
            materialProperties.lookup<word>("thermodynamics")
        );
        matDict.add
        (
            "kappaModel",
            materialProperties.lookup<word>("kappaModel")
        );
        matDict.add
        (
            "energy",
            materialProperties.lookupOrDefault<word>
            (
                "energy", "sensibleEnthalpy"
            )
        );

        matDict.add
        (
            "molWeight",
            materialProperties.lookup<scalar>("molWeight")
        );
        if (materialProperties.found("alphahModel"))
        {
            matDict.add
            (
                "alphahModel",
                materialProperties.lookup<word>("alphahModel")
            );
            matDict.add
            (
                "alphahModelCoeffs",
                materialProperties.subOrEmptyDict("alphahModelCoeffs")
            );
        }
        matDict.add("thermodynamicsCoeffs", materialProperties.subDict("thermodynamicsCoeffs"));
        matDict.add("kappaModelCoeffs", materialProperties.subDict("kappaModelCoeffs"));
        matDict.add
        (
            "equationOfStateCoeffs",
            materialProperties.subOrEmptyDict("equationOfStateCoeffs")
        );
    }
    else
    {
        constant().add("thermophysicalProperties", dictionary(), false);
        dictionary& thermDict(constant().subDict("thermophysicalProperties"));
        thermDict.add("materialName", materialName_);

        thermDict.add
        (
            "thermoType",
            materialProperties.subDict("thermoType")
        );

        //add mixture properties (pure mixture)
        thermDict.add("mixture", dictionary());
        dictionary& mixture(thermDict.subDict("mixture"));

        mixture.add("specie", materialProperties.subDict("specie"));
        mixture.add("thermodynamics", materialProperties.subDict("thermodynamics"));
        mixture.add("transport", materialProperties.subDict("transport"));
        mixture.add
        (
            "equationOfState",
            materialProperties.subOrEmptyDict("equationOfState")
        );
    }

    regionState::initialise();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace stateFunction
} // End namespace Foam

// ************************************************************************* //
