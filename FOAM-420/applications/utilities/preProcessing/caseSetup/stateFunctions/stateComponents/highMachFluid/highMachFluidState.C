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
    (c) 2016-2022 Esi Ltd

\*---------------------------------------------------------------------------*/

#include "highMachFluidState.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
namespace stateFunctions
{

defineTypeNameAndDebug(highMachFluidState, 0);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

highMachFluidState::highMachFluidState
(
    word region,
    const dictionary& input,
    const dictionary& defaults,
    const stateFunction& master,
    const stateIndex& index,
    word meshName
)
:
    multiComponentFluidState
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


void highMachFluidState::initialise()
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
            isMaterials ? word(materialName()) + "Material" : word(materialName())
        )
    );
    if (isMaterials)
    {
        //set default thermoType to sensibleInternalEnergy
        materialProperties.add("materialType","fluid",true); // not defalut psi thermo anymore
        materialProperties.add("energy","sensibleInternalEnergy",true);

        materialProperties.merge(defaults().subDict("materialProperties").subDict("referenceValues"));
        materialProperties.merge
        (
            input().subOrEmptyDict("materialProperties").subOrEmptyDict(materialName())
        );
        if (compressibility() == ctComp)
        {
            system().add("materialProperties", dictionary(), false);
            system().subDict("materialProperties").merge
            (
                setupMaterialProperties
                (
                    materialName(),
                    materialProperties
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
    else
    {
        //set default thermoType to sensibleInternalEnergy
        materialProperties.add("thermoType", dictionary(), false);
        materialProperties.subDict("thermoType").add("type","hePsiThermo",true);
        materialProperties.subDict("thermoType").add("energy","sensibleInternalEnergy",true);

        if (input().found("materialProperties"))
        {
            materialProperties.merge
            (
                input().subDict("materialProperties").subOrEmptyDict(materialName())
            );
        }
        if (compressibility() == ctComp)
        {
            constant().add("thermophysicalProperties", dictionary(), false);

            constant().subDict("thermophysicalProperties").merge
            (
                singlePhaseFluidState::compressibleProperties
                (
                    materialName(),
                    materialProperties
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

    turbulenceModelState::initialise();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace stateFunction
} // End namespace Foam

// ************************************************************************* //
