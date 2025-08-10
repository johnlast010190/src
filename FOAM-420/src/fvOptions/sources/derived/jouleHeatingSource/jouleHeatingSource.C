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
    (c) 2016-2017 OpenCFD Ltd.
    (c) 2019-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "jouleHeatingSource.H"
#include "fvMatrices/fvMatrices.H"
#include "finiteVolume/fvm/fvmLaplacian.H"
#include "finiteVolume/fvc/fvcGrad.H"
#include "fields/fvPatchFields/basic/zeroGradient/zeroGradientFvPatchField.H"
#include "basicThermo/basicThermo.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(jouleHeatingSource, 0);

    addToRunTimeSelectionTable
    (
        option,
        jouleHeatingSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::jouleHeatingSource::jouleHeatingSource
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const objectRegistry& obr
)
:
    option(sourceName, modelType, dict, obr)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::jouleHeatingSource::~jouleHeatingSource()
{}


// * * * * * * * * * * * * * * Member functions  * * * * * * * * * * * * * * //

bool Foam::fv::jouleHeatingSource::initialise()
{
    // Set the field name to that of the energy field from which the temperature
    // is obtained

    const basicThermo& thermo =
        obr_.lookupObject<basicThermo>(basicThermo::dictName);

    fieldNames_.setSize(1, thermo.he().name());

    applied_.setSize(fieldNames_.size(), false);

    return true;
}


void Foam::fv::jouleHeatingSource::addSourceDependencies
(
    SolveTable<solveList>& dependencies
)
{
    const word thisRegionName
    (
        obr_.name() == word::null ? mesh_.name() : obr_.name()
    );
    solveID key(fieldNames_[0], thisRegionName);
    dependencies.insert(key, solveList());
    dependencies[key].append
    (
        solveID("electrical_V", thisRegionName)
    );
}

void Foam::fv::jouleHeatingSource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    DebugInformation<< name() << ": applying source to " << eqn.psi().name() << endl;

    // Add the Joule heating source term

    const volScalarField& V =
        obr_.lookupObject<volScalarField>("electrical_V");
    const volVectorField gradV(fvc::grad(V));

    if (mesh().foundObject<volSymmTensorField>("electrical_sigma"))
    {
        // Anisotropic conductivity
        const volSymmTensorField& sigma =
            obr_.lookupObject<volSymmTensorField>("electrical_sigma");
        eqn += (sigma & gradV) & gradV;
    }
    else
    {
        const volScalarField& sigma =
            obr_.lookupObject<volScalarField>("electrical_sigma");
        eqn += (sigma*gradV) & gradV;
    }
}


// ************************************************************************* //
