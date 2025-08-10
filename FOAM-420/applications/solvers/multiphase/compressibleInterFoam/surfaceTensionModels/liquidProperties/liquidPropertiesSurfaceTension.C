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
    (c) 2017 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "liquidPropertiesSurfaceTension.H"
#include "rhoThermo/liquidThermo.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceTensionModels
{
    defineTypeNameAndDebug(liquidProperties, 0);
    addToRunTimeSelectionTable
    (
        surfaceTensionModel,
        liquidProperties,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceTensionModels::liquidProperties::liquidProperties
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    surfaceTensionModel(mesh),
    phaseName_(dict.lookup("phase"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfaceTensionModels::liquidProperties::~liquidProperties()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::surfaceTensionModels::liquidProperties::sigma() const
{
    const heRhoThermopureMixtureliquidProperties& thermo =
        mesh_.lookupObject<heRhoThermopureMixtureliquidProperties>
        (
             IOobject::groupName(basicThermo::dictName, phaseName_)
        );

    const Foam::liquidProperties& liquid = thermo.mixture().properties();

    tmp<volScalarField> tsigma
    (
        new volScalarField
        (
            IOobject
            (
                "sigma",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimSigma
        )
    );
    volScalarField& sigma = tsigma.ref();

    const volScalarField& T = thermo.T();
    const volScalarField& p = thermo.p();

    volScalarField::Internal& sigmai = sigma;
    const volScalarField::Internal& pi = p;
    const volScalarField::Internal& Ti = T;

    forAll(sigmai, celli)
    {
        sigmai[celli] = liquid.sigma(pi[celli], Ti[celli]);
    }

    volScalarField::Boundary& sigmaBf = sigma.boundaryFieldRef();
    const volScalarField::Boundary& pBf = p.boundaryField();
    const volScalarField::Boundary& TBf = T.boundaryField();

    forAll(sigmaBf, patchi)
    {
        scalarField& sigmaPf = sigmaBf[patchi];
        const scalarField& pPf = pBf[patchi];
        const scalarField& TPf = TBf[patchi];

        forAll(sigmaPf, facei)
        {
            sigmaPf[facei] = liquid.sigma(pPf[facei], TPf[facei]);
        }
    }

    return tsigma;
}


bool Foam::surfaceTensionModels::liquidProperties::readDict
(
    const dictionary& dict
)
{
    return true;
}


bool Foam::surfaceTensionModels::liquidProperties::writeData
(
    Ostream& os
) const
{
    if (surfaceTensionModel::writeData(os))
    {
        return os.good();
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
