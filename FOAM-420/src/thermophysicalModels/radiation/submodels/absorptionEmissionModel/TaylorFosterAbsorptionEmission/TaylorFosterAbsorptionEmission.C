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
    (c) 2017 Esi Ltd.
    (c) 2011-2017 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "submodels/absorptionEmissionModel/TaylorFosterAbsorptionEmission/TaylorFosterAbsorptionEmission.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(TaylorFosterAbsorptionEmission, 0);

        addToRunTimeSelectionTable
        (
            constantAbsorptionEmission,
            TaylorFosterAbsorptionEmission,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::TaylorFosterAbsorptionEmission::TaylorFosterAbsorptionEmission
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    constantAbsorptionEmission(dict, mesh, typeName),
    rhoName_("rho"),
    sootName_(coeffsDict_.lookupOrDefault<word>("sootName", "smoke")),
    TName_(coeffsDict_.lookupOrDefault<word>("TName", "T"))
{
    const surfaceScalarField& phi
    = mesh_.lookupObject<surfaceScalarField>("phi");

    if (phi.dimensions() == dimVelocity * dimArea)
    {
        rhoName_="rhoEff";
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::TaylorFosterAbsorptionEmission::~TaylorFosterAbsorptionEmission()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiation::TaylorFosterAbsorptionEmission::aCont(const label bandI) const
{
    const volScalarField& rho = mesh_.lookupObject<volScalarField>(rhoName_);
    const volScalarField& smoke = mesh_.lookupObject<volScalarField>(sootName_);
    const volScalarField& T = mesh_.lookupObject<volScalarField>(TName_);

    dimensionedScalar b1("b1", dimensionSet(-1,2,0,0,0,0,0), 1232.4);
    dimensionedScalar bT("bT", dimensionSet(0,0,0,-1,0,0,0), 4.8E-4);
    dimensionedScalar constT("constT", dimensionSet(0,0,0,1,0,0,0), 2000);

    // gas + smoke absorption coefficient
    tmp<volScalarField> tags
        (
            new volScalarField("tags", a_ + b1*rho*smoke*(1+bT*(T-constT)))
    );

    // clippers and warning
    volScalarField& ags = tags.ref();

    bool upperbound = false;
    bool lowerbound = false;

    //cells
    forAll(ags,cellI)
    {
    if (ags[cellI] > 1.0)
    {
        ags[cellI] = 1.0;
        upperbound = true;
    }
    if (ags[cellI] < 0.0)
    {
        ags[cellI] = 0.0;
        lowerbound = true;
    }
    }

    //boundary
    forAll(ags.boundaryField(), patchI)
    {
    scalarField& absBoundary = ags.boundaryFieldRef()[patchI];

    forAll(absBoundary,faceI)
    {
             if (absBoundary[faceI] > 1.0)
             {
                 absBoundary[faceI] = 1.0;
        upperbound = true;
             }
             if (absBoundary[faceI] < 0.0)
             {
                 absBoundary[faceI] = 0.0;
        lowerbound = true;
             }
    }
    }

    if (upperbound)
    {
    WarningInFunction
        <<"absorption coeff. clipped to 1.0 "<<endl;
    }
    if (lowerbound)
    {
    WarningInFunction
        <<"absorption coeff. clipped to 0.0 "<<endl;
    }

    return tags;

}

// ************************************************************************* //
