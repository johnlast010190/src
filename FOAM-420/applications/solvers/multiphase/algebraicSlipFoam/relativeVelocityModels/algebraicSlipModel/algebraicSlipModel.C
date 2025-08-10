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
    (c) 2014-2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "algebraicSlipModel.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace relativeVelocityModels
{
    defineTypeNameAndDebug(algebraicSlipModel, 0);
    addToRunTimeSelectionTable(relativeVelocityModel, algebraicSlipModel, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::relativeVelocityModels::algebraicSlipModel::algebraicSlipModel
(
    const dictionary& dict,
    const immiscibleIncompressibleTwoPhaseMixture& mixture
)
:
    relativeVelocityModel(dict, mixture),
    dd_("d", dimLength, dict.lookupOrDefault("d", 0.)),
    Cd_("Cd", dimless, dict.lookupOrDefault("Cd", 0.44)),
    U_(mixture.U()),
    rhoPhi_(U_.db().lookupObject<surfaceScalarField>("rhoPhi")),
    g_(
        IOobject
        (
            "g",
            U_.mesh().time().constant(),
            U_.mesh(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    residualAlpha_("residualAlpha", alphad_.dimensions(), dict.lookupOrDefault("residualAlpha", 1e-06))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::relativeVelocityModels::algebraicSlipModel::~algebraicSlipModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::relativeVelocityModels::algebraicSlipModel::correct()
{
    const volScalarField limitedAlphad
    (
        "limitedAlphad",
        min(max(alphad_, scalar(0)), scalar(1))
    );

    // get mixture density
    volScalarField rhom(mixture_.rho()());

//Info<< "diameter: " << dd_ << endl;
//Info<< "drag coeff: " << Cd_ << endl;
//Info<< "rho continuous: " << rhoc_ << endl;
//Info<< "rho disperse: " << rhod_ << endl;
//Info<< "gravity: " << g_ << endl;
//Info<< "residualAlpha: " << residualAlpha_ << endl;

    // compute relative velocity based on algebraic slip model (compute mag(Ur)*Ur)
    volVectorField UrelSqr( -(4./3.)*dd_/(rhoc_*Cd_)*(rhod_-rhom)*((fvc::ddt(rhom,U_) + fvc::div(rhoPhi_,U_))/rhom - g_) );

    // scale by square root
    Urel_ = UrelSqr/(sqrt(mag(UrelSqr)) + dimensionedScalar("U", dimVelocity, VSMALL));

    //scalar beta = 0.75;
    //Urel_ *= pos(beta - limitedAlphad); //min(1., 1./beta*(1.-limitedAlphad));

    // set Urel boundaryField to zeroGradient
    Urel_.correctBoundaryConditions();

    // set Urel to zero at walls
    forAll(Urel_.boundaryField(), patchI)
    {
        //Info<<"patch " << patchI << " boundary type: " << U_.boundaryField()[patchI].type() << endl;
        //Info<<"patch " << patchI << " patch type: " << U_.boundaryField()[patchI].patch().type() << endl;
        if (U_.boundaryField()[patchI].patch().type() == "wall")
        {
            forAll(Urel_.boundaryField()[patchI], faceI)
            {
                Urel_.boundaryFieldRef()[patchI][faceI] = vector::zero; //U_.boundaryField()[patchI][faceI];
            }
        }
        //Info<<"patch " << patchI << " patch name: " << U_.boundaryField()[patchI].patch().name() << endl;
        // set to a constant at the inlet:
        if (U_.boundaryField()[patchI].patch().name() == "inlet")
        {
            forAll(Urel_.boundaryField()[patchI], faceI)
            {
                Urel_.boundaryFieldRef()[patchI][faceI] = Foam::vector(0,0.1,0); //U_.boundaryField()[patchI][faceI];
            }
        }/*
        // set Urel to zero at outlet !cells!
        if (U_.boundaryField()[patchI].patch().name() == "outlet")
        {
            forAll(Urel_.boundaryField()[patchI], faceI)
            {
                Urel_.boundaryField()[patchI][faceI] = Foam::vector(0,0,0);

                label fIstart = Urel_.mesh().boundaryMesh()[patchI].start();
                const label own = Urel_.mesh().owner()[faceI + fIstart];
                Urel_.internalField()[own] = Foam::vector(0,0,0);
            }
        }*/
    }
    Info<< "min/max(Urel): " << min(mag(Urel_)) << " , " << max(mag(Urel_)) << endl;
}


// ************************************************************************* //
