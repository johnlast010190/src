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
    (c) 2015 OpenCFD Ltd.
    (c) 2015 OpenFOAM Foundation
    (c) 2017-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "CourantNumberLimit.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fvMatrices/fvMatrix/fvMatrix.H"
#include "fields/volFields/volFields.H"
#include "finiteVolume/fvc/fvcSurfaceIntegrate.H"
#include "fields/fvPatchFields/basic/zeroGradient/zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(CourantNumberLimit, 0);
    addToRunTimeSelectionTable
    (
        option,
        CourantNumberLimit,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::CourantNumberLimit::CourantNumberLimit
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const objectRegistry& obr
)
:
    cellSetOption(name, modelType, dict, obr)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::CourantNumberLimit::correct(volVectorField& U)
{
    //calculate CFL number
    //based on mean cell thickness in flow direction
    const fvMesh& mesh(U.mesh());
    const Time& runTime(mesh.time());

    Info<< "phiName_: " <<phiName_<<endl;


    const surfaceScalarField phi
    (
        mesh.lookupObject<surfaceScalarField>(phiName_)
    );


    volScalarField Co
    (
        IOobject
        (
            "CourantNumber",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dimensionedScalar("0", dimless, 0),
        zeroGradientFvPatchScalarField::typeName
    );


    if (phi.dimensions() == dimensionSet(1, 0, -1, 0, 0))
    {
        // compressible
        const volScalarField& rho =
            mesh.lookupObject<volScalarField>(rhoName_);

        Co.ref() =
            (0.5*runTime.deltaT())
           *fvc::surfaceSum(mag(phi))().internalField()
           /(rho*mesh.V());
        Co.correctBoundaryConditions();
    }
    else if (phi.dimensions() == dimensionSet(0, 3, -1, 0, 0))
    {
        Co.ref() =
            (0.5*runTime.deltaT())
           *fvc::surfaceSum(mag(phi))().internalField()
           /mesh.V();
        Co.correctBoundaryConditions();
    }
    else
    {
        FatalErrorInFunction
            << "Incorrect dimensions of phi: " << phi.dimensions()
            << abort(FatalError);
    }



    // constrain to interface if phase field is present
    if (mesh.foundObject<volScalarField>(phaseName_))
    {
        const volScalarField& alpha
        (
            mesh.lookupObject<volScalarField>(phaseName_)
        );

        dimensionedScalar offset("offset", dimless, 0.01);
        const DimensionedField<scalar, volMesh>& aif
        (
            alpha.internalField()
        );

        Co.ref()
            *= pos0(aif - offset)*pos0(1 - offset - aif);
    }

    //Clip velocity
    label nCellsClipped = 0;

    forAll(Co, i)
    {
        if (Co[i] > cflMax_)
        {
            nCellsClipped++;
            vector Uold = U[i];
            U[i] *= cflMax_/Co[i];

            U[i] += (Uold - U[i])*residualScale_;
        }
    }

    // handle boundaries in the case of 'all'
    if (selectionMode_ == smAll)
    {
        volVectorField::Boundary& bf = U.boundaryFieldRef();

        forAll(bf, patchi)
        {
            fvPatchVectorField& Up = bf[patchi];
            const fvPatchScalarField& Cop
            (
                Co.boundaryField()[patchi]
            );

            if (!Up.fixesValue())
            {
                forAll(Up, facei)
                {
                    if (Cop[facei] > cflMax_)
                    {
                        vector Uold = Up[facei];
                        Up[facei] *= cflMax_/Co[facei];

                        U[facei] += (Uold - Up[facei])*residualScale_;
                    }
                }
            }
        }
    }


    reduce(nCellsClipped, sumOp<label>());

    if (nCellsClipped > 0)
    {
        Info<< "Clipping velocity in "
             << nCellsClipped << " cells." << endl;
    }
}


bool Foam::fv::CourantNumberLimit::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        cflMax_ = readScalar(coeffs_.lookup("cflMax"));

        residualScale_
            = coeffs_.lookupOrDefault<scalar>("residualScale", 0.0);

        phaseName_ = coeffs_.lookupOrDefault<word>("phaseName", "none");

        rhoName_ = coeffs_.lookupOrDefault<word>("rhoName", "none");

        phiName_ = coeffs_.lookupOrDefault<word>("phiName", "phi");

        if (coeffs_.found("UName"))
        {
            word UName(coeffs_.lookup("UName"));
            fieldNames_ = wordList(1, UName);
        }
        else
        {
            fieldNames_ = wordList(1, "U");
        }

        applied_.setSize(fieldNames_.size(), false);

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
