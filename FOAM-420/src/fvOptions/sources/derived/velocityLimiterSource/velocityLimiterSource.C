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
    (c) 2012-2013 OpenFOAM Foundation
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "velocityLimiterSource.H"
#include "fvMesh/fvMesh.H"
#include "fvMatrices/fvMatrices.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(velocityLimiterSource, 0);
    addToRunTimeSelectionTable
    (
        option,
        velocityLimiterSource,
        dictionary
    );
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::velocityLimiterSource::addDamping(fvMatrix<vector>& eqn)
{
    //copied from OF+ velocityDampingConstraint
    const scalarField& vol = mesh_.V();
    const volVectorField& U = eqn.psi();
    scalarField& diag = eqn.diag();

    label nDamped = 0;

    scalarField damping(U.size(), 0.0);

    forAll(U, cellI)
    {
        scalar magU = mag(U[cellI]);
        if (magU > Ulimit_)
        {
            scalar scale = sqr(Foam::cbrt(vol[cellI]));

            damping[cellI] = scale*(magU-Ulimit_);
            diag[cellI] += damping[cellI];

            nDamped++;
        }
    }

    reduce(nDamped, sumOp<label>());

    printDampingInformation(nDamped, U, damping);
}


void Foam::fv::velocityLimiterSource::addDamping(fvBlockMatrix<vector>& eqn)
{
    //copied from OF+ velocityDampingConstraint
    const scalarField& vol = mesh_.V();
    const volVectorField& U = eqn.psi();

    label nCoeffs = pTraits<vector>::nComponents;
    CoeffField<vector>& diagCoeff = eqn.diag();
    tensorField& bDiag = diagCoeff.asSquare();

    label nDamped = 0;

    scalarField damping(U.size(), 0.0);

    forAll(U, cellI)
    {
        scalar magU = mag(U[cellI]);
        if (magU > Ulimit_)
        {
            scalar scale = sqr(Foam::cbrt(vol[cellI]));

            damping[cellI] = scale*(magU-Ulimit_);
            for (label iC=0; iC<nCoeffs; iC++)
            {
                bDiag[cellI][iC*nCoeffs+iC] += damping[cellI];
            }

            nDamped++;
        }
    }

    reduce(nDamped, sumOp<label>());

    printDampingInformation(nDamped, U, damping);
}


void Foam::fv::velocityLimiterSource::printDampingInformation
(
    const label& nDamped,
    const volVectorField& U,
    const scalarField& damping
) const
{
    if (verbose_)
    {
        Info<< type() << " " << name_ << " damped "
            << nDamped << " ("
            << 100*scalar(nDamped)/mesh_.globalData().nTotalCells()
            << "%) of cells" << endl;
    }
    if (outputDiagnostics_ && U.mesh().time().outputTime())
    {
        volScalarField alphaLimit
        (
            IOobject
            (
                "alphaLimit",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("al", dimVolume/dimTime, 0),
            calculatedFvPatchField<scalar>::typeName
        );

        alphaLimit.primitiveFieldRef() = damping;

        alphaLimit.write();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::velocityLimiterSource::velocityLimiterSource
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

void Foam::fv::velocityLimiterSource::constrain
(
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    if (velocityDamp_)
    {
        addDamping(eqn);
    }
}



void Foam::fv::velocityLimiterSource::constrain
(
    fvBlockMatrix<vector>& eqn,
    const label fieldI
)
{
    if (velocityDamp_)
    {
        addDamping(eqn);
    }
}


void Foam::fv::velocityLimiterSource::correct(volVectorField& U)
{
    //clip U magnitude
    if (velocityClip_)
    {
        label nCellsClipped = 0;
        vectorField Uctmp = U.internalField();

        forAll(U, i)
        {
            scalar magUi(mag(U[i]));

            if (magUi > Ulimit_)
            {
                U[i] *= Ulimit_/magUi;
                nCellsClipped++;
            }
        }
        Uctmp -= U.internalField();

        reduce(nCellsClipped, sumOp<label>());

        if (verbose_)
        {
            Info<< type() << " " << name_ << " clipped "
                << nCellsClipped << " ("
                << 100*scalar(nCellsClipped)/mesh_.globalData().nTotalCells()
                << "%) of cells" << endl;
        }
        if (outputDiagnostics_ && U.mesh().time().outputTime())
        {
            volVectorField Uclip
            (
                IOobject
                (
                    "Uclip",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                U,
                calculatedFvPatchField<vector>::typeName
            );

            Uclip.primitiveFieldRef() = Uctmp;

            Uclip.write();
        }
    }
}


void Foam::fv::velocityLimiterSource::writeData(Ostream& os) const
{
    os  << indent << name_ << endl;
    dict_.write(os);
}


bool Foam::fv::velocityLimiterSource::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        Ulimit_ = readScalar(coeffs_.lookup("Ulimit"));

        if (coeffs_.found("UNames"))
        {
            coeffs_.lookup("UNames") >> fieldNames_;
        }
        else if (coeffs_.found("UName"))
        {
            word UName(coeffs_.lookup("UName"));
            fieldNames_ = wordList(1, UName);
        }
        else
        {
            fieldNames_ = wordList(1, "U");
        }

        applied_.setSize(fieldNames_.size(), false);

        velocityClip_ = coeffs_.lookupOrDefault<Switch>("velocityClip", true);
        velocityDamp_ = coeffs_.lookupOrDefault<Switch>("velocityDamp", false);
        verbose_ = coeffs_.lookupOrDefault<Switch>("verbose", false);
        outputDiagnostics_
            = coeffs_.lookupOrDefault<Switch>("outputDiagnostics", false);

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
