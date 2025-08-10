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
    (c) 2021 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "scalarFieldLimiterSource.H"
#include "fvMesh/fvMesh.H"
#include "fvMatrices/fvMatrices.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(scalarFieldLimiterSource, 0);
    addToRunTimeSelectionTable
    (
        option,
        scalarFieldLimiterSource,
        dictionary
    );
}
}

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::fv::scalarFieldLimiterSource::addDamping(fvMatrix<scalar>& eqn)
{
    //copied from OF+ scalarFieldDampingConstraint
    const scalarField& vol = mesh_.V();
    const volScalarField& vf = eqn.psi();
    scalarField& diag = eqn.diag();

    label nDamped = 0;

    scalarField damping(vf.size(), 0.0);

    forAll(vf, cellI)
    {
        if (vf[cellI] > limit_)
        {
            scalar scale = sqr(Foam::cbrt(vol[cellI]));

            damping[cellI] = scale*(vf[cellI]-limit_);
            diag[cellI] += damping[cellI];

            nDamped++;
        }
    }

    reduce(nDamped, sumOp<label>());

    if (verbose_)
    {
        Info<< type() << " " << name_ << " damped "
            << nDamped << " ("
            << 100*scalar(nDamped)/mesh_.globalData().nTotalCells()
            << "%) of cells" << endl;
    }
    if (outputDiagnostics_ && vf.mesh().time().outputTime())
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

Foam::fv::scalarFieldLimiterSource::scalarFieldLimiterSource
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

void Foam::fv::scalarFieldLimiterSource::constrain
(
    fvMatrix<scalar>& eqn,
    const label fieldI
)
{
    if (scalarFieldDamp_)
    {
        addDamping(eqn);
    }
}


void Foam::fv::scalarFieldLimiterSource::correct(volScalarField& vf)
{
    //clip U magnitude
    if (scalarFieldClip_)
    {
        label nCellsClipped = 0;
        scalarField vfctmp = vf.internalField();

        forAll(vf, i)
        {
            if (vf[i] > limit_)
            {
                vf[i] *= limit_/vf[i];
                nCellsClipped++;
            }
        }
        vfctmp -= vf.internalField();

        reduce(nCellsClipped, sumOp<label>());

        if (verbose_)
        {
            Info<< type() << " " << name_ << " clipped "
                << nCellsClipped << " ("
                << 100*scalar(nCellsClipped)/mesh_.globalData().nTotalCells()
                << "%) of cells" << endl;
        }
        if (outputDiagnostics_ && vf.mesh().time().outputTime())
        {
            volScalarField vfclip
            (
                IOobject
                (
                    "vfclip",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                vf,
                calculatedFvPatchField<scalar>::typeName
            );

            vfclip.primitiveFieldRef() = vfctmp;

            vfclip.write();
        }
    }
}


void Foam::fv::scalarFieldLimiterSource::writeData(Ostream& os) const
{
    os  << indent << name_ << endl;
    dict_.write(os);
}


bool Foam::fv::scalarFieldLimiterSource::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        limit_ = readScalar(coeffs_.lookup("limit"));

        if (coeffs_.found("fieldNames"))
        {
            coeffs_.lookup("fieldNames") >> fieldNames_;
        }
        else if (coeffs_.found("fieldName"))
        {
            word fieldName(coeffs_.lookup("fieldName"));
            fieldNames_ = wordList(1, fieldName);
        }
        else
        {
            //FoamFatalError
        }

        applied_.setSize(fieldNames_.size(), false);

        scalarFieldClip_ = coeffs_.lookupOrDefault<Switch>("scalarFieldClip", true);
        scalarFieldDamp_ = coeffs_.lookupOrDefault<Switch>("scalarFieldDamp", false);
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
