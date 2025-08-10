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
    (c) 2018 Esi Ltd.
    (c) 2011-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "rotatingVelocityForce.H"
#include "fvMatrices/fvMatrices.H"
#include "fields/DimensionedFields/DimensionedField/DimensionedField.H"
#include "db/IOstreams/Fstreams/IFstream.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(rotatingVelocityForce, 0);

    addToRunTimeSelectionTable
    (
        option,
        rotatingVelocityForce,
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fv::rotatingVelocityForce::writeProps
(
    const scalar Comega
) const
{
    // Only write on output time
    if (mesh_.time().writeTime())
    {
        IOdictionary propsDict
        (
            IOobject
            (
                name_ + "Properties",
                mesh_.time().timeName(),
                "uniform",
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            )
        );
        propsDict.add("Comega", Comega);
        propsDict.regIOobject::write();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::rotatingVelocityForce::rotatingVelocityForce
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const objectRegistry& obr
)
:
    cellSetOption(sourceName, modelType, dict, obr),
    omega_(Function1<scalar>::New("omega", dict)),
    axis_(dict.lookup("axis")),
    origin_(dict.lookup("origin")),
    Comega0_(0.0),
    dComega_(0.0),
    relaxation_(coeffs_.lookupOrDefault<scalar>("relaxation", 1.0)),
    rAPtr_(nullptr)
{
    coeffs_.lookup("fields") >> fieldNames_;

    if (fieldNames_.size() != 1)
    {
        FatalErrorInFunction
            << "settings are:" << fieldNames_ << exit(FatalError);
    }

    applied_.setSize(fieldNames_.size(), false);

    // Read the initial pressure gradient from file if it exists
    IFstream propsFile
    (
        mesh_.time().timePath()/"uniform"/(name_ + "Properties")
    );

    if (propsFile.good())
    {
        Info<< "    Reading force coefficient from file" << endl;
        dictionary propsDict(dictionary::null, propsFile);
        propsDict.lookup("Comega") >> Comega0_;
    }

    Info<< "    Initial force coefficient = " << Comega0_ << nl << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::fv::rotatingVelocityForce::omegaAve
(
    const volVectorField& U
) const
{
    scalar omegaAve = 0.0;

    const scalarField& cv = mesh_.V();
    const vectorField& cc = mesh_.C();

    //omega * (axis ^ r ) = Ut : omega unknown
    //omega * (axis ^ r ).(axis ^ r )/magSqr(axis ^ r ) = Ut.(axis ^ r )/magSqr(axis ^ r )
    //omega = Ut.(axis ^ r )/magSqr(axis ^ r )

    forAll(cells_, i)
    {
        label celli = cells_[i];
        scalar volCell = cv[celli];
        vector R = cc[celli] - origin_;
        R -= axis_*(R & axis_);

        vector axR = (axis_ ^ R);

        omegaAve += (U[celli] & axR) / magSqr(axR)*volCell;
    }

    reduce(omegaAve, sumOp<scalar>());

    omegaAve /= V_;

    return omegaAve;
}

Foam::tmp<Foam::vectorField>
Foam::fv::rotatingVelocityForce::targetVelocity() const
{
    tmp<vectorField> Utmp(new vectorField(cells_.size(), Zero));

    const vectorField& cc = mesh_.C();

    scalar omega = omega_->value(mesh().time().value());

    forAll(cells_, i)
    {
        label celli = cells_[i];
        vector R = cc[celli] - origin_;
        R -= axis_*(R & axis_);

        Utmp->operator[](i) = omega*(axis_ ^ R);
    }

    return Utmp;
}


void Foam::fv::rotatingVelocityForce::correct(volVectorField& U)
{
    const scalarField& rAU = rAPtr_();

    // Integrate flow variables over cell set
    scalar rAUave = 0.0;
    const scalarField& cv = mesh_.V();
    forAll(cells_, i)
    {
        label celli = cells_[i];
        scalar volCell = cv[celli];
        rAUave += rAU[celli]*volCell;
    }

    // Collect across all processors
    reduce(rAUave, sumOp<scalar>());

    // Volume averages
    rAUave /= V_;

    scalar omegaAve = this->omegaAve(U);

    // Calculate the force coefficient increment needed to adjust the average
    // flow-rate to the desired value
    dComega_ = relaxation_*(omega_->value(mesh().time().value()) - omegaAve)/rAUave;

    tmp<vectorField> Utarget(targetVelocity());

    // Apply correction to velocity field
    forAll(cells_, i)
    {
        label celli = cells_[i];
        U[celli] += rAU[celli]*dComega_*Utarget()[i];
    }

    scalar Comega = Comega0_ + dComega_;

    Info<< "Swirl source: uncorrected mean rotation, omegaAve = " << omegaAve
        << ", force coefficient = " << Comega << endl;

    writeProps(Comega);
}


void Foam::fv::rotatingVelocityForce::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    volVectorField::Internal Su
    (
        IOobject
        (
            name_ + fieldNames_[fieldi] + "Sup",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("zero", eqn.dimensions()/dimVolume, Zero)
    );

    scalar Comega = Comega0_ + dComega_;

    UIndirectList<vector>(Su, cells_) = Comega*targetVelocity();

    eqn += Su;
}


void Foam::fv::rotatingVelocityForce::addSup
(
    fvBlockMatrix<vector>& eqn,
    const label fieldi
)
{
    volVectorField::Internal Su
    (
        IOobject
        (
            name_ + fieldNames_[fieldi] + "Sup",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("zero", eqn.dimensionSets()[0]/dimVolume, Zero)
    );

    scalar Comega = Comega0_ + dComega_;

    UIndirectList<vector>(Su, cells_) = Comega*targetVelocity();

    eqn += Su;
}


void Foam::fv::rotatingVelocityForce::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    this->addSup(eqn, fieldi);
}


void Foam::fv::rotatingVelocityForce::addSup
(
    const volScalarField& rho,
    fvBlockMatrix<vector>& eqn,
    const label fieldi
)
{
    this->addSup(eqn, fieldi);
}


void Foam::fv::rotatingVelocityForce::constrain
(
    fvMatrix<vector>& eqn,
    const label
)
{
    if (rAPtr_.empty())
    {
        rAPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    name_ + ":rA",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                1.0/eqn.A()
            )
        );
    }
    else
    {
        rAPtr_() = 1.0/eqn.A();
    }

    Comega0_ += dComega_;
    dComega_ = 0.0;
}


void Foam::fv::rotatingVelocityForce::constrain
(
    fvBlockMatrix<vector>& eqn,
    const label
)
{
    const Field<tensor>& d = eqn.diag().asSquare();
    if (rAPtr_.empty())
    {
        rAPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    name_ + ":rA",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                mesh_,
                dimensionedScalar
                (
                    "zero",
                    eqn.dimensionSets()[0],
                    0.0
                ),
                zeroGradientFvPatchField<tensor>::typeName
            )
        );
    }

    forAll(d, cI)
    {
        rAPtr_()[cI] = 3*mesh_.V()[cI]/
            ((d[cI].component(0)+d[cI].component(4)+d[cI].component(8)));
    }

    Comega0_ += dComega_;
    dComega_ = 0.0;
}


// ************************************************************************* //
