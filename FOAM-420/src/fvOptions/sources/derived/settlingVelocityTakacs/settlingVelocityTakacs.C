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
    (c) 2020 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "settlingVelocityTakacs.H"
#include "fvMatrices/fvMatrices.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "interpolation/surfaceInterpolation/surfaceInterpolation/surfaceInterpolate.H"
#include "finiteVolume/fvc/fvcReconstruct.H"
#include "interpolation/surfaceInterpolation/limitedSchemes/upwind/upwind.H"
#include "interpolation/surfaceInterpolation/schemes/downwind/downwind.H"
#include "fields/volFields/volFields.H"
#include "fields/UniformDimensionedFields/uniformDimensionedFields.H"
#include "finiteVolume/fvc/fvcFlux.H"
#include "finiteVolume/fvm/fvmDiv.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(settlingVelocityTakacs, 0);

    addToRunTimeSelectionTable
    (
        option,
        settlingVelocityTakacs,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::settlingVelocityTakacs::computeFlux(const fvMesh& mesh)
{
    if (constant_)
    {
        *UPtr_ = V0_;
        limitCellSet();
        UPtr_->correctBoundaryConditions();
        *phiPtr_ = mesh.Sf() & fvc::interpolate(*UPtr_);
    }
    else
    {
        // get access to concentration tracer
        const volScalarField& X = mesh.lookupObject<volScalarField>(XName_);

        // compute blending
        tmp<volScalarField> blending = mag(X)*0.0;
        forAll(X, fI)
        {
            if (X[fI] > X34_.value())
            {
                blending.ref()[fI] = 1.0;
            }
        }
        forAll(X.boundaryField(), patchI)
        {
            forAll(X.boundaryField()[patchI], fI)
            {
                if (X.boundaryField()[patchI][fI] > X34_.value())
                {
                    blending.ref().boundaryFieldRef()[patchI][fI] = 1.0;
                }
            }
        }

        // make sure exponent doesnt have dimensions but value is in mg/L
        dimensionedScalar fac("fac", Xns_.dimensions()/X.dimensions(), fac_);

        // computing settling velocity U
        *UPtr_ = blending() * V0t_ * exp(-rt_*(X*fac-Xns_))
            + (1.-blending()) * V0_
            * max((exp(-rh_*(X*fac-Xns_)) - exp(-rp_*(X*fac-Xns_))), scalar(0));
        blending.clear();

        // limit to selected cells
        limitCellSet();

        // compute settling flux direction for interpolation scheme
        tmp<surfaceScalarField> fluxDir;
        fluxDir = V0_ & mesh.Sf();

        // compute face-interpolated velocity (=flux)
        tmp<surfaceVectorField> Uf =
            downwind<vector>(mesh, fluxDir()).interpolate
            (
                UPtr_
            );
        tmp<surfaceVectorField> Uupw =
            upwind<vector>(mesh, fluxDir()).interpolate
            (
                UPtr_
            );
        fluxDir.clear();
        tmp<surfaceVectorField> Ulin = linearInterpolate(*UPtr_);
        Uf.ref() *= blend_;
        Uf.ref() += (1.-blend_) * Ulin();

        // blend downwind-upwind interpolation based on Xf
        scalar X23
        (
            (Foam::log(rp_.value()) - Foam::log(rh_.value()))
          / (rp_.value() - rh_.value())
          + Xns_.value()
        );

        tmp<surfaceScalarField> Xf = fvc::interpolate(X);
        forAll(Xf(), fI)
        {
            if (Xf()[fI] < X23)
            {
                Uf.ref()[fI] = blend_*Uupw()[fI] + (1.-blend_)*Ulin()[fI];
            }
        }
        forAll(Xf().boundaryField(), patchI)
        {
            forAll(Xf().boundaryField()[patchI], fI)
            {
                if (X.boundaryField()[patchI][fI] < X23)
                {
                    Uf.ref().boundaryFieldRef()[patchI][fI] =
                        blend_*Uupw().boundaryField()[patchI][fI]
                      + (1.-blend_)*Ulin().boundaryField()[patchI][fI];
                }
            }
        }
        Xf.clear();
        Uupw.clear();
        Ulin.clear();

        Info<< "Computing solid particle flux " << phiName_<< endl;
        *phiPtr_ = mesh.Sf() & Uf;
    }

    if (phaseName_ != "none")
    {
        if (mesh.foundObject<volScalarField>(phaseName_))
        {
            const volScalarField& alpha = mesh.lookupObject<volScalarField>(phaseName_);
            volScalarField alphaLim(alpha);

            if (limitBulk_)
            {
                forAll(alpha, cI)
                {
                    // limit based on alpha threshold
                    if (alpha[cI] < threshold_)
                    {
                        const labelList& nbr = mesh.cellCells()[cI];
                        bool limit(true);
                        forAll(nbr, idx)
                        {
                            if (alpha[nbr[idx]] > threshold_ && nbr[idx] != cI)
                            {
                                limit = false;
                            }
                        }

                        if (limit)
                        {
                            alphaLim[cI] = 0;
                        }
                    }

                    // limit based on position
                    vector unitDir = -V0_.value()/(mag(V0_).value()+VSMALL);
                    if ((mesh.C()[cI] & unitDir) < limitHeight_)
                    {
                        alphaLim[cI] = 0;
                    }
                }
                alphaLim.correctBoundaryConditions();
            }

            *phiPtr_ *= fvc::interpolate(alphaLim);
            *UPtr_ *= alphaLim;
        }
        else
        {
            WarningInFunction
                << "Phase " << phaseName_ << " specified but not found."
                << nl << "Continue without flux scaling." << nl << endl;
        }
    }
}


void Foam::fv::settlingVelocityTakacs::limitCellSet()
{
    if (debug)
    {
        Info<<"limit drift fluxes to " << cells_.size() << " cells" << endl;
    }
    volVectorField Ulim(*UPtr_);
    Ulim *= 0.0;
    forAll(cells_, i)
    {
        label celli = cells_[i];
        Ulim[celli] = (*UPtr_)[celli];
    }
    *UPtr_ = Ulim;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::settlingVelocityTakacs::settlingVelocityTakacs
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const objectRegistry& obr
)
:
    cellSetOption(sourceName, modelType, dict, obr),
    XName_(),
    phiName_(dict.lookupOrDefault<word>("phiName", "phiSolid")),
    UName_(dict.lookupOrDefault<word>("UName", "USolid")),
    phaseName_(dict.lookupOrDefault<word>("phaseName", "none")),
    V0_(
        dimensionedVector
        (
            "V0",
            dimVelocity,
            dict.lookup("V0")
        )
    ),
    V0t_(dimensionedVector("V0t", dimVelocity, dict.lookup("V0t"))),
    rh_(dimensionedScalar("rh", dimless/dimDensity, readScalar(dict.lookup("rh")))),
    rp_(dimensionedScalar("rp", dimless/dimDensity, readScalar(dict.lookup("rp")))),
    rt_(dimensionedScalar("rt", dimless/dimDensity, readScalar(dict.lookup("rt")))),
    Xns_(dimensionedScalar("Xns", dimDensity, readScalar(dict.lookup("Xns")))),
    X34_(dimensionedScalar("X34", dimDensity, readScalar(dict.lookup("X34")))),
    fac_(),
    blend_(dict.lookupOrDefault<scalar>("blend", 0.99)),
    constant_(dict.lookupOrDefault<Switch>("constant", false)),
    limitBulk_(dict.lookupOrDefault<Switch>("limitBulk", false)),
    threshold_(dict.lookupOrDefault<scalar>("alphaThreshold", 0.2)),
    limitHeight_(dict.lookupOrDefault<scalar>("limitHeight", -GREAT)),
    phiPtr_(nullptr),
    UPtr_(nullptr)
{
    // from options: add field name to call source for
    if (coeffs_.found("fieldNames"))
    {
        coeffs_.lookup("fieldNames") >> fieldNames_;

        if (fieldNames_.size() != 1)
        {
            FatalErrorInFunction
                << "settings are:" << fieldNames_ << exit(FatalError);
        }
    }
    else if (coeffs_.found("fieldName"))
    {
        word fieldName(coeffs_.lookup("fieldName"));
        fieldNames_ = wordList(1, fieldName);
    }

    applied_.setSize(fieldNames_.size(), false);
}


bool Foam::fv::settlingVelocityTakacs::initialise()
{
    UPtr_ = new volVectorField
    (
        IOobject
        (
            UName_,
            mesh().time().timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh()
    );

    if (mesh().foundObject<IOdictionary>("transportProperties"))
    {
        const IOdictionary& transpDict = obr_.lookupObject<IOdictionary>("transportProperties");
        if (transpDict.found("BokilBewtraCoeffs"))
        {
            XName_ = word(transpDict.subDict("BokilBewtraCoeffs").lookup("XName"));
            fac_ = transpDict.subDict("BokilBewtraCoeffs").lookupOrDefault<scalar>("Xscaling", 1.0);
        }
        else
        {
            XName_ = coeffs().lookupOrDefault<word>("XName", fieldNames_[0]);
            fac_ = coeffs().lookupOrDefault<scalar>("Xscaling", 1.0);
        }
    }
    else
    {
        XName_ = coeffs().lookupOrDefault<word>("XName", fieldNames_[0]);
        fac_ = coeffs().lookupOrDefault<scalar>("Xscaling", 1.0);
    }

    phiPtr_ = new surfaceScalarField
    (
        IOobject
        (
            phiName_,
            mesh().time().timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::flux(*UPtr_)
    );

    return true;
}


void Foam::fv::settlingVelocityTakacs::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    computeFlux(eqn.psi().mesh());
    eqn -= fvm::div(*phiPtr_,eqn.psi(),"div(" +  phiName_+ "," + Foam::word(eqn.psi().name()) + ")");
}


void Foam::fv::settlingVelocityTakacs::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    computeFlux(eqn.psi().mesh());
    eqn -= fvm::div(*phiPtr_*fvc::interpolate(rho),eqn.psi(),"div(" +  phiName_+ "," + Foam::word(eqn.psi().name()) + ")");
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::settlingVelocityTakacs::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
