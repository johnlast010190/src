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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2010-2017, 2020 Esi Ltd

\*---------------------------------------------------------------------------*/

#include "radiationModels/fvDOM/fvDOM/fvDOM.H"
#include "submodels/absorptionEmissionModel/absorptionEmissionModel/absorptionEmissionModel.H"
#include "submodels/scatterModel/scatterModel/scatterModel.H"
#include "global/constants/constants.H"
#include "finiteVolume/fvm/fvm.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "submodels/boundaryRadiationProperties/boundaryRadiationProperties.H"

using namespace Foam::constant;
using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(fvDOM, 0);
        addToRadiationRunTimeSelectionTables(fvDOM);
    }

    template<>
    const char* NamedEnum<radiation::fvDOM::solarLoadType, 3>::names[] =
    {
        "viewFactor",
        "DOM",
        "none"
    };

    const NamedEnum<radiation::fvDOM::solarLoadType, 3>
        radiation::fvDOM::solarLoadTypeNames_(0);

}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::radiation::fvDOM::initialise()
{
    // 3D
    if (mesh_.nSolutionD() == 3)
    {
        if (!coeffs_.found("initialMaxIter"))
        {
            if (participating())
            {
                label nCells = mesh_.nCells();
                reduce(nCells, sumOp<label>());
                initMaxIter_ = int(pow(nCells, 1.0/3.0));
            }
            else
            {
                initMaxIter_ = max(initMaxIter_, Pstream::nProcs());
            }
        }

        if (nRayInput_)
        {
            // nTheta and nPhi(Theta)  calculated based on user specified
            // value for nRay_ - See 6.12 and 6.13 FDS Technical Ref.

            nTheta_ = 2.0*round(0.5*1.17*pow(nRay_,(1.0/2.26)));

            labelList nPhiOfTheta(nTheta_);

            label nRayUser = nRay_;
            nRay_ = 0; // reset and compute actual value

            for (label i = 1; i<= nTheta_; i++)
            {
                scalar lowerTheta = pi*(i-1)/nTheta_;
                scalar upperTheta = pi*i/nTheta_;

                nPhi_ =
                    4.0*round
                    (
                        0.25*max
                        (
                            4.0,
                            round(0.5*nRayUser*(Foam::cos(lowerTheta)-Foam::cos(upperTheta)))
                        )
                    );

                nPhiOfTheta[i-1] = nPhi_;

                nRay_ +=nPhi_;
            }

            IRay_.setSize(nRay_);
            scalar deltaTheta = pi/nTheta_;
            label i = 0;
            for (label n = 1; n <= nTheta_; n++)
            {
                scalar thetai = (2.0*n - 1.0)*deltaTheta/2.0;

                for (label m = 1; m <= nPhiOfTheta[n-1]; m++)
                {
                    scalar deltaPhi = 2*pi/nPhiOfTheta[n-1];
                    scalar phii = (2.0*m - 1.0)*deltaPhi/2.0;

                    IRay_.set
                    (
                        i,
                        new radiativeIntensityRay
                        (
                            *this,
                            mesh_,
                            phii,
                            thetai,
                            deltaPhi,
                            deltaTheta,
                            nLambda_,
                            absorptionEmission_,
                            blackBody_,
                            i
                       )
                   );
                   i++;

                }
            }

        }
        else
        {
            nRay_ = 4*nPhi_*nTheta_;
            IRay_.setSize(nRay_);
            scalar deltaPhi = pi/(2.0*nPhi_);
            scalar deltaTheta = pi/nTheta_;
            label i = 0;
            for (label n = 1; n <= nTheta_; n++)
            {
                for (label m = 1; m <= 4*nPhi_; m++)
                {
                    scalar thetai = (2.0*n - 1.0)*deltaTheta/2.0;
                    scalar phii = (2.0*m - 1.0)*deltaPhi/2.0;
                    IRay_.set
                    (
                        i,
                        new radiativeIntensityRay
                        (
                            *this,
                            mesh_,
                            phii,
                            thetai,
                            deltaPhi,
                            deltaTheta,
                            nLambda_,
                            absorptionEmission_,
                            blackBody_,
                            i
                       )
                   );
                   i++;
                }
            }
        }
    }
    // 2D
    else if (mesh_.nSolutionD() == 2)
    {
        scalar thetai = piByTwo;
        scalar deltaTheta = pi;
        if (!coeffs_.found("initialMaxIter"))
        {
            label nCells = mesh_.nCells();
            reduce(nCells, sumOp<label>());
            initMaxIter_ = int(Foam::sqrt(scalar(nCells)));
        }

        if (nRayInput_)
        {
            label nRayUser = nRay_;

            // Ref. FDS code
            nTheta_ = 1;
            nPhi_ = 4.0*round(0.25*nRayUser);

            nRay_ = nTheta_*nPhi_;
            IRay_.setSize(nRay_);

            scalar deltaPhi = 2.0*pi/nPhi_;
            label i = 0;
            for (label m = 1; m <= nPhi_; m++)
            {
                scalar phii = (2.0*m - 1.0)*deltaPhi/2.0;
                IRay_.set
                (
                    i,
                    new radiativeIntensityRay
                    (
                        *this,
                        mesh_,
                        phii,
                        thetai,
                        deltaPhi,
                        deltaTheta,
                        nLambda_,
                        absorptionEmission_,
                        blackBody_,
                        i
                    )
                );
                i++;
            }
        }
        else
        {
            nRay_ = 4*nPhi_;
            IRay_.setSize(nRay_);
            scalar deltaPhi = pi/(2.0*nPhi_);
            label i = 0;
            for (label m = 1; m <= 4*nPhi_; m++)
            {
                scalar phii = (2.0*m - 1.0)*deltaPhi/2.0;
                IRay_.set
                (
                    i,
                    new radiativeIntensityRay
                    (
                        *this,
                        mesh_,
                        phii,
                        thetai,
                        deltaPhi,
                        deltaTheta,
                        nLambda_,
                        absorptionEmission_,
                        blackBody_,
                        i
                    )
                );
                i++;
            }
        }
    }
    // 1D
    else
    {
        scalar thetai = piByTwo;
        scalar deltaTheta = pi;
        if (!coeffs_.found("initialMaxIter"))
        {
            label nCells = mesh_.nCells();
            reduce(nCells, sumOp<label>());
            initMaxIter_ = nCells;
        }

        if (nRayInput_)
        {
            NotImplemented;
        }
        else
        {
            nRay_ = 2;
            IRay_.setSize(nRay_);
            scalar deltaPhi = pi;
            label i = 0;
            for (label m = 1; m <= 2; m++)
            {
                scalar phii = (2.0*m - 1.0)*deltaPhi/2.0;
                IRay_.set
                (
                    i,
                    new radiativeIntensityRay
                    (
                        *this,
                        mesh_,
                        phii,
                        thetai,
                        deltaPhi,
                        deltaTheta,
                        nLambda_,
                        absorptionEmission_,
                        blackBody_,
                        i
                     )
                );
                i++;
            }
        }
    }


    // Construct absorption field for each wavelength
    forAll(aLambda_, lambdaI)
    {
        aLambda_.set
        (
            lambdaI,
            new volScalarField
            (
                IOobject
                (
                    "aLambda_" + Foam::name(lambdaI) ,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                a_
            )
        );
    }

    Info<< "fvDOM : Allocated " << IRay_.size()
        << " rays with average orientation:" << nl;

    forAll(IRay_, rayId)
    {
        if (omegaMax_ <  IRay_[rayId].omega())
        {
            omegaMax_ = IRay_[rayId].omega();
        }
        Info<< '\t' << IRay_[rayId].rayName() << " : " << "dAve : "
            << '\t' << IRay_[rayId].dAve() << nl;
    }

    Info<< endl;

    solarLoadMode_ = solarLoadTypeNames_
        [this->lookupOrDefault<word>("solarLoadMode", "none")];

    switch (solarLoadMode_)
    {
        case sltViewFactor :
        {
            const dictionary& solarDict = this->subDict("solarLoadCoeffs");
            vfSolarLoad_.reset
            (
                new solarLoad(solarDict, T_, externalRadHeatFieldName_)
            );

            if (vfSolarLoad_->nBands() > 1)
            {
                FatalErrorInFunction
                    << "Requested solar radiation with fvDOM. Using "
                    << "more than one band for the solar load is not allowed"
                    << exit(FatalError);
            }

            Info<< "Creating Viewfactor Solar Load Model " << nl;

        } break;

        case sltDOM :
        {
            const dictionary& solarDict = this->subDict("domSolarCoeffs");

            domSolarLoad_.reset
            (
                new domSolar(solarDict, T_, externalRadHeatFieldName_)
            );

            Info<< "Creating DOM Solar Load Model " << nl;
        } break;

        case sltnone :
        {
        } break;

        default:
        {
            FatalErrorInFunction
                << "Unsupported solarLoadMode "
                << solarLoadMode_
                << exit(FatalError);
        }
    }

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::fvDOM::fvDOM(const volScalarField& T)
:
    radiationModel(typeName, T),
    G_
    (
        IOobject
        (
            "G",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("G", dimMass/pow3(dimTime), 0.0)
    ),
    qr_
    (
        IOobject
        (
            "Qr",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("Qr", dimMass/pow3(dimTime), 0.0)
    ),
    qin_
    (
        IOobject
        (
            "qin",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("qin", dimMass/pow3(dimTime), 0.0)
    ),
    /*
    qem_
    (
        IOobject
        (
            "qem",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("qem", dimMass/pow3(dimTime), 0.0)
    ),
    */
    qg_
    (
        IOobject
        (
            "QrG",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("QrG", dimMass/pow3(dimTime), 0.0)
    ),
    a_
    (
        IOobject
        (
            "a",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
            // IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("a", dimless/dimLength, 0.0)
    ),
    // nTheta_(readLabel(coeffs_.lookup("nTheta"))),
    // nPhi_(readLabel(coeffs_.lookup("nPhi"))),
    nTheta_(coeffs_.lookupOrDefault<label>("nTheta", 4)),
    nPhi_(coeffs_.lookupOrDefault<label>("nPhi", 2)),
    nRayInput_(coeffs_.lookupOrDefault<Switch>("nRayInput", false)),
    nRay_(coeffs_.lookupOrDefault<label>("nRay", 0)),
    nLambda_(absorptionEmission_->nBands()),
    aLambda_(nLambda_),
    blackBody_(nLambda_, T),
    IRay_(0),
    convergence_(coeffs_.lookupOrDefault<scalar>("convergence", 0.0)),
    maxIter_(coeffs_.lookupOrDefault<label>("maxIter", 50)),
    initConvergence_
    (
        coeffs_.lookupOrDefault<scalar>("initialConvergence", convergence_)
    ),
    initMaxIter_(coeffs_.lookupOrDefault<label>("initialMaxIter", maxIter_)),
    omegaMax_(0),
    solarLoadMode_
    (
        solarLoadTypeNames_
        [this->lookupOrDefault<word>("solarLoadMode", "none")]
    ),
    vfSolarLoad_(),
    domSolarLoad_(),
    meshOrientation_
    (
        coeffs_.lookupOrDefault<vector>("meshOrientation", Zero)
    )
{
    initialise();
}
Foam::radiation::fvDOM::fvDOM
(
    const dictionary& dict,
    const volScalarField& T
)
:
    radiationModel(typeName, dict, T),
    G_
    (
        IOobject
        (
            "G",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("G", dimMass/pow3(dimTime), 0.0)
    ),
    qr_
    (
        IOobject
        (
            "Qr",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("Qr", dimMass/pow3(dimTime), 0.0)
    ),
    qin_
    (
        IOobject
        (
            "qin",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("qin", dimMass/pow3(dimTime), 0.0)
    ),
    /*
    qem_
    (
        IOobject
        (
            "qem",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("qem", dimMass/pow3(dimTime), 0.0)
    ),
    */
    qg_
    (
        IOobject
        (
            "QrG",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("QrG", dimMass/pow3(dimTime), 0.0)
    ),
    a_
    (
        IOobject
        (
            "a",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("a", dimless/dimLength, 0.0)
    ),
    // nTheta_(readLabel(coeffs_.lookup("nTheta"))),
    // nPhi_(readLabel(coeffs_.lookup("nPhi"))),
    nTheta_(coeffs_.lookupOrDefault<label>("nTheta", 4)),
    nPhi_(coeffs_.lookupOrDefault<label>("nPhi", 2)),
    nRayInput_(coeffs_.lookupOrDefault<Switch>("nRayInput", false)),
    nRay_(coeffs_.lookupOrDefault<label>("nRay", 0)),
    nLambda_(absorptionEmission_->nBands()),
    aLambda_(nLambda_),
    blackBody_(nLambda_, T),
    IRay_(0),
    convergence_(coeffs_.lookupOrDefault<scalar>("convergence", 0.0)),
    maxIter_(coeffs_.lookupOrDefault<label>("maxIter", 50)),
    initConvergence_
    (
        coeffs_.lookupOrDefault<scalar>("initialConvergence", convergence_)
    ),
    initMaxIter_(coeffs_.lookupOrDefault<label>("initialMaxIter", maxIter_)),
    omegaMax_(0),
    solarLoadMode_
    (
        solarLoadTypeNames_
        [this->lookupOrDefault<word>("solarLoadMode", "none")]
    ),
    vfSolarLoad_(),
    domSolarLoad_(),
    meshOrientation_
    (
        coeffs_.lookupOrDefault<vector>("meshOrientation", Zero)
    )
{
    initialise();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::fvDOM::~fvDOM()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::radiation::fvDOM::read()
{
    if (radiationModel::read())
    {
        // Only reading solution parameters - not changing ray geometry
        coeffs_.readIfPresent("convergence", convergence_);
        coeffs_.readIfPresent("maxIter", maxIter_);

        return true;
    }
    else
    {
        return false;
    }
}


void Foam::radiation::fvDOM::calculate()
{
    absorptionEmission_->correct(a_, aLambda_);

    updateBlackBodyEmission();

    switch (solarLoadMode_)
    {
        case sltViewFactor :
        {
            vfSolarLoad_->calculate();
        } break;

        case sltDOM :
        {
            domSolarLoad_->calculate();
        } break;

        case sltnone :
        {
        } break;

        default:
        {
            FatalErrorInFunction
                << "Unsupported solarLoadMode "
                << solarLoadMode_
                << exit(FatalError);
        }
    }

    // Set rays convergence false
    List<bool> rayIdConv(nRay_, false);

    scalar maxResidual = 0.0;
    label radIter = 0;

    scalar conv = convergence_;
    scalar maxi = maxIter_;

    if (time_.timeIndex() == 0)
    {
        conv = initConvergence_;
        maxi = initMaxIter_;
        Info<< nl <<  "Initial fvDOM solution " << maxi
             << " iterations." << endl;
    }

    do
    {
        Info<< "Radiation solver iter: " << radIter << endl;

        radIter++;
        maxResidual = 0.0;
        forAll(IRay_, rayI)
        {
            if (!rayIdConv[rayI])
            {
                scalar maxBandResidual = IRay_[rayI].correct();
                maxResidual = max(maxBandResidual, maxResidual);

                reduce(maxBandResidual, maxOp<scalar>());

                if (maxBandResidual < convergence_)
                {
                    rayIdConv[rayI] = true;
                }
            }
        }

        reduce(maxResidual, maxOp<scalar>());
        Info<< "Global Max initial residual: " <<maxResidual << endl;

    } while (maxResidual > conv && radIter < maxi);

    updateG();
}


Foam::tmp<Foam::volScalarField> Foam::radiation::fvDOM::Rp() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "Rp",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            4.0*a_*physicoChemical::sigma //absorptionEmission_->a()
        )
    );
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::radiation::fvDOM::Ru() const
{

    const volScalarField::Internal& G = G_();
    const volScalarField::Internal E = absorptionEmission_->ECont()()();
    const volScalarField::Internal a = a_.internalField();

    return a*G - E;
}


void Foam::radiation::fvDOM::updateBlackBodyEmission()
{
    for (label j=0; j < nLambda_; j++)
    {
        blackBody_.correct(j, absorptionEmission_->bands(j));
    }
}


void Foam::radiation::fvDOM::updateG()
{
    G_ = dimensionedScalar("zero",dimMass/pow3(dimTime), 0.0);
    qr_ = dimensionedScalar("zero",dimMass/pow3(dimTime), 0.0);
    //qem_ = dimensionedScalar("zero", dimMass/pow3(dimTime), 0.0);
    qg_ = dimensionedScalar("zero", dimMass/pow3(dimTime), 0.0);
    qin_ = dimensionedScalar("zero", dimMass/pow3(dimTime), 0.0);

    forAll(IRay_, rayI)
    {
        G_ += IRay_[rayI].I()*IRay_[rayI].omega();
        qr_.boundaryFieldRef() += IRay_[rayI].qr();
        //qem_.boundaryFieldRef() += IRay_[rayI].qem().boundaryField();
        //qem_.boundaryFieldRef() += IRay_[rayI].qem();
        qg_.boundaryFieldRef() += IRay_[rayI].qg();
        qin_.boundaryFieldRef() += IRay_[rayI].qin();
    }

    if
    (
        mesh_.foundObject<volScalarField>(externalRadHeatFieldName_)
     && solarLoadMode_ != sltnone
    )
    {
        const volScalarField& Qext =
            mesh_.lookupObject<volScalarField>(externalRadHeatFieldName_);

        forAll(mesh_.boundaryMesh(), patchI)
        {
            if (!boundaryProperties().radBoundaryProperties()[patchI].empty())
            {
                qin_.boundaryFieldRef()[patchI] -=
                    boundaryProperties().emissivity(patchI)*
                        Qext.boundaryField()[patchI];

                qr_.boundaryFieldRef()[patchI] +=
                    (1-boundaryProperties().transmissivity(patchI))*
                        Qext.boundaryField()[patchI];

                qg_.boundaryFieldRef()[patchI] -= Qext.boundaryField()[patchI];
            }
        }
    }
}


void Foam::radiation::fvDOM::setRayIdLambdaId
(
    const word& name,
    label& rayId,
    label& lambdaId
) const
{
    // assuming name is in the form: CHARS_rayId_lambdaId
    size_type i1 = name.find_first_of("_");
    size_type i2 = name.find_last_of("_");

    rayId = readLabel(IStringStream(name.substr(i1+1, i2-1))());
    lambdaId = readLabel(IStringStream(name.substr(i2+1, name.size()-1))());
}

Foam::tmp<Foam::scalarField> Foam::radiation::fvDOM::emittedRadiantIntensity
(
    label patchI,
    const scalarField& Tp
) const
{
    tmp<scalarField> efPtr(new scalarField(Tp.size(), 0.0));
    scalarField& ef = efPtr.ref();

    forAll(IRay_, rayI)
    {
        ef += IRay_[rayI].emittedRadiantIntensity(patchI, Tp);
    }

    return efPtr;
}

const Foam::radiation::domSolar& Foam::radiation::fvDOM::getDomSolarObj() const
{
    return domSolarLoad_;
}

// ************************************************************************* //
