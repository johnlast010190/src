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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2016 OpenCFD Ltd.
    (c) 2010-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "SpalartAllmarasDES.H"
#include "fvMesh/fvPatches/derived/wall/wallFvPatch.H"
#include "fields/fvPatchFields/basic/fixedValue/fixedValueFvPatchFields.H"
#include "finiteVolume/snGradSchemes/uncorrectedSnGrad/uncorrectedSnGrad.H"
#include "cfdTools/general/fvOptions/fvOptions.H"
#include "fluidThermo/fluidThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool SpalartAllmarasDES<BasicTurbulenceModel>::isIsochoric() const
{
    bool isIsoc = true;

    if
    (
        this->mesh_.objectRegistry::template foundObject<fluidThermo>
        (
            basicThermo::dictName
        )
    )
    {
        const fluidThermo* thermoPtr =
            this->db_.objectRegistry::template lookupObjectPtr<fluidThermo>
            (
                basicThermo::dictName
            );
        if (!thermoPtr->isochoric())
        {
            isIsoc = false;
        }
    }

    return isIsoc;
}


template<class BasicTurbulenceModel>
tmp<volTensorField> SpalartAllmarasDES<BasicTurbulenceModel>::gradUnwc() const
{
    if (gradUstd_)
    {
        return(fvc::grad(this->U_));
    }

    //add wall shear velocity to correct wall-adjacent cell centre gradient
    volVectorField Usa
    (
        IOobject
        (
            "Usa",
            this->mesh_.time().timeName(),
            this->db_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedVector("vector", dimVelocity, vector(0,0,0)),
        fixedValueFvPatchVectorField::typeName
    );

    Usa.forceAssign(this->U_);


    const tmp<volScalarField> tmuEff = this->muEff();

    volVectorField::Boundary& Usabf = Usa.boundaryFieldRef();

    forAll(Usa.boundaryField(),pI)
    {
        if (isA<wallFvPatch>(this->mesh_.boundary()[pI]))
        {
            const labelUList& faceCells =
                this->mesh_.boundary()[pI].faceCells();
            tmp<vectorField> n(this->mesh_.boundary()[pI].nf());

            forAll(faceCells,pfI)
            {
                label fcI = faceCells[pfI];

                Usabf[pI][pfI] = this->U_[fcI] -
                    (tmuEff->boundaryField()[pI][pfI]/tmuEff()[fcI])
                    *(this->U_[fcI] - this->U_.boundaryField()[pI][pfI]);

                //set wall normal component equal to cell value
                //for zero wall normal normal stress

                Usabf[pI][pfI] += n()[pfI]*(n()[pfI] & this->U_[fcI])
                    - n()[pfI]*(n()[pfI] & Usa.boundaryField()[pI][pfI]);
            }
        }
    }

    word gradName = "grad(" + this->U_.name() + ')';

    if (this->mesh().solution().cache(gradName))
    {
        return(fvc::grad(Usa));
    }
    else
    {
        return(fvc::grad(Usa, gradName));
    }
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasDES<BasicTurbulenceModel>::Cb1
(
    const volScalarField& SStar
) const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "Cb1",
                this->mesh_.time().timeName(),
                this->db_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh_,
            Cb1_,
            zeroGradientFvPatchScalarField::typeName
        )
    );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasDES<BasicTurbulenceModel>::Cw1
(
    const volScalarField& SStar
) const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "Cw1",
                this->mesh_.time().timeName(),
                this->db_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            Cb1(SStar)/sqr(kappa_) + (1.0 + Cb2_)/sigmaNut_
        )
    );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasDES<BasicTurbulenceModel>::chi() const
{
    return nuTilda_/this->nu();
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasDES<BasicTurbulenceModel>::fv1
(
    const volScalarField& chi
) const
{
    const volScalarField chi3("chi3", pow3(chi));
    return chi3/(chi3 + pow3(Cv1_));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasDES<BasicTurbulenceModel>::fv2
(
    const volScalarField& chi,
    const volScalarField& fv1
) const
{
    return 1.0 - chi/(1.0 + chi*fv1);
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasDES<BasicTurbulenceModel>::ft2
(
    const volScalarField& chi
) const
{
    if (ft2Term_)
    {
        return Ct3_*exp(-Ct4_*sqr(chi));
    }
    else
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    type() + ":ft2",
                    this->time().timeName(),
                    this->mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                this->mesh(),
                dimensionedScalar("zero", dimless, 0.0)
            )
        );
    }
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasDES<BasicTurbulenceModel>::Omega
(
    const volTensorField& gradU
) const
{
    return sqrt(2.0)*mag(skew(gradU));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasDES<BasicTurbulenceModel>::Stilda
(
    const volScalarField& chi,
    const volScalarField& fv1,
    const volScalarField& Omega,
    const volScalarField& dTilda
) const
{
    return
    (
        max
        (
            Omega
          + fv2(chi, fv1)*nuTilda_/sqr(kappa_*dTilda),
            Cs_*Omega
        )
    );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasDES<BasicTurbulenceModel>::r
(
    const volScalarField& nur,
    const volScalarField& Stilda,
    const volScalarField& dTilda
) const
{
    tmp<volScalarField> tr
    (
        min
        (
            nur
           /(
                max
                (
                    Stilda,
                    dimensionedScalar("SMALL", Stilda.dimensions(), SMALL)
                )
               *sqr(kappa_*dTilda)
            ),
            scalar(10)
        )
    );
    tr.ref().boundaryFieldRef().forceAssign(0.0);

    return tr;
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasDES<BasicTurbulenceModel>::fw
(
    const volScalarField& Stilda,
    const volScalarField& dTilda
) const
{
    const volScalarField r(this->r(nuTilda_, Stilda, dTilda));
    const volScalarField g(r + Cw2_*(pow6(r) - r));

    return g*pow((1 + pow6(Cw3_))/(pow6(g) + pow6(Cw3_)), 1.0/6.0);
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasDES<BasicTurbulenceModel>::psi
(
    const volScalarField& chi,
    const volScalarField& fv1,
    const volTensorField& gradU
) const
{
    tmp<volScalarField> tpsi
    (
        new volScalarField
        (
            IOobject
            (
                type() + ":psi",
                this->time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedScalar("one", dimless, 1)
        )
    );

    if (lowReCorrection_)
    {
        volScalarField& psi = tpsi.ref();

        const volScalarField fv2(this->fv2(chi, fv1));
        const volScalarField ft2(this->ft2(chi));

        const volScalarField Omega(this->Omega(gradU));

        psi =
            sqrt
            (
                min
                (
                    scalar(100),
                    (
                        1 - Cb1(Omega)
                        /(Cw1(Omega)*sqr(kappa_)*fwStar_)
                        *(ft2 + (1 - ft2)*fv2)
                    )
                    /max(SMALL, (fv1 * max(scalar(1e-10), 1 - ft2)))
                )
            );
    }

    return tpsi;
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasDES<BasicTurbulenceModel>::dTilda
(
    const volScalarField& chi,
    const volScalarField& fv1,
    const volTensorField& gradU
) const
{
    tmp<volScalarField> tdTilda(psi(chi, fv1, gradU)*CDES_*this->delta());
    const volScalarField& y(wallDist::New(this->mesh_).y());
    min(tdTilda.ref().ref(), tdTilda(), y);

    return tdTilda;
}


template<class BasicTurbulenceModel>
void SpalartAllmarasDES<BasicTurbulenceModel>::correctNut
(
    const volScalarField& fv1
)
{
    this->nut_ = nuTilda_*fv1;
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_, this->nut_.db()).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


template<class BasicTurbulenceModel>
void SpalartAllmarasDES<BasicTurbulenceModel>::correctNut()
{
    correctNut(fv1(this->chi()));
}


template<class BasicTurbulenceModel>
geometricZeroField
SpalartAllmarasDES<BasicTurbulenceModel>::gradNuTildaDotGradRho
(
    const volVectorField& gradNuTilda,
    const geometricOneField& rho
) const
{
    return geometricZeroField();
}


template<class BasicTurbulenceModel>
tmp<volScalarField>
SpalartAllmarasDES<BasicTurbulenceModel>::gradNuTildaDotGradRho
(
    const volVectorField& gradNuTilda,
    const volScalarField& rho
) const
{
    return
    (
        (gradNuTilda & fvc::grad(rho))
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
SpalartAllmarasDES<BasicTurbulenceModel>::SpalartAllmarasDES
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    DESModel<BasicTurbulenceModel>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    sigmaNut_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaNut",
            this->coeffDict_,
            0.66666
        )
    ),
    sigmaNu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaNu",
            this->coeffDict_,
            0.66666
        )
    ),
    kappa_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "kappa",
            this->coeffDict_,
            0.41
        )
    ),
    Cb1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cb1",
            this->coeffDict_,
            0.1355
        )
    ),
    Cb2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cb2",
            this->coeffDict_,
            0.622
        )
    ),
    Cw2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cw2",
            this->coeffDict_,
            0.3
        )
    ),
    Cw3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cw3",
            this->coeffDict_,
            2.0
        )
    ),
    Cv1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cv1",
            this->coeffDict_,
            7.1
        )
    ),
    Cs_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cs",
            this->coeffDict_,
            0.3
        )
    ),
    CDES_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CDES",
            this->coeffDict_,
            0.65
        )
    ),
    ck_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ck",
            this->coeffDict_,
            0.07
        )
    ),
    ft2Term_
    (
        Switch::lookupOrAddToDict
        (
            "ft2Term",
            this->coeffDict_,
            false
        )
    ),
    gradUstd_
    (
        Switch::lookupOrAddToDict
        (
            "gradUstd",
            this->coeffDict_,
            false
        )
    ),
    lowReCorrection_
    (
        Switch::lookupOrAddToDict
        (
            "lowReCorrection",
            this->coeffDict_,
            true
        )
    ),
    Ct3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ct3",
            this->coeffDict_,
            1.2
        )
    ),
    Ct4_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ct4",
            this->coeffDict_,
            0.5
        )
    ),
    fwStar_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "fwStar",
            this->coeffDict_,
            0.4241
        )
    ),

    nuTilda_
    (
        IOobject
        (
            "nuTilda",
            this->runTime_.timeName(),
            this->db_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    )
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool SpalartAllmarasDES<BasicTurbulenceModel>::read()
{
    if (DESModel<BasicTurbulenceModel>::read())
    {
        sigmaNut_.readIfPresent(this->coeffDict());
        sigmaNu_.readIfPresent(this->coeffDict());
        kappa_.readIfPresent(*this);

        Cb1_.readIfPresent(this->coeffDict());
        Cb2_.readIfPresent(this->coeffDict());
        Cw2_.readIfPresent(this->coeffDict());
        Cw3_.readIfPresent(this->coeffDict());
        Cv1_.readIfPresent(this->coeffDict());
        Cs_.readIfPresent(this->coeffDict());

        CDES_.readIfPresent(this->coeffDict());
        ck_.readIfPresent(this->coeffDict());
        ft2Term_.readIfPresent("ft2Term", this->coeffDict());
        gradUstd_.readIfPresent("gradUstd", this->coeffDict());

        lowReCorrection_.readIfPresent("lowReCorrection", this->coeffDict());
        Ct3_.readIfPresent(this->coeffDict());
        Ct4_.readIfPresent(this->coeffDict());
        fwStar_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasDES<BasicTurbulenceModel>::
DnuTildaEff() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            "DnuTildaEff",
            nuTilda_/sigmaNut_ + this->nu()/sigmaNu_
        )
    );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasDES<BasicTurbulenceModel>::k() const
{
    const volScalarField chi(this->chi());
    const volScalarField fv1(this->fv1(chi));

    return sqr(this->nut()/ck_/dTilda(chi, fv1, gradUnwc()));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmarasDES<BasicTurbulenceModel>::LESRegion() const
{
    const volScalarField chi(this->chi());
    const volScalarField fv1(this->fv1(chi));
    const volScalarField& y(wallDist::New(this->mesh_).y());

    tmp<volScalarField> tLESRegion
    (
        new volScalarField
        (
            IOobject
            (
                "DES::LESRegion",
                this->mesh_.time().timeName(),
                this->db_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            neg(dTilda(chi, fv1, gradUnwc()) - y)
        )
    );

    return tLESRegion;
}


template<class BasicTurbulenceModel>
void SpalartAllmarasDES<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    fv::options& fvOptions(fv::options::New(this->mesh_, nuTilda_.db()));

    DESModel<BasicTurbulenceModel>::correct();

    const volScalarField chi(this->chi());
    const volScalarField fv1(this->fv1(chi));
    const volScalarField ft2(this->ft2(chi));

    tmp<volTensorField> tgradU = gradUnwc();
    const volScalarField Omega(this->Omega(tgradU()));
    const volScalarField dTilda("dTilda", this->dTilda(chi, fv1, tgradU()));
    if (this->mesh_.time().outputTime())
    {
        dTilda.write();
    }
    const volScalarField Stilda(this->Stilda(chi, fv1, Omega, dTilda));
    const volVectorField gradNuTilda(fvc::grad(nuTilda_));

    tmp<fvScalarMatrix> nuTildaEqn
    (
        fvm::ddt(alpha, rho, nuTilda_)
      + fvm::div(alphaRhoPhi, nuTilda_)
      - fvm::laplacian(alpha*rho*DnuTildaEff(), nuTilda_)
      - Cb2_/sigmaNut_*alpha*rho*magSqr(gradNuTilda)
      // - Cb2_/sigmaNut_*alpha*(fvc::grad(nuTilda_) & fvc::grad(rho()*nuTilda_))
     ==
        Cb1(Omega)*alpha*rho*Stilda*nuTilda_*(scalar(1)-ft2)
      - fvm::Sp
        (
            (Cw1(Omega)*fw(Stilda, dTilda)-Cb1(Omega)/sqr(kappa_)*ft2)
            *alpha*rho*nuTilda_/sqr(dTilda),
            nuTilda_
        )
      + fvOptions(alpha, rho, nuTilda_)
    );

    if (!isIsochoric())
    {
        nuTildaEqn.ref() +=
            (nuTilda_/sigmaNut_ + this->nu()/sigmaNu_)*alpha
            *gradNuTildaDotGradRho(gradNuTilda, rho);
    }

    nuTildaEqn.ref().relax();
    fvOptions.constrain(nuTildaEqn.ref());
    solve(nuTildaEqn);
    fvOptions.correct(nuTilda_);
    bound(nuTilda_, dimensionedScalar("0", nuTilda_.dimensions(), 0.0));
    nuTilda_.correctBoundaryConditions();

    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
