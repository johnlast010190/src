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
    (c) 2010-2016 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "SpalartAllmaras.H"
#include "cfdTools/general/fvOptions/fvOptions.H"
#include "cfdTools/general/bound/bound.H"
#include "fvMesh/wallDist/wallDist/wallDist.H"
#include "fvMesh/fvPatches/derived/wall/wallFvPatch.H"
#include "fields/fvPatchFields/basic/fixedValue/fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volTensorField> SpalartAllmaras<BasicTurbulenceModel>::gradUnwc() const
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

    const tmp<volScalarField> tnuEff = this->nuEff();

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
                    (tnuEff->boundaryField()[pI][pfI]/tnuEff()[fcI])
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
tmp<volScalarField> SpalartAllmaras<BasicTurbulenceModel>::chi() const
{
    return nuTilda_/this->nu();
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmaras<BasicTurbulenceModel>::fv1
(
    const volScalarField& chi
) const
{
    const volScalarField chi3(pow3(chi));
    return chi3/(chi3 + pow3(Cv1_));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmaras<BasicTurbulenceModel>::fv2
(
    const volScalarField& chi,
    const volScalarField& fv1
) const
{
    return 1.0 - chi/(1.0 + chi*fv1);
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmaras<BasicTurbulenceModel>::ft2
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
tmp<volScalarField> SpalartAllmaras<BasicTurbulenceModel>::Stilda
(
    const volScalarField& chi,
    const volScalarField& fv1
) const
{
    volScalarField Omega(::sqrt(2.0)*mag(skew(gradUnwc())));
    const volScalarField& y(wallDist::New(this->mesh_).y());

    return
    (
        max
        (
            Omega
          + fv2(chi, fv1)*nuTilda_/sqr(kappa_*y),
            Cs_*Omega
        )
    );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmaras<BasicTurbulenceModel>::fw
(
    const volScalarField& Stilda
) const
{
    const volScalarField& y(wallDist::New(this->mesh_).y());
    volScalarField r
    (
        min
        (
            nuTilda_
           /(
               max
               (
                   Stilda,
                   dimensionedScalar("SMALL", Stilda.dimensions(), SMALL)
               )
              *sqr(kappa_*y)
            ),
            scalar(10.0)
        )
    );
    r.boundaryFieldRef().forceAssign(0.0);

    const volScalarField g(r + Cw2_*(pow6(r) - r));

    return g*pow((1.0 + pow6(Cw3_))/(pow6(g) + pow6(Cw3_)), 1.0/6.0);
}


template<class BasicTurbulenceModel>
void SpalartAllmaras<BasicTurbulenceModel>::correctNut
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
void SpalartAllmaras<BasicTurbulenceModel>::correctNut()
{
    correctNut(fv1(this->chi()));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
SpalartAllmaras<BasicTurbulenceModel>::SpalartAllmaras
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
    eddyViscosity<RASModel<BasicTurbulenceModel>>
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
    Cw1_(Cb1_/sqr(kappa_) + (1.0 + Cb2_)/sigmaNut_),
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
    ft2Term_
    (
        Switch::lookupOrAddToDict
        (
            "ft2Term",
            this->coeffDict_,
            true
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
bool SpalartAllmaras<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
        sigmaNut_.readIfPresent(this->coeffDict());
        kappa_.readIfPresent(this->coeffDict());

        Cb1_.readIfPresent(this->coeffDict());
        Cb2_.readIfPresent(this->coeffDict());
        Cw1_ = Cb1_/sqr(kappa_) + (1.0 + Cb2_)/sigmaNut_;
        Cw2_.readIfPresent(this->coeffDict());
        Cw3_.readIfPresent(this->coeffDict());
        Cv1_.readIfPresent(this->coeffDict());
        Cs_.readIfPresent(this->coeffDict());

        ft2Term_.readIfPresent("ft2Term", this->coeffDict());
        gradUstd_.readIfPresent("gradUstd", this->coeffDict());

        Ct3_.readIfPresent(this->coeffDict());
        Ct4_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmaras<BasicTurbulenceModel>::DnuTildaEff() const
{
    return tmp<volScalarField>
    (
        new volScalarField("DnuTildaEff", (nuTilda_ + this->nu())/sigmaNut_)
    );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmaras<BasicTurbulenceModel>::k() const
{
    WarningInFunction
        << "Turbulence kinetic energy not defined for "
        << "Spalart-Allmaras model. Returning zero field"
        << endl;

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "k",
                this->runTime_.timeName(),
                this->db_
            ),
            this->mesh_,
            dimensionedScalar("0", dimensionSet(0, 2, -2, 0, 0), 0),
            zeroGradientFvPatchField<vector>::typeName
        )
    );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> SpalartAllmaras<BasicTurbulenceModel>::epsilon() const
{
    WarningInFunction
        << "Turbulence kinetic energy dissipation rate not defined for "
        << "Spalart-Allmaras model. Returning zero field"
        << endl;

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "epsilon",
                this->runTime_.timeName(),
                this->db_
            ),
            this->mesh_,
            dimensionedScalar("0", dimensionSet(0, 2, -3, 0, 0), 0),
            zeroGradientFvPatchField<vector>::typeName
        )
    );
}


template<class BasicTurbulenceModel>
tmp<scalarField> SpalartAllmaras<BasicTurbulenceModel>::yPlus
(
    const label patchi,
    const scalar Cmu
) const
{
    vectorField nf(this->mesh_.boundary()[patchi].nf()());

    const volScalarField& yvsf(wallDist::New(this->mesh_).y());
    const scalarField y(yvsf.boundaryField()[patchi].patchInternalField());
    const fvPatchVectorField& Uw = this->U().boundaryField()[patchi];
    const scalarField nuw(this->nu(patchi));
    const scalarField nuEffw(this->nuEff(patchi));

    scalarField magGradU(this->size(), 0.0);
    {
        vectorField gradU(Uw.snGrad());
        gradU -= nf*(gradU & nf);
        magGradU = mag(gradU);
        magGradU = max(magGradU, SMALL);
    }

    return y*sqrt( max(SMALL,nuEffw) * magGradU)/nuw;
}


template<class BasicTurbulenceModel>
void SpalartAllmaras<BasicTurbulenceModel>::correct()
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

    eddyViscosity<RASModel<BasicTurbulenceModel>>::correct();

    const volScalarField chi(this->chi());
    const volScalarField fv1(this->fv1(chi));
    const volScalarField ft2(this->ft2(chi));

    const volScalarField Stilda(this->Stilda(chi, fv1));

    const volScalarField& y(wallDist::New(this->mesh_).y());

    tmp<fvScalarMatrix> nuTildaEqn
    (
        fvm::ddt(alpha, rho, nuTilda_)
      + fvm::div(alphaRhoPhi, nuTilda_)
      - fvm::laplacian(alpha*rho*DnuTildaEff(), nuTilda_)
      - Cb2_/sigmaNut_*alpha*rho*magSqr(fvc::grad(nuTilda_))
     ==
        Cb1_*alpha*rho*Stilda*nuTilda_*(scalar(1)-ft2)
      - fvm::Sp
        (
            (Cw1_*fw(Stilda)-Cb1_/sqr(kappa_)*ft2)*alpha*rho*nuTilda_/sqr(y),
            nuTilda_
        )
        + fvOptions(alpha, rho, nuTilda_)
    );

    nuTildaEqn.ref().relax();
    fvOptions.constrain(nuTildaEqn.ref());
    solve(nuTildaEqn);
    fvOptions.correct(nuTilda_);
    bound(nuTilda_, dimensionedScalar("0", nuTilda_.dimensions(), 0.0));
    nuTilda_.correctBoundaryConditions();

    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
