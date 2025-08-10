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
    (c) 2016 Esi Ltd.
    (c) 2011-2015 OpenFOAM Foundation
    (c) 2019 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/fvOptions/fvOptionList.H"
#include "fvMatrices/fvMatrix/fvMatrix.H"
#include "fields/fvsPatchFields/basic/fixedValue/fixedValueFvsPatchFields.H"

// * * * * * * * * * * * * * * *  Frame Functions  * * * * * * * * * * * * * //



void Foam::fv::optionList::addAcceleration
(
    const volVectorField& U,
    volVectorField& ddtU
) const
{
    forAll(*this, i)
    {
        if (this->operator[](i).isMRF())
        {
            operator[](i).addAcceleration(U, ddtU);
        }
    }
}


void Foam::fv::optionList::addAcceleration
(
    fvVectorMatrix& UEqn,
    bool rhs
) const
{
    forAll(*this, i)
    {
        if (this->operator[](i).isMRF())
        {
            operator[](i).addAcceleration(UEqn, rhs);
        }
    }
}


void Foam::fv::optionList::addAcceleration
(
    const volScalarField& rho,
    fvVectorMatrix& UEqn,
    bool rhs
) const
{
    forAll(*this, i)
    {
        if (this->operator[](i).isMRF())
        {
            operator[](i).addAcceleration(rho, UEqn, rhs);
        }
    }
}


void Foam::fv::optionList::addAcceleration(fvBlockMatrix<vector>& UEqn) const
{
    forAll(*this, i)
    {
        if (this->operator[](i).isMRF())
        {
            operator[](i).addAcceleration(UEqn);
        }
    }
}


void Foam::fv::optionList::addAcceleration
(
    const volScalarField& rho,
    fvBlockMatrix<vector>& UEqn
) const
{
    forAll(*this, i)
    {
        if (this->operator[](i).isMRF())
        {
            operator[](i).addAcceleration(rho, UEqn);
        }
    }
}


Foam::tmp<Foam::volVectorField>  Foam::fv::optionList::MRFDDt
(
    const volVectorField& U
)
{
    tmp<volVectorField> tacceleration
    (
        new volVectorField
        (
            IOobject
            (
                "MRFSourceList:acceleration",
                U.mesh().time().timeName(),
                U.mesh()
            ),
            U.mesh(),
            dimensionedVector("0", U.dimensions()/dimTime, Zero)
        )
    );
    volVectorField& acceleration = tacceleration.ref();

    forAll(*this, i)
    {
        option& source = this->operator[](i);

        if (source.isMRF())
        {
            if (source.isActive())
            {
                if (debug)
                {
                    Info<< "Applying source " << source.name() << " to field "
                        << U.name() << endl;
                }
                source.addAcceleration(U, acceleration);
            }
        }
    }

    return tacceleration;
}


Foam::tmp<Foam::volVectorField> Foam::fv::optionList::MRFDDt
(
    const volScalarField& rho,
    const volVectorField& U
)
{
    return rho*MRFDDt(U);
}


void Foam::fv::optionList::makeRelative(volVectorField& U) const
{
    forAll(*this, i)
    {
        if (this->operator[](i).isMRF())
        {
            operator[](i).makeRelative(U);
        }
    }
}


void Foam::fv::optionList::makeRelative(surfaceScalarField& phi) const
{
    forAll(*this, i)
    {
        if (this->operator[](i).isMRF())
        {
            operator[](i).makeRelative(phi);
        }
    }
}


Foam::tmp<Foam::surfaceScalarField> Foam::fv::optionList::relative
(
    const tmp<surfaceScalarField>& phi
) const
{
    tmp<surfaceScalarField> rphi(phi.ptr());
    makeRelative(rphi.ref());
    return rphi;
}


Foam::tmp<Foam::FieldField<Foam::fvsPatchField, Foam::scalar>>
Foam::fv::optionList::relative
(
    const tmp<FieldField<fvsPatchField, scalar>>& phi
) const
{
    tmp<FieldField<fvsPatchField, scalar>> rphi(phi.ptr());

    forAll(*this, i)
    {
        if (this->operator[](i).isMRF())
        {
            operator[](i).makeRelative(rphi.ref());
        }
    }

    return rphi;
}


Foam::tmp<Foam::Field<Foam::scalar>>
Foam::fv::optionList::relative
(
    const tmp<Field<scalar>>& tphi,
    const label patchi
) const
{
    if (size())
    {
        tmp<Field<scalar>> rphi(New(tphi, true));

        forAll(*this, i)
        {
            if (this->operator[](i).isMRF())
            {
                operator[](i).makeRelative(rphi.ref(), patchi);
            }
        }

        tphi.clear();

        return rphi;
    }
    else
    {
        return tmp<Field<scalar>>(tphi, true);
    }
}


void Foam::fv::optionList::makeRelative
(
    const surfaceScalarField& rho,
    surfaceScalarField& phi
) const
{
    forAll(*this, i)
    {
        if (this->operator[](i).isMRF())
        {
            operator[](i).makeRelative(rho, phi);
        }
    }
}


void Foam::fv::optionList::makeAbsolute(volVectorField& U) const
{
    forAll(*this, i)
    {
        if (this->operator[](i).isMRF())
        {
            operator[](i).makeAbsolute(U);
        }
    }
}


void Foam::fv::optionList::makeAbsolute(surfaceScalarField& phi) const
{
    forAll(*this, i)
    {
        if (this->operator[](i).isMRF())
        {
            operator[](i).makeAbsolute(phi);
        }
    }
}


Foam::tmp<Foam::surfaceScalarField> Foam::fv::optionList::absolute
(
    const tmp<surfaceScalarField>& phi
) const
{
    tmp<surfaceScalarField> rphi(phi.ptr());
    makeAbsolute(rphi.ref());
    return rphi;
}


void Foam::fv::optionList::makeAbsolute
(
    const surfaceScalarField& rho,
    surfaceScalarField& phi
) const
{
    forAll(*this, i)
    {
        if (this->operator[](i).isMRF())
        {
            operator[](i).makeAbsolute(rho, phi);
        }
    }
}


void Foam::fv::optionList::correctBoundaryFlux
(
    const volVectorField& U,
    surfaceScalarField& phi
) const
{

    FieldField<fvsPatchField, scalar> phibf
    (
        relative(mesh_.Sf().boundaryField() & U.boundaryField())
    );

    forAll(mesh_.boundary(), patchi)
    {
        if
        (
            isA<fixedValueFvsPatchScalarField>(phi.boundaryField()[patchi])
        )
        {
            phi.boundaryFieldRef()[patchi].forceAssign(phibf[patchi]);
        }
    }

}

Foam::tmp<Foam::surfaceScalarField> Foam::fv::optionList::zeroFilter
(
    const tmp<surfaceScalarField>& tphi
) const
{
    if (size())
    {
        tmp<surfaceScalarField> zphi
        (
            New
            (
                tphi,
                "zeroFilter(" + tphi().name() + ')',
                tphi().dimensions(),
                true
            )
        );

        forAll(*this, i)
        {
            if (this->operator[](i).isMRF())
            {
                operator[](i).zero(zphi.ref());
            }
        }

        return zphi;
    }
    else
    {
        return tmp<surfaceScalarField>(tphi, true);
    }
}
// ************************************************************************* //
