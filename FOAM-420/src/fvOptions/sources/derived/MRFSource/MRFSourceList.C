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
    (c) 2012-2016 OpenFOAM Foundation
    (c) 2017-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "MRFSourceList.H"
#include "fields/fvsPatchFields/basic/fixedValue/fixedValueFvsPatchFields.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //
void Foam::fv::MRFSourceList::initialise()
{
    forAll(*this, i)
    {
        MRFSource& pm = this->operator[](i);
        pm.initialise();
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::MRFSourceList::MRFSourceList
(
    const objectRegistry& obr,
    const dictionary& dict
)
:
    PtrList<MRFSource>(),
    obr_(obr)
{
    reset(dict);

    active(true);
    initialise();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::MRFSourceList::~MRFSourceList()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::MRFSourceList::active(const bool warn) const
{
    bool a = false;
    forAll(*this, i)
    {
        a = a || this->operator[](i).active();
    }

    if (warn && this->size() && !a)
    {
        Info<< "    No MRF sources active" << endl;
    }

    return a;
}


void Foam::fv::MRFSourceList::reset(const dictionary& dict)
{
    label count = 0;
    forAllConstIter(dictionary, dict, iter)
    {
        if (iter().isDict())
        {
            word modelType = iter().dict().lookup("type");

            if (modelType == "MRFSource")
            {
                count++;
            }
        }
    }

    this->setSize(count);

    label i = 0;
    forAllConstIter(dictionary, dict, iter)
    {
        if (iter().isDict())
        {
            const word& name = iter().keyword();

            word modelType = iter().dict().lookup("type");
            if (modelType == "MRFSource")
            {

                const dictionary& modelDict =
                    iter().dict().subDict(modelType + "Coeffs");

                bool active =
                    iter().dict().lookupOrDefault<bool>("active", true);

                if
                (
                  active
                )
                {
                    this->set
                    (
                        i++,
                        new MRFSource(name, modelType, modelDict, obr_)
                    );
                }
            }
        }
    }
}


const Foam::fv::MRFSource& Foam::fv::MRFSourceList::getFromName
(
    const word& name
) const
{
    DynamicList<word> names;
    for (const auto& mrf: *this)
    {
        if (mrf.name() == name)
        {
            return mrf;
        }

        names.append(mrf.name());
    }

    FatalErrorInFunction
        << "Unable to find MRFSource" << name
        << ". Available zones are: " << names
        << exit(FatalError);

    return first();
}


bool Foam::fv::MRFSourceList::read(const dictionary& dict)
{
    bool allOk = true;
    forAll(*this, i)
    {
        MRFSource& pm = this->operator[](i);
        bool ok = pm.read(dict.subDict(pm.name()));
        allOk = (allOk && ok);
    }
    return allOk;
}


bool Foam::fv::MRFSourceList::writeData(Ostream& os) const
{
    forAll(*this, i)
    {
        os  << nl;
        this->operator[](i).writeData(os);
    }

    return os.good();
}


void Foam::fv::MRFSourceList::addAcceleration
(
    const volVectorField& U,
    volVectorField& ddtU
) const
{
    forAll(*this, i)
    {
        operator[](i).addAcceleration(U, ddtU);
    }
}


void Foam::fv::MRFSourceList::addAcceleration
(
    fvVectorMatrix& UEqn,
    bool rhs
) const
{
    forAll(*this, i)
    {
        operator[](i).addAcceleration(UEqn, rhs);
    }
}


void Foam::fv::MRFSourceList::addAcceleration
(
    const volScalarField& rho,
    fvVectorMatrix& UEqn,
    bool rhs
) const
{
    forAll(*this, i)
    {
        operator[](i).addAcceleration(rho, UEqn, rhs);
    }
}


Foam::tmp<Foam::volVectorField> Foam::fv::MRFSourceList::DDt
(
    const volVectorField& U
) const
{
    tmp<volVectorField> tacceleration
    (
        new volVectorField
        (
            IOobject
            (
                "MRFSourceList:acceleration",
                U.mesh().time().timeName(),
                obr_
            ),
            U.mesh(),
            dimensionedVector("0", U.dimensions()/dimTime, Zero)
        )
    );
    volVectorField& acceleration = tacceleration.ref();

    forAll(*this, i)
    {
        operator[](i).addAcceleration(U, acceleration);
    }

    return tacceleration;
}


Foam::tmp<Foam::volVectorField> Foam::fv::MRFSourceList::DDt
(
    const volScalarField& rho,
    const volVectorField& U
) const
{
    return rho*DDt(U);
}


void Foam::fv::MRFSourceList::makeRelative(volVectorField& U) const
{
    forAll(*this, i)
    {
        operator[](i).makeRelative(U);
    }
}


void Foam::fv::MRFSourceList::makeRelative(surfaceScalarField& phi) const
{
    forAll(*this, i)
    {
        operator[](i).makeRelative(phi);
    }
}


Foam::tmp<Foam::surfaceScalarField> Foam::fv::MRFSourceList::relative
(
    const tmp<surfaceScalarField>& tphi
) const
{
    if (size())
    {
        tmp<surfaceScalarField> rphi
        (
            New
            (
                tphi,
                "relative(" + tphi().name() + ')',
                tphi().dimensions(),
                true
            )
        );

        makeRelative(rphi.ref());

        tphi.clear();

        return rphi;
    }
    else
    {
        return tmp<surfaceScalarField>(tphi, true);
    }
}


Foam::tmp<Foam::FieldField<Foam::fvsPatchField, Foam::scalar>>
Foam::fv::MRFSourceList::relative
(
    const tmp<FieldField<fvsPatchField, scalar>>& tphi
) const
{
    if (size())
    {
        tmp<FieldField<fvsPatchField, scalar>> rphi(New(tphi, true));

        forAll(*this, i)
        {
            operator[](i).makeRelative(rphi.ref());
        }

        tphi.clear();

        return rphi;
    }
    else
    {
        return tmp<FieldField<fvsPatchField, scalar>>(tphi, true);
    }
}


Foam::tmp<Foam::Field<Foam::scalar>>
Foam::fv::MRFSourceList::relative
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
            operator[](i).makeRelative(rphi.ref(), patchi);
        }

        tphi.clear();

        return rphi;
    }
    else
    {
        return tmp<Field<scalar>>(tphi, true);
    }
}


void Foam::fv::MRFSourceList::makeRelative
(
    const surfaceScalarField& rho,
    surfaceScalarField& phi
) const
{
    forAll(*this, i)
    {
        operator[](i).makeRelative(rho, phi);
    }
}


void Foam::fv::MRFSourceList::makeAbsolute(volVectorField& U) const
{
    forAll(*this, i)
    {
        operator[](i).makeAbsolute(U);
    }
}


void Foam::fv::MRFSourceList::makeAbsolute(surfaceScalarField& phi) const
{
    forAll(*this, i)
    {
        operator[](i).makeAbsolute(phi);
    }
}


Foam::tmp<Foam::surfaceScalarField> Foam::fv::MRFSourceList::absolute
(
    const tmp<surfaceScalarField>& tphi
) const
{
    if (size())
    {
        tmp<surfaceScalarField> rphi
        (
            New
            (
                tphi,
                "absolute(" + tphi().name() + ')',
                tphi().dimensions(),
                true
            )
        );

        makeAbsolute(rphi.ref());

        tphi.clear();

        return rphi;
    }
    else
    {
        return tmp<surfaceScalarField>(tphi, true);
    }
}


void Foam::fv::MRFSourceList::makeAbsolute
(
    const surfaceScalarField& rho,
    surfaceScalarField& phi
) const
{
    forAll(*this, i)
    {
        operator[](i).makeAbsolute(rho, phi);
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::fv::operator<<
(
    Ostream& os,
    const MRFSourceList& models
)
{
    models.writeData(os);
    return os;
}


// ************************************************************************* //
