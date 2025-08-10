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
    (c) 2017 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/ABLProfile/tabulatedProfile.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
namespace ABLProfiles
{

defineTypeNameAndDebug(tabulatedProfile, 0);
addToRunTimeSelectionTable(ABLProfile, tabulatedProfile, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void tabulatedProfile::applyFunction1ToField
(
    const Function1<Type>& data,
    const scalarField& Z,
    Field<Type>& field
)
{
    field.setSize(Z.size());

    forAll(field, i)
    {
        field[i] = data.value(Z[i]);
    }
}

tmp<tensorField> tabulatedProfile::tableGradient
(
    const Function1<vector>& Uf1,
    const patchDistanceFunction& Z

) const
{
    scalarField Lf ( sqrt(patch_.magSf()) );

    tmp<tensorField> gradU(new tensorField(patch_.size(), tensor::zero));
    tmp<vectorField> direction(Z.dir());

    forAll(Lf, fI)
    {
        scalar Zp = Z.value()[fI] + Lf[fI];
        scalar Zm = Z.value()[fI] - Lf[fI];
        gradU.ref()[fI] += (direction()[fI]*(Uf1.value(Zp) - Uf1.value(Zm)));
    }

    return gradU;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

tabulatedProfile::tabulatedProfile
(
    const fvPatch& patch
)
:
    ABLProfile(patch),
    Cmu_(0.0),
    Uf1_(),
    kf1_(),
    if1_(),
    epsilonf1_(),
    Lf1_(),
    Lxyzf1_(),
    Rf1_(),
    fLxyz_(vector::one)
{}


tabulatedProfile::tabulatedProfile
(
    const fvPatch& patch,
    const dictionary& dict
)
:
    ABLProfile(patch, dict),
    Cmu_
    (
        dict.subDict(this->typeName + "Coeffs")
        .lookupOrDefault<scalar>("Cmu", 0.09)
    ),
    Uf1_
    (
        Function1<vector>::New("U", dict.subDict(this->typeName + "Coeffs"))
    ),
    kf1_(),
    if1_(),
    epsilonf1_(),
    Lf1_(),
    Lxyzf1_(),
    Rf1_(),
    fLxyz_
    (
        dict.subDict(this->typeName + "Coeffs")
        .lookupOrDefault<vector>("fLxyz", vector::one)
    )
{
    //read function tables
    const dictionary& coeffDict(dict.subDict(this->typeName + "Coeffs"));

    if (coeffDict.found("k"))
    {
        kf1_ = Function1<scalar>::New("k", coeffDict);
    }
    if (coeffDict.found("i"))
    {
        if1_ = Function1<scalar>::New("i", coeffDict);
    }
    if (coeffDict.found("epsilon"))
    {
        epsilonf1_ = Function1<scalar>::New("epsilon", coeffDict);
    }
    if (coeffDict.found("L"))
    {
        Lf1_ = Function1<scalar>::New("L", coeffDict);
    }
    if (coeffDict.found("Lxyz"))
    {
        Lxyzf1_ = Function1<vector>::New("Lxyz", coeffDict);
    }
    if (coeffDict.found("R"))
    {
        Rf1_ = Function1<symmTensor>::New("R", coeffDict);
    }


    //check combination
    //valid turbulent kinetic energy
    if
    (
        (Rf1_.valid() || if1_.valid() || kf1_.valid()
         || (epsilonf1_.valid() && (Lf1_.valid() || Lxyzf1_.valid())))
        && (epsilonf1_.valid() || (Lf1_.valid() || Lxyzf1_.valid()))
    )
    {}
    else
    {
        FatalError << "Unsupported combination of user defined data." << nl
            << "Supported combinations are: " << nl
            << tab << "U, R, epsilon" << nl
            << tab << "U, R, L" << nl
            << tab << "U, k, epsilon" << nl
            << tab << "U, k, L" << nl
            << tab << "U, i, epsilon" << nl
            << tab << "U, i, L" << nl
            << tab << "U, epsilon, L" << exit(FatalError);
    }

}

tabulatedProfile::tabulatedProfile
(
    const tabulatedProfile& abl
)
:
    ABLProfile(abl),
    Cmu_(abl.Cmu_),
    Uf1_(abl.Uf1_, false),
    kf1_(abl.kf1_.valid() ? abl.kf1_->clone().ptr() : nullptr),
    if1_(abl.if1_.valid() ? abl.if1_->clone().ptr() : nullptr),
    epsilonf1_
    (
        abl.epsilonf1_.valid() ? abl.epsilonf1_->clone().ptr() : nullptr
    ),
    Lf1_(abl.Lf1_.valid() ? abl.Lf1_->clone().ptr() : nullptr),
    Lxyzf1_(abl.Lxyzf1_.valid() ? abl.Lxyzf1_->clone().ptr() : nullptr),
    Rf1_(abl.Rf1_.valid() ? abl.Rf1_->clone().ptr() : nullptr),
    fLxyz_(abl.fLxyz_)
{}

// * * * * * * * * * * * * * * * * * Functions* * * * * * * * * * * * * * * //

void tabulatedProfile::setProfiles
(
    vectorField& U,
    scalarField& L,
    symmTensorField& R
)
{
    // set U
    applyFunction1ToField(Uf1_(), Z_->value(), U);


    if (Rf1_.valid())
    {
        applyFunction1ToField(Rf1_(), Z_->value(), R);

        //only need  L
        if (!Lf1_.valid())
        {
            //epsilon must be valid
            scalarField k ( 0.5*tr(R) );
            scalarField epsilon(k.size(), Zero);
            applyFunction1ToField(epsilonf1_(), Z_->value(), epsilon);

            L = k*sqrt(k)/epsilon * pow(Cmu_,0.75);
        }
        else
        {
            applyFunction1ToField(Lf1_(), Z_->value(), L);
        }
    }
    else //need to calculate R, need gradU, k and nut
    {
        scalarField k(patch_.size(), Zero);
        scalarField nut(patch_.size(), 0.0);

        //calculate k
        if (kf1_.valid())
        {
            applyFunction1ToField(kf1_(), Z_->value(), k);

        }
        else if (if1_.valid())
        {
            scalarField i(patch_.size(), 0.0);
            applyFunction1ToField(if1_(), Z_->value(), i);
            k = 3.0/2.0*sqr(i*mag(U));

        }
        else if (epsilonf1_.valid() && Lf1_.valid())
        {
            scalarField epsilon(patch_.size(), Zero);

            applyFunction1ToField(epsilonf1_(), Z_->value(), epsilon);
            applyFunction1ToField(Lf1_(), Z_->value(), L);

            k = pow(pow(Cmu_,-0.75)*L*epsilon, (2.0/3.0));
        }

        //calculate nut
        if (epsilonf1_.valid())
        {
            scalarField epsilon(patch_.size(), Zero);
            applyFunction1ToField(epsilonf1_(), Z_->value(), epsilon);

            nut = Cmu_*sqr(k)/epsilon;
            L = k* sqrt(k)/epsilon *pow(Cmu_, 0.75);
        }
        else if (Lf1_.valid())
        {
            applyFunction1ToField(Lf1_(), Z_->value(), L);

            nut = pow(Cmu_, 0.25)*sqrt(k)*L;
        }

        //generate grad(U) tensor
        {
            tmp<tensorField> gradU = tableGradient(Uf1_(), Z_());

            R = 2.0/3.0*I*k - nut*symm(gradU);
        }
    }
}

void tabulatedProfile::setProfiles
(
    vectorField& U,
    vectorField& L,
    symmTensorField& R
)
{
    // set U
    applyFunction1ToField(Uf1_(), Z_->value(), U);


    if (Rf1_.valid())
    {
        applyFunction1ToField(Rf1_(), Z_->value(), R);

        //only need  L
        if (!Lxyzf1_.valid())
        {
            if (Lf1_.valid())
            {
                FatalErrorInFunction
                    << "anisotropic Lxyz should be specified instead of L: "
                    << abort(FatalError);
            }

            //epsilon must be valid
            scalarField k ( 0.5*tr(R) );
            scalarField epsilon(k.size(), Zero);
            applyFunction1ToField(epsilonf1_(), Z_->value(), epsilon);

            L = fLxyz_*k*sqrt(k)/epsilon * pow(Cmu_,0.75);
        }
        else
        {
            applyFunction1ToField(Lxyzf1_(), Z_->value(), L);
        }
    }
    else //need to calculate R, need gradU, k and nut
    {
        scalarField k(patch_.size(), Zero);
        scalarField nut(patch_.size(), 0.0);

        //calculate k
        if (kf1_.valid())
        {
            applyFunction1ToField(kf1_(), Z_->value(), k);

        }
        else if (if1_.valid())
        {
            scalarField i(patch_.size(), 0.0);
            applyFunction1ToField(if1_(), Z_->value(), i);
            k = 3.0/2.0*sqr(i*mag(U));

        }
        else if (epsilonf1_.valid() && Lf1_.valid())
        {
            scalarField epsilon(patch_.size(), Zero);

            applyFunction1ToField(epsilonf1_(), Z_->value(), epsilon);
            applyFunction1ToField(Lxyzf1_(), Z_->value(), L);

            scalarField lscale(mag(L));
            k = pow(pow(Cmu_,-0.75)*lscale*epsilon, (2.0/3.0));
        }

        //calculate nut
        if (epsilonf1_.valid())
        {
            scalarField epsilon(patch_.size(), Zero);
            applyFunction1ToField(epsilonf1_(), Z_->value(), epsilon);

            nut = Cmu_*sqr(k)/epsilon;

            L = fLxyz_*k* sqrt(k)/epsilon *pow(Cmu_, 0.75);
        }
        else if (Lxyzf1_.valid())
        {
            applyFunction1ToField(Lxyzf1_(), Z_->value(), L);

            scalarField lscale(mag(L));
            nut = pow(Cmu_, 0.25)*sqrt(k)*lscale;
        }

        //generate grad(U) tensor
        {
            tmp<tensorField> gradU = tableGradient(Uf1_(), Z_());

            R = 2.0/3.0*I*k - nut*symm(gradU);
        }
    }
}

void tabulatedProfile::write(Ostream& os) const
{
    ABLProfile::write(os);

    os.writeEntry("profileType", typeName);

    os.beginBlock(word(typeName + word("Coeffs")));
    if (Cmu_ != 0.09)
    {
        os.writeEntry("Cmu", Cmu_);
    }
    if (kf1_.valid())
    {
        kf1_->writeData(os);
    }
    if (if1_.valid())
    {
        if1_->writeData(os);
    }
    if (epsilonf1_.valid())
    {
        epsilonf1_->writeData(os);
    }
    if (Uf1_.valid())
    {
        Uf1_->writeData(os);
    }
    if (Lf1_.valid())
    {
        Lf1_->writeData(os);
    }
    if (Lxyzf1_.valid())
    {
        Lxyzf1_->writeData(os);
    }
    if (Rf1_.valid())
    {
        Rf1_->writeData(os);
    }

    if (fLxyz_ != vector::one)
    {
        os.writeEntry("fLxyz", fLxyz_);
    }

    os.endBlock() << endl;
}

// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //



} // End namespace ABLProfiles
} // End namespace Foam

// ************************************************************************* //
