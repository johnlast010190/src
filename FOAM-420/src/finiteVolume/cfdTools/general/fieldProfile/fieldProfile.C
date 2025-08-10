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

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/fieldProfile/fieldProfile.H"
#include "primitives/functions/Function1/Constant/Constant.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fieldProfile<Type>::fieldProfile(const fvPatch& patch)
:
    patch_(patch),
    dict_(dictionary()),
    distribution_(new Function1Types::Constant<Type>("profile", Zero)),
    distance_(new patchDistanceFunction(patch)),
    xscale_(new Function1Types::Constant<scalar>("xscale", 1.0)),
    yscale_(new Function1Types::Constant<scalar>("yscale", 1.0)),
    xoffset_(Zero),
    yoffset_(new Function1Types::Constant<Type>("yoffset", Zero))
{}


template<class Type>
Foam::fieldProfile<Type>::fieldProfile
(
    const fvPatch& patch,
    const dictionary& dict
)
:
    patch_(patch),
    dict_(dict),
    distribution_(Function1<Type>::New("profile", dict)),
    distance_(new patchDistanceFunction(patch, dict)),
    xscale_
    (
        dict.found("xscale")
      ? Function1<scalar>::New("xscale", dict)
      : autoPtr<Function1<scalar>>
        (
            new Function1Types::Constant<scalar>("xscale", 1.0)
        )
    ),
    yscale_
    (
        dict.found("yscale")
      ? Function1<scalar>::New("yscale", dict)
      : autoPtr<Function1<scalar>>
        (
            new Function1Types::Constant<scalar>("yscale", 1.0)
        )
    ),
    xoffset_(dict.lookupOrDefault<scalar>("xoffset", 0.0)),
    yoffset_
    (
        dict.found("yoffset")
      ? Function1<Type>::New("yoffset", dict)
      : autoPtr<Function1<Type>>
        (
            new Function1Types::Constant<Type>("yoffset", Zero)
        )
    )
{}


template<class Type>
Foam::fieldProfile<Type>::fieldProfile(const fieldProfile& cpy)
:
    tmp<fieldProfile<Type>>::refCount(),
    patch_(cpy.patch_),
    dict_(cpy.dict_),
    distribution_(cpy.distribution_, false),
    distance_(cpy.distance_, false),
    xscale_(cpy.xscale_, false),
    yscale_(cpy.yscale_, false),
    xoffset_(cpy.xoffset_),
    yoffset_(cpy.yoffset_, false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::fieldProfile<Type>::update()
{
    return distance_->update();
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::fieldProfile<Type>::value() const
{
    scalar t = patch_.boundaryMesh().mesh().time().timeOutputValue();

    tmp<Field<Type>> f
    (
        distribution_->value((Z() - xoffset_)/xscale_->value(t))
      * yscale_->value(t) + yoffset_->value(t)
    );
    return f;
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::fieldProfile<Type>::gradient() const
{
    const scalar t = patch_.boundaryMesh().mesh().time().timeOutputValue();

    tmp<scalarField> deltaPtr(distance_->delta());

    tmp<Field<Type>> minVal
    (
        distribution_->value
        (
            (Z() - xoffset_ - deltaPtr())
           /xscale_->value(t)
        )
       *yscale_->value(t)
    );
    tmp<Field<Type>> maxVal
    (
        distribution_->value
        (
            (Z() - xoffset_ + deltaPtr())
           /xscale_->value(t)
        )
       *yscale_->value(t)
    );

    return tmp<Field<Type>>(0.5*(maxVal - minVal)/deltaPtr());
}


template<class Type>
void Foam::fieldProfile<Type>::write(Ostream& os, bool writeFrame) const
{
    distribution_->writeData(os);
    distance_->write(os, writeFrame);
    xscale_->writeData(os);
    yscale_->writeData(os);
    yoffset_->writeData(os);
    os.writeEntry("xoffset",  xoffset_);
}


// ************************************************************************* //
