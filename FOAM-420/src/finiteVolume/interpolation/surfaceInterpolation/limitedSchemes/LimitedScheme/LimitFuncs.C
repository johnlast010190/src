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
    (c) 2011-2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
inline Foam::tmp<Foam::volScalarField>
Foam::limitFuncs::magSqr<Type>::operator()
(
    const GeometricField<Type, fvPatchField, volMesh>& phi
) const
{
    return Foam::magSqr(phi);
}


template<>
inline Foam::tmp<Foam::volScalarField>
Foam::limitFuncs::magSqr<Foam::scalar>::operator()
(
    const volScalarField& phi
) const
{
    return phi;
}


template<>
inline Foam::tmp<Foam::volScalarField>
Foam::limitFuncs::magSqr<Foam::symmTensor>::operator()
(
    const volSymmTensorField& phi
) const
{
    return Foam::tr(phi);
}


template<>
inline Foam::tmp<Foam::volScalarField>
Foam::limitFuncs::magSqr<Foam::tensor>::operator()
(
    const volTensorField& phi
) const
{
    return Foam::tr(phi);
}


template<class Type>
inline Foam::tmp<Foam::volScalarField>
Foam::limitFuncs::rhoMagSqr<Type>::operator()
(
    const GeometricField<Type, fvPatchField, volMesh>& phi
) const
{
    const volScalarField& rho =
        phi.db().objectRegistry::template lookupObject<volScalarField>("rho");
    return Foam::magSqr(phi/rho);
}


template<>
inline Foam::tmp<Foam::volScalarField>
Foam::limitFuncs::rhoMagSqr<Foam::scalar>::operator()
(
    const volScalarField& phi
) const
{
    const volScalarField& rho =
        phi.db().objectRegistry::lookupObject<volScalarField>("rho");
    return phi/rho;
}


// ************************************************************************* //
