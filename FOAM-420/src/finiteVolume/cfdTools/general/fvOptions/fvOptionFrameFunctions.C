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

#include "cfdTools/general/fvOptions/fvOption.H"

// * * * * * * * * * * * * * * *  Frame Functions  * * * * * * * * * * * * * //

void Foam::fv::option::frameError(word classname) const
{
    word errmsg("void Foam::fv::option::" + classname + "(...)");
    FatalErrorIn(errmsg)
        << "Frame call not supported by non-MRF class: " << type() <<"."
        << exit(FatalError);
}

void Foam::fv::option::addAcceleration
(
    const volVectorField& U,
    volVectorField& ddtU
) const
{
    frameError("addAcceleration");
}


void Foam::fv::option::addAcceleration
(
    fvVectorMatrix& UEqn,
    bool rhs
) const
{
    frameError("addAcceleration");
}


void Foam::fv::option::addAcceleration
(
    const volScalarField& rho,
    fvVectorMatrix& UEqn,
    bool rhs
) const
{
    frameError("addAcceleration");
}


void Foam::fv::option::addAcceleration
(
    fvBlockMatrix<vector>& UEqn
) const
{
    frameError("addAcceleration");
}



void Foam::fv::option::addAcceleration
(
    const volScalarField& rho,
    fvBlockMatrix<vector>& UEqn
) const
{
    frameError("addAcceleration");
}


void Foam::fv::option::makeRelative(surfaceScalarField& phi) const
{
    frameError("makeRelative");
}


void Foam::fv::option::makeRelative
(
    FieldField<fvsPatchField, scalar>& phi
) const
{
    frameError("makeRelative");
}


void Foam::fv::option::makeRelative
(
    Field<scalar>& phi,
    const label patchi
) const
{
    frameError("makeRelative");
}


void Foam::fv::option::makeRelative
(
    const surfaceScalarField& rho,
    surfaceScalarField& phi
) const
{
    frameError("makeRelative");
}

void Foam::fv::option::makeRelative
(
    volVectorField& U
) const
{
    frameError("makeRelative");
}


void Foam::fv::option::makeAbsolute(surfaceScalarField& phi) const
{
    frameError("makeAbsolute");
}


void Foam::fv::option::makeAbsolute
(
    const surfaceScalarField& rho,
    surfaceScalarField& phi
) const
{
    frameError("makeAbsolute");
}

void Foam::fv::option::makeAbsolute
(
    volVectorField& U
) const
{
    frameError("makeAbsolute");
}

void Foam::fv::option::zero
(
	surfaceScalarField& tphi
) const
{
	frameError("zero");
}
// ************************************************************************* //
