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
    (c) 2016 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::injectedParticle::injectedParticle
(
    const polyMesh& mesh,
    const vector& position,
    const label tag,
    const scalar soi,
    const scalar d,
    const vector& U
)
:
    particle(mesh, position, -1, false),
    tag_(tag),
    soi_(soi),
    d_(d),
    U_(U)
{}


inline Foam::label Foam::injectedParticle::tag() const
{
    return tag_;
}


inline Foam::scalar Foam::injectedParticle::soi() const
{
    return soi_;
}


inline Foam::scalar Foam::injectedParticle::d() const
{
    return d_;
}


inline const Foam::vector& Foam::injectedParticle::U() const
{
    return U_;
}


inline Foam::label& Foam::injectedParticle::tag()
{
    return tag_;
}


inline Foam::scalar& Foam::injectedParticle::soi()
{
    return soi_;
}


inline Foam::scalar& Foam::injectedParticle::d()
{
    return d_;
}


inline Foam::vector& Foam::injectedParticle::U()
{
    return U_;
}


// ************************************************************************* //
