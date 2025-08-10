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
    (c) 2015-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "primitives/ops/flipOp.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<>
Foam::scalar Foam::flipOp::operator()(const scalar& v) const
{
    return -v;
}


template<> Foam::vector Foam::flipOp::operator()(const vector& v) const
{
    return -v;
}


template<>Foam::sphericalTensor Foam::flipOp::operator()
(
    const sphericalTensor& v
) const
{
    return -v;
}


template<> Foam::symmTensor Foam::flipOp::operator()
(
    const symmTensor& v
) const
{
    return -v;
}


template<> Foam::tensor Foam::flipOp::operator()(const tensor& v) const
{
    return -v;
}


template<> Foam::triad Foam::flipOp::operator()
(
    const triad& v
) const
{
    return -v;
}


// ************************************************************************* //
