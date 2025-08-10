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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::regionCoupledBase::interpolate
(
    const Field<Type>& fld
) const
{
    if (owner())
    {
        return AMI().interpolateToSource(fld);
    }
    else
    {
        return nbrPatch().AMI().interpolateToTarget(fld);
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::regionCoupledBase::interpolate
(
    const tmp<Field<Type>>& tFld
) const
{
    if (owner())
    {
        return AMI().interpolateToSource(tFld);
    }
    else
    {
        return nbrPatch().AMI().interpolateToTarget(tFld);
    }
}


template<class Type, class BinaryOp>
void Foam::regionCoupledBase::interpolate
(
    const UList<Type>& fld,
    const BinaryOp& bop,
    List<Type>& result
) const
{
    if (owner())
    {
        AMI().interpolateToSource(fld, bop, result);
    }
    else
    {
        nbrPatch().AMI().interpolateToTarget(fld, bop, result);
    }
}


// ************************************************************************* //
