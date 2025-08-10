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

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::fv::interRegionHeatTransferModel::interpolate
(
    const interRegionHeatTransferModel& nbrModel,
    const Field<Type>& field
) const
{
    if (master_)
    {
        return meshInterp().mapTgtToSrc(field);
    }
    else
    {
        return (nbrModel.meshInterp().mapSrcToTgt(field));
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::fv::interRegionHeatTransferModel::interpolate
(
    const Field<Type>& field
) const
{
    return interpolate(nbrModel(), field);
}


template<class Type>
void Foam::fv::interRegionHeatTransferModel::interpolate
(
    const interRegionHeatTransferModel& nbrModel,
    const Field<Type>& field,
    Field<Type>& result
) const
{
    if (master_)
    {
        meshInterp().mapTgtToSrc(field, plusEqOp<scalar>(), result);
    }
    else
    {
        nbrModel.meshInterp().mapSrcToTgt(field, plusEqOp<scalar>(), result);
    }
}


template<class Type>
void Foam::fv::interRegionHeatTransferModel::interpolate
(
    const Field<Type>& field,
    Field<Type>& result
) const
{
    return interpolate(nbrModel(), field, result);
}


// ************************************************************************* //
