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
void Foam::mappedPatchBase::distribute(List<Type>& lst) const
{
    switch (mode_)
    {
        case NEARESIPATCHFACEAMI:
        {
            lst = AMI().interpolateToSource(Field<Type>(lst.xfer()));
            break;
        }
        default:
        {
            map().distribute(lst);
        }
    }
}


template<class Type>
void Foam::mappedPatchBase::distribute
(
    List<Type>& lst,
    const UList<Type>& defaultValues
) const
{
    switch (mode_)
    {
        case NEARESIPATCHFACEAMI:
        {
            lst =
                AMI().interpolateToSource
                (
                    Field<Type>(lst.xfer()), defaultValues
                );
            break;
        }
        default:
        {
            map().distribute(lst);
        }
    }
}


template<class Type, class CombineOp>
void Foam::mappedPatchBase::distribute
(
    List<Type>& lst,
    const CombineOp& cop
) const
{
    switch (mode_)
    {
        case NEARESIPATCHFACEAMI:
        {
            lst = AMI().interpolateToSource
                (
                    Field<Type>(lst.xfer()),
                    cop
                );
            break;
        }
        default:
        {
            mapDistributeBase::distribute
            (
                Pstream::defaultCommsType,
                map().schedule(),
                map().constructSize(),
                map().subMap(),
                false,
                map().constructMap(),
                false,
                lst,
                cop,
                flipOp(),
                Type(Zero)
            );
        }
    }
}


template<class Type>
void Foam::mappedPatchBase::reverseDistribute(List<Type>& lst) const
{
    switch (mode_)
    {
        case NEARESIPATCHFACEAMI:
        {
            lst = AMI().interpolateToTarget(Field<Type>(lst.xfer()));
            break;
        }
        default:
        {
            map().reverseDistribute(sampleSize(), lst);
            break;
        }
    }
}


template<class Type>
void Foam::mappedPatchBase::reverseDistribute
(
    List<Type>& lst,
    const UList<Type>& defaultValues
) const
{
    switch (mode_)
    {
        case NEARESIPATCHFACEAMI:
        {
            lst =
                AMI().interpolateToTarget
                (
                    Field<Type>(lst.xfer()), defaultValues
                );
            break;
        }
        default:
        {
            map().reverseDistribute(sampleSize(), lst);
        }
    }
}


template<class Type, class CombineOp>
void Foam::mappedPatchBase::reverseDistribute
(
    List<Type>& lst,
    const CombineOp& cop
) const
{
    switch (mode_)
    {
        case NEARESIPATCHFACEAMI:
        {
            lst = AMI().interpolateToTarget
                (
                    Field<Type>(lst.xfer()),
                    cop
                );
            break;
        }
        default:
        {
            label cSize = sampleSize();
            mapDistributeBase::distribute
            (
                Pstream::defaultCommsType,
                map().schedule(),
                cSize,
                map().constructMap(),
                false,
                map().subMap(),
                false,
                lst,
                cop,
                flipOp(),
                Type(Zero)
            );
            break;
        }
    }
}


// ************************************************************************* //
