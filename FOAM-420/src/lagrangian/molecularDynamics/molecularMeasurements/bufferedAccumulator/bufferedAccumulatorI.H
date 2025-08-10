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
    (c) 2011 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
inline Field<Type>& bufferedAccumulator<Type>::accumulationBuffer()
{
    return (*this)[nBuffers()];
}


template<class Type>
inline const Field<Type>& bufferedAccumulator<Type>::accumulationBuffer() const
{
    return (*this)[nBuffers()];
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
inline label bufferedAccumulator<Type>::averagesTaken() const
{
    return averagesTaken_;
}


template<class Type>
inline label bufferedAccumulator<Type>::nBuffers() const
{
    return bufferOffsets_.size();
}


template<class Type>
inline label bufferedAccumulator<Type>::bufferLength() const
{
    return (*this)[0].size();
}


template<class Type>
inline const List<label>& bufferedAccumulator<Type>::bufferOffsets() const
{
    return bufferOffsets_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
