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
    (c) 2021-2022 OpenFOAM Foundation
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * Private Static Member Functions * * * * * * * * * * //

template<class SubListA, class SubListB>
inline void Foam::patchToPatch::transferListList
(
    List<SubListA>& a,
    List<SubListB>& b
)
{
    a.resize(b.size());
    forAll(a, i)
    {
        a[i].transfer(b[i]);
    }
}


template<class Type>
inline void Foam::patchToPatch::rDistributeListList
(
    const label size,
    const mapDistribute& map,
    List<List<Type>>& data
)
{
    mapDistributeBase::distribute
    (
        Pstream::commsTypes::nonBlocking,
        List<labelPair>(),
        size,
        map.constructMap(),
        false,
        map.subMap(),
        false,
        data,
        ListAppendEqOp<Type>(),
        flipOp(),
        List<Type>()
    );
}


template<class Type>
inline void Foam::patchToPatch::rDistributeListList
(
    const label size,
    const mapDistribute& map,
    List<DynamicList<Type>>& data
)
{
    List<List<Type>> tData;
    transferListList(tData, data);
    rDistributeListList(size, map, tData);
    transferListList(data, tData);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::patchToPatch::reverse() const
{
    return reverse_;
}


inline Foam::label Foam::patchToPatch::singleProcess() const
{
    return singleProcess_;
}


inline bool Foam::patchToPatch::isSingleProcess() const
{
    return singleProcess_ != -1;
}


inline Foam::List<Foam::List<Foam::patchToPatch::procFace>>
Foam::patchToPatch::srcTgtProcFaces() const
{
    return
        isSingleProcess()
      ? localFacesToProcFaces(srcLocalTgtFaces_)
      : localFacesToProcFaces(srcLocalTgtFaces_, localTgtProcFacesPtr_());
}


inline Foam::List<Foam::List<Foam::patchToPatch::procFace>>
Foam::patchToPatch::tgtSrcProcFaces() const
{
    return
        isSingleProcess()
      ? localFacesToProcFaces(tgtLocalSrcFaces_)
      : localFacesToProcFaces(tgtLocalSrcFaces_, localSrcProcFacesPtr_());
}


// ************************************************************************* //
