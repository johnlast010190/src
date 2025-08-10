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
    (c) 2011-2013 OpenFOAM Foundation

InClass
    Foam

Description
    Combination-Reduction operation for a parallel run.  The
    information from all nodes is collected on the master node,
    combined using the given combination function and the result is
    broadcast to all nodes


\*---------------------------------------------------------------------------*/

#ifndef PstreamCombineReduceOps_H
#define PstreamCombineReduceOps_H

#include "db/IOstreams/Pstreams/UPstream.H"
#include "db/IOstreams/Pstreams/Pstream.H"
#include "primitives/ops/ops.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class T, class CombineOp>
void combineReduce
(
    const List<UPstream::commsStruct>& comms,
    T& Value,
    const CombineOp& cop,
    const int tag,
    const label comm
)
{
    Pstream::combineGather(comms, Value, cop, tag, comm);
    Pstream::combineScatter(comms, Value, tag, comm);
}


template<class T, class CombineOp>
void combineReduce
(
    T& Value,
    const CombineOp& cop,
    const int tag = Pstream::msgType(),
    const label comm = Pstream::worldComm
)
{
    if (UPstream::nProcs(comm) < UPstream::nProcsSimpleSum)
    {
        Pstream::combineGather
        (
            UPstream::linearCommunication(comm),
            Value,
            cop,
            tag,
            comm
        );
        Pstream::combineScatter
        (
            UPstream::linearCommunication(comm),
            Value,
            tag,
            comm
        );
    }
    else
    {
        Pstream::combineGather
        (
            UPstream::treeCommunication(comm),
            Value,
            cop,
            tag,
            comm
        );
        Pstream::combineScatter
        (
            UPstream::treeCommunication(comm),
            Value,
            tag,
            comm
        );
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
