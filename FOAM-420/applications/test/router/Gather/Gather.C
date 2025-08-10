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
    (c) 2011-2017 OpenFOAM Foundation

Description

\*---------------------------------------------------------------------------*/

#include "Gather.H"
#include "db/IOstreams/Pstreams/IPstream.H"
#include "db/IOstreams/Pstreams/OPstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from component
template<class T0>
Gather<T0>::Gather(const T0& localData, const bool redistribute)
:
    List<T0>(0),
    nProcs_(max(1, Pstream::nProcs()))
{
    this->setSize(nProcs_);

    //
    // Collect sizes on all processor
    //

    if (Pstream::parRun())
    {
        if (Pstream::master())
        {
            this->operator[](0) = localData;

            // Receive data
            for
            (
                int slave = Pstream::firstSlave(), procIndex = 1;
                slave <= Pstream::lastSlave();
                slave++, procIndex++
            )
            {
                IPstream fromSlave(Pstream::commsTypes::scheduled, slave);
                fromSlave >> this->operator[](procIndex);
            }

            // Send data
            for
            (
                int slave = Pstream::firstSlave(), procIndex = 1;
                slave <= Pstream::lastSlave();
                slave++, procIndex++
            )
            {
                OPstream toSlave(Pstream::commsTypes::scheduled, slave);

                if (redistribute)
                {
                    toSlave << *this;
                }
                else
                {
                    // Dummy send just to balance sends/receives
                    toSlave << 0;
                }
            }
        }
        else
        {
            // Slave: send my local data to master
            {
                OPstream toMaster
                (
                    Pstream::commsTypes::scheduled,
                    Pstream::masterNo()
                );
                toMaster << localData;
            }

            // Receive data from master
            {
                IPstream fromMaster
                (
                    Pstream::commsTypes::scheduled,
                    Pstream::masterNo()
                );
                if (redistribute)
                {
                    fromMaster >> *this;
                }
                else
                {
                    label dummy;
                    fromMaster >> dummy;
                }
            }
        }
    }
    else
    {
        this->operator[](0) = localData;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
