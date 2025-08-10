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
    (c) 2012 Esi Ltd.

\*---------------------------------------------------------------------------*/


#include "switchParallel.H"

namespace Foam
{
    switchParallel::switchParallel(bool newValue)
    {
        oldValue = Pstream::parRun();
        if (oldValue)
        {
            // Wait for all threads to get here.
            Pstream::barrier();
        }
        Pstream::parRun()=newValue;
    }

    switchParallel::~switchParallel()
    {
        Pstream::parRun()=oldValue;
        if (Pstream::parRun())
        {
            // Wait for all threads to get here.
            Pstream::barrier();
        }
    }

    void switchParallel::switchParallelBack()
    {
        Pstream::parRun()=oldValue;
        if (Pstream::parRun())
        {
            // Wait for all threads to get here.
            Pstream::barrier();
        };
    }
}
