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
    (c) 2016-2017 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "ensightOutputCloud.H"
#include "ensight/type/ensightPTraits.H"

#include "db/IOobjects/IOField/IOField.H"
#include "db/Time/Time.H"
#include "meshes/polyMesh/globalMeshData/globalIndex.H"


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::ensightCloud::writeCloudField
(
    const Foam::IOField<Type>& field,
    Foam::ensightFile& os
)
{
    const bool exists = (returnReduce(field.size(), sumOp<label>()) > 0);

    if (exists)
    {
        if (Pstream::master())
        {
            // 6 values per line
            label count = 0;

            // Master
            forAll(field, i)
            {
                Type val = field[i];

                if (mag(val) < 1e-90)
                {
                    val = Zero;
                }

                for (direction d=0; d < pTraits<Type>::nComponents; ++d)
                {
                    label cmpt = ensightPTraits<Type>::componentOrder[d];
                    os.write(component(val, cmpt));

                    if (++count % 6 == 0)
                    {
                        os.newline();
                    }
                }
            }

            // Slaves
            for (int slave=1; slave<Pstream::nProcs(); ++slave)
            {
                IPstream fromSlave(Pstream::commsTypes::scheduled, slave);
                Field<Type> slaveData(fromSlave);

                forAll(slaveData, i)
                {
                    Type val = slaveData[i];

                    if (mag(val) < 1e-90)
                    {
                        val = Zero;
                    }

                    for (direction d=0; d < pTraits<Type>::nComponents; ++d)
                    {
                        label cmpt = ensightPTraits<Type>::componentOrder[d];
                        os.write(component(val, cmpt));

                        if (++count % 6 == 0)
                        {
                            os.newline();
                        }
                    }
                }
            }

            // add final newline if required
            if (count % 6)
            {
                os.newline();
            }
        }
        else
        {
            OPstream toMaster
            (
                Pstream::commsTypes::scheduled,
                Pstream::masterNo()
            );

            toMaster
                << field;
        }
    }

    return exists;
}


template<class Type>
bool Foam::ensightCloud::writeCloudField
(
    const Foam::IOobject& fieldObject,
    const bool exists,
    Foam::autoPtr<Foam::ensightFile>& output
)
{
    if (exists)
    {
        IOField<Type> field(fieldObject);
        writeCloudField(field, output.rawRef());
    }

    return true;
}


// ************************************************************************* //
