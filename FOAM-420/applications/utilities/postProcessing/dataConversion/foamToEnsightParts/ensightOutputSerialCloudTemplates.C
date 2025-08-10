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
    (c) 2016 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "ensightOutputSerialCloud.H"
#include "ensight/output/ensightSerialOutput.H"
#include "ensight/type/ensightPTraits.H"

#include "passiveParticle/passiveParticle.H"
#include "db/IOobjects/IOField/IOField.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::ensightSerialCloud::writeCloudField
(
    const IOField<Type>& field,
    ensightFile& os
)
{
    // 6 values per line
    label count = 0;

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

    // add final newline if required
    if (count % 6)
    {
        os.newline();
    }

    return true;
}


template<class Type>
bool Foam::ensightSerialCloud::writeCloudField
(
    const IOobject& fieldObject,
    autoPtr<ensightFile> output
)
{
    IOField<Type> field(fieldObject);
    return writeCloudField(field, output.rawRef());
}


// ************************************************************************* //
