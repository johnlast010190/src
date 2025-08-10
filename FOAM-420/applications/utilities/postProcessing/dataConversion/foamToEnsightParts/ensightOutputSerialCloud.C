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
#include "ensight/type/ensightPTraits.H"

#include "passiveParticle/passiveParticle.H"
#include "db/IOobjects/IOField/IOField.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

void Foam::ensightSerialCloud::writePositions
(
    const polyMesh& mesh,
    const word& cloudName,
    autoPtr<ensightFile> output
)
{
    label nTotParcels = 0;
    autoPtr<Cloud<passiveParticle>> cloudPtr;

    cloudPtr.reset(new Cloud<passiveParticle>(mesh, cloudName, false));
    nTotParcels = cloudPtr().size();

    Cloud<passiveParticle> parcels(mesh, cloudName, false);

    if (Pstream::master())
    {
        ensightFile& os = output();
        os.beginParticleCoordinates(nTotParcels);

        // binary write is Ensight6 - first ids, then positions
        if (os.format() == IOstream::BINARY)
        {
            // 1-index
            for (label parcelId = 0; parcelId < nTotParcels; ++parcelId)
            {
                os.write(parcelId+1);
            }


            forAllConstIter(Cloud<passiveParticle>, cloudPtr(), elmnt)
            {
                const vector& p = elmnt().position();

                os.write(p.x());
                os.write(p.y());
                os.write(p.z());
            }
        }
        else
        {
            // ASCII id + position together

            label parcelId = 0;
            forAllConstIter(Cloud<passiveParticle>, cloudPtr(), elmnt)
            {
                const vector& p = elmnt().position();

                os.write(++parcelId, 8); // unusual width
                os.write(p.x());
                os.write(p.y());
                os.write(p.z());
                os.newline();
            }
        }
    }
}


// ************************************************************************* //
