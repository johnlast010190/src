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
    (c) 2017 OpenCFD Ltd.

Class
    Foam::waveMethod

Description
    Direct (one-to-one cell correspondence) mesh-to-mesh interpolation class

SourceFiles
    waveMethod.C

\*---------------------------------------------------------------------------*/

#ifndef waveMethod_H
#define waveMethod_H

#include "meshToMesh/calcMethod/meshToMeshMethod/meshToMeshMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class waveMethod Declaration
\*---------------------------------------------------------------------------*/

class waveMethod
:
    public meshToMeshMethod
{
    // Private member functions

        //- Disallow default bitwise copy construct
        waveMethod(const waveMethod&);

        //- Disallow default bitwise assignment
        void operator=(const waveMethod&);


public:

    //- Run-time type information
    TypeName("wave");

    //- Construct from source and target meshes
    waveMethod
    (
        const polyMesh& src,
        const polyMesh& tgt,
        const dictionary& dict
    );

    //- Destructor
    virtual ~waveMethod();


    // Member Functions

        // Evaluate

            //- Calculate addressing
            static void calculate
            (
                const polyMesh& src,
                const polyMesh& tgt,
                labelList& srcToTgtAddr
            );

            //- Calculate addressing and weights and optionally offset vectors
            virtual void calculate
            (
                labelListList& srcToTgtAddr,
                scalarListList& srcToTgtWght,
                pointListList& srcToTgtVec,
                labelListList& tgtToSrcAddr,
                scalarListList& tgtToSrcWght,
                pointListList& tgtToSrcVec
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
