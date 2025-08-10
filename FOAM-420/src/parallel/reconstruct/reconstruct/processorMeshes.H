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
    (c) 2011-2014 OpenFOAM Foundation
    (c) 2016 OpenCFD Ltd.
    (c) 2022 Esi Ltd.

Class
    Foam::processorMeshes

Description
    Container for processor mesh addressing.

SourceFiles
    processorMeshes.C

\*---------------------------------------------------------------------------*/

#ifndef processorMeshes_H
#define processorMeshes_H

#include "containers/Lists/PtrList/PtrList.H"
#include "fvMesh/fvMesh.H"
#include "db/IOobjectList/IOobjectList.H"
#include "primitives/ints/lists/labelIOList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                     Class processorMeshes Declaration
\*---------------------------------------------------------------------------*/

class processorMeshes
{
    // Private data

        const word meshName_;

        //- Processor databases
        PtrList<Time>& databases_;

        //- List of processor meshes
        PtrList<fvMesh> meshes_;

        //- List of processor point addressing lists
        PtrList<labelIOList> pointProcAddressing_;

        //- List of processor face addressing lists
        PtrList<labelIOList> faceProcAddressing_;

        //- List of processor cell addressing lists
        PtrList<labelIOList> cellProcAddressing_;


    // Private Member Functions

        //- Read all meshes
        void read();

        //- Disallow default bitwise copy construct
        processorMeshes(const processorMeshes&);

        //- Disallow default bitwise assignment
        void operator=(const processorMeshes&);


public:

    //- Runtime type information
    ClassName("processorMeshes");


    // Constructors

        //- Construct from components
        processorMeshes(PtrList<Time>& databases, const word& meshName);


    // Member Functions

        //- Update the meshes based on the mesh files saved in time directories
        fvMesh::readUpdateState readUpdate();

        //- Reconstruct point position after motion in parallel
        void reconstructPoints(fvMesh&);

        const PtrList<fvMesh>& meshes() const
        {
            return meshes_;
        }

        PtrList<fvMesh>& meshes()
        {
            return meshes_;
        }

        const PtrList<labelIOList>& pointProcAddressing() const
        {
            return pointProcAddressing_;
        }

        PtrList<labelIOList>& faceProcAddressing()
        {
            return faceProcAddressing_;
        }

        const PtrList<labelIOList>& cellProcAddressing() const
        {
            return cellProcAddressing_;
        }

        //- Helper: remove all procAddressing files from mesh instance
        static void removeFiles(const polyMesh& mesh);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
