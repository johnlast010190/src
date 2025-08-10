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
    (c) Creative Fields, Ltd.

Authors
    Franjo Juretic (franjo.juretic@c-fields.com)

Description

\*---------------------------------------------------------------------------*/

#include "include/demandDrivenData.H"
#include "utilities/smoothers/topology/topologicalCleaner/topologicalCleaner.H"
#include "utilities/faceDecomposition/decomposeFaces.H"

//#define DEBUGCleaner

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void topologicalCleaner::checkNonConsecutiveBoundaryVertices()
{
    Info<< "Checking for invalid face connections" << endl;

    const faceListPMG& faces = mesh_.faces();
    const labelList& owner = mesh_.owner();
    const cellListPMG& cells = mesh_.cells();

    const label nIntFaces = mesh_.nInternalFaces();

    boolList decomposeFace(faces.size(), false);
    bool changed(false);

    const PtrList<boundaryPatch>& boundaries = mesh_.boundaries();
    forAll(boundaries, patchI)
    {
        const label start = boundaries[patchI].patchStart();
        const label end = start + boundaries[patchI].patchSize();

        for (label faceI=start;faceI<end;++faceI)
        {
            # ifdef USE_OMP
            # pragma omp task shared(decomposeFace,faces,cells,owner,faceI)
            # endif
            {
                const face& bf = faces[faceI];

                # ifdef DEBUGCleaner
                Info<< "Checking boundary face " << faceI << " with vertices "
                    << bf << endl;
                # endif

                const cell& c = cells[owner[faceI]];

                forAll(c, fI)
                    if (
                        (c[fI] < nIntFaces) ||
                        (mesh_.faceIsInProcPatch(c[fI]) != -1)
                    )
                    {
                        const face& f = faces[c[fI]];

                        DynList<label> shN;

                        forAll(bf, pI)
                            forAll(f, pJ)
                                if (bf[pI] == f[pJ])
                                {
                                    shN.append(pI);
                                }

                        # ifdef DEBUGCleaner
                        Info<< "Shared vertices with internal face " << f
                            << " are " << shN << endl;
                        # endif

                        if (shN.size() > 2)
                        {
                            # ifdef DEBUGCleaner
                            Info<< "1. Face has to be split" << endl;
                            # endif

                            decomposeFace[faceI] = true;
                            decomposeCell_[owner[faceI]] = true;
                            changed = true;
                        }
                        else if (shN.size() == 2)
                        {
                            if (!(
                                    (shN[0] == bf.fcIndex(shN[1])) ||
                                    (shN[0] == bf.rcIndex(shN[1]))
                                )
                            )
                            {
                                # ifdef DEBUGCleaner
                                Info<< "2. Face has to be split" << endl;
                                # endif

                                decomposeFace[faceI] = true;
                                decomposeCell_[owner[faceI]] = true;
                                changed = true;
                            }
                        }
                    }
            }
        }
    }

    if (Pstream::parRun())
        reduce(changed, maxOp<bool>());

    if (changed)
    {
        changed_ = true;
        decomposeFaces(mesh_).decomposeMeshFaces(decomposeFace);
    }

    Info<< "Finished checking for invalid face connections" << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
