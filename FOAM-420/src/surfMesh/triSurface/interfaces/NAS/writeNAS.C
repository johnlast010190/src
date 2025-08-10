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
    (c) 2010-2011, Esi Ltd

\*---------------------------------------------------------------------------*/

#include "triSurface/triSurface.H"
#include "meshes/primitiveMesh/primitivePatch/primitivePatch.H"
#include "db/IOstreams/IOstreams/IOmanip.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void triSurface::writeNASASCII(Ostream& os) const
{
    //Write out header
    os << "$" <<endl;
    os << "$ OpenFOAM generated Nastran file" <<endl;
    os << "$" <<endl;

    labelList faceMap;

    surfacePatchList myPatches(calcPatches(faceMap));

    label faceIndex = 0;

    const pointField& ps = points();


    os.setf(ios_base::scientific);
    // Write vertex coords
    forAll(ps, pointi)
    {
        os  << setw(8) << "GRID*   "
            << setw(16) << pointi
            << setw(16) << "        "
            << setw(16) << ps[pointi].x()
            << setw(16) << ps[pointi].y()
            << endl;
        os  << setw(8) << "*       "  << setw(16) << ps[pointi].z()
            << endl;
    }

    forAll(myPatches, patchI)
    {
        // Print all faces belonging to this region
        const surfacePatch& patch = myPatches[patchI];
        for
        (
            label patchFaceI = 0;
            patchFaceI < patch.size();
            patchFaceI++
        )
        {
            const label faceI = faceMap[faceIndex++];
            const labelledTri& f = (*this)[faceI];

            os  << setw(8) << "CTRIA3  "
                << setw(8) << faceIndex
                << setw(8) << patchI
                << setw(8) << f[0]
                << setw(8) << f[1]
                << setw(8) << f[2]
                << endl;
        }
    }

    forAll(myPatches, patchI)
    {
        // Print all faces belonging to this region
        const surfacePatch& patch = myPatches[patchI];

        os << "$ANSA_NAME;" << patchI << ";PSHELL;~" << endl;
        os << "$" << patch.name() << endl;
    }

}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
