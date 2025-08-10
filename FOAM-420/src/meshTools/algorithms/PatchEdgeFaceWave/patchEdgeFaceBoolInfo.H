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
    (c) 2018 Esi Ltd.
    (c) 2011 OpenFOAM Foundation

Class
    Foam::patchEdgeFaceBoolInfo

Description

SourceFiles
    patchEdgeFaceBoolInfoI.H
    patchEdgeFaceBoolInfo.C

\*---------------------------------------------------------------------------*/

#ifndef patchEdgeFaceBoolInfo_H
#define patchEdgeFaceBoolInfo_H

#include "meshes/primitiveShapes/point/point.H"
#include "primitives/ints/label/label.H"
#include "primitives/Scalar/scalar/scalar.H"
#include "primitives/Tensor/tensor/tensor.H"
#include "primitives/pTraits/pTraits.H"
#include "meshes/primitiveMesh/primitivePatch/primitivePatch.H"
#include "meshes/primitiveMesh/primitivePatch/indirectPrimitivePatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyPatch;
class polyMesh;
class patchEdgeFaceBoolInfo;
Ostream& operator<<(Ostream&, const patchEdgeFaceBoolInfo&);
Istream& operator>>(Istream&, patchEdgeFaceBoolInfo&);

/*---------------------------------------------------------------------------*\
                           Class patchEdgeFaceBoolInfo Declaration
\*---------------------------------------------------------------------------*/

class patchEdgeFaceBoolInfo
{
    // Private data

        //- position of nearest wall center
        bool origin_;



public:

    // Constructors

        //- Construct null
        inline patchEdgeFaceBoolInfo();

        //- Construct from origin, distance
        inline patchEdgeFaceBoolInfo(const bool&);

    // Member Functions

        // Access

            inline const bool& origin() const;


        // Needed by meshWave

            //- Check whether origin has been changed at all or
            //  still contains original (invalid) value.
            template<class TrackingData>
            inline bool valid(TrackingData& td) const;

            //- Apply rotation matrix
            template<class TrackingData>
            inline void transform
            (
                const polyMesh& mesh,
                const indirectPrimitivePatch& patch,
                const tensor& rotTensor,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of face on edge
            template<class TrackingData>
            inline bool updateEdge
            (
                const polyMesh& mesh,
                const indirectPrimitivePatch& patch,
                const label edgeI,
                const label faceI,
                const patchEdgeFaceBoolInfo& faceInfo,
                const scalar tol,
                TrackingData& td
            );

            //- New information for edge (from e.g. coupled edge)
            template<class TrackingData>
            inline bool updateEdge
            (
                const polyMesh& mesh,
                const indirectPrimitivePatch& patch,
                const patchEdgeFaceBoolInfo& edgeInfo,
                const bool sameOrientation,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of edge on face.
            template<class TrackingData>
            inline bool updateFace
            (
                const polyMesh& mesh,
                const indirectPrimitivePatch& patch,
                const label faceI,
                const label edgeI,
                const patchEdgeFaceBoolInfo& edgeInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Same (like operator==)
            template<class TrackingData>
            inline bool equal(const patchEdgeFaceBoolInfo&, TrackingData& td) const;


    // Member Operators

        // Needed for List IO
        inline bool operator==(const patchEdgeFaceBoolInfo&) const;
        inline bool operator!=(const patchEdgeFaceBoolInfo&) const;


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const patchEdgeFaceBoolInfo&);
        friend Istream& operator>>(Istream&, patchEdgeFaceBoolInfo&);

};


//- Data associated with patchEdgeFaceBoolInfo type are contiguous
template<>
constexpr bool contiguous<patchEdgeFaceBoolInfo>()
{
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "algorithms/PatchEdgeFaceWave/patchEdgeFaceBoolInfoI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
