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
    (c) 2017 Esi Ltd.

Class
    Foam::cellFeatures

Description
    Cell analysis class.

    Constructs feature edges and feature points, which are edges/points with
    and angle > given specification.
    Can be asked for 'superFaces' which can be used to see if a cell is a
    'splitHex'.

SourceFiles
    cellFeatures.C

\*---------------------------------------------------------------------------*/

#ifndef cellFeatures_H
#define cellFeatures_H

#include "meshes/meshShapes/face/faceList.H"
#include "primitives/ints/lists/labelList.H"
#include "primitives/bools/lists/boolList.H"
#include "containers/HashTables/HashSet/HashSet.H"
#include "containers/HashTables/Map/Map.H"
#include "containers/Lists/DynamicList/DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class primitiveMesh;

/*---------------------------------------------------------------------------*\
                           Class cellFeatures Declaration
\*---------------------------------------------------------------------------*/

class cellFeatures
{
    // Private data

        const primitiveMesh& mesh_;

        //- Cos of angle between two connected faces or two connected edges on
        //  same face before edge/point is 'feature'.
        scalar minCos_;

        label celli_;

        //- Feature edges
        labelHashSet featureEdge_;

        //- (demand driven) Faces after removing internal points&edges
        mutable faceList* facesPtr_;

        //- New to old face mapping
        mutable List<DynamicList<label>> faceMap_;


    // Private Member Functions

        bool faceAlignedEdge(const label, const label) const;

        label nextEdge
        (
            const Map<label>& toSuperFace,
            const label superFacei,
            const label thisEdgeI,
            const label thisVertI
        ) const;

        bool isCellFeatureEdge(const scalar, const label) const;

        void walkSuperFace
        (
            const label facei,
            const label superFacei,
            Map<label>& toSuperFace
        ) const;

        void calcSuperFaces() const;


        //- Disallow default bitwise copy construct
        cellFeatures(const cellFeatures&);

        //- Disallow default bitwise assignment
        void operator=(const cellFeatures&);

public:

    // Constructors

        //- Construct from cell in mesh
        cellFeatures
        (
            const primitiveMesh&,
            const scalar minCos,    // angle to use for feature recognition.
            const label celli
        );


    //- Destructor
    ~cellFeatures();


    // Member Functions

        // Access

            const labelHashSet& featureEdge() const
            {
                return featureEdge_;
            }

            const faceList& faces() const
            {
                if (!facesPtr_)
                {
                    calcSuperFaces();
                }
                return *facesPtr_;
            }

            //- New to old faceMap. Guaranteed to be shrunk.
            const List<DynamicList<label>>& faceMap() const
            {
                if (!facesPtr_)
                {
                    calcSuperFaces();
                }
                return faceMap_;
            }

            //- Return labels of super cell vertices
            labelList labels() const;


        // Check

            //- Is edge a feature edge (uniquely determined since on cell
            //  only two faces sharing edge)
            bool isFeatureEdge(const label edgeI) const
            {
                return featureEdge().found(edgeI);
            }

            //- Are two edges connected at feature point?
            //  Is local to face since point might be seen as feature point
            //  from one face but not from another.
            bool isFeaturePoint(const label edge0, const label edge1) const;

            //- Is vertexI on facei used by two edges that form feature
            //  point
            bool isFeatureVertex(const label facei, const label vertI) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
