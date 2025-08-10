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
    (c) 2016 Esi Ltd.

Class
    Foam::twistCorrectedGaussData

Description
    Twist corrected gauss gradient scheme - calculate skew correction tensors
    tensors (face area vector times skew correction vector summed over
    each triangle in triangle-decomposed face), and the cell volumes
    consistent with triangle-decomposed faces.

SourceFiles
    twistCorrectedGaussData.C

\*---------------------------------------------------------------------------*/

#ifndef twistCorrectedGaussData_H
#define twistCorrectedGaussData_H

#include "meshes/MeshObject/MeshObject.H"
#include "fvMesh/fvMesh.H"
#include "fields/surfaceFields/surfaceFieldsFwd.H"
#include "primitives/Pair/labelPair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class twistCorrectedGaussData Declaration
\*---------------------------------------------------------------------------*/

class twistCorrectedGaussData
:
    public MeshObject<fvMesh, MoveableMeshObject, twistCorrectedGaussData>
{
    // Private data

        //- Skew-correction tensors
        mutable surfaceTensorField* skewCorrTensorsPtr_;

        //- Consistent volumes
        mutable DimensionedField<scalar, volMesh>* consistentVolsPtr_;

    // Private Member Functions

        //- Construct skew correction tensors and consistent volumes
        void makeData() const;


public:

    // Declare name of the class and its debug switch
    TypeName("twistCorrectedGaussData");


    // Constructors

        //- Construct given an fvMesh
        explicit twistCorrectedGaussData(const fvMesh&);


    //- Destructor
    virtual ~twistCorrectedGaussData();


    // Member functions

        //- Return reference to effective edge surface normals
        const surfaceTensorField& skewCorrTensors() const;

        const DimensionedField<scalar, volMesh>& consistentVols() const;

        //- Delete the pseudo volume vectors when the mesh moves
        virtual bool movePoints();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
