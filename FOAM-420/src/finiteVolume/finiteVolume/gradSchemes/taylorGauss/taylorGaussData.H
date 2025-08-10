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
    Foam::taylorGaussData

Description
    TaylorGauss gradient scheme - compute inverse 'pseudo volume' tensors
    derived from Taylor expansion.

SourceFiles
    taylorGaussData.C

\*---------------------------------------------------------------------------*/

#ifndef taylorGaussData_H
#define taylorGaussData_H

#include "meshes/MeshObject/MeshObject.H"
#include "fvMesh/fvMesh.H"
#include "fields/surfaceFields/surfaceFieldsFwd.H"
#include "primitives/Pair/labelPair.H"
#include "containers/HashTables/HashPtrTable/HashPtrTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class taylorGaussData Declaration
\*---------------------------------------------------------------------------*/

class taylorGaussData
:
    public MeshObject<fvMesh, MoveableMeshObject, taylorGaussData>
{
    // Private data

        //- geometric inverse pseudo-volume tensor
        mutable tensorField* invPseudoVolPtr_;

        //- field specific inverse pseudo-volume tensors
        mutable HashPtrTable<Field<tensor>, word, string::hash>*
            fieldInvPseudoVolPtrs_;


    // Private Member Functions

        //- constrain
        void constrainSolution(tensorField&) const;

        //- Construct pseudo-volume tensors
        tmp<tensorField> pseudoVolumes() const;

        //- Calculate field specific pseudo-volume modifications
        template<class Type>
        tmp<tensorField> pseudoVolumes
        (
            GeometricField<Type, fvPatchField, volMesh>
        ) const;

public:

    // Declare name of the class and its debug switch
    TypeName("taylorGaussData");


    // Constructors

        //- Construct given an fvMesh
        explicit taylorGaussData(const fvMesh&);


    //- Destructor
    virtual ~taylorGaussData();


    // Member functions

        //- Return reference to inverse geometric pseudo volumes
        const tensorField& invPseudoVolumes() const;

        //- Return reference to inverse field plus geometric pseudo volumes
        const tensorField& invPseudoVolumes(const volScalarField&) const;

        //- Return reference to inverse  field plus geometric pseudo volumes
        const tensorField& invPseudoVolumes(const volVectorField&) const;

        //- Delete the pseudo volume vectors when the mesh moves
        virtual bool movePoints();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
