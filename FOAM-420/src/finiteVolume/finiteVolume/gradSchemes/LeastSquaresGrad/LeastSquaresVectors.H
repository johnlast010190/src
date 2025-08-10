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
    (c) 2013-2016 OpenFOAM Foundation

Class
    Foam::fv::LeastSquaresVectors

Description
    Least-squares gradient scheme vectors

See also
    Foam::fv::LeastSquaresGrad

SourceFiles
    LeastSquaresVectors.C

\*---------------------------------------------------------------------------*/

#ifndef LeastSquaresVectors_H
#define LeastSquaresVectors_H

#include "fvMesh/extendedStencil/cellToCell/extendedCentredCellToCellStencil.H"
#include "meshes/MeshObject/MeshObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                     Class LeastSquaresVectors Declaration
\*---------------------------------------------------------------------------*/

template<class Stencil>
class LeastSquaresVectors
:
    public MeshObject<fvMesh, MoveableMeshObject, LeastSquaresVectors<Stencil>>
{
    // Private data

        //- Least-squares gradient vectors
        List<List<vector>> vectors_;


    // Private Member Functions

        //- Calculate Least-squares gradient vectors
        void calcLeastSquaresVectors();


public:

    // Declare name of the class and its debug switch
    TypeName("LeastSquaresVectors");


    // Constructors

        //- Construct given an fvMesh and the minimum determinant criterion
        LeastSquaresVectors
        (
            const fvMesh&
        );


    //- Destructor
    virtual ~LeastSquaresVectors();


    // Member functions

        //- Return reference to the stencil
        const extendedCentredCellToCellStencil& stencil() const
        {
            return Stencil::New(this->mesh_);
        }

        //- Return reference to the least square vectors
        const List<List<vector>>& vectors() const
        {
            return vectors_;
        }

        //- Update the least square vectors when the mesh moves
        virtual bool movePoints();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "finiteVolume/gradSchemes/LeastSquaresGrad/LeastSquaresVectors.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
