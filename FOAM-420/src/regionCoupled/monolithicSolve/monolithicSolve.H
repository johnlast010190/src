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
    (c) 2015-2020 Esi Ltd.
    (c) 2011-2013 OpenFOAM Foundation

Class
    monolithicSolve

Description
    Solve a region-coupled collection of matrices as a monolithic system.

SourceFiles
    monolithicSolve.C
    monolithicSolveTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef monolithicSolve_H
#define monolithicSolve_H

#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "matrices/lduMatrix/lduMatrix/lduMatrix.H"
#include "fvMatrices/fvMatrices.H"
#include "memory/tmp/tmp.H"
#include "memory/autoPtr/autoPtr.H"
#include "dimensionedTypes/dimensionedTypes.H"
#include "primitives/zero/zero.H"
#include "db/typeInfo/className.H"

#include "matrices/LduMatrix/LduMatrix/LduMatrix.H"
#include "fields/Fields/diagTensorField/diagTensorField.H"

#include "offsetLduInterface/offsetLduInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                           Function Declarations
\*---------------------------------------------------------------------------*/

template<class Type>
class monolithicSolve
{

protected:

    static bool getCoupledMatrixAndOffset
    (
        const lduInterface& interface,
        const UPtrList<fvMatrix<Type>>& matrices,
        const label matrixI,
        const label patchI,
        label& foreignMatrix, // Warns and returns -1 if not found
        label& foreignOffset
    );

    static void makeSuperMatrix
    (
        autoPtr<lduMesh>& superMesh,
        autoPtr<lduMatrix>& superMatrix,
        UPtrList<fvMatrix<Type>>& matrices,
        Time& obr // object registry to register new mesh
    );

    // Update necessary MeshObjects for a mesh move.
    static void superMeshMoved(objectRegistry& obr);

    // Update necessary MeshObjects for a mesh topo change.
    static void superMeshTopoChanged(objectRegistry& obr);

public:

    static solverPerformance solve
    (
        UPtrList<fvMatrix<Type>>& matrices,
        const dictionary& solverControls,
        const bool meshMoved, // has there been a mesh move since last solve?
        const bool meshTopoChanged // topo change since last solve?
    );

    static word makeSuperName(const stringList& names);

    static solverPerformance solve
    (
        UPtrList<fvMatrix<Type>>& matrices,
        const solution& solutionControls,
        const bool finalIter,
        const bool meshMoved, // has there been a mesh move since last solve?
        const bool meshTopoChanged // topo change since last solve?
    );

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "monolithicSolveTemplates.C"
#endif

namespace Foam
{

// Specialisation for scalars
template<>
solverPerformance monolithicSolve<scalar>::solve
(
    UPtrList<fvMatrix<scalar>>& matrices,
    const dictionary& solverControls,
    const bool meshMoved,
    const bool meshTopoChanged
);

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
