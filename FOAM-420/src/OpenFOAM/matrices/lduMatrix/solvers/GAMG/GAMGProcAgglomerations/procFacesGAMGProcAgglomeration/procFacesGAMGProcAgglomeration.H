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
    (c) 2013 OpenFOAM Foundation

Class
    Foam::procFacesGAMGProcAgglomeration

Description
    Processor agglomeration of GAMGAgglomerations. Needs nAgglomeratingCells
    which is when to start agglomerating processors. Processors get agglomerated
    by constructing a single cell mesh for each processor with each
    processor interface a face. This then gets agglomerated using
    the pairGAMGAgglomeration algorithm with the number of faces
    on the original processor interface as face weight.

SourceFiles
    procFacesGAMGProcAgglomeration.C

\*---------------------------------------------------------------------------*/

#ifndef procFacesGAMGProcAgglomeration_H
#define procFacesGAMGProcAgglomeration_H

#include "matrices/lduMatrix/solvers/GAMG/GAMGProcAgglomerations/GAMGProcAgglomeration/GAMGProcAgglomeration.H"
#include "containers/Lists/DynamicList/DynamicList.H"
#include "fields/Fields/labelField/labelField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class GAMGAgglomeration;
class lduMesh;
class lduPrimitiveMesh;

/*---------------------------------------------------------------------------*\
                    Class procFacesGAMGProcAgglomeration Declaration
\*---------------------------------------------------------------------------*/

class procFacesGAMGProcAgglomeration
:
    public GAMGProcAgglomeration
{
    // Private data

        //- When to processor agglomerate
        const label nAgglomeratingCells_;

        //- Allocated communicators
        DynamicList<label> comms_;

    // Private Member Functions

        //- Return (on master) all single-cell meshes collected. single-cell
        //  meshes are just one cell with all proc faces intact.
        autoPtr<lduPrimitiveMesh> singleCellMesh
        (
            const label singleCellMeshComm,
            const lduMesh& mesh,
            scalarField& faceWeights
        ) const;

        //- Construct processor agglomeration: for every processor the
        //  coarse processor-cluster it agglomerates onto
        tmp<labelField> processorAgglomeration(const lduMesh&) const;

        //- Do we need to agglomerate across processors?
        bool doProcessorAgglomeration(const lduMesh&) const;

        //- Disallow default bitwise copy construct
        procFacesGAMGProcAgglomeration
        (
            const procFacesGAMGProcAgglomeration&
        );

        //- Disallow default bitwise assignment
        void operator=(const procFacesGAMGProcAgglomeration&);


public:

    //- Runtime type information
    TypeName("procFaces");


    // Constructors

        //- Construct given agglomerator and controls
        procFacesGAMGProcAgglomeration
        (
            GAMGAgglomeration& agglom,
            const dictionary& controlDict
        );


    //- Destructor
    virtual ~procFacesGAMGProcAgglomeration();


    // Member Functions

       //- Modify agglomeration. Return true if modified
        virtual bool agglomerate();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
