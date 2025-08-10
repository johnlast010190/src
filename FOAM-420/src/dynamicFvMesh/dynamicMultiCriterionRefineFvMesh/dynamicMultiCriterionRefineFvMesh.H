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

Class
    Foam::dynamicMultiCriterionRefineFvMesh

Description
    A fvMesh with built-in multi-criterion refinement.

    Determines which cells to refine/unrefine and does all in update().
    Uses functionality of dynamicRefineFvMesh.

SourceFiles
    dynamicRefineFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicMultiCriterionRefineFvMesh_H
#define dynamicMultiCriterionRefineFvMesh_H

#include "dynamicRefineFvMesh/dynamicRefineFvMesh.H"
//#include "fields/volFields/volFields.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dynamicRefineFvMesh Declaration
\*---------------------------------------------------------------------------*/

class dynamicMultiCriterionRefineFvMesh
:
    public dynamicRefineFvMesh
{
private:
        //- Get access to refinement level and refinement field
        volScalarField* internalRefinementFieldPtr_;
        volScalarField* targetLevelPtr_;
        volScalarField* isLevelPtr_;

        //- Refinement based on field min/max
        HashTable< dictionary > fields_;

        //- Refinement based on mag(gradFields) min/max
        HashTable< List<scalar>> gradFields_;

        //- Refinement based on curl(field) min/max
        HashTable< List<scalar>> curlFields_;

        //- Refinement based on interface (multiphase solver)
        HashTable< dictionary > interface_;

        //- Refinement based on prescribed region
        PtrList<entry> refinedRegions_;

        //- Switch for standard refinement vs multi-criterion
        Switch enableRefinementControl_;

        //- Switch for multi-criterion refinement based on fixed target field
        Switch fixedTargetField_;

        //- Switch for target level smoothing
        Switch targetLevelSmoothing_;

        //- Switch to (dis-)allow target level growth
        Switch allowOnlyLevelGrowth_;

        //- non-dimensionalized diffusion coefficient for elliptical relaxation
        scalar diffCoeff_;

        //- smooth target refinement field only in vicinity of refinement jumps
        //  value controls how many neighbour cells are affected
        label nLayerSmooth_;

        //- (Re-)Calculate the refinement field
        void updateRefinementField();

        //- read dynamicMeshDict
        void readDict();

        //- Disallow default bitwise copy construct
        dynamicMultiCriterionRefineFvMesh(const dynamicMultiCriterionRefineFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const dynamicMultiCriterionRefineFvMesh&);

        label nBufferLayers_;
        label refineInterval_;
        label maxRefinement_;

    // virtual functions to be overwritten from base
    // to prevent consecutive refinement and unrefnement of cells
        //- Get point min of connected cell
        virtual scalarField minCellField(const volScalarField&) const;

        // overwrite protected member function of dynamicRefineFvMesh
        virtual void selectRefineCandidates
        (
            const scalar lowerRefineLevel,
            const scalar upperRefineLevel,
            const scalarField& vFld,
            PackedBoolList& candidateCell
        ) const;

    // smoothing algorithm for target field
       void smoothTargetLevel
       (
           const volScalarField& targetLevelField,
           labelList& targetLevel
       );

public:

    //- Runtime type information
    TypeName("dynamicMultiCriterionRefineFvMesh");


    // Constructors

        //- Construct from IOobject
        explicit dynamicMultiCriterionRefineFvMesh(const IOobject& io);


    //- Destructor
    virtual ~dynamicMultiCriterionRefineFvMesh();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
