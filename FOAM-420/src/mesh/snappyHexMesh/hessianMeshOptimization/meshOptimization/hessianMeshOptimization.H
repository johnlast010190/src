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
    (c) 2014 Esi Ltd.

\*---------------------------------------------------------------------------*/

#ifndef hessianMeshOptimization_H
#define hessianMeshOptimization_H

#include "fvMesh/fvMesh.H"
#include "cfdTools/general/include/fvCFD.H"
#include "hessianMeshOptimization/meshMetric/meshMetric.H"
#include "hessianMeshOptimization/meshOptimizationMethod/meshOptimizationMethod.H"
#include "hessianMeshOptimization/smoother/generalMeshSmoother.H"
#include "memory/tmp/tmp.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

namespace Foam
{

class hessianMeshOptimization

{
        // Private data

        //- Description of data_
        dictionary dict_;

        fvMesh& mesh_;

        vectorField state_;

        meshMetric metric_;

        autoPtr<meshOptimizationMethod> methodPtr_;

        label numberOfIterations_;

        generalMeshSmoother smoother_;

        label currentIteration_;

        label smoothingIterationThreshold_;

        // Private Member Functions
        void meshOptimization();

        void optimizePoint(const label&);

        void optimizeProcessorPoints();

        bool converge();

        void satisfyWolfeConditions
        (
            const label& pI,
            const scalar& initObjective,
            const vector& disp
        );

        void satisfyWolfeConditions
        (
            const labelList& pList,
            const scalarList& initObjective,
            const vectorField& disp,
            const vectorField& derivatives
        );

public:

    // Constructors
        hessianMeshOptimization
        (
            fvMesh& mesh,
            const pointField& state,
            const dictionary& dict,
            const bool local
        );

        hessianMeshOptimization
        (
            fvMesh& mesh,
            const dictionary& dict,
            const bool local
        );

        hessianMeshOptimization
        (
            fvMesh& mesh,
            const dictionary& dict
        );

        static autoPtr<hessianMeshOptimization> fullMeshOptimization
        (
            fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    ~hessianMeshOptimization()
    {}

    // Member Functions
    tmp<pointField> newPoints();

    tmp<pointField> newPoints(const labelList& pointSet);

};

}//end namespace Foam

    #endif
