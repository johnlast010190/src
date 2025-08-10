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

#ifndef meshMetric_H
#define meshMetric_H

#include "hessianMeshOptimization/meshGeometry/meshGeometry.H"
#include "primitives/Tensor/modifiedNewtonsMethod/modifiedNewtonsMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

namespace Foam
{
class meshMetric : public meshGeometry
{
    // Private data

            //- Description of data_

            const fvMesh& mesh_;

            scalarField spher_;

            scalarField nonSphericity_;

            scalar objective_;

            vector derivative_;

            tensor hessian_;

            label minSphericityCell_;

            label pow_;

            void sphericityCalculation(const label& pI);

            scalar sphericityValue(const label& cL);

            vector sphericityDer
            (
                const scalar& Vol, const scalar& Surf,
                const vector& d_Vol, const vector& d_Surf
            );

            tensor sphericity2ndDer
            (
                const scalar& Vol, const scalar& Surf,
                const vector& d_Vol, const vector& d_Surf,
                const tensor& d_2Sur
            );

            void nonSphericityHessian
            (
                tensor& t,
                const vector& der,
                const label& cI
            );

            meshMetric(const meshMetric&);

            void operator=(const meshMetric&);

public:

    meshMetric
    (
        const fvMesh& mesh,
        pointField& state,
        const dictionary& dict
    );

    ~meshMetric();

     // Member Functions

    //-Calculate and return the sphericity based point value
    scalar pointValue(const label& pI);

    //-Calculate synchronize and return the
    //-sphericity based point value for processor points
    tmp<scalarField> pointValueGlobal(const labelList& pList);

    //-Calculate first and second gradients
    void calculateGradients(const label& pI);

    //-Calculate and return the displacement based on the the gradients
    vector displacement(const label& pI);

    //-Update the geometry after point movement
    void update(const label& pI);

    //-Update all the shared(processor points)
    void updateGlobal(const labelList& pList);

    //-Point which is currently being optimize
    const label& currentLabel() const;

    //-Get point's derivative
    const vector& getDerivative(const label& pI) const;

    //-Get point's hessian matrix
    const tensor& getHessian(const label& pI) const;

    void printInfo();


};

} /* namespace Foam */

#endif /* meshMetric_H*/
