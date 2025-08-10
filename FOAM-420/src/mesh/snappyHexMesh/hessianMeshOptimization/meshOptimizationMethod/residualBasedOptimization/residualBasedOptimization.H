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

Class
    Foam::residualBasedOptimization

Description
    Optimizing the points that are actually moving during the optimization

SourceFiles
    qualityBasedOptimization.C

\*---------------------------------------------------------------------------*/

#ifndef residualBasedOptimization_H
#define qualityBasedOptimization_H

#include "fvMesh/fvMesh.H"
#include "hessianMeshOptimization/meshOptimizationMethod/meshOptimizationMethod.H"
#include "hessianMeshOptimization/meshQualityStatus/meshQualityStatus.H"

namespace Foam
{
class residualBasedOptimization
:
    public meshOptimizationMethod
{
    // Private Data

    const fvMesh& mesh_;

    const pointField& points_;

    pointField oldPoints_;

    scalarField maxDisp_;

    labelField counter_;

    label pL_;
    // Private Member Functions

    //- Disallow default bitwise copy construct
    residualBasedOptimization(const residualBasedOptimization&);

    //- Disallow default bitwise assignment
    void operator=(const residualBasedOptimization&);

    void initialActiveSet();

    void calculateActiveSet();

public    :

    //- Runtime type information
    TypeName("residualBasedOptimization");

    //Constructors
        residualBasedOptimization
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const pointField& points
        );

    // Member Functions
    bool updateActiveSet();

    bool converge();

    void merge(const labelList& pointSet);

    //Destructor
    virtual ~residualBasedOptimization();
};

} /* namespace Foam */

#endif /* residualBasedOptimization_H */
