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
    (c) Esi Ltd

InNamespace
    Foam::explicitScalarWaveSolve

Description
    Perform a 1st order fixed direction convection sweep of a variable using the
    provided cell order

SourceFiles
    .C

\*---------------------------------------------------------------------------*/

#ifndef explicitScalarWaveSolve_H
#define explicitScalarWaveSolve_H

#include "dimensionedTypes/dimensionedScalar/dimensionedScalar.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFieldsFwd.H"
#include "primitives/Tuple2/Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //


// waveSolverCoeffs  class
class waveSolverCoeffs
:
    public labelList //order of solution
{
    protected:

        //- diagonal contribution
        scalarField A_;

        //- off-diagonal coefficient and neighbour address
        List<List <Tuple2 <scalar, label>>> HCoeffs_;

        //- boundary contributions
        volScalarField::Boundary Jppos_;

    public:
    // Constructors

        //- from components
        waveSolverCoeffs
        (
            const fvMesh& mesh,
            const vector& waveDirection
        );


    //- perform fixed direction wave propagation.
    friend scalar explicitScalarWaveSolve
    (
        volScalarField& vsf,
        surfaceScalarField& ssf,
        const surfaceScalarField& flux,
        const List<label>& order
    );

    //- perform fixed direction wave propagation with precalculated coefficients.
    friend scalar explicitScalarWaveSolve
    (
        volScalarField& vsf,
        const waveSolverCoeffs& wsc
    );

};

//- perform fixed direction wave propagation.
scalar explicitScalarWaveSolve
(
    volScalarField& vsf,
    surfaceScalarField& ssf,
    const surfaceScalarField& flux,
    const List<label>& order
);

//- perform fixed direction wave propagation with precalculated coefficients.
scalar explicitScalarWaveSolve
(
    volScalarField& vsf,
    const waveSolverCoeffs& wsc
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
