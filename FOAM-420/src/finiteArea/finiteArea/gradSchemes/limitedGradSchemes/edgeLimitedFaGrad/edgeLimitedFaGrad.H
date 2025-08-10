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
    (c) 1991-2005 OpenCFD Ltd.
    (c) 2016 Esi Ltd.

Class
    Foam::fa::edgeLimitedGrad

Description
    edgeLimitedGrad gradient scheme applied to a runTime selected base gradient
    scheme.

    The scalar limiter based on limiting the extrapolated edge values
    between the edge-neighbour cell values and is applied to all components
    of the gradient.

SourceFiles
    edgeLimitedFaGrads.C

\*---------------------------------------------------------------------------*/

#ifndef edgeLimitedGrad_H
#define edgeLimitedGrad_H

#include "finiteArea/gradSchemes/faGradScheme/faGradScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fa
{

/*---------------------------------------------------------------------------*\
                       Class edgeLimitedGrad Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class edgeLimitedGrad
:
    public fa::gradScheme<Type>
{
    // Private Data

        tmp<fa::gradScheme<Type>> basicGradScheme_;

        //- Limiter coefficient
        scalar k_;


    // Private Member Functions

        inline void limitEdge
        (
            scalar& limiter,
            const scalar maxDelta,
            const scalar minDelta,
            const scalar extrapolate
        ) const;


        //- Disallow default bitwise copy construct
        edgeLimitedGrad(const edgeLimitedGrad&);

        //- Disallow default bitwise assignment
        void operator=(const edgeLimitedGrad&);


public:

    //- RunTime type information
    TypeName("edgeLimited");


    // Constructors

        //- Construct from mesh and schemeData
        edgeLimitedGrad(const faMesh& mesh, Istream& schemeData)
        :
            gradScheme<Type>(mesh),
            basicGradScheme_(fa::gradScheme<Type>::New(mesh, schemeData)),
            k_(readScalar(schemeData))
        {
            if (k_ < 0 || k_ > 1)
            {
                FatalIOErrorInFunction
                (
                    schemeData
                )   << "coefficient = " << k_
                    << " should be >= 0 and <= 1"
                    << exit(FatalIOError);
            }
        }

        //Destructor
        virtual ~edgeLimitedGrad()
        {}

    // Member Functions

        tmp
        <
            GeometricField
            <typename outerProduct<vector, Type>::type, faPatchField, areaMesh>
        > grad
        (
            const GeometricField<Type, faPatchField, areaMesh>&
        ) const;

};

// * * * * * * * * * * * * Inline Member Function  * * * * * * * * * * * * * //

template<class Type>
inline void edgeLimitedGrad<Type>::limitEdge
(
    scalar& limiter,
    const scalar maxDelta,
    const scalar minDelta,
    const scalar extrapolate
) const
{
    if (extrapolate > maxDelta + VSMALL)
    {
        limiter = min(limiter, maxDelta/extrapolate);
    }
    else if (extrapolate < minDelta - VSMALL)
    {
        limiter = min(limiter, minDelta/extrapolate);
    }
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
