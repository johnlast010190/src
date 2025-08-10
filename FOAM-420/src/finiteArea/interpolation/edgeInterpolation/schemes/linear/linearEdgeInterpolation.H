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
    linearEdgeInterpolation

Description
    Central-differencing interpolation scheme class

SourceFiles
    linearEdgeInterpolation.C

\*---------------------------------------------------------------------------*/

#ifndef linearEdgeInterpolation_H
#define linearEdgeInterpolation_H

#include "interpolation/edgeInterpolation/edgeInterpolationScheme/edgeInterpolationScheme.H"
#include "fields/areaFields/areaFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class linear Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class linearEdgeInterpolation
:
    virtual public edgeInterpolationScheme<Type>
{
    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const linearEdgeInterpolation&);

public:

    //- Runtime type information
    TypeName("linear");


    // Constructors

        //- Construct from mesh
        linearEdgeInterpolation(const faMesh& mesh)
        :
            edgeInterpolationScheme<Type>(mesh)
        {}

        //- Construct from Istream
        linearEdgeInterpolation(const faMesh& mesh, Istream&)
        :
            edgeInterpolationScheme<Type>(mesh)
        {}

        //- Construct from faceFlux and Istream
        linearEdgeInterpolation
        (
            const faMesh& mesh,
            const edgeScalarField&,
            Istream&
        )
        :
            edgeInterpolationScheme<Type>(mesh)
        {}

        // Destructor

        virtual ~linearEdgeInterpolation()
        {}

    // Member Functions

        //- Return the interpolation weighting factors
        tmp<edgeScalarField> weights
        (
            const GeometricField<Type, faPatchField, areaMesh>&
        ) const
        {
            return this->mesh().edgeInterpolation::weights();
        }
};


template<class Type>
tmp<GeometricField<Type, faePatchField, faEdgeMesh>>
linearEdgeInterpolate(const GeometricField<Type, faPatchField, areaMesh>& vf)
{
    return edgeInterpolationScheme<Type>::interpolate
    (
        vf,
        vf.mesh().edgeInterpolation::weights()
    );
}


template<class Type>
tmp<GeometricField<Type, faePatchField, faEdgeMesh>>
linearEdgeInterpolate
(
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tvf
)
{
    tmp<GeometricField<Type, faePatchField, faEdgeMesh>> tinterp =
        linearEdgeInterpolate(tvf());
    tvf.clear();
    return tinterp;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
