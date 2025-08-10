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
    (c) 2016-2017 Wikki Ltd

Class
    Foam::boundedBackwardFaDdtScheme

Description
    Second-order backward-differencing implicit ddt using the current and
    two previous time-step values.

SourceFiles
    boundedBackwardFaDdtScheme.C

\*---------------------------------------------------------------------------*/

#ifndef boundedBackwardFaDdtScheme_H
#define boundedBackwardFaDdtScheme_H

#include "finiteArea/ddtSchemes/faDdtScheme/faDdtScheme.H"
#include "faMatrices/faMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fa
{

/*---------------------------------------------------------------------------*\
                       Class boundedBackwardFaDdtScheme Declaration
\*---------------------------------------------------------------------------*/

class boundedBackwardFaDdtScheme
:
    public fa::faDdtScheme<scalar>
{
    // Private Member Functions

        //- Return the current time-step
        scalar deltaT_() const;

        //- Return the previous time-step
        scalar deltaT0_() const;

        //- Return the previous time-step or GREAT if the old timestep field
        //  wasn't available in which case Euler ddt is used
        template<class GeoField>
        scalar deltaT0_(const GeoField& vf) const
        {
            if (vf.oldTime().timeIndex() == vf.oldTime().oldTime().timeIndex())
            {
                return GREAT;
            }
            else
            {
                return deltaT0_();
            }
        }


        //- Disallow default bitwise copy construct
        boundedBackwardFaDdtScheme(const boundedBackwardFaDdtScheme&);

        //- Disallow default bitwise assignment
        void operator=(const boundedBackwardFaDdtScheme&);


public:

    //- Runtime type information
    TypeName("boundedBackward");


    // Constructors

        //- Construct from mesh
        boundedBackwardFaDdtScheme(const faMesh& mesh)
        :
            faDdtScheme<scalar>(mesh)
        {}

        //- Construct from mesh and Istream
        boundedBackwardFaDdtScheme(const faMesh& mesh, Istream& is)
        :
            faDdtScheme<scalar>(mesh, is)
        {}


    // Member Functions

        //- Return mesh reference
        const faMesh& mesh() const
        {
            return fa::faDdtScheme<scalar>::mesh();
        }

        tmp<areaScalarField> facDdt
        (
            const dimensionedScalar
        );

        tmp<areaScalarField> facDdt0
        (
            const dimensionedScalar
        );

        tmp<areaScalarField> facDdt
        (
            const areaScalarField&
        );

        tmp<areaScalarField> facDdt0
        (
            const areaScalarField&
        );

        tmp<areaScalarField> facDdt
        (
            const dimensionedScalar&,
            const areaScalarField&
        );

        tmp<areaScalarField> facDdt0
        (
            const dimensionedScalar&,
            const areaScalarField&
        );

        tmp<areaScalarField> facDdt
        (
            const areaScalarField&,
            const areaScalarField&
        );

        tmp<areaScalarField> facDdt0
        (
            const areaScalarField&,
            const areaScalarField&
        );

        tmp<faScalarMatrix> famDdt
        (
            const areaScalarField&
        );

        tmp<faScalarMatrix> famDdt
        (
            const dimensionedScalar&,
            const areaScalarField&
        );

        tmp<faScalarMatrix> famDdt
        (
            const areaScalarField&,
            const areaScalarField&
        );

        tmp<edgeScalarField> meshPhi
        (
            const areaScalarField&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
