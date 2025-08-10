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
    (c) 2017 Wikki Ltd.
    (c) 2016 Esi Ltd.

Class
    Foam::fa::limitedLnGrad

Description
    Central-difference lnGrad scheme with limited non-orthogonal correction.

    The limiter is controlled by a coefficient with a value between 0 and 1
    which when 0 switches the correction off and the scheme behaves as
    uncorrectedSnGrad, when set to 1 the full correction is applied and the
    scheme behaves as correctedSnGrad and when set to 0.5 the limiter is
    calculated such that the non-orthogonal contribution does not exceed the
    orthogonal part.

SourceFiles
    limitedLnGrad.C

\*---------------------------------------------------------------------------*/

#ifndef limitedLnGrad_H
#define limitedLnGrad_H

#include "finiteArea/lnGradSchemes/lnGradScheme/lnGradScheme.H"
#include "finiteArea/lnGradSchemes/correctedLnGrad/correctedLnGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fa
{

/*---------------------------------------------------------------------------*\
                           Class limitedLnGrad Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class limitedLnGrad
:
    public lnGradScheme<Type>
{
    // Private data

        tmp<lnGradScheme<Type>> correctedScheme_;

        //- Limiter.  0 = no limiting, 1 = full limiting
        scalar limitCoeff_;


    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const limitedLnGrad&);

        //- Lookup function for the corrected to support backward compatibility
        //  of dictionary specification
        tmp<lnGradScheme<Type>> lookupCorrectedScheme(Istream& schemeData)
        {
            token nextToken(schemeData);

            if (nextToken.isNumber())
            {
                limitCoeff_ = nextToken.number();
                return tmp<lnGradScheme<Type>>
                (
                    new correctedLnGrad<Type>(this->mesh())
                );
            }
            else
            {
                schemeData.putBack(nextToken);
                tmp<lnGradScheme<Type>> tcorrectedScheme
                (
                    fa::lnGradScheme<Type>::New(this->mesh(), schemeData)
                );

                schemeData >> limitCoeff_;

                return tcorrectedScheme;
            }
        }


public:

    //- Runtime type information
    TypeName("limited");


    // Constructors

        //- Construct from mesh
        limitedLnGrad(const faMesh& mesh)
        :
            lnGradScheme<Type>(mesh),
            correctedScheme_(new correctedLnGrad<Type>(this->mesh())),
            limitCoeff_(1)
        {}


        //- Construct from mesh and data stream
        limitedLnGrad(const faMesh& mesh, Istream& schemeData)
        :
            lnGradScheme<Type>(mesh),
            correctedScheme_(lookupCorrectedScheme(schemeData))
        {
            if (limitCoeff_ < 0 || limitCoeff_ > 1)
            {
                FatalIOErrorInFunction
                (
                    schemeData
                )   << "limitCoeff is specified as " << limitCoeff_
                    << " but should be >= 0 && <= 1"
                    << exit(FatalIOError);
            }
        }

    // Destructor
        virtual ~limitedLnGrad();


    // Member Functions

        //- Return the interpolation weighting factors for the given field
        virtual tmp<edgeScalarField> deltaCoeffs
        (
            const GeometricField<Type, faPatchField, areaMesh>&
        ) const
        {
            return this->mesh().nonOrthDeltaCoeffs();
        }

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const
        {
            return true;
        }

        //- Return the explicit correction to the limitedLnGrad
        //  for the given field
        virtual tmp<GeometricField<Type, faePatchField, faEdgeMesh>>
        correction(const GeometricField<Type, faPatchField, areaMesh>&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "finiteArea/lnGradSchemes/limitedLnGrad/limitedLnGrad.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
