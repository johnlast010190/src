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
    (c) 2011 OpenFOAM Foundation
    (c) 2010-2016 Esi Ltd.

Class
    Foam::wallCorrectedLinearUpwind

Description
    wallCorrectedLinearUpwind interpolation scheme class derived from upwind and returns
    upwind weighting factors and also applies a gradient-based explicit
    correction.

SourceFiles
    wallCorrectedLinearUpwind.C

\*---------------------------------------------------------------------------*/

#ifndef wallCorrectedLinearUpwind_H
#define wallCorrectedLinearUpwind_H

#include "interpolation/surfaceInterpolation/limitedSchemes/upwind/upwind.H"
#include "finiteVolume/gradSchemes/gaussGrad/gaussGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class wallCorrectedLinearUpwind Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class wallCorrectedLinearUpwind
:
    public upwind<Type>
{
    // Private Data

        word gradSchemeName_;
        tmp<fv::gradScheme<Type>> gradScheme_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        wallCorrectedLinearUpwind(const wallCorrectedLinearUpwind&);

        //- Disallow default bitwise assignment
        void operator=(const wallCorrectedLinearUpwind&);


public:

    //- Runtime type information
    TypeName("wallCorrectedLinearUpwind");


    // Constructors

        //- Construct from faceFlux
        wallCorrectedLinearUpwind
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux
        )
        :
            upwind<Type>(mesh, faceFlux),
            gradSchemeName_("grad"),
            gradScheme_
            (
                new fv::gaussGrad<Type>(mesh)
            )
        {}

        //- Construct from Istream.
        //  The name of the flux field is read from the Istream and looked-up
        //  from the mesh objectRegistry
        wallCorrectedLinearUpwind
        (
            const fvMesh& mesh,
            const objectRegistry& db,
            Istream& schemeData
        )
        :
            upwind<Type>(mesh, db, schemeData),
            gradSchemeName_(schemeData),
            gradScheme_
            (
                fv::gradScheme<Type>::New
                (
                    mesh,
                    db,
                    mesh.schemes().gradScheme(gradSchemeName_)
                )
            )
        {
            if (!schemeData.eof())
            {
                IOWarningInFunction(schemeData)
                    << "unexpected additional entries in stream." << nl
                    << "    Only the name of the gradient scheme in the"
                       " 'gradSchemes' dictionary should be specified."
                    << endl;
            }
        }

        //- Construct from faceFlux and Istream
        wallCorrectedLinearUpwind
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& schemeData
        )
        :
            upwind<Type>(mesh, faceFlux, schemeData),
            gradSchemeName_(schemeData),
            gradScheme_
            (
                fv::gradScheme<Type>::New
                (
                    mesh,
                    faceFlux.db(),
                    mesh.schemes().gradScheme(gradSchemeName_)
                )
            )
        {

            if (!schemeData.eof())
            {
                IOWarningInFunction(schemeData)
                    << "unexpected additional entries in stream." << nl
                    << "    Only the name of the gradient scheme in the"
                       " 'gradSchemes' dictionary should be specified."
                    << endl;
            }
        }


    // Member Functions

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const
        {
            return true;
        }

        //- Return the explicit correction to the face-interpolate
        virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        correction
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
