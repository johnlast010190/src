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
    (c) 2012-2016 OpenFOAM Foundation
    (c) 2022 Esi Ltd.

Class
    Foam::fv::smartBoundedDdtScheme

Group
    grpFvDdtSchemes

Description
    Bounded form of the selected ddt scheme.
    Universal bounded scheme for incompressible and compressible temporal
    schemes. "Smart" in a sense that understands if the flow is incompressible
    or compressible from the thermo library.

    Check bounded scheme

SourceFiles
    smartBoundedDdtScheme.C

\*---------------------------------------------------------------------------*/

#ifndef smartBoundedDdtScheme_H
#define smartBoundedDdtScheme_H

#include "finiteVolume/ddtSchemes/boundedDdtScheme/boundedDdtScheme.H"
#include "basicThermo/basicThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

class options;

/*---------------------------------------------------------------------------*\
                       Class smartBoundedDdtScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class smartBoundedDdtScheme
:
    public boundedDdtScheme<Type>
{
    // Private data

        mutable fv::options* fvOptionsPtr_;

    // Private Member Functions

        //- Find/construct and get the thermo
        const basicThermo& thermo(const objectRegistry& obr) const
        {
            return basicThermo::lookupOrCreate(obr);
        }

        //- Construct and return fvOptions source
        tmp<volScalarField> fvOptSrc(const volScalarField& rho) const;

        //- Disallow default bitwise copy construct
        smartBoundedDdtScheme(const smartBoundedDdtScheme&);

        //- Disallow default bitwise assignment
        void operator=(const smartBoundedDdtScheme&);


public:

    //- Runtime type information
    TypeName("smartBounded");


    // Constructors

        //- Construct from mesh and Istream
        smartBoundedDdtScheme(const fvMesh& mesh, Istream& is)
        :
            boundedDdtScheme<Type>(mesh, is),
            fvOptionsPtr_(nullptr)
        {}


    // Member Functions

        tmp<GeometricField<Type, fvPatchField, volMesh>> fvcDdt
        (
            const volScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<GeometricField<Type, fvPatchField, volMesh>> fvcDdt
        (
            const volScalarField& alpha,
            const volScalarField& rho,
            const GeometricField<Type, fvPatchField, volMesh>& psi
        );

        tmp<fvMatrix<Type>> fvmDdt
        (
            const volScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        tmp<fvMatrix<Type>> fvmDdt
        (
            const volScalarField& alpha,
            const volScalarField& rho,
            const GeometricField<Type, fvPatchField, volMesh>& psi
        );


        typedef typename ddtScheme<Type>::fluxFieldType fluxFieldType;

        virtual tmp<fluxFieldType> fvcDdtUfCorr
        (
            const GeometricField<Type, fvPatchField, volMesh>& U,
            const GeometricField<Type, fvsPatchField, surfaceMesh>& Uf,
            const word interName
        )
        {
            return this->boundedDdtScheme<Type>::fvcDdtUfCorr(U, Uf, interName);
        }


        virtual tmp<fluxFieldType> fvcDdtPhiCorr
        (
            const GeometricField<Type, fvPatchField, volMesh>& U,
            const fluxFieldType& phi,
            const word interName
        )
        {
            return this->boundedDdtScheme<Type>::fvcDdtPhiCorr
            (
                U, phi, interName
            );
        }

        virtual tmp<fluxFieldType> fvcDdtUfCorr
        (
            const volScalarField& rho,
            const GeometricField<Type, fvPatchField, volMesh>& U,
            const GeometricField<Type, fvsPatchField, surfaceMesh>& Uf,
            const word interName
        )
        {
            return this->boundedDdtScheme<Type>::fvcDdtUfCorr
            (
                rho, U, Uf, interName
            );
        }

        virtual tmp<fluxFieldType> fvcDdtPhiCorr
        (
            const volScalarField& rho,
            const GeometricField<Type, fvPatchField, volMesh>& U,
            const fluxFieldType& phi,
            const word interName
        )
        {
            return this->boundedDdtScheme<Type>::fvcDdtPhiCorr
            (
                rho, U, phi, interName
            );
        }

};

template<>
tmp<surfaceScalarField> smartBoundedDdtScheme<scalar>::fvcDdtUfCorr
(
    const GeometricField<scalar, fvPatchField, volMesh>& U,
    const GeometricField<scalar, fvsPatchField, surfaceMesh>& Uf,
    const word
);

template<>
tmp<surfaceScalarField> smartBoundedDdtScheme<scalar>::fvcDdtPhiCorr
(
    const volScalarField& U,
    const surfaceScalarField& phi,
    const word
);

template<>
tmp<surfaceScalarField> smartBoundedDdtScheme<scalar>::fvcDdtUfCorr
(
    const volScalarField& rho,
    const volScalarField& U,
    const surfaceScalarField& Uf,
    const word
);

template<>
tmp<surfaceScalarField> smartBoundedDdtScheme<scalar>::fvcDdtPhiCorr
(
    const volScalarField& rho,
    const volScalarField& U,
    const surfaceScalarField& phi,
    const word
);



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "finiteVolume/ddtSchemes/smartBoundedDdtScheme/smartBoundedDdtScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
