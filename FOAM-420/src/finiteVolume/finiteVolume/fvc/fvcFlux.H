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
    (c) 2011-2016 OpenFOAM Foundation

InNamespace
    Foam::fvc

Description
    Calculate the face-flux of the given field.

SourceFiles
    fvcFlux.C

\*---------------------------------------------------------------------------*/


#ifndef fvcFlux_H
#define fvcFlux_H

#include "fields/volFields/volFieldsFwd.H"
#include "fields/surfaceFields/surfaceFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Namespace fvc functions Declaration
\*---------------------------------------------------------------------------*/

namespace fvc
{
    //- Return the face-flux field obtained from the given volVectorField
    tmp<surfaceScalarField> flux(const volVectorField& vvf);

    //- Return the face-flux field obtained from the given tmp volVectorField
    tmp<surfaceScalarField> flux(const tmp<volVectorField>& tvvf);

    template<class Type>
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> flux
    (
        const surfaceScalarField&,
        const tmp<GeometricField<Type, fvPatchField, volMesh>>&
    );

    template<class Type>
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> flux
    (
        const surfaceScalarField&,
        const GeometricField<Type, fvPatchField, volMesh>&,
        const word& name
    );

    template<class Type>
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> flux
    (
        const tmp<surfaceScalarField>&,
        const GeometricField<Type, fvPatchField, volMesh>&,
        const word& name
    );

    template<class Type>
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> flux
    (
        const surfaceScalarField&,
        const tmp<GeometricField<Type, fvPatchField, volMesh>>&,
        const word& name
    );

    template<class Type>
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> flux
    (
        const tmp<surfaceScalarField>&,
        const tmp<GeometricField<Type, fvPatchField, volMesh>>&,
        const word& name
    );


    template<class Type>
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> flux
    (
        const surfaceScalarField&,
        const GeometricField<Type, fvPatchField, volMesh>&
    );

    template<class Type>
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> flux
    (
        const tmp<surfaceScalarField>&,
        const GeometricField<Type, fvPatchField, volMesh>&
    );

    template<class Type>
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> flux
    (
        const surfaceScalarField&,
        const tmp<GeometricField<Type, fvPatchField, volMesh>>&
    );

    template<class Type>
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> flux
    (
        const tmp<surfaceScalarField>&,
        const tmp<GeometricField<Type, fvPatchField, volMesh>>&
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "finiteVolume/fvc/fvcFluxTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
