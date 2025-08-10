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
    (c) 2011-2014 OpenFOAM Foundation
    (c) 2022 Esi Ltd

InNamespace
    Foam::fvc

Description
    Calculate the mesh motion flux and convert fluxes from absolute to relative
    and back.

SourceFiles
    fvcMeshPhi.C

\*---------------------------------------------------------------------------*/

#ifndef fvcMeshPhi_H
#define fvcMeshPhi_H

#include "fields/volFields/volFieldsFwd.H"
#include "fields/surfaceFields/surfaceFieldsFwd.H"
#include "dimensionedTypes/dimensionedTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Namespace fvc functions Declaration
\*---------------------------------------------------------------------------*/

namespace fvc
{
    template<class Type>
    tmp<surfaceScalarField> meshPhi
    (
        const GeometricField<Type, fvPatchField, volMesh>& vf
    );

    template<class Type>
    tmp<surfaceScalarField> meshPhi
    (
        const dimensionedScalar& rho,
        const GeometricField<Type, fvPatchField, volMesh>& vf
    );

    template<class Type>
    tmp<surfaceScalarField> meshPhi
    (
        const volScalarField& rho,
        const GeometricField<Type, fvPatchField, volMesh>& vf
    );


    //- Make the given flux relative
    template<class Type>
    void makeRelative
    (
        surfaceScalarField& phi,
        const GeometricField<Type, fvPatchField, volMesh>& vf
    );

    //- Make the given flux relative
    template<class Type>
    void makeRelative
    (
        surfaceScalarField& phi,
        const dimensionedScalar& rho,
        const GeometricField<Type, fvPatchField, volMesh>& vf
    );

    //- Make the given flux relative
    template<class Type>
    void makeRelative
    (
        surfaceScalarField& phi,
        const volScalarField& rho,
        const GeometricField<Type, fvPatchField, volMesh>& vf
    );


    //- Make the given flux relative
    template<class Type>
    void makeRelative
    (
        surfaceScalarField& phi,
        const volScalarField& rho,
        const surfaceScalarField& rhof,
        const GeometricField<Type, fvPatchField, volMesh>& vf
    );


    //- Make the given flux absolute
    template<class Type>
    void makeAbsolute
    (
        surfaceScalarField& phi,
        const GeometricField<Type, fvPatchField, volMesh>& vf
    );

    //- Make the given flux absolute
    template<class Type>
    void makeAbsolute
    (
        surfaceScalarField& phi,
        const dimensionedScalar& rho,
        const GeometricField<Type, fvPatchField, volMesh>& vf
    );

    //- Make the given flux absolute
    template<class Type>
    void makeAbsolute
    (
        surfaceScalarField& phi,
        const volScalarField& rho,
        const GeometricField<Type, fvPatchField, volMesh>& vf
    );


    //- Return the given absolute flux in relative form
    template<class Type>
    tmp<surfaceScalarField> relative
    (
        const tmp<surfaceScalarField>& tphi,
        const GeometricField<Type, fvPatchField, volMesh>& vf
    );

    //- Return the given absolute flux in relative form
    template<class Type>
    tmp<surfaceScalarField> relative
    (
        const tmp<surfaceScalarField>& tphi,
        const volScalarField& rho,
        const GeometricField<Type, fvPatchField, volMesh>& vf
    );


    //- Return the given relative flux in absolute form
    template<class Type>
    tmp<surfaceScalarField> absolute
    (
        const tmp<surfaceScalarField>& tphi,
        const GeometricField<Type, fvPatchField, volMesh>& vf
    );

    //- Return the given relative flux in absolute form
    template<class Type>
    tmp<surfaceScalarField> absolute
    (
        const tmp<surfaceScalarField>& tphi,
        const volScalarField& rho,
        const GeometricField<Type, fvPatchField, volMesh>& vf
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "finiteVolume/fvc/fvcMeshPhiTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
