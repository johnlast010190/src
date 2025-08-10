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
    (c) 2017 Esi Ltd.

InNamespace
    Foam::fvc

Description
    Reconstruct volField from a face flux field.

SourceFiles
    fvcReconstruct.C

\*---------------------------------------------------------------------------*/


#ifndef fvcReconstruct_H
#define fvcReconstruct_H

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
    template<class Type>
    tmp
    <
        GeometricField
        <typename outerProduct<vector, Type>::type, fvPatchField, volMesh>
    > reconstruct
    (
        const GeometricField<Type, fvsPatchField, surfaceMesh>&
    );

    template<class Type>
    tmp
    <
        GeometricField
        <typename outerProduct<vector, Type>::type, fvPatchField, volMesh>
    > reconstruct
    (
        const tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>&
    );

    tmp<volScalarField> reconstructMag(const surfaceScalarField&);
    tmp<volScalarField> reconstructMag(const tmp<surfaceScalarField>&);

    template<class Type>
    tmp
    <
        GeometricField
        <typename outerProduct<vector, Type>::type, fvPatchField, volMesh>
    > reconstructV
    (
        const GeometricField<Type, fvsPatchField, surfaceMesh>&
    );

    template<class Type>
    tmp
    <
        GeometricField
        <typename outerProduct<vector, Type>::type, fvPatchField, volMesh>
    > reconstructV
    (
        const tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>&
    );


    template<class Type>
    tmp<GeometricField <Type, fvPatchField, volMesh>> reconstructSmooth
    (
        const GeometricField<Type, fvPatchField, volMesh> &,
        const GeometricField<Type, fvsPatchField, surfaceMesh> &
    );

    template<class Type>
    tmp<GeometricField <Type, fvPatchField, volMesh>> reconstructSmooth
    (
        const GeometricField<Type, fvPatchField, volMesh> &,
        const tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>&
    );

    template<class Type>
    tmp<GeometricField <Type, fvPatchField, volMesh>> reconstructSmooth
    (
        const GeometricField<Type, fvPatchField, volMesh> &
    );

    template<class Type>
    tmp<GeometricField <Type, fvPatchField, volMesh>> reconstructSmooth
    (
        const tmp<GeometricField<Type, fvPatchField, volMesh>>&
    );


}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "finiteVolume/fvc/fvcReconstruct.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
