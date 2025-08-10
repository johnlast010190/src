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

InNamespace
    Foam::fac

Description
    Calculate the divergence of the given field.

SourceFiles
    facDiv.C

\*---------------------------------------------------------------------------*/


#ifndef facNDiv_H
#define facNDiv_H

#include "fields/areaFields/areaFieldsFwd.H"
#include "fields/edgeFields/edgeFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Namespace fac functions Declaration
\*---------------------------------------------------------------------------*/

namespace fac
{
    template<class Type>
    tmp<GeometricField<Type, faPatchField, areaMesh>> ndiv
    (
        const GeometricField<Type, faePatchField, faEdgeMesh>&
    );

    template<class Type>
    tmp<GeometricField<Type, faPatchField, areaMesh>> ndiv
    (
        const tmp<GeometricField<Type, faePatchField, faEdgeMesh>>&
    );


    template<class Type>
    tmp
    <
        GeometricField
        <typename innerProduct<vector, Type>::type, faPatchField, areaMesh>
    > ndiv
    (
        const GeometricField<Type, faPatchField, areaMesh>&,
        const word& name
    );

    template<class Type>
    tmp
    <
        GeometricField
        <typename innerProduct<vector, Type>::type, faPatchField, areaMesh>
    > ndiv
    (
        const tmp<GeometricField<Type, faPatchField, areaMesh>>&,
        const word& name
    );


    template<class Type>
    tmp
    <
        GeometricField
        <typename innerProduct<vector, Type>::type, faPatchField, areaMesh>
    > ndiv
    (
        const GeometricField<Type, faPatchField, areaMesh>&
    );

    template<class Type>
    tmp
    <
        GeometricField
        <typename innerProduct<vector, Type>::type, faPatchField, areaMesh>
    > ndiv
    (
        const tmp<GeometricField<Type, faPatchField, areaMesh>>&
    );


    template<class Type>
    tmp<GeometricField<Type, faPatchField, areaMesh>> ndiv
    (
        const edgeScalarField&,
        const GeometricField<Type, faPatchField, areaMesh>&,
        const word& name
    );

    template<class Type>
    tmp<GeometricField<Type, faPatchField, areaMesh>> ndiv
    (
        const tmp<edgeScalarField>&,
        const GeometricField<Type, faPatchField, areaMesh>&,
        const word& name
    );

    template<class Type>
    tmp<GeometricField<Type, faPatchField, areaMesh>> ndiv
    (
        const edgeScalarField&,
        const tmp<GeometricField<Type, faPatchField, areaMesh>>&,
        const word& name
    );

    template<class Type>
    tmp<GeometricField<Type, faPatchField, areaMesh>> ndiv
    (
        const tmp<edgeScalarField>&,
        const tmp<GeometricField<Type, faPatchField, areaMesh>>&,
        const word& name
    );


    template<class Type>
    tmp<GeometricField<Type, faPatchField, areaMesh>> ndiv
    (
        const edgeScalarField&,
        const GeometricField<Type, faPatchField, areaMesh>&
    );

    template<class Type>
    tmp<GeometricField<Type, faPatchField, areaMesh>> ndiv
    (
        const tmp<edgeScalarField>&,
        const GeometricField<Type, faPatchField, areaMesh>&
    );

    template<class Type>
    tmp<GeometricField<Type, faPatchField, areaMesh>> ndiv
    (
        const edgeScalarField&,
        const tmp<GeometricField<Type, faPatchField, areaMesh>>&
    );

    template<class Type>
    tmp<GeometricField<Type, faPatchField, areaMesh>> ndiv
    (
        const tmp<edgeScalarField>&,
        const tmp<GeometricField<Type, faPatchField, areaMesh>>&
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "finiteArea/fac/facNDiv.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
