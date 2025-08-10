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

Namespace
    fac

Description

SourceFiles
    edgeInterpolate.C

\*---------------------------------------------------------------------------*/

#ifndef edgeInterpolate_H
#define edgeInterpolate_H

#include "memory/tmp/tmp.H"
#include "fields/areaFields/areaFieldsFwd.H"
#include "fields/edgeFields/edgeFieldsFwd.H"
#include "interpolation/edgeInterpolation/edgeInterpolationScheme/edgeInterpolationScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Namespace fac functions Declaration
\*---------------------------------------------------------------------------*/

namespace fac
{
    //- Return weighting factors for scheme given from Istream
    template<class Type>
    static tmp<edgeInterpolationScheme<Type>> scheme
    (
        const edgeScalarField& faceFlux,
        Istream& schemeData
    );

    //- Return weighting factors for scheme given by name in dictionary
    template<class Type>
    static tmp<edgeInterpolationScheme<Type>> scheme
    (
        const edgeScalarField& faceFlux,
        const word& name
    );


    //- Return weighting factors for scheme given from Istream
    template<class Type>
    static tmp<edgeInterpolationScheme<Type>> scheme
    (
        const faMesh& mesh,
        Istream& schemeData
    );

    //- Return weighting factors for scheme given by name in dictionary
    template<class Type>
    static tmp<edgeInterpolationScheme<Type>> scheme
    (
        const faMesh& mesh,
        const word& name
    );


    //- Interpolate field onto faces using scheme given by Istream
    template<class Type>
    static tmp<GeometricField<Type, faePatchField, faEdgeMesh>> interpolate
    (
        const GeometricField<Type, faPatchField, areaMesh>& tvf,
        const edgeScalarField& faceFlux,
        Istream& schemeData
    );

    //- Interpolate field onto faces using scheme given by name in faSchemes
    template<class Type>
    static tmp<GeometricField<Type, faePatchField, faEdgeMesh>> interpolate
    (
        const GeometricField<Type, faPatchField, areaMesh>& tvf,
        const edgeScalarField& faceFlux,
        const word& name
    );

    //- Interpolate field onto faces using scheme given by name in faSchemes
    template<class Type>
    static tmp<GeometricField<Type, faePatchField, faEdgeMesh>> interpolate
    (
        const tmp<GeometricField<Type, faPatchField, areaMesh>>& tvf,
        const edgeScalarField& faceFlux,
        const word& name
    );

    //- Interpolate field onto faces using scheme given by name in faSchemes
    template<class Type>
    static tmp<GeometricField<Type, faePatchField, faEdgeMesh>> interpolate
    (
        const GeometricField<Type, faPatchField, areaMesh>& tvf,
        const tmp<edgeScalarField>& faceFlux,
        const word& name
    );

    //- Interpolate field onto faces using scheme given by name in faSchemes
    template<class Type>
    static tmp<GeometricField<Type, faePatchField, faEdgeMesh>> interpolate
    (
        const tmp<GeometricField<Type, faPatchField, areaMesh>>& tvf,
        const tmp<edgeScalarField>& faceFlux,
        const word& name
    );


    //- Interpolate field onto faces using scheme given by Istream
    template<class Type>
    static tmp<GeometricField<Type, faePatchField, faEdgeMesh>> interpolate
    (
        const GeometricField<Type, faPatchField, areaMesh>& tvf,
        Istream& schemeData
    );

    //- Interpolate field onto faces using scheme given by name in faSchemes
    template<class Type>
    static tmp<GeometricField<Type, faePatchField, faEdgeMesh>> interpolate
    (
        const GeometricField<Type, faPatchField, areaMesh>& tvf,
        const word& name
    );

    //- Interpolate field onto faces using scheme given by name in faSchemes
    template<class Type>
    static tmp<GeometricField<Type, faePatchField, faEdgeMesh>> interpolate
    (
        const tmp<GeometricField<Type, faPatchField, areaMesh>>& tvf,
        const word& name
    );


    //- Interpolate tmp field onto faces using central differencing
    template<class Type>
    static tmp<GeometricField<Type, faePatchField, faEdgeMesh>> interpolate
    (
        const tmp<GeometricField<Type, faPatchField, areaMesh>>& tvf
    );

    //- Interpolate field onto faces using central differencing
    template<class Type>
    static tmp<GeometricField<Type, faePatchField, faEdgeMesh>> interpolate
    (
        const GeometricField<Type, faPatchField, areaMesh>& tvf
    );

    //- Interpolate boundary field onto faces (simply a type conversion)
    template<class Type>
    static tmp<FieldField<faePatchField, Type>> interpolate
    (
        const FieldField<faPatchField, Type>& fvpff
    );

    //- Interpolate boundary field onto faces (simply a type conversion)
    template<class Type>
    static tmp<FieldField<faePatchField, Type>> interpolate
    (
        const tmp<FieldField<faPatchField, Type>>& tfvpff
    );

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "interpolation/edgeInterpolation/edgeInterpolate.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
