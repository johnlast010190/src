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

Description

\*---------------------------------------------------------------------------*/

#include "interpolation/edgeInterpolation/edgeInterpolate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fac
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Return weighting factors for scheme given by name in dictionary
template<class Type>
tmp<edgeInterpolationScheme<Type>> scheme
(
    const edgeScalarField& faceFlux,
    Istream& streamData
)
{
    return edgeInterpolationScheme<Type>::New
    (
        faceFlux.mesh(),
        faceFlux,
        streamData
    );
}


// Return weighting factors for scheme given by name in dictionary
template<class Type>
tmp<edgeInterpolationScheme<Type>> scheme
(
    const edgeScalarField& faceFlux,
    const word& name
)
{
    return edgeInterpolationScheme<Type>::New
    (
        faceFlux.mesh(),
        faceFlux,
        faceFlux.mesh().schemes().interpolationScheme(name)
    );
}


// Return weighting factors for scheme given by name in dictionary
template<class Type>
tmp<edgeInterpolationScheme<Type>> scheme
(
    const faMesh& mesh,
    Istream& streamData
)
{
    return edgeInterpolationScheme<Type>::New
    (
        mesh,
        streamData
    );
}


// Return weighting factors for scheme given by name in dictionary
template<class Type>
tmp<edgeInterpolationScheme<Type>> scheme
(
    const faMesh& mesh,
    const word& name
)
{
    return edgeInterpolationScheme<Type>::New
    (
        mesh,
        mesh.schemes().interpolationScheme(name)
    );
}


// Interpolate field onto faces using scheme given by name in dictionary
template<class Type>
tmp<GeometricField<Type, faePatchField, faEdgeMesh>>
interpolate
(
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    const edgeScalarField& faceFlux,
    Istream& schemeData
)
{
#   ifdef DEBUGInterpolations
    if (edgeInterpolation::debug)
    {
        Info<< "interpolate"
            << "(const GeometricField<Type, faPatchField, areaMesh>&, "
            << "const edgeScalarField&, Istream&) : "
            << "interpolating GeometricField<Type, faPatchField, areaMesh> "
            << endl;
    }
#   endif

    return scheme<Type>(faceFlux, schemeData)().interpolate(vf);
}


// Interpolate field onto faces using scheme given by name in dictionary
template<class Type>
tmp<GeometricField<Type, faePatchField, faEdgeMesh>>
interpolate
(
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    const edgeScalarField& faceFlux,
    const word& name
)
{
#   ifdef DEBUGInterpolations
    if (edgeInterpolation::debug)
    {
        Info<< "interpolate"
            << "(const GeometricField<Type, faPatchField, areaMesh>&, "
            << "const edgeScalarField&, const word&) : "
            << "interpolating GeometricField<Type, faPatchField, areaMesh> "
            << "using " << name
            << endl;
    }
#   endif

    return scheme<Type>(faceFlux, name)().interpolate(vf);
}

// Interpolate field onto faces using scheme given by name in dictionary
template<class Type>
tmp<GeometricField<Type, faePatchField, faEdgeMesh>>
interpolate
(
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tvf,
    const edgeScalarField& faceFlux,
    const word& name
)
{
    tmp<GeometricField<Type, faePatchField, faEdgeMesh>> tsf =
        interpolate(tvf(), faceFlux, name);

    tvf.clear();

    return tsf;
}

// Interpolate field onto faces using scheme given by name in dictionary
template<class Type>
tmp<GeometricField<Type, faePatchField, faEdgeMesh>>
interpolate
(
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    const tmp<edgeScalarField>& tFaceFlux,
    const word& name
)
{
    tmp<GeometricField<Type, faePatchField, faEdgeMesh>> tsf =
        interpolate(vf, tFaceFlux(), name);

    tFaceFlux.clear();

    return tsf;
}

// Interpolate field onto faces using scheme given by name in dictionary
template<class Type>
tmp<GeometricField<Type, faePatchField, faEdgeMesh>>
interpolate
(
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tvf,
    const tmp<edgeScalarField>& tFaceFlux,
    const word& name
)
{
    tmp<GeometricField<Type, faePatchField, faEdgeMesh>> tsf =
        interpolate(tvf(), tFaceFlux(), name);

    tvf.clear();
    tFaceFlux.clear();

    return tsf;
}


// Interpolate field onto faces using scheme given by name in dictionary
template<class Type>
tmp<GeometricField<Type, faePatchField, faEdgeMesh>>
interpolate
(
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    Istream& schemeData
)
{
#   ifdef DEBUGInterpolations
    if (edgeInterpolation::debug)
    {
        Info<< "interpolate"
            << "(const GeometricField<Type, faPatchField, areaMesh>&, "
            << "Istream&) : "
            << "interpolating GeometricField<Type, faPatchField, areaMesh> "
            << endl;
    }
#   endif

    return scheme<Type>(vf.mesh(), schemeData)().interpolate(vf);
}

// Interpolate field onto faces using scheme given by name in dictionary
template<class Type>
tmp<GeometricField<Type, faePatchField, faEdgeMesh>>
interpolate
(
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    const word& name
)
{
#   ifdef DEBUGInterpolations
    if (edgeInterpolation::debug)
    {
        Info<< "interpolate"
            << "(const GeometricField<Type, faPatchField, areaMesh>&, "
            << "const word&) : "
            << "interpolating GeometricField<Type, faPatchField, areaMesh> "
            << "using " << name
            << endl;
    }
#   endif

    return scheme<Type>(vf.mesh(), name)().interpolate(vf);
}

// Interpolate field onto faces using scheme given by name in dictionary
template<class Type>
tmp<GeometricField<Type, faePatchField, faEdgeMesh>>
interpolate
(
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tvf,
    const word& name
)
{
    tmp<GeometricField<Type, faePatchField, faEdgeMesh>> tsf =
        interpolate(tvf(), name);

    tvf.clear();

    return tsf;
}


// Interpolate field onto faces using central differencing
template<class Type>
tmp<GeometricField<Type, faePatchField, faEdgeMesh>>
interpolate
(
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{

#   ifdef DEBUGInterpolations
    if (edgeInterpolation::debug)
    {
        Info<< "interpolate"
            << "(const GeometricField<Type, faPatchField, areaMesh>&) : "
            << "interpolating GeometricField<Type, faPatchField, areaMesh> "
            << "using run-time selected scheme"
            << endl;
    }
#   endif

    return interpolate(vf, "interpolate(" + vf.name() + ')');
}


// Interpolate field onto faces using central differencing
template<class Type>
tmp<GeometricField<Type, faePatchField, faEdgeMesh>>
interpolate
(
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tvf
)
{
    tmp<GeometricField<Type, faePatchField, faEdgeMesh>> tsf =
        interpolate(tvf());
    tvf.clear();
    return tsf;
}

template<class Type>
tmp<FieldField<faePatchField, Type>> interpolate
(
    const FieldField<faPatchField, Type>& fvpff
)
{
    FieldField<faePatchField, Type>* fvspffPtr
    (
        new FieldField<faePatchField, Type>(fvpff.size())
    );

    forAll(*fvspffPtr, patchi)
    {
        fvspffPtr->set
        (
            patchi,
            faePatchField<Type>::NewCalculatedType(fvpff[patchi].patch()).ptr()
        );
        (*fvspffPtr)[patchi] = fvpff[patchi];
    }

    return tmp<FieldField<faePatchField, Type>>(fvspffPtr);
}


template<class Type>
tmp<FieldField<faePatchField, Type>> interpolate
(
    const tmp<FieldField<faPatchField, Type>>& tfvpff
)
{
    tmp<FieldField<faePatchField, Type>> tfvspff = interpolate(tfvpff());
    tfvpff.clear();
    return tfvspff;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fac

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
