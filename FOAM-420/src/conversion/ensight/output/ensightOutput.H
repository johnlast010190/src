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
    (c) 2016 OpenCFD Ltd.

Class
    Foam::ensightOutput

Description
    A collection of functions for writing ensight file content in parallel.

SourceFiles
    ensightOutput.C
    ensightOutputTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef ensightOutput_H
#define ensightOutput_H

#include "ensight/file/ensightFile.H"
#include "ensight/mesh/ensightMesh.H"
#include "memory/autoPtr/autoPtr.H"
#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class ensightOutput Declaration
\*---------------------------------------------------------------------------*/

class ensightOutput
{
    // Private Methods

    //- Write field content (component-wise) for the given ensight element type
    template<template<typename> class FieldContainer, class Type>
    static void writeFieldContent
    (
        const char* key,
        const FieldContainer<Type>& fld,
        ensightFile& os
    );

    //- Write a field of faces values as an indirect list,
    //  using the face ids from ensightFaces
    template<class Type>
    static bool writeFaceField
    (
        const Field<Type>& fld,
        const ensightFaces&,
        ensightFile& os
    );

    //- Write a field of faces values as a sublist,
    //  using the sublist sizes ensightFaces
    template<class Type>
    static bool writeFaceSubField
    (
        const Field<Type>& fld,
        const ensightFaces&,
        ensightFile& os
    );

    //- Write a field of cell values as an indirect list,
    //  using the cell ids from ensightCells
    template<class Type>
    static bool writeCellField
    (
        const Field<Type>& fld,
        const ensightCells&,
        ensightFile& os,
        const bool deprecatedOrder = false
    );


    //- Write volume field component-wise
    template<class Type>
    static bool writeField
    (
        const GeometricField<Type, fvPatchField, volMesh>& vf,
        const ensightMesh& ensMesh,
        ensightFile& os
    );

    //- Write point field component-wise
    template<class Type>
    static bool ensightPointField
    (
        const GeometricField<Type, pointPatchField, pointMesh>& pf,
        const ensightMesh& ensMesh,
        ensightFile& os
    );


    //- Disallow null constructor
    ensightOutput() = delete;

public:

    // Public Methods


    //- Write volume field component-wise
    template<class Type>
    static bool writeField
    (
        const GeometricField<Type, fvPatchField, volMesh>&,
        const ensightMesh& ensMesh,
        ensightFile& os,
        const bool nodeValues
    );


    //- Write volume field component-wise
    template<class Type>
    static inline bool writeField
    (
        const GeometricField<Type, fvPatchField, volMesh>& vf,
        const ensightMesh& ensMesh,
        autoPtr<ensightFile>& output,
        const bool nodeValues = false
    )
    {
        return writeField(vf, ensMesh, output.rawRef(), nodeValues);
    }

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ensight/output/ensightOutputTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
