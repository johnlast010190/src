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
    (c) 2016-2019 OpenCFD Ltd.

Namespace
    Foam::vtk

Description
    Namespace for handling VTK output.
    Contains classes and functions for writing VTK file content.

Namespace
    Foam::vtk::legacy

Description
    Namespace for legacy VTK output constants and functions.

SourceFiles
    foamVtkOutput.C
    foamVtkOutputTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef foamVtkOutput_H
#define foamVtkOutput_H

#include "memory/autoPtr/autoPtr.H"
#include "containers/Lists/PackedList/PackedBoolList.H"
#include "primitives/enums/Enum.H"
#include "vtk/core/foamVtkCore.H"
#include "vtk/format/foamVtkFormatter.H"
#include "primitives/Scalar/floatScalar/floatScalar.H"
#include "primitives/SymmTensor/symmTensor/symmTensor.H"
#include "db/IOstreams/IOstreams/IOstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class instant;
class globalIndex;

namespace vtk
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// General Functions

    //- Return a default asciiFormatter
    autoPtr<vtk::formatter> newFormatter
    (
        std::ostream& os,
        unsigned prec = IOstream::defaultPrecision()
    );

    //- Return a new formatter based on the specified format type
    autoPtr<vtk::formatter> newFormatter
    (
        std::ostream& os,
        const enum formatType fmtType,
        unsigned prec = IOstream::defaultPrecision()
    );


    //- Write a list of uint8_t values.
    //  The output does not include the payload size.
    void writeList
    (
        vtk::formatter& fmt,
        const UList<uint8_t>& values
    );

    //- Write a list of uint8_t values.
    //  The output does not include the payload size.
    void writeListParallel
    (
        vtk::formatter& fmt,
        const UList<uint8_t>& values
    );

    //- Write a value component-wise.
    template<class Type>
    inline void write
    (
        vtk::formatter& fmt,
        const Type& val
    );


    //- Write a list of values.
    //  The output does not include the payload size.
    template<class Type>
    void writeList
    (
        vtk::formatter& fmt,
        const UList<Type>& values
    );

    //- Write a list of values.
    //  The output does not include the payload size.
    template<class Type, unsigned N>
    void writeList
    (
        vtk::formatter& fmt,
        const FixedList<Type, N>& values
    );


    //- Write a list of values via indirect addressing.
    //  The output does not include the payload size.
    template<class Type>
    void writeList
    (
        vtk::formatter& fmt,
        const UList<Type>& values,
        const labelUList& addressing
    );

    //- Write a list of values via indirect addressing.
    //  The output does not include the payload size.
    template<class Type>
    void writeList
    (
        vtk::formatter& fmt,
        const UList<Type>& values,
        const PackedBoolList& selected
    );

    //- Write a list of values and a list of values via indirect addressing.
    //  The output does not include the payload size.
    template<class Type>
    void writeLists
    (
        vtk::formatter& fmt,
        const UList<Type>& values1,
        const UList<Type>& values2,
        const labelUList& addressing
    );


    //- Write a list of values.
    //  The output does not include the payload size.
    template<class Type>
    void writeListParallel
    (
        vtk::formatter& fmt,
        const UList<Type>& values
    );

    //- Write a list of values, with constant per-processor offset
    //  The output does not include the payload size.
    void writeListParallel
    (
        vtk::formatter& fmt,
        const UList<label>& values,
        const globalIndex& procOffset
    );

    //- Write a list of values via indirect addressing.
    //  The output does not include the payload size.
    template<class Type>
    void writeListParallel
    (
        vtk::formatter& fmt,
        const UList<Type>& values,
        const labelUList& addressing
    );

    //- Write a list of values via indirect addressing.
    //  The output does not include the payload size.
    template<class Type>
    void writeListParallel
    (
        vtk::formatter& fmt,
        const UList<Type>& values,
        const PackedBoolList& selected
    );

    //- Write a list of values and another list of values.
    //  The output does not include the payload size.
    template<class Type>
    void writeListsParallel
    (
        vtk::formatter& fmt,
        const UList<Type>& values1,
        const UList<Type>& values2
    );

    //- Write a list of values and a list of values via indirect addressing.
    //  The output does not include the payload size.
    template<class Type>
    void writeListsParallel
    (
        vtk::formatter& fmt,
        const UList<Type>& values1,
        const UList<Type>& values2,
        const labelUList& addressing
    );


/*---------------------------------------------------------------------------*\
                               Namespace legacy
\*---------------------------------------------------------------------------*/

namespace legacy
{

// Functions

    //- Emit header for legacy file.
    //  Writes "ASCII" or "BINARY" depending on specified type.
    void fileHeader(std::ostream& os, const std::string& title, bool binary);

    //- Emit header for legacy file, with "ASCII" or "BINARY" depending on
    //- the formatter type.
    //  If the contentType is non-empty, it is used for "DATASET" line.
    void fileHeader
    (
        vtk::formatter& fmt,
        const std::string& title,
        const std::string& contentType
    );

    //- Emit header for legacy file, with "ASCII" or "BINARY" depending on
    //- the formatter type.
    //  Includes "DATASET" with the specified dataset type.
    inline void fileHeader
    (
        vtk::formatter& fmt,
        const std::string& title,
        vtk::fileTag contentType
    );

    //- Emit header for legacy file, with "ASCII" or "BINARY" depending on
    //- the formatter type.
    //  Includes "DATASET" of the templated dataset type.
    template<vtk::fileTag ContentType>
    inline void fileHeader(vtk::formatter& fmt, const std::string& title);

    //- Emit header for POINTS (with trailing newline).
    inline void beginPoints(std::ostream& os, label nPoints);

    //- Emit header for POLYGONS (with trailing newline).
    //  The nConnectivity is the sum of all connectivity points used,
    //  but \b without additional space for the size prefixes.
    //  The additional prefix sizes are added internally.
    inline void beginPolys(std::ostream& os, label nPolys, label nConnectivity);


    //- Emit "FIELD FieldData <n>"
    inline void fieldData(vtk::formatter& fmt, label nFields);

    //- Emit legacy FIELD FieldData nFields.
    inline void beginFieldData(vtk::formatter& fmt, label nFields);

    //- Emit legacy CELL_DATA nCells, FIELD FieldData nFields.
    inline void beginCellData
    (
        vtk::formatter& fmt,
        label nCells,
        label nFields
    );

    //- Emit legacy POINT_DATA nPoints, FIELD FieldData nFields.
    inline void beginPointData
    (
        vtk::formatter& fmt,
        label nPoints,
        label nFields
    );


    // TODO:  This is included for backwards-compatability with FOAM, and needs runtime testing
    //- Use the enumerations vtk::fileTag::CELL_DATA, vtk::fileTag::POINT_DATA,
    //  to emit a legacy CELL_DATA, POINT_DATA element.
    //  The nEntries corresponds similarly to the number of cells or points,
    //  respectively.
    inline void dataHeader
    (
        std::ostream& os,
        const vtk::fileTag& dataTypeTag,
        const label nEntries,
        const label nFields
    );


    //- Emit "TimeValue" for a FIELD entry (name as per Catalyst output)
    inline void writeTimeValue(vtk::formatter& fmt, scalar timeValue);

    inline void floatField
    (
        std::ostream& os,
        const word& fieldName,
        const int nCmpt,
        const label nEntries
    )
    {
        os  << fieldName << ' ' << nCmpt << ' ' << nEntries << " float" << nl;
    }

    //- Start output of float field with the specified name.
    template<direction nComp>
    inline void floatField
    (
        vtk::formatter& fmt,
        const word& name,
        const label nEntries
    );

    //- Start output of double field with the specified name.
    template<direction nComp>
    inline void doubleField
    (
        vtk::formatter& fmt,
        const word& name,
        const label nEntries
    );

    //- Start output of int field with the specified name.
    template<direction nComp>
    inline void intField
    (
        vtk::formatter& fmt,
        const word& name,
        const label nEntries
    );

    inline void intField
    (
        std::ostream& os,
        const word& fieldName,
        const int nCmpt,
        const label nEntries
    )
    {
        os  << fieldName << ' ' << nCmpt << ' ' << nEntries << " int" << nl;
    }

} // End namespace legacy


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Template specialization for symmTensor ordering
template<>
inline void write(vtk::formatter& fmt, const symmTensor& val)
{
    // symmTensor ( XX, XY, XZ, YY, YZ, ZZ )
    // VTK order  ( XX, YY, ZZ, XY, YZ, XZ ) -> (0, 3, 5, 1, 4, 2)

    fmt.write(component(val, 0)); // XX
    fmt.write(component(val, 3)); // YY
    fmt.write(component(val, 5)); // ZZ
    fmt.write(component(val, 1)); // XY
    fmt.write(component(val, 4)); // YZ
    fmt.write(component(val, 2)); // XZ
}


} // End namespace vtk
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "vtk/output/foamVtkOutputI.H"

#ifdef NoRepository
    #include "vtk/output/foamVtkOutputTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
