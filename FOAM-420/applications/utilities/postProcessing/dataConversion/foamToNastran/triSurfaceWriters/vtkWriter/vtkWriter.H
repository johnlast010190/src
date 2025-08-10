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
    (c) 2010-2012 Esi Ltd.
    (c) 1991-2007 OpenCFD Ltd.

Class
    vtkWriter

Description

SourceFiles
    vtkWriter.C

\*---------------------------------------------------------------------------*/

#ifndef vtkWriter_H
#define vtkWriter_H

#include "triSurfaceWriter/triSurfaceWriter.H"
#include "primitives/Scalar/floatScalar/floatScalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class triSurface;

/*---------------------------------------------------------------------------*\
                           Class vtkWriter Declaration
\*---------------------------------------------------------------------------*/

class vtkWriter
:
    public triSurfaceWriter
{

    // Private Member Functions

        static void swapWord(label& word32);

        static void swapWords(const label nWords, label* words32);

        static void write
        (
            std::ostream& os,
            const bool binary,
            List<floatScalar>& fField
        );

        static void write
        (
            std::ostream& os,
            const bool binary,
            labelList& elems
        );

        static void insert(const point&, DynamicList<floatScalar>&);
        static void insert(const labelList&, DynamicList<label>&);
        static void insert(const scalarField&, DynamicList<floatScalar>&);
        static void insert
        (
            const List<point>&,
            DynamicList<floatScalar>&
        );

        static void writeGeometry
        (
            const bool binary,
            const triSurface& s,
            std::ostream& os
        );

public:

    //- Runtime type information
    TypeName("vtk");


    // Constructors

        //- Construct null
        vtkWriter();

    // Destructor

        virtual ~vtkWriter();


    // Member Functions

        //- Writes surface + data to file.
        virtual void write
        (
            const fileName& path,           // <root>/<case>/sampleSurfaces
            const fileName& timeDir,        // time name
            const fileName& surfaceName,    // name of surface
            const triSurface& surface,
            PtrList<triSurfaceScalarField>& scalarFields,
            PtrList<triSurfaceVectorField>& vectorFields,
            PtrList<triSurfaceTensorField>& tensorFields
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
