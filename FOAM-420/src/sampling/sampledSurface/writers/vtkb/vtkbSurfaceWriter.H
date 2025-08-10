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
    (c) 2011 OpenFOAM Foundation
    (c) 2015-2016 OpenCFD Ltd.

Class
    Foam::vtkbSurfaceWriter

Description
    A surfaceWriter for VTK legacy binary format.

SourceFiles
    vtkbSurfaceWriter.C

\*---------------------------------------------------------------------------*/

#ifndef vtkbSurfaceWriter_H
#define vtkbSurfaceWriter_H

#include "sampledSurface/writers/surfaceWriter.H"
#include "vtk/output/foamVtkOutputOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class vtkbSurfaceWriter Declaration
\*---------------------------------------------------------------------------*/

class vtkbSurfaceWriter
:
    public surfaceWriter
{
    // Private data

    // Private Member Functions

        static void writeGeometry
        (
            std::ofstream&,
            autoPtr<vtk::formatter>&,
            const meshedSurf&
        );


        //- Templated write operation
        template<class Type>
        fileName writeTemplate
        (
            const fileName& outputDir,
            const fileName& surfaceName,
            const meshedSurf& surf,
            const word& fieldName,
            const Field<Type>& values,
            const bool isNodeValues,
            const bool verbose
        ) const;


public:

    //- Runtime type information
    TypeName("vtkb");


    // Constructors

        //- Construct null
        vtkbSurfaceWriter();

        //- Construct with some output options
        vtkbSurfaceWriter(const dictionary& options);


    //- Destructor
    virtual ~vtkbSurfaceWriter();


    // Member Functions

        //- Write single surface geometry to file.
        virtual fileName write
        (
            const fileName& outputDir,
            const fileName& surfaceName,
            const meshedSurf& surf,
            const bool verbose = false
        ) const; // override


        //- Write scalarField for a single surface to file.
        //  One value per face or vertex (isNodeValues = true)
        virtual fileName write
        (
            const fileName& outputDir,      // <case>/surface/TIME
            const fileName& surfaceName,    // name of surface
            const meshedSurf& surf,
            const word& fieldName,          // name of field
            const Field<scalar>& values,
            const bool isNodeValues,
            const bool verbose = false
        ) const; // override

        //- Write vectorField for a single surface to file.
        //  One value per face or vertex (isNodeValues = true)
        virtual fileName write
        (
            const fileName& outputDir,      // <case>/surface/TIME
            const fileName& surfaceName,    // name of surface
            const meshedSurf& surf,
            const word& fieldName,          // name of field
            const Field<vector>& values,
            const bool isNodeValues,
            const bool verbose = false
        ) const; // override

        //- Write sphericalTensorField for a single surface to file.
        //  One value per face or vertex (isNodeValues = true)
        virtual fileName write
        (
            const fileName& outputDir,      // <case>/surface/TIME
            const fileName& surfaceName,    // name of surface
            const meshedSurf& surf,
            const word& fieldName,          // name of field
            const Field<sphericalTensor>& values,
            const bool isNodeValues,
            const bool verbose = false
        ) const; // override

        //- Write symmTensorField for a single surface to file.
        //  One value per face or vertex (isNodeValues = true)
        virtual fileName write
        (
            const fileName& outputDir,      // <case>/surface/TIME
            const fileName& surfaceName,    // name of surface
            const meshedSurf& surf,
            const word& fieldName,          // name of field
            const Field<symmTensor>& values,
            const bool isNodeValues,
            const bool verbose = false
        ) const; // override

        //- Write tensorField for a single surface to file.
        //  One value per face or vertex (isNodeValues = true)
        virtual fileName write
        (
            const fileName& outputDir,      // <case>/surface/TIME
            const fileName& surfaceName,    // name of surface
            const meshedSurf& surf,
            const word& fieldName,          // name of field
            const Field<tensor>& values,
            const bool isNodeValues,
            const bool verbose = false
        ) const; // override

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "sampledSurface/writers/vtkb/vtkbSurfaceWriterTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
