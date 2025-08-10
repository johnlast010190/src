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
    (c) 2020 Esi Ltd.

Class
    Foam::CADSurface

Description
    Generate a triSurfaceMesh of a CAD file (IGES or STEP)
    with OpenCASCADE libraries

Note

    Example usage:
    \verbatim
    diffuser.igs
    {
        type CADSurface;

        // Accuracy of triangulation
        linDeflection 0.01;
        angDeflection 0.5;

        // Force triangulation even if a triangulated file
        // with required linDeflection has been previously saved
        forceTriangulation false; //default false

        // shapes used to identify regions
        nameRegionsBy   solid;//shell;
    }
    \endverbatim


SourceFiles
    CADSurface.C

\*---------------------------------------------------------------------------*/

#ifndef CADSurface_H
#define CADSurface_H

#include "searchableSurfaces/triSurfaceMesh/triSurfaceMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class CADSurface Declaration
\*---------------------------------------------------------------------------*/

class CADSurface
:
    public triSurfaceMesh
{

public:

    //- Runtime type information
    TypeName("CADSurface");

    // Constructors

        //- Construct from IO and dictionary.
        CADSurface
        (
            const IOobject& io,
            const dictionary& dict
        );


    //- Destructor
    virtual ~CADSurface();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
