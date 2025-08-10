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
    (c) 2010-2016 Esi Ltd.

Class
    Foam::surfaceBoundarySetup

Description
    Interpret a structured input dictionary to create a boundary for
    a surfaceField. Note this is not the standard boundaryField format dictionary.

SourceFiles
    surfaceBoundarySetup.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceBoundarySetup_H
#define surfaceBoundarySetup_H

#include "surfaceMesh/surfaceMesh.H"
#include "fields/DimensionedFields/DimensionedField/DimensionedField.H"
#include "fields/fvsPatchFields/fvsPatchField/fvsPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class surfaceBoundarySetup Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class surfaceBoundarySetup
:
    public PtrList<fvsPatchField<Type>>
{

    //- default fvPatch-type wide boundaries
    void defaultBoundaries
        (const surfaceMesh&, const DimensionedField<Type, surfaceMesh>&);

    //- fvPatch-type wide boundaries
    void typeBoundaries
    (
        const surfaceMesh&,
        const DimensionedField<Type, surfaceMesh>&,
        const dictionary&
    );

    //- boundaries for named patches
    void namedBoundaries
    (
        const surfaceMesh&,
        const DimensionedField<Type, surfaceMesh>&,
        const dictionary&,
        bool strictPatchNameChecking
    );

public:

    // Constructors

        //- Construct from components.
        surfaceBoundarySetup
        (
            const surfaceMesh&,
            const DimensionedField<Type, surfaceMesh>&,
            const dictionary&,
            bool strictPatchNameChecking = true
        );

    // Destructor

        ~surfaceBoundarySetup();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#include "cfdTools/general/surfaceBoundarySetup/surfaceBoundarySetup.C"
#endif

#endif

// ************************************************************************* //
