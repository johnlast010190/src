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

Application
    slicedFieldTest

Description

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/include/fvCFD.H"
#include "fields/GeometricFields/SlicedGeometricField/SlicedGeometricField.H"
#include "fields/fvPatchFields/basic/sliced/slicedFvPatchFields.H"
#include "fields/surfaceFields/slicedSurfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"

    Info<< "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    //Info<< min(p, p);

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    SlicedGeometricField<vector, fvPatchField, slicedFvPatchField, volMesh>
    C
    (
        IOobject
        (
            "C2",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimLength,
        mesh.cellCentres(),
        mesh.faceCentres()
    );

    Info<< C << endl;
    Info<< (C & U) << endl;

    SlicedGeometricField
    <
         vector,
         fvsPatchField,
         slicedFvsPatchField,
         surfaceMesh
    >
    Sf
    (
        IOobject
        (
            "Sf2",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimArea,
        mesh.faceAreas()
    );

    //Info<< Sf << endl;

    return 0;
}


// ************************************************************************* //
