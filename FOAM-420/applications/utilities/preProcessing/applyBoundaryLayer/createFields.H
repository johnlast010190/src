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
    (c) 2011-2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

Info<< "Calculating wall distance field" << endl;
volScalarField y
(
    IOobject
    (
        "y",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimLength, 0.0),
    zeroGradientFvPatchScalarField::typeName
);
y.primitiveFieldRef() = wallDist::New(mesh).y().primitiveField();
y.correctBoundaryConditions();


// Set the mean boundary-layer thickness
dimensionedScalar ybl("ybl", dimLength, 0);

if (args.optionFound("ybl"))
{
    // If the boundary-layer thickness is provided use it
    ybl.value() = args.optionRead<scalar>("ybl");
}
else if (args.optionFound("Cbl"))
{
    // Calculate boundary layer thickness as Cbl*mean distance to wall
    ybl.value() = gAverage(y)*args.optionRead<scalar>("Cbl");
}

Info<< "\nCreating boundary-layer for U of thickness "
    << ybl.value() << " m" << nl << endl;

Info<< "Creating mask field" << endl;
volScalarField mask
(
    IOobject
    (
        "mask",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0),
    zeroGradientFvPatchScalarField::typeName
);


// ************************************************************************* //
