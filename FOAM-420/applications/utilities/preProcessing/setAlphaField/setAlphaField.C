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
    (c) 2016-2017 DHI
    (c) 2017 OpenCFD Ltd.

Application
    setAlphaField

Description
    Uses isoCutCell to create a volume fraction field from either a cylinder,
    a sphere or a plane.

    Original code supplied by Johan Roenby, DHI (2016)

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/include/fvCFD.H"
#include "fvMatrices/solvers/isoAdvection/isoCutFace/isoCutFace.H"
#include "fvMatrices/solvers/isoAdvection/isoCutCell/isoCutCell.H"
#include "primitives/enums/Enum.H"
#include "global/constants/mathematical/mathematicalConstants.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class shapeSelector
{
    public:

        enum class shapeType
        {
            PLANE,
            SPHERE,
            CYLINDER,
            SIN
        };

    static const Foam::Enum<shapeType> shapeTypeNames;
};

const Foam::Enum<shapeSelector::shapeType> shapeSelector::shapeTypeNames
{
    { shapeSelector::shapeType::PLANE, "plane" },
    { shapeSelector::shapeType::SPHERE, "sphere" },
    { shapeSelector::shapeType::CYLINDER, "cylinder" },
    { shapeSelector::shapeType::SIN, "sin" },
};


int main(int argc, char *argv[])
{
    #include "include/setRootCase.H"
    #include "include/createTime.H"
    #include "include/createMesh.H"

    Info<< "Reading setAlphaFieldDict\n" << endl;

    IOdictionary dict
    (
        IOobject
        (
            "setAlphaFieldDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const shapeSelector::shapeType surfType
    (
        shapeSelector::shapeTypeNames.read(dict.lookup("type"))
    );
    const vector centre(dict.lookup("centre"));
    const word fieldName(dict.lookup("field"));

    Info<< "Reading field " << fieldName << "\n" << endl;
    volScalarField alpha1
    (
        IOobject
        (
            fieldName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    scalar f0 = 0.0;
    scalarField f(mesh.points().size());

    Info<< "Processing type '" << shapeSelector::shapeTypeNames[surfType]
        << "'" << endl;

    switch (surfType)
    {
        case shapeSelector::shapeType::PLANE:
        {
            const vector direction(dict.lookup("direction"));

            f = -(mesh.points() - centre) & (direction/mag(direction));
            f0 = 0.0;
            break;
        }
        case shapeSelector::shapeType::SPHERE:
        {
            const scalar radius(readScalar(dict.lookup("radius")));

            f = -mag(mesh.points() - centre);
            f0 = -radius;
            break;
        }
        case shapeSelector::shapeType::CYLINDER:
        {
            const scalar radius(readScalar(dict.lookup("radius")));
            const vector direction(dict.lookup("direction"));

            f = -sqrt
            (
                sqr(mag(mesh.points() - centre))
              - sqr(mag((mesh.points() - centre) & direction))
            );
            f0 = -radius;
            break;
        }
        case shapeSelector::shapeType::SIN:
        {
            const scalar period(readScalar(dict.lookup("period")));
            const scalar amplitude(readScalar(dict.lookup("amplitude")));
            const vector up(dict.lookup("up"));
            const vector direction(dict.lookup("direction"));

            const scalarField xx
            (
                (mesh.points() - centre) & direction/mag(direction)
            );
            const scalarField zz((mesh.points() - centre) & up/mag(up));

            f = amplitude*Foam::sin(2*mathematical::pi*xx/period) - zz;
            f0 = 0;
            break;
        }
    }


    // Define function on mesh points and isovalue

    // Calculating alpha1 volScalarField from f = f0 isosurface
    isoCutCell icc(mesh, f);
    icc.volumeOfFluid(alpha1, f0);

    // Writing volScalarField alpha1
    ISstream::defaultPrecision(18);
    alpha1.write();

    Info<< nl << "Phase-1 volume fraction = "
        << alpha1.weightedAverage(mesh.Vsc()).value()
        << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
        << "  Max(" << alpha1.name() << ") - 1 = " << max(alpha1).value() - 1
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
