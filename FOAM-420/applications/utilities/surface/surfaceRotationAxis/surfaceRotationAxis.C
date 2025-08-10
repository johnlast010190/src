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
    (c) 2015 OpenCFD Ltd.

Application
    surfaceRotationAxis

Group
    grpSurfaceUtilities

Description
    Determines whether a surface is rotational to within a tolerance and
    provides the axis of rotation if it is.

    The tolerance is the mean difference in radius between the principle
    axes in the plane of rotation. The tolerance is calculated by default
    as 0.1% of the diagonal dimension of the bounding box encapsulating
    the input surface. In the context of AMI sliding interfaces, the
    tolerance should be a fraction (0.5 or less) of the minimum cell size
    adjacent to the AMI.

\*---------------------------------------------------------------------------*/

#include "global/argList/argList.H"
#include "containers/Lists/ListOps/ListOps.H"
#include "triSurface/triSurface.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "meshTools/meshTools.H"
#include "primitives/random/Random/Random.H"
#include "primitives/transform/transform.H"
#include "db/IOstreams/IOstreams/IOmanip.H"
#include "primitives/Pair/Pair.H"
#include "momentOfInertia/momentOfInertia.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

tensor triSurfInertia
(
    const point& a,
    const point& b,
    const point& c,
    const point& cM
)
{
    vector arm = (a + b + c)/3.0 - cM;

    scalar area = triPointRef(a,b,c).mag();

    return (area*(sqr(arm)));
}


void surfaceInertia
(
    const triSurface& surf,
    vector& cM,
    tensor& moI
)
{
    bool doReduce = false;

    triFaceList triFaces(surf.size());

    forAll(surf, i)
    {
        triFaces[i] = triFace(surf[i]);
    }

    const pointField& pts(surf.points());

    // Reset properties for accumulation

    scalar weightSum = 0.0;
    cM = Zero;
    moI = Zero;

    // Find centre of mass

    forAll(triFaces, i)
    {
        const triFace& tri(triFaces[i]);

        triPointRef t
        (
            pts[tri[0]],
            pts[tri[1]],
            pts[tri[2]]
        );

        scalar triMag = t.mag();

        cM +=  triMag*t.centre();

        weightSum += triMag;
    }

    if (doReduce)
    {
        reduce(cM, sumOp<vector>());
    }

    cM /= weightSum;

    // Find inertia around centre of mass (cM)

    forAll(triFaces, i)
    {
        const triFace& tri(triFaces[i]);

        moI += triSurfInertia
        (
            pts[tri[0]],
            pts[tri[1]],
            pts[tri[2]],
            cM
        );

    }

    moI /= weightSum;

    if (doReduce)
    {
        reduce(moI, sumOp<tensor>());
    }

}

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Calculates origin and axis of rotation for a rotating body."
    );

    argList::noParallel();
    argList::validArgs.append("surfaceFile");
    argList::addOption
    (
        "rotationalTolerance",
        "scalar",
        "provide the tolerance (in [m]) used to determine if two Eigen values"
        "of the moment of inertia tensor are equivalent"
    );
    argList::addBoolOption
    (
        "verbose",
        "Output more information"
    );

    argList args(argc, argv);

    const fileName surfFileName = args[1];

    triSurface surf(surfFileName);

    scalar rotTol = 0.001 * boundBox(surf.points()).mag();
    args.optionReadIfPresent("rotationalTolerance", rotTol);
    const bool verbose = args.optionFound("verbose");


    vector cM = Zero;
    tensor moI = Zero;

    //Comment: it might be a good idea to normalize everything to reduce
    //rounding errors at extreme scales

    surfaceInertia(surf, cM, moI);

    vector eVal = eigenValues(moI);

    tensor eVec = eigenVectors(moI);

    //Two identical Eigen values indicate a rotational body
    //Characterize the Eigen properties and extract axis

    //output axis
    vector rotAxis = Zero;

    scalarList dEig(3, Zero);
    dEig[0] = mag(Foam::sqrt(mag(eVal.y())) - Foam::sqrt(mag(eVal.z())));
    dEig[1] = mag(Foam::sqrt(mag(eVal.z())) - Foam::sqrt(mag(eVal.x())));
    dEig[2] = mag(Foam::sqrt(mag(eVal.x())) - Foam::sqrt(mag(eVal.y())));

    bool rotational = false;

    if (max(dEig) < rotTol)
    {
        //sperical
        //need to check this initially as any axis is valid
        Info<< "The input surface is spherical and has no prefered rotational"
             << " axis." << endl;

        return 0;
    }
    else
    {
        List<vector> principal(3);

        principal[0] = eVec.x();
        principal[1] = eVec.y();
        principal[2] = eVec.z();

        forAll(dEig, i)
        {
            if (dEig[i] < rotTol)
            {
                rotAxis = principal[i];
                rotAxis /= mag(rotAxis);
                rotational = true;
                break;
            }
        }
    }

    if (verbose)
    {
        Info<< "Surface moment: " << nl
             << "(" << nl
             << tab << moI.x() << nl
             << tab << moI.y() << nl
             << tab << moI.z() << nl
             << ")" << endl;

        Info<< "Eigen values/vector : " << eVal[0] << " " <<  eVec.x() << nl
             << "Eigen values/vector : " << eVal[1] << " " <<  eVec.y() << nl
             << "Eigen values/vector : " << eVal[2] << " " <<  eVec.z() << nl
             << endl;

        Info<< "Rotational tolerance: " << rotTol << " [m]" << nl
             << "Mean rotational deviations around principle axes: " << dEig
             << " [m]" << endl;


    }

    if (rotational)
    {
        Info<< nl << setprecision(12)
            << "Centre of mass: " << cM << nl
            << "Rotation axis: " << rotAxis << endl;

        scalar scale = 0.5 * boundBox(surf.points()).mag();

        OFstream str("rotationAxis.obj");
        meshTools::writeOBJ(str, cM);
        meshTools::writeOBJ(str, scale*rotAxis);
    }
    else
    {
        Info<< nl
             << "The input surface is not rotationally symmetric to within the"
             << " specified tolerance:" << nl
             << tab << "- Tolerance: " << rotTol << " [m]" << nl
             << tab << "- Mean deviation from circular: " << min(dEig)
             << " [m]" << endl;
    }


    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
