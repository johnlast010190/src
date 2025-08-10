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
    (c) 1991-2009 OpenCFD Ltd.
    (c) 2017 Esi Ltd.

Application
    transformPoints

Description
    Stretches mesh between start and end point along specified vector.
    Automatic smoothing of endpoints using supplied expansion ratio.
    Currently stretches entire mesh in stretch zone - extension using
    sets to localise stretching would be useful.

Usage
    system/stretchMeshDict

\*---------------------------------------------------------------------------*/

#include "global/argList/argList.H"
#include "db/Time/Time.H"
#include "fvMesh/fvMesh.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "fields/ReadFields/ReadFields.H"
#include "fields/GeometricFields/pointFields/pointFields.H"
#include "fields/Fields/transformField/transformField.H"
#include "fields/GeometricFields/transformGeometricField/transformGeometricField.H"
#include "db/IOstreams/StringStreams/IStringStream.H"

using namespace Foam;
using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class GeoField>
void readAndRotateFields
(
    PtrList<GeoField>& flds,
    const fvMesh& mesh,
    const tensor& T,
    const IOobjectList& objects
)
{
    ReadFields(mesh, objects, flds);
    forAll(flds, i)
    {
        Info<< "Transforming " << flds[i].name() << endl;
        dimensionedTensor dimT("t", flds[i].dimensions(), T);
        transform(flds[i], dimT, flds[i]);
    }
}


//  Main program:

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Perform stretching of mesh according to stretchMeshDict"
    );

#include "include/addRegionOption.H"
    argList::addOption
    (
        "regions",
        "list",
        "List of regions to stretch"
    );

#include "include/setRootCase.H"
#include "include/createTime.H"

    wordList regionNames;
    wordList regionDirs;
    if (args.optionFound("regions"))
    {
        regionNames = args.optionRead<wordList>("regions");
        regionDirs = regionNames;
    }
    else if (args.optionFound("region"))
    {
        word regionName = args.optionRead<word>("region");
        regionNames = wordList(1, regionName);
        regionDirs = regionNames;
    }
    else
    {
        regionNames = wordList(1, fvMesh::defaultRegion);
        regionDirs = wordList(1, word::null);
    }

    forAll(regionNames, regionI)
    {
        const word& regionName = regionNames[regionI];
        const word& regionDir = regionDirs[regionI];

        Info<< "\n\nStretching for mesh " << regionName << nl
            << endl;

        pointIOField points
        (
            IOobject
            (
                "points",
                runTime.findInstance(regionDir/polyMesh::meshSubDir, "points"),
                regionDir/polyMesh::meshSubDir,
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );


        IOdictionary smDict
        (
            IOobject
            (
                "stretchMeshDict",
                runTime.system(),
                regionDir,
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        // Allow multiples transforms on a single mesh
        PtrList<dictionary> transformDicts(smDict.lookup("transforms"));

        forAll(transformDicts, dictI)
        {
            const dictionary& tDict = transformDicts[dictI];

            //read input

            vector basePoint = tDict.lookup("basePoint");
            vector stretchV = tDict.lookup("stretchDirection");
            stretchV /= mag(stretchV);
            scalar initL = readScalar(tDict.lookup("initialLength"));
            scalar finL = readScalar(tDict.lookup("stretchLength"));
            scalar xr = readScalar(tDict.lookup("expansionRatio"));
            if (xr <= 1)
            {
                FatalError << "Expansion ratio " << xr << " out of range."
                           << " Expansion ratio must be > 1."
                           << " In dictionary " << nl << tDict
                           << exit(FatalError);
            }
            scalar rlnr = 1/::log(xr);
            scalar delta = readScalar(tDict.lookup("delta"));

            Switch sym = tDict.lookupOrDefault<Switch>("symmetric", true);
            scalar Xf2 = finL/delta;
            scalar Xi2 = initL/delta;
            if (sym)
            {
                Xf2 /= 2;
                Xi2 /= 2;
            }

            //check expansion ratio is large enough to fill the space
            scalar maxXf2 = (::pow(xr, Xi2) - 1)*rlnr;
            if (maxXf2 < Xf2)
            {
                FatalError << "Expansion ratio " << xr << " is too small"
                           << " to stretch the mesh by the specified amount."
                           << " In dictionary " << nl << tDict
                           << exit(FatalError);
            }

            //calculate the length of the expanding region
            scalar Xi1 = 1;
            {
                scalar maxErr = 1e-10;
                scalar err = GREAT;
                label maxIter = 1000;
                label iter = 0;

                do
                {
                    scalar Xi1old = Xi1;

                    Xi1 = rlnr*::log((Xf2+rlnr)/(rlnr+Xi2-Xi1));

                    err = 2*mag(Xi1 - Xi1old)/(Xi1 + Xi1old);
                }
                while (err > maxErr && ++iter < maxIter);

                if (iter == maxIter)
                {
                    FatalError << "Could not calculate expansion region length."
                               << " Using dictionary " << nl << tDict
                               << exit(FatalError);
                }

                Info<< "Expansion section length = " << Xi1*delta << endl;
            }

            scalar xrPowXi1 = ::pow(xr, Xi1);
            scalar Xf1 = rlnr*(xrPowXi1 - 1);

            //convert mesh points into relative frame
            scalarField x( ((points - basePoint) & stretchV)/delta );
            scalarField xi = x;


            forAll(x, j)
            {
                //expansion section
                if (x[j] > 0)
                {
                    if (x[j] <= Xi1)
                    {
                        x[j] = rlnr*(::pow(xr, x[j]) -1);
                    }
                    else if (x[j] <= Xi2)
                    {
                        x[j] = Xf1 + xrPowXi1*(x[j] - Xi1);
                    }
                    else if (sym)
                    {
                        //transform coord to start from other side
                        scalar x_ = 2*Xi2 - x[j];

                        if (x_ >= Xi1)
                        {
                            x_ = Xf1 + xrPowXi1*(x_ - Xi1);
                        }
                        else if (x_ >= 0)
                        {
                            x_ = rlnr*(::pow(xr, x_) -1);
                        }

                        //transform back and assign
                        x[j] = 2*Xf2 - x_;
                    }
                    //displace any beyond the stretch region to make space
                    else
                    {
                        x[j] += Xf2 - Xi2;
                    }
                }
            }

            //modify points
            points += (x - xi)*delta*stretchV;

        }

        // Set the precision of the points data to 10
        IOstream::defaultPrecision(10);

        Info<< "Writing points into directory " << points.path() << nl << endl;
        points.write();
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
