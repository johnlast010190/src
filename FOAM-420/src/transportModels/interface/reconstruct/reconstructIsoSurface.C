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
    (c) held by original author
    (c) 2022 Esi Ltd.

Class
    reconstruct

SourceFiles
    reconstruct.C

Authors
    Daniel Rettenmaier < rettenmaier@gsc.tu-darmstadt.de>
    Daniel Deising     < deising@mma.tu-darmstadt.de>
    All rights reserved.

Description

    You may refer to this software as :
    //- full bibliographic data to be provided

    This code has been developed by :
        Daniel Rettenmaier < rettenmaier@gsc.tu-darmstadt.de> (main developer).

    Method Development and Intellectual Property :
        Daniel Rettenmaier <rettenmaier@gsc.tu-darmstadt.de>
        Daniel Deising <deising@mma.tu-darmstadt.de>
        Holger Marschall <marschall@csi.tu-darmstadt.de>
        Dieter Bothe <bothe@csi.tu-darmstadt.de>
        Cameron Tropea <ctropea@sla.tu-darmstadt.de>

        Mathematical Modeling and Analysis
        Institute for Fluid Mechanics and Aerodynamics
        Center of Smart Interfaces
        Technische Universitaet Darmstadt

    If you use this software for your scientific work or your publications,
    please don't forget to acknowledge explicitly the use of it.

\*---------------------------------------------------------------------------*/


#include "reconstruct/reconstructIsoSurface.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug( reconstructIsoSurface, 0 );
addToRunTimeSelectionTable( reconstruct, reconstructIsoSurface, dictionary );

//- Calculates the 0.5 isoSurface using the alphaP field
//  Sets the distance of the 0.5 isoSurface to the cell center of cells
//  containing the interface
void reconstructIsoSurface::interfacePosition()
{
    interfacePoints_.clear();

    // ensure large enough space to prevent copying
    interfacePoints_.resize(mesh_.nEdges());
    label idx = 0; // index to fill interfacePoints_ list
    label listSize = 0; // to resize interfacePoints_ list

    //- Determine interface position and calculate normals
    forAll(mesh_.cells(), iCell)
    {
        cell faces = mesh_.cells()[iCell];
        int startFace = -1;

        //- count the number of cuts of faces with interface
        int cutFaces = 0;
        forAll(faces, j)
        {
            label iFace = faces[j];

            const face& f = mesh_.faces()[iFace];
            labelList points(f.size(), Zero);
            points = f;

            const label pointsSize(points.size());

            label own = mesh_.faceOwner()[iFace];
            if (own != iCell)
            {
                labelList pointsTmp(points);

                for (int i = 0; i < pointsSize; i++)
                {
                    points[i] = pointsTmp[pointsSize-1-i];
                }
            }
            bool cutFound = 0;
            int i = 0;
            while ((!cutFound) && (i < pointsSize))
            {
                label ip1 = (i+1) % pointsSize;
                scalar di = alphaP_[points[i]];
                scalar dip1 = alphaP_[points[ip1]];

                //- avoid points having alpha == 0.5, leads to a non-unique isosurface
                if (mag(di-0.5) < 1e-9)
                {
                    di += 1e-8;                 //Interface moves slightly
                }
                if (mag(dip1-0.5) < 1e-9)
                {
                    dip1 += 1e-8;
                }

                if ((di < 0.5) && (dip1 > 0.5))
                {
                    cutFound = 1;
                    if (startFace == -1)
                    {
                        startFace = j;
                    }
                }
                i += 1;
            }
            if (cutFound)
            {
                cutFaces += 1;
            }
        }
        //- if there are cuts found, find the according points
        DynamicList<vector> xCuttingPlane(0);
        xCuttingPlane.reserve(faces.size());

        label nextFace = startFace;
        int N = 0;
        while (N < cutFaces)
        {
            label iFace = faces[nextFace];
            const face& f = mesh_.faces()[iFace];
            labelList points(f.size(), Zero);
            points = f;

            const label pointsSize(points.size());

            label own = mesh_.faceOwner()[iFace];
            if (own != iCell)
            {
                labelList pointsTmp(points);
                for (int i = 0; i < pointsSize; i++)
                {
                    points[i] = pointsTmp[pointsSize-1-i];
                }
            }
            for (int i = 0; i < pointsSize; i++)
            {
                label ip1 = (i+1) % pointsSize;
                vector xi = mesh_.points()[points[i]];
                vector xip1 = mesh_.points()[points[ip1]];
                scalar di = alphaP_[points[i]];
                scalar dip1 = alphaP_[points[ip1]];

                //- avoid points having alpha == 0.5, leads to a non-unique isosurface
                if (mag(di-0.5) < 1e-9)
                {
                    di += 1e-8;
                }
                if (mag(dip1-0.5) < 1e-9)
                {
                    dip1 += 1e-8;
                }

                if ((di < 0.5) && (dip1 > 0.5))
                {
                    N += 1;
                    vector xNew = xi+(0.5-di)/(dip1-di)*(xip1-xi);
                    bool exists = 0;
                    const label xCuttingPlaneSize(xCuttingPlane.size());
                    for (int n = 0; n < xCuttingPlaneSize; n++)
                    {
                        if (xNew == xCuttingPlane[n])
                        {
                            exists = 1;
                        }
                    }
                    if (!exists)
                    {
                        xCuttingPlane.append(xNew);
                        interfacePoints_[idx] = xNew;
                        idx++;
                    }
                    //- find the next face to be visited, i.e. the one sharing
                    //  the edge at which the previously found cut is
                    bool found = 0;
                    int j2 = -1;
                    while (!found)
                    {
                        j2 += 1;
                        label iFace2 = faces[j2];

                        const face& f2 = mesh_.faces()[iFace2];
                        labelList points2(f2.size(), Zero);
                        points2 = f2;

                        const label points2Size(points2.size());


                        label own2 = mesh_.faceOwner()[iFace2];
                        if (own2 != iCell)
                        {
                            labelList pointsTmp(points2);
                            for (int i = 0; i < points2Size; i++)
                            {
                                points2[i] = pointsTmp[points2Size-1-i];
                            }
                        }
                        for (int i2 = 0; i2 < points2Size; i2++)
                        {
                            label i2p1 = (i2+1) % points2Size;
                            if ((points2[i2] == points[ip1]) && (points2[i2p1] == points[i]))
                            {
                                nextFace = j2;
                                found = 1;
                            }
                        }
                    }
                }
            }

        }

        // add to global list size
        listSize += xCuttingPlane.size();

        //- calculate area vector of reconstructed surface
        vector S (0, 0, 0);
        interfaceDensity_[iCell] = 0.0;
        nHatv_[iCell] = S;
        interfaceCenters_[iCell] = vector(0,0,0);

        //- create a face object to use its capability of calculating the centre
        face pseudoFace (xCuttingPlane.size());
        pointField pseudoPoints (xCuttingPlane.size(), vector(0,0,0));

        for (int i = 0; i < xCuttingPlane.size(); i++)
        {
            label ip1 = i+1;
            if (ip1 == xCuttingPlane.size())
            {
                ip1 = 0;
            }
            vector xi = xCuttingPlane[i];
            vector xip1 = xCuttingPlane[ip1];
            S += 0.5*(xi^xip1);
            pseudoFace[i] = i;
            pseudoPoints[i] = xCuttingPlane[i];
        }
        scalar magS = mag(S);
        if (magS/mesh_.V()[iCell] > SMALL)
        {
            nHatv_[iCell] = S/magS;
            interfaceCenters_[iCell] = pseudoFace.centre(pseudoPoints);
            interfaceDensity_[iCell] = magS/mesh_.V()[iCell];
        }
//      interfaceDensity_[iCell] = magS/mesh_.V()[iCell];

        //- calculate plane offset and distance
        scalar lambda = 0.0;
        if (xCuttingPlane.size() > 0)
        {
            for (int i = 0; i < xCuttingPlane.size(); i++)
            {
                lambda += nHatv_[iCell] & xCuttingPlane[i];
            }
            lambda /= xCuttingPlane.size();
            distance_[iCell] = (nHatv_[iCell] & mesh_.C()[iCell]) - lambda;
        }
    }

    // resize list back to correct size
    interfacePoints_.resize(listSize);
}

/*
 *  Propagates nHat and if selcected, the distance of cellPoints to the interface
    is calculated by passing the interface information, namely the interface
    normal and the current distance to neighborcells, within the isInterface cells
    or the user selected minimum distance to the interface
 */
void reconstructIsoSurface::propagateNHatv()
{
    boolList set(mesh_.nCells(), false);
    boolList isInterfaceList(mesh_.nCells(), false);

    forAll(set, iCell)
    {
        if (interfaceDensity_[iCell] > VSMALL)
        {
            set[iCell] = true;
        }
        isInterfaceList[iCell] = isInterface_[iCell];
    }

    volScalarField sharpInterface("sharpInterface", isInterface_*0);
    forAll(sharpInterface, cellI)
    {
        sharpInterface[cellI] = pos0(interfaceDensity_[cellI] - VSMALL);
    }
    volScalarField distributeRegion(isInterface_);

    if (!calcSignedDistance_)
    {
        distributeField_.distributeVolField
        (
                nHatv_,          //field to distribute
                sharpInterface,  //list numbered as field with current distribution status
                isInterface_,    //list numbered as field with width of the destination field
                40,              //maxLoops
                3,               //minLoops
                true             //normalizeField
        );

/*
        distributeField_.distributeVolField
        (
                nHatv_,         //field to distribute
                set,            //list numbered as field with current distribution status
                isInterfaceList,//list numbered as field with width of the destination field
                20,             //maxLoops
                3,              //minLoops
                true            //normalizeField
        );
*/

    } else {

        //- initialize distance field in the non-interface cells
        //  in cells not containing interface
        //     +1 --> liquid side
        //     -1 --> vapor side
        forAll(distance_, iCell)
        {
            if (interfaceDensity_[iCell] < VSMALL)
            {
                distance_[iCell] = 2.0*alpha1_[iCell]-1.0;
            }
        }
        volScalarField::Boundary& distancebf = distance_.boundaryFieldRef();
        forAll(distancebf, iPatch)
        {
            forAll(distancebf[iPatch], iFace)
            {
                label own = mesh_.faceOwner()[iFace + mesh_.boundaryMesh()[iPatch].start()];
                distancebf[iPatch][iFace] = distance_[own];
            }
        }

        distributeField_.calculateDistance
        (
            set,                //list numbered as field holding distance
            distance_,          //information distance to interface
            nHatv_,             //nHat used for distance calculation
            isInterfaceList,    //minimum field to cover
            distanceThreshold_, //up to how far from the interface
            20,                //maxLoops
            3                   //minLoop
        );
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

reconstructIsoSurface::reconstructIsoSurface
(
    const word& name,
    const volScalarField& alpha,
    const dictionary& transpProp,
    const List<bool>& isWallPatch,
    const volScalarField& isInterface
)
:
    reconstruct( name, alpha, transpProp, isWallPatch, isInterface ),

    mesh_(alpha.mesh()),

    distance_
    (
        IOobject
        (
            "distance",
            alpha.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("distance", dimensionSet(0,1,0,0,0,0,0), 0.0)
    ),

    interfacePoints_
    (
        IOobject
        (
            "interfacePoints",
            alpha.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_.boundaryMesh().size()
    ),

    interfaceCenters_
    (
        IOobject
        (
            "interfaceCenters",
            alpha.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("interfaceCenters", dimensionSet(0,1,0,0,0,0,0), vector(0,0,0))
    ),

    distanceThreshold_(readScalar(
        transpProp.subDict("reconstruct").lookup("distanceThreshold"))),

    calcSignedDistance_(transpProp.subDict("reconstruct").lookup("calcSignedDistance")),
/*
    alphaInterpolator_
    (
        alphaInterpolation::New

        (   transpProp.subDict("reconstruct").lookup("alphaInterpolationMethod"),
            alpha1_,
            isWallPatch_
        )
    ),
*/
    //alphaP_(alphaInterpolator_->alphaP())
    alphaP_(mesh_.nPoints(), 0.0)

{
    if (calcSignedDistance_)
    {
        distance_.writeOpt() = IOobject::AUTO_WRITE;
    }

    distributeField_ = distributeField();

    //DR
    //This line magically solves a bug in line ~260: distance_[iCell] = (nHatv_[iCell] & mesh_.C()[iCell]) - lambda;
    //Which only occurs in some parallel configurations with dynamic refined mesh.
    //Not all processors reach this line, but when they do they get stuck calling mesh_.C()
    //Processors that skip the line will wait at the next mpi-barrier which is misleading while debugging.
    //If you understand why please write me a mail :)
    mesh_.C();

    reconstructInterface();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void reconstructIsoSurface::reconstructInterface()
{
    //original    const surfaceScalarField alphaSnGrad = fvc::snGrad(alpha1_);
    //original #include correctBoundaries.H

    //- interpolate alpha from cells to points
    //alphaInterpolator_->interpolate();

    const volPointInterpolation& pInterp = volPointInterpolation::New(mesh_);
    pointScalarField pAlpha(pInterp.interpolate(alpha1_));

    forAll(pAlpha, pointI)
    {
        alphaP_[pointI] = pAlpha[pointI];
    }
//Info<< "alphaP_.size(): " << alphaP_.size() << endl;
//Info<< "pAlpha.size(): " << pAlpha.size() << endl;
/*
//GeometricField<Type, pointPatchField, pointMesh>& pf
//const polyBoundaryMesh& bm = mesh.boundaryMesh();
//forAll(bm, patchi)
//{
//    const labelList& bp = bm[patchi].boundaryPoints();
//    const labelList& meshPoints = bm[patchi].meshPoints();
//
//    forAll(bp, pointi)
//    {
//        label ppp = meshPoints[bp[pointi]];
//        Info<<boundary value of pf is: " << pf[ppp] << endl;;
//    }
//}
    forAll(pAlpha.boundaryField(), patchI)
    {
        //const scalarField pAlphaP = pAlpha.boundaryField()[patchI];
        label startPoint = mesh_.boundary()[patchI].start();

        forAll(pAlpha.boundaryField()[patchI], pointI)
        {
            alphaP_[pointI+startPoint] = pAlpha.boundaryField()[patchI][pointI];
        }
    }
*/

    //original alphaInterpolator_->correctWallFaces(alphaSnGrad);

    //- calculate nHatv interfaceDensity and interfacePoints
    //  alphaP -> edgeCutting at isoValue 0.5 of interpolated alphaP
    //  -> cuttingPolygon -> surfaceArea and nHatv
    if (debug)
    {
        Info<< "computing interface position" << endl;
    }
    interfacePosition();

    //- transport nHatv away from the interface and calculates the distance
    //  of the cell points to the interface if selected by the user
    if (debug)
    {
        Info<< "propagating nHat" << endl;
    }
    propagateNHatv();

    //- interpolate the transported nHatv_ field
    nHatfv_ = fvc::interpolate(nHatv_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
