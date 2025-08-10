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
    curvature

SourceFiles
    curvature.C

Authors
    Daniel Rettenmaier < rettenmaier@gsc.tu-darmstadt.de>
    Daniel Deising     < deising@mma.tu-darmstadt.de>
    All rights reserved.

Description

    You may refer to this software as :
    //- full bibliographic data to be provided

    This code has been developed by :
        Daniel Rettenmaier <rettenmaier@gsc.tu-darmstadt.de> (main developer).

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


#include "curvature/curvature.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(curvature, 0);
defineRunTimeSelectionTable(curvature, dictionary);

// * * * * * * * * * * * *  Protected Member Function * * * * * * * * * * * //

void curvature::initialize()
{

    //- determine cells attached to symmetryPlanes --> curvature will not be
    //  smoothed there /depending of the setting of smoothSym_)
    if (!smoothSym_)
    {
        forAll(mesh_.boundaryMesh(), iPatch)
        {
            if ((alpha1_.boundaryField()[iPatch].type() == "symmetryPlane"))
            {
                forAll(alpha1_.boundaryField()[iPatch], iFace)
                {
                    label own = mesh_.faceOwner()
                                [iFace + mesh_.boundaryMesh()[iPatch].start()];
                    attachedToSym_[own] = true;
                }
            }
        }

        //- check also for rotational symmetry (no faces there!)
        forAll(mesh_.cells(), iCell)
        {
            cell faces = mesh_.cells()[iCell];
            //- not valid for unstructured wedge meshes!!!
            if (faces.size() == 5)
            {
                int nWedges = 0;
                forAll(faces, i)
                {
                    label iPatch = mesh_.boundaryMesh().whichPatch(faces[i]);
                    if (iPatch != -1)
                    {
                        if (mesh_.boundaryMesh()[iPatch].type() == "wedge")
                        {
                            nWedges++;
                        }
                    }
                }

                if (nWedges == 2)
                {
                    attachedToSym_[iCell] = true;
                }
            }
        }
    }

    if (!smoothWall_)
    {
        forAll(mesh_.boundaryMesh(), iPatch)
        {
            if (isWallPatch_[iPatch])
            {
                forAll(alpha1_.boundaryField()[iPatch], iFace)
                {
                    label own = mesh_.faceOwner()
                                [iFace + mesh_.boundaryMesh()[iPatch].start()];
                    attachedToWall_[own] = true;

                    /*-DR
                                        //also prevent smoothing on neighbors
                                        labelList cellCells = mesh_.cellCells()[own];
                                        forAll(cellCells, i)
                                        {
                                            label neigh = cellCells[i];
                                            attachedToWall_[neigh] = true;
                                        }
                    */
                }
            }
        }
    }
}

/*
    Smoothens the kappa field by averaging the curvature in two loops hereby
    the neighbors of the cells with the variable set[iCell] == true are
    taken into account.
 */
void curvature::smooth(const boolList& set, volScalarField& kappa,  int nLoops)
{
    //- smooth curvature field
    for (int loop = 0; loop < nLoops; loop++)
    {
        scalarList neighK (mesh_.nCells(), 0.0);
        scalarList count (mesh_.nCells(), 0.0);

        //- sum up kappa of neighboring cells
        forAll(kappa, iCell)
        {
            if (set[iCell] == true)
            {
                labelList cellCells = mesh_.cellCells()[iCell];
                forAll(cellCells, i)
                {
                    label neigh = cellCells[i];
                    if (set[neigh] == true)
                    {
                        neighK[iCell] += kappa[neigh];
                        count[iCell] += 1.0;
                    }
                }
            }
        }

        //- take into account coupled boundaries
        scalarList KProc (mesh_.nFaces(), 0.0);
        scalarList countProc (mesh_.nFaces(), 0.0);
        forAll(kappa.boundaryField(), iPatch)
        {
            if (
                (alpha1_.boundaryField()[iPatch].type() == "processor")
                || (alpha1_.boundaryField()[iPatch].type() == "cyclic")
            )
            {
                forAll(kappa.boundaryField()[iPatch], iFace)
                {
                    label iFaceGlobal = iFace + mesh_.boundaryMesh()[iPatch].start();
                    label own = mesh_.faceOwner()[iFaceGlobal];
                    if (set[own] == true)
                    {
                        KProc[iFaceGlobal] = kappa[own];
                        countProc[iFaceGlobal] += 1.0;
                    }
                }
            }
        }

        syncTools::syncFaceList(mesh_, KProc, plusEqOp<scalar>());
        syncTools::syncFaceList(mesh_, countProc, plusEqOp<scalar>());

        forAll(kappa.boundaryField(), iPatch)
        {
            if (
                (alpha1_.boundaryField()[iPatch].type() == "processor")
                || (alpha1_.boundaryField()[iPatch].type() == "cyclic")
            )
            {
                forAll(kappa.boundaryField()[iPatch], iFace)
                {
                    label iFaceGlobal = iFace + mesh_.boundaryMesh()[iPatch].start();
                    label own = mesh_.faceOwner()[iFaceGlobal];
                    if ((set[own] == true) && (countProc[iFaceGlobal] == 2.0))
                    {
                        //- the value of the owner has to be subtracted to get
                        //  only the value of the neighbor
                        neighK[own] += KProc[iFaceGlobal] - kappa[own];
                        count[own] += 1.0;
                    }
                }
            }
        }

        //- average kappa using neighbor values
        forAll(kappa, iCell)
        {
            if ((set[iCell]) && (!attachedToSym_[iCell]) && (!attachedToWall_[iCell])) //&& (interfaceDensity_[iCell] != 0.0)
            {
                kappa[iCell] =
                    (neighK[iCell] + kappa[iCell]) / (1.0 + count[iCell]);
            }
        }
    }
}

//- Calculate the curvature with an optional correction term
void curvature::calcK()
{
    kappa_ = -fvc::div(nHatfv_ & mesh_.Sf());

    // Complex expression for curvature.
    // Correction is formally zero but numerically non-zero.
    if (complexDivCorrection_)
    {
        kappa_ += (nHatv_ & fvc::grad(nHatfv_) & nHatv_);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

curvature::curvature
(
    const word& name,
    const volScalarField& alpha,
    const volVectorField& nHatv,
    const surfaceVectorField& nHatfv,
    const volScalarField& interfaceDensity,
    const List<bool>& isWallPatch,
    const dictionary& transpProp,
    const volScalarField& isInterface
)
    :
    kappa_ //volScalarField
    (
       IOobject
       (
           "kappa",
           alpha.mesh().time().timeName(),
           alpha.mesh(),
           IOobject::NO_READ,
           IOobject::AUTO_WRITE
       ),
       alpha.mesh(),
       dimensionedScalar("kappa", dimensionSet(0, -1, 0, 0, 0, 0, 0), 0)
    ),
    name_(name),
    mesh_(alpha.mesh()),
    alpha1_(alpha),
    nHatv_(nHatv),
    nHatfv_(nHatfv),
    interfaceDensity_(interfaceDensity),
    isWallPatch_(isWallPatch),
    isInterface_(isInterface),

    smoothSym_(transpProp.subDict("curvature").lookup("smoothSym")),

    attachedToSym_(mesh_.nCells(), false),

    smoothWall_(transpProp.subDict("curvature").lookup("smoothWall")),

    attachedToWall_(mesh_.nCells(), false),

    preSmoothCycles_(readScalar(transpProp.subDict("curvature").lookup("preSmoothCycles"))),

    postSmoothCycles_(readScalar(transpProp.subDict("curvature").lookup("postSmoothCycles"))),

    complexDivCorrection_(transpProp.subDict("curvature").lookupOrDefault<Switch>("complexDivCorrection", false))
{

    initialize();

    distributeField distributeField_;
}

/*
    Reduces kappa to the interfaceDensity > 0 field.
    Smoothes the reduced field. This smothened field is then propagated
    to all isInterface cells and smoothened again.
    Use only recommended for reconstructed fields!
*/
void curvature::prePostSmooth
(
    volScalarField& kappa
)
{
    //- mesh size might have changed, adapt the size of the marker-fields
    if (mesh_.changing())
    {
        attachedToSym_ = List<bool>(mesh_.nCells(), false);
        attachedToWall_= List<bool>(mesh_.nCells(), false);

        initialize();
    }

    //- prepare smoothing list
    boolList isInterfaceList(mesh_.nCells(), false);
    boolList set(mesh_.nCells(), true);

    forAll(set, iCell)
    {
        isInterfaceList[iCell] = isInterface_[iCell];
        if (interfaceDensity_[iCell] == 0.0)
        {
            kappa[iCell] = 0.0; //keep kappa only in interface cells
            set[iCell] = false;
        }
    }

    volScalarField sharpInterface(isInterface_*0);
    forAll(sharpInterface, cellI)
    {
        sharpInterface[cellI] = pos0(interfaceDensity_[cellI] - VSMALL);
    }

    //- pre-smooth the curvature
    smooth(set, kappa, preSmoothCycles_);              //TODO leads to errors at wall patches
                                        //     however necessary for good behavior in static drop test

    //- distributes kappaSet and kappa on isInterface stencil
    distributeField_.distributeVolField
    (
        kappa,           //field to distribute
        sharpInterface,  //list numbered as field with current distribution status
        isInterface_,    //list numbered as field with width of the destination field
        40,              //maxLoops
        3,               //minLoops
        false            //normalize field
    );
    //- distributes kappaSet and kappa on isInterface stencil
    //distributeField_.distributeVolField
    //(
    //    kappa,          //field to distribute
    //    set,            //current distribution state (bool)
    //    isInterfaceList,//target marker field
    //    10,             //maximum loops
    //    0,              //minimum loops
    //    false           //normalize field
    //);

    //- post-smooth the curvature
    smooth(set, kappa, postSmoothCycles_);

    //- update the boundaryField of kappa
    forAll(kappa.boundaryField(), iPatch)
    {
        if (
            (kappa.boundaryField()[iPatch].type() != "processor")
            && (kappa.boundaryField()[iPatch].type() != "cyclic")
        )
        {
            forAll(kappa.boundaryField()[iPatch], iFace)
            {
                label iFaceGlobal = iFace + mesh_.boundaryMesh()[iPatch].start();
                label own = mesh_.faceOwner()[iFaceGlobal];
                kappa.boundaryFieldRef()[iPatch][iFace] = kappa[own];
            }
        }

        if (isWallPatch_[iPatch])
        {
            forAll(kappa.boundaryField()[iPatch], iFace)
            {
                kappa.boundaryFieldRef()[iPatch][iFace] = 0.0;    //WHY?
            }
        }
    }
    kappa.correctBoundaryConditions();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Ostream& operator<<(Ostream& os, curvature& curv)
{
    curv.write(os);
    return os;
}

void curvature::write(Ostream& os) const
{
    volScalarField test = kappa_; //*this;
    os << test;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
