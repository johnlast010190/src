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
    (c) 2017 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "dynamicGIBFvMesh/dynamicGIBFvMesh/dynamicGIBFvMesh.H"
#include "meshTools/meshTools.H"
#include "motionSmoother/motionSmoother.H"
#include "sets/topoSets/faceSet.H"
#include "meshes/meshTools/simpleVTKWriter.H"
#include "fields/fvsPatchFields/fvsPatchField/fvsPatchField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::dynamicGIBFvMesh::writeProblematicCells
(
    label cI
) const
{
    if (debugMode_)
    {
        Pout<< "In cell" << tab << cI << tab << "at location"
             << tab << this->C()[cI] << endl;

        const labelList& cellsI = this->cells()[cI];

        Pout<< "cell faces:" << endl;

        forAll(cellsI, fI)
        {
            const label& faceI = cellsI[fI];
            if (faceI<this->nInternalFaces())
            {
                Pout<< fI << tab << faceI << tab << this->Cf()[faceI] << endl;
            }
            else
            {
                label patchI = this->boundaryMesh().whichPatch(faceI);
                label lfI = faceI - this->boundaryMesh()[patchI].start();
                Pout<< fI << tab << faceI << tab
                     << patchI << tab << this->boundary()[patchI].name() << tab
                     << this->Cf().boundaryField()[patchI][lfI] << endl;
            }
        }

        OFstream str
        (
            this->time().path()/
            "cells_zeroOldVol_"+this->time().timeName()+".obj"
        );

        meshTools::writeOBJ
        (
            str,
            this->cells(),
            this->faces(),
            this->points(),
            labelList(1, cI)
        );
    }
}

void Foam::dynamicGIBFvMesh::writeProblematicCells
(
    label cI,
    const surfaceScalarField& ssf
) const
{
    if (debugMode_)
    {
        Pout<< "In cell" << tab << cI << tab << "at location"
             << tab << this->C()[cI] << endl;

        const labelList& cellsI = this->cells()[cI];

        Pout<< "cell faces:" << endl;

        forAll(cellsI, fI)
        {
            const label& faceI = cellsI[fI];
            if (faceI<this->nInternalFaces())
            {
                Pout<< fI << tab << faceI << tab
                     << this->Cf()[faceI] << tab
                     << ssf[faceI] << endl;
            }
            else
            {
                label patchI = this->boundaryMesh().whichPatch(faceI);
                Pout<< fI << tab << faceI << tab
                     << patchI << tab << this->boundary()[patchI].name() <<endl;
            }
        }

        OFstream str
        (
            this->time().path()/
            "cells_zeroOldVol_"+this->time().timeName()+".obj"
        );

        meshTools::writeOBJ
        (
            str,
            this->cells(),
            this->faces(),
            this->points(),
            labelList(1, cI)
        );
    }
}



void Foam::dynamicGIBFvMesh::writeProblematicCells() const
{
    if (debugMode_)
    {
        DynamicList<label> dlV0(cells().size());
        DynamicList<label> dlV(cells().size());
        const scalarField& cV0 = this->V0();
        const scalarField& cV = this->V();
        forAll(this->cells(), cI)
        {
            if (cV0[cI]<=0)
            {
                dlV0.append(cI);
                Pout<< "cell " << cI << tab << "negative V0 found" << endl;
            }
            if (cV[cI]<=0)
            {
                dlV.append(cI);
                Pout<< "cell " << cI << tab << "negative V found" << endl;
            }
        }
        dlV0.shrink();
        dlV.shrink();

        const labelList lV0 = labelList(dlV0);
        const labelList lV = labelList(dlV);

        label sizeV0 = lV0.size();
        label sizeV = lV.size();

        Foam::reduce(
            std::tie(sizeV0, sizeV),
            ParallelOp<sumOp<label>, sumOp<label>>{},
            comm()
        );

        if (sizeV0 != 0)
        {
            OFstream str
            (
                this->time().path()/
                "probCells_V0_"+this->time().timeName()+".obj"
            );

            meshTools::writeOBJ
            (
                str,
                this->cells(),
                this->faces(),
                this->points(),
                lV0
            );
        }
        if (sizeV != 0)
        {
            OFstream str
            (
                this->time().path()/
                "probCells_V_"+this->time().timeName()+".obj"
            );

            meshTools::writeOBJ
            (
                str,
                this->cells(),
                this->faces(),
                this->points(),
                lV
            );
        }
    }
}


void Foam::dynamicGIBFvMesh::writeProblematicCells
(
    const boolList& interPoints
) const
{
    if (debugMode_)
    {
        const pointField& basePoints = *basePoints_;
//        const vectorField& baseCc = *baseCC_;
        DynamicList<label> dliP(cells().size());
        forAll(this->cells(), cI)
        {
            const labelList& cP = this->cellPoints()[cI];
            bool foundUnsnapped = false;
            forAll(cP, pI)
            {
                const label& cPI = cP[pI];
                if (!interPoints[cPI])
                {
                    foundUnsnapped = true;
                }
            }
            if (!foundUnsnapped)
            {
                dliP.append(cI);
           //     Pout<< "cell " << cI << tab << baseCc[cI] <<tab
          //           << "all points are snapped" << endl;
            }
        }
        dliP.shrink();

        const labelList liP = labelList(dliP);

        label sizeiP = liP.size();
        reduce(sizeiP, sumOp<label>());

        if (sizeiP != 0)
        {
            OFstream str
            (
                this->time().path()/
                "probCells_snapPoints_"+this->time().timeName()+".obj"
            );

            meshTools::writeOBJ
            (
                str,
                this->cells(),
                this->faces(),
                basePoints,
                liP
            );
        }
    }
}


void Foam::dynamicGIBFvMesh::checkRegion(const labelList& cR) const
{
    if (debugMode_)
    {
        bool regFound = false;
        forAll(cR, cI)
        {
            if
            (
                (cR[cI]!=0) &&
                (regFound==false)
            )
            {
                regFound = true;
            }
        }
        reduceToMaster(regFound, orOp<bool>());
        if (!regFound)
        {
            WarningInFunction
                << "in dynamicGIBFvMesh::modifyRegionLabels"
                << "(labelList& cellIndi): "
                << " different regions not found" << endl;
        }
    }
}


void Foam::dynamicGIBFvMesh::writeRegionInterface
(
    word name,
    const boolList& regFace
) const
{
    DynamicList<label> dregFaceL(regFace.size());
    const pointField& basePoints = *basePoints_;
    forAll(regFace, fI)
    {
        if (regFace[fI])
        {
            dregFaceL.append(fI);
        }
    }
    dregFaceL.shrink();
    labelList regFaceL(dregFaceL);
    simpleVTKWriter
    (
        this->faces(),
        regFaceL,
        basePoints
    ).write
    (
        name+"Base_"+this->time().timeName()+".vtk"
    );
    simpleVTKWriter
    (
        this->faces(),
        regFaceL,
        this->points()
    ).write
    (
        name+this->time().timeName()+".vtk"
    );
}

void Foam::dynamicGIBFvMesh::writeRegionInterface
(
    word name,
    const boolList& regFace,
    label cI
) const
{
    DynamicList<label> dregFaceL(regFace.size());
    const pointField& basePoints = *basePoints_;

    if (Pstream::myProcNo()==6)
    {
        const labelList& cFaces = this->cells()[cI];
        forAll(cFaces, cFI)
        {
            const label& cFacesI = cFaces[cFI];
            if (regFace[cFacesI])
            {
                dregFaceL.append(cFacesI);
            }
        }
    }

    dregFaceL.shrink();
    labelList regFaceL(dregFaceL);
    simpleVTKWriter
    (
        this->faces(),
        regFaceL,
        basePoints
    ).write
    (
        name+"Base_"+this->time().timeName()+".vtk"
    );
}



void Foam::dynamicGIBFvMesh::qualityControl() const
{
    const dictionary& dmCoeffs = dynamicMeshCoeffs_.
        subDict("meshQualityControls");

    faceSet errorFaces
    (
        *this,
        "errorFaces",
        this->nFaces()
    );

    bool hasErrors = motionSmoother::checkMesh
    (
        false,  // report
        *this,
        dmCoeffs,
        errorFaces
    );
    if (hasErrors)
    {
        Info<< "Mesh has errors" << endl;
    }
    errorFaces.write();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// ************************************************************************* //
