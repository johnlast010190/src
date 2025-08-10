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

Application
    meshUpdate

Description
    executes mesh.update() based on constant/dynamicMeshDict.
    Specialized for reading/executing AMR (dynamic(-MultiCriterion-)RefineFvMesh)

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/include/fvCFD.H"
#include "dynamicFvMesh/dynamicFvMesh.H"
#include "dynamicRefineFvMesh/dynamicRefineFvMesh.H"

#include "singlePhaseTransportModel/singlePhaseTransportModel.H"
#include "turbulentTransportModels/turbulentTransportModel.H"

#include "fields/ReadFields/ReadFields.H"
#include "db/IOobjectList/IOobjectList.H"

#include "cfdTools/general/CorrectPhi/CorrectPhi.H"

#include <iostream>
#include <fstream>

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // add arguments to function call:
    timeSelector::addOptions(true, false);
#include "include/addOverwriteOption.H"
#include "include/addRegionOption.H"

    #include "include/setRootCase.H"
    #include "include/createTime.H"
    #include "include/createDynamicFvMesh.H"

    IOobjectList objects(mesh, runTime.timeName());

    // Read vol fields.
    PtrList<volScalarField> vsFlds;
    ReadFields(mesh, objects, vsFlds);

    PtrList<volVectorField> vvFlds;
    ReadFields(mesh, objects, vvFlds);

    // Read surface fields.
    PtrList<surfaceScalarField> ssFlds;
    ReadFields(mesh, objects, ssFlds);

    PtrList<surfaceVectorField> svFlds;
    ReadFields(mesh, objects, svFlds);

    // get reference to velocity field
    volVectorField& U(const_cast<volVectorField&>(mesh.lookupObject<volVectorField>("U")));

    // create flux field
    surfaceScalarField& phi(const_cast<surfaceScalarField&>(mesh.lookupObject<surfaceScalarField>("phi")));
    //#include "cfdTools/incompressible/createPhi.H"

    phi = fvc::flux(U);

    singlePhaseTransportModel laminarTransport(U, phi);
    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, laminarTransport)
    );

    IOdictionary refineDict
    (
        IOobject
        (
            "dynamicMeshDict",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    label refineInterval = refineDict.subDict("dynamicRefineFvMeshCoeffs").lookupOrDefault<label>("refineInterval", 1);
    label maxRefinement = refineDict.subDict("dynamicRefineFvMeshCoeffs").lookupOrDefault<label>("maxRefinement", 1);

    #include "correctPhi.H"

    for (label i=0; i < (refineInterval * maxRefinement) + ceil(maxRefinement/2); i++)
    {
        runTime++;
        Info<< "updating mesh for time " << runTime.timeName() << endl;
        mesh.update();
        Info<< "finished mesh update.  " << endl;
    }

    runTime.writeNow();
    return(0);
}


// ************************************************************************* //
