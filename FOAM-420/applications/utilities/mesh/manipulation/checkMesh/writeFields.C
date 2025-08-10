#include "writeFields.H"
#include "fields/volFields/volFields.H"
#include "meshes/polyMesh/polyMeshCheck/polyMeshTools.H"
#include "fields/fvPatchFields/basic/zeroGradient/zeroGradientFvPatchFields.H"
#include "meshes/polyMesh/syncTools/syncTools.H"

using namespace Foam;

void maxFaceToCell
(
    const scalarField& faceData,
    volScalarField& cellData
)
{
    const cellList& cells = cellData.mesh().cells();

    scalarField& cellFld = cellData.ref();

    cellFld = -GREAT;
    forAll(cells, cellI)
    {
        const cell& cFaces = cells[cellI];
        forAll(cFaces, i)
        {
            cellFld[cellI] = max(cellFld[cellI], faceData[cFaces[i]]);
        }
    }

    forAll(cellData.boundaryField(), patchI)
    {
        fvPatchScalarField& fvp = cellData.boundaryFieldRef()[patchI];

        fvp = fvp.patch().patchSlice(faceData);
    }
    cellData.correctBoundaryConditions();
}


void minFaceToCell
(
    const scalarField& faceData,
    volScalarField& cellData
)
{
    const cellList& cells = cellData.mesh().cells();

    scalarField& cellFld = cellData.ref();

    cellFld = GREAT;
    forAll(cells, cellI)
    {
        const cell& cFaces = cells[cellI];
        forAll(cFaces, i)
        {
            cellFld[cellI] = min(cellFld[cellI], faceData[cFaces[i]]);
        }
    }

    forAll(cellData.boundaryField(), patchI)
    {
        fvPatchScalarField& fvp = cellData.boundaryFieldRef()[patchI];

        fvp = fvp.patch().patchSlice(faceData);
    }
    cellData.correctBoundaryConditions();
}


void Foam::writeFields
(
    const fvMesh& mesh,
    const HashSet<word>& selectedFields
)
{
    if (selectedFields.empty())
    {
        return;
    }

    Info<< "Writing fields " << endl;

    // cell type
    // ~~~~~~~~~
    if (selectedFields.found("cellShapes"))
    {
        volScalarField shape
        (
            IOobject
            (
                "cellShapes",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            mesh,
            dimensionedScalar("cellShapes", dimless, 0),
            zeroGradientFvPatchScalarField::typeName
        );
        const cellShapeList& cellShapes = mesh.cellShapes();
        forAll(cellShapes, cellI)
        {
            const cellModel& model = cellShapes[cellI].model();
            shape[cellI] = model.index();
        }
        shape.correctBoundaryConditions();
        Info<< "    Writing cell shape (hex, tet etc.) to " << shape.name()
            << endl;
        shape.write();
    }

    if (selectedFields.found("cellVolume"))
    {
        volScalarField V
        (
            IOobject
            (
                "cellVolume",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            mesh,
            dimensionedScalar("cellVolume", dimVolume, 0),
            calculatedFvPatchScalarField::typeName
        );
        V.ref() = mesh.V();
        Info<< "    Writing cell volume to " << V.name() << endl;
        V.write();
    }

    if (selectedFields.found("cellCentres"))
    {
        volVectorField C
        (
            IOobject
            (
                "cellCentres",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            mesh,
            dimensionedVector("cellCentres", dimLength, vector::zero),
            calculatedFvPatchVectorField::typeName
        );
        C.ref() = mesh.C();
        Info<< "    Writing cell centres to " << C.name() << endl;
        C.write();
    }

    Info<< endl;
}
