/*---------------------------------------------------------------------------*\
| Modified 2010-2016 Copyright (C) Esi Ltd                                  |
\*---------------------------------------------------------------------------*/

// Check that walls are tangential to flow direction
// count number of wall faces

label nWallFaces(0);
volScalarField y(wallDist::New(mesh).y());

const fvPatchList& patches = mesh.boundary();

// Some modifications to allow calculatiuon of 2D ducts
forAll(patches, patchi)
{
    const fvPatch& currPatch = patches[patchi];

    if (isA<wallFvPatch>(currPatch))
    {
        const vectorField nf(currPatch.nf());
        nWallFaces += currPatch.size();

        if (gMax(mag(nf & flowDirection)) > 0.01)
        {
            FatalErrorInFunction
                << "Detected misaligned wall face normal in patch "
                << currPatch.name()
                << ", normal must be orthogonal to flow direction."
                << nl << exit(FatalError);
        }
    }
}

if (nWallFaces == 0)
{
    FatalErrorInFunction
        << exit(FatalError);
}


// if only 1 or 2 wall faces are detected keep original behaviour
bool graphData = false;
if (nWallFaces == 1 || nWallFaces == 2 || !Pstream::parRun())
{
    graphData = true;
}

label faceId(-1);
label patchId(-1);
vector wallNormal12(Zero);

if (graphData)
{
    forAll(patches, patchi)
    {
        const fvPatch& currPatch = patches[patchi];

        if (isA<wallFvPatch>(currPatch))
        {
            const vectorField nf(currPatch.nf());

            wallNormal12 = -nf[0];
            faceId = 0;
            patchId = patchi;
        }
    }

    // store local id of near-walll cell to process
    label cellId = patches[patchId].faceCells()[faceId];

    Info<< "    Height to first cell centre y0 = " << y[cellId] << endl;
}
else
{
    scalar y0(0);
    scalar vsum(0);

    forAll(patches, patchi)
    {
        const fvPatch& currPatch = patches[patchi];

        if (isA<wallFvPatch>(currPatch))
        {
            const labelList& faceCells(currPatch.faceCells());

            forAll(faceCells, fci)
            {
                y0 += y[faceCells[fci]] * mesh.V()[faceCells[fci]];
                vsum += mesh.V()[faceCells[fci]];
            }
        }
    }

    reduceToMaster(std::tie(y0, vsum), UniformParallelOp<sumOp<scalar>, 2>{});
    y0 /= vsum;


    Info<< "   Mean first cell centre height, y0_mean: "
         << y0 << " [m]." << endl;

}
