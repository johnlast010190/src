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
    (c) 1991-2005 OpenCFD Ltd.
    (c) 2016 Esi Ltd.

Description
    Face to edge interpolation scheme. Included in faMesh.

\*---------------------------------------------------------------------------*/

#include "faMesh/faMesh.H"
#include "fields/areaFields/areaFields.H"
#include "fields/edgeFields/edgeFields.H"
#include "include/demandDrivenData.H"
#include "fields/faPatchFields/faPatchField/faPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(edgeInterpolation, 0);


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void edgeInterpolation::clearOut()
{
    deleteDemandDrivenData(lPN_);
    deleteDemandDrivenData(weights_);
    deleteDemandDrivenData(deltaCoeffs_);
    deleteDemandDrivenData(nonOrthDeltaCoeffs_);
    deleteDemandDrivenData(nonOrthCorrectionVectors_);
    deleteDemandDrivenData(skewCorrectionVectors_);
//     deleteDemandDrivenData(leastSquarePvectors_);
//     deleteDemandDrivenData(leastSquareNvectors_);
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

edgeInterpolation::edgeInterpolation(const faMesh& fam)
:
    faMesh_(fam),
    lPN_(nullptr),
    weights_(nullptr),
    deltaCoeffs_(nullptr),
    nonOrthDeltaCoeffs_(nullptr),
    nonOrthCorrectionVectors_(nullptr),
    skew_(true),
    skewCorrectionVectors_(nullptr)
//     leastSquarePvectors_(nullptr),
//     leastSquareNvectors_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

edgeInterpolation::~edgeInterpolation()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const edgeScalarField& edgeInterpolation::lPN() const
{
    if (!lPN_)
    {
        makeLPN();
    }

    return (*lPN_);
}


const edgeScalarField& edgeInterpolation::weights() const
{
    if (!weights_)
    {
        makeWeights();
    }

    return (*weights_);
}


const edgeScalarField& edgeInterpolation::deltaCoeffs() const
{
    if (!deltaCoeffs_)
    {
        makeDeltaCoeffs();
    }

    return (*deltaCoeffs_);
}

const edgeScalarField& edgeInterpolation::nonOrthDeltaCoeffs() const
{
    if (!nonOrthDeltaCoeffs_)
    {
        makeNonOrthDeltaCoeffs();
    }

    return (*nonOrthDeltaCoeffs_);
}

const edgeVectorField& edgeInterpolation::nonOrthCorrectionVectors() const
{
    if (!nonOrthCorrectionVectors_)
    {
        makeNonOrthCorrectionVectors();
    }

    return (*nonOrthCorrectionVectors_);
}


bool edgeInterpolation::skew() const
{
    if (skew_ == true && !skewCorrectionVectors_)
    {
        makeSkewCorrectionVectors();
    }

    return skew_;
}


const edgeVectorField& edgeInterpolation::skewCorrectionVectors() const
{
    if (!skew())
    {
        FatalErrorInFunction
            << "cannot return skewCorrectionVectors; mesh is now skewed"
            << abort(FatalError);
    }

    return (*skewCorrectionVectors_);
}


// const edgeVectorField& edgeInterpolation::leastSquarePvectors() const
// {
//     if (!leastSquarePvectors_)
//     {
//         makeLeastSquareVectors();
//     }

//     return (*leastSquarePvectors_);
// }


// const edgeVectorField& edgeInterpolation::leastSquareNvectors() const
// {
//     if (!leastSquareNvectors_)
//     {
//         makeLeastSquareVectors();
//     }

//     return (*leastSquareNvectors_);
// }


// Do what is neccessary if the mesh has moved
bool edgeInterpolation::movePoints() const
{
    deleteDemandDrivenData(lPN_);
    deleteDemandDrivenData(weights_);
    deleteDemandDrivenData(deltaCoeffs_);
    deleteDemandDrivenData(nonOrthDeltaCoeffs_);
    deleteDemandDrivenData(nonOrthCorrectionVectors_);

    skew_ = true;
    deleteDemandDrivenData(skewCorrectionVectors_);

//     deleteDemandDrivenData(leastSquarePvectors_);
//     deleteDemandDrivenData(leastSquareNvectors_);

    return true;
}


void edgeInterpolation::makeLPN() const
{
    if (debug)
    {
        Info<< "edgeInterpolation::makeLPN() : "
            << "Constructing geodesic distance between points P and N"
            << endl;
    }

    lPN_ = new edgeScalarField
    (
        IOobject
        (
            "lPN",
            faMesh_.time().constant(),
            mesh().thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh(),
        dimLength
    );
    edgeScalarField& lPN = *lPN_;


    // Set local references to mesh data
    const edgeVectorField& edgeCentres = mesh().edgeCentres();
    const areaVectorField& faceCentres = mesh().areaCentres();
    const labelUList& owner = mesh().owner();
    const labelUList& neighbour = mesh().neighbour();

    forAll(owner, edgeI)
    {
        vector curSkewCorrVec = vector::zero;

        if (skew())
        {
            curSkewCorrVec = skewCorrectionVectors()[edgeI];
        }

        scalar lPE =
            mag
            (
                edgeCentres[edgeI]
              - curSkewCorrVec
              - faceCentres[owner[edgeI]]
            );

        scalar lEN =
            mag
            (
                faceCentres[neighbour[edgeI]]
              - edgeCentres[edgeI]
              + curSkewCorrVec
            );

        lPN.primitiveFieldRef()[edgeI] = (lPE + lEN);
    }


    forAll(lPN.boundaryField(), patchI)
    {
        mesh().boundary()[patchI].makeDeltaCoeffs
        (
            lPN.boundaryFieldRef()[patchI]
        );

        lPN.boundaryFieldRef()[patchI] = 1.0/lPN.boundaryField()[patchI];
    }


    if (debug)
    {
        Info<< "edgeInterpolation::makeLPN() : "
            << "Finished constructing geodesic distance PN"
            << endl;
    }
}


void edgeInterpolation::makeWeights() const
{
    if (debug)
    {
        Info<< "edgeInterpolation::makeWeights() : "
            << "Constructing weighting factors for edge interpolation"
            << endl;
    }


    weights_ = new edgeScalarField
    (
        IOobject
        (
            "weightingFactors",
            mesh()().pointsInstance(),
            mesh()(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh(),
        dimless
    );
    edgeScalarField& weights = *weights_;


    // Set local references to mesh data
    const edgeVectorField& edgeCentres = mesh().edgeCentres();
    const areaVectorField& faceCentres = mesh().areaCentres();
    const labelUList& owner = mesh().owner();
    const labelUList& neighbour = mesh().neighbour();

    forAll(owner, edgeI)
    {
        vector curSkewCorrVec = vector::zero;

        if (skew())
        {
            curSkewCorrVec = skewCorrectionVectors()[edgeI];
        }

        scalar lPE =
            mag
            (
                edgeCentres[edgeI]
              - curSkewCorrVec
              - faceCentres[owner[edgeI]]
            );

        scalar lEN =
            mag
            (
                faceCentres[neighbour[edgeI]]
              - edgeCentres[edgeI]
              + curSkewCorrVec
            );

        weights.primitiveFieldRef()[edgeI] =
            lEN
            /(
                lPE
#               ifdef BAD_MESH_STABILISATION
              + VSMALL
#               endif
              + lEN
            );
    }

    forAll(mesh().boundary(), patchI)
    {
        mesh().boundary()[patchI].makeWeights
        (
            weights.boundaryFieldRef()[patchI]
        );
    }

    if (debug)
    {
        Info<< "edgeInterpolation::makeWeights() : "
            << "Finished constructing weighting factors for face interpolation"
            << endl;
    }
}


void edgeInterpolation::makeDeltaCoeffs() const
{
    if (debug)
    {
        Info<< "edgeInterpolation::makeDeltaCoeffs() : "
            << "Constructing differencing factors array for edge gradient"
            << endl;
    }

    // Force the construction of the weighting factors
    // needed to make sure deltaCoeffs are calculated for parallel runs.
    weights();

    deltaCoeffs_ = new edgeScalarField
    (
        IOobject
        (
            "differenceFactors_",
            mesh()().pointsInstance(),
            mesh()(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        faMesh_,
        dimless/dimLength
    );
    edgeScalarField& DeltaCoeffs = *deltaCoeffs_;
    scalarField& dc = DeltaCoeffs.primitiveFieldRef();


    // Set local references to mesh data
    const areaVectorField& faceCentres = mesh().areaCentres();
    const edgeVectorField& edgeCentres = mesh().edgeCentres();
    const labelUList& owner = mesh().owner();
    const labelUList& neighbour = mesh().neighbour();

    forAll(owner, edgeI)
    {
        // Calc PN arc length
                vector curSkewCorrVec = vector::zero;

                if (skew())
                {
                    curSkewCorrVec = skewCorrectionVectors()[edgeI];
                }

                scalar lPE =
                    mag
                    (
                        edgeCentres[edgeI]
                      - curSkewCorrVec
                      - faceCentres[owner[edgeI]]
                    );

                scalar lEN =
                    mag
                    (
                        faceCentres[neighbour[edgeI]]
                      - edgeCentres[edgeI]
                      + curSkewCorrVec
                    );

                scalar lPN = lPE + lEN;
        // PA
        dc[edgeI] = 1.0/mag(lPN);
    }


    forAll(DeltaCoeffs.boundaryField(), patchi)
    {
        DeltaCoeffs.boundaryFieldRef()[patchi] =
            1.0/mag(mesh().boundary()[patchi].delta());
    }

//    forAll(DeltaCoeffs.boundaryField(), patchI)
//    {
//        mesh().boundary()[patchI].makeDeltaCoeffs
//        (
//            DeltaCoeffs.boundaryField()[patchI]
//        );
//    }
}


void edgeInterpolation::makeNonOrthDeltaCoeffs() const
{
    if (debug)
    {
        Pout<< "edgeInterpolation::makeNonOrthDeltaCoeffs() : "
            << "Constructing differencing factors array for face gradient"
            << endl;
    }

    // Force the construction of the weighting factors
    // needed to make sure deltaCoeffs are calculated for parallel runs.
    weights();

    nonOrthDeltaCoeffs_ = new edgeScalarField
    (
        IOobject
        (
            "nonOrthDeltaCoeffs",
            mesh()().pointsInstance(),
            mesh()(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh(),
        dimless/dimLength
    );
    edgeScalarField& nonOrthDeltaCoeffs = *nonOrthDeltaCoeffs_;
    scalarField& dc = nonOrthDeltaCoeffs.primitiveFieldRef();


    // Set local references to mesh data
    const edgeVectorField& edgeCentres = mesh().edgeCentres();
    const areaVectorField& faceCentres = mesh().areaCentres();
    const labelUList& owner = mesh().owner();
    const labelUList& neighbour = mesh().neighbour();
    const edgeVectorField& lengths = mesh().Le();

    const edgeList& edges = mesh().edges();
    const pointField& points = mesh().points();

    forAll(owner, edgeI)
    {
        // Edge normal - area normal
        vector edgeNormal = lengths[edgeI]^edges[edgeI].vec(points);

        // Unit delta vector
        vector unitDelta =
            faceCentres[neighbour[edgeI]]
          - faceCentres[owner[edgeI]];

        scalar magEdgeNormal = mag(edgeNormal);
        if (magEdgeNormal!=0)
        {
            edgeNormal /= magEdgeNormal;
        }
        else
        {
            edgeNormal = Zero;
        }

        unitDelta -=
            edgeNormal*(edgeNormal&unitDelta);

        if (mag(unitDelta)<VSMALL)
        {
            unitDelta = vector::zero;
        }
        else
        {
            unitDelta /= mag(unitDelta);
        }

        // Calc PN arc length
        vector curSkewCorrVec = vector::zero;

        if (skew())
        {
            curSkewCorrVec = skewCorrectionVectors()[edgeI];
        }

        scalar lPE =
            mag
            (
                edgeCentres[edgeI]
              - curSkewCorrVec
              - faceCentres[owner[edgeI]]
            );

        scalar lEN =
            mag
            (
                faceCentres[neighbour[edgeI]]
              - edgeCentres[edgeI]
              + curSkewCorrVec
            );

        scalar lPN = lPE + lEN;


        // Edge normal - area tangent
        scalar magL = mag(lengths[edgeI]);
        if (magL != 0)
        {
            edgeNormal = lengths[edgeI]/magL;
        }
        else
        {
            edgeNormal = Zero;
        }

        // Stabilised form for bad meshes.  HJ, 23/Jul/2009
        dc[edgeI] = 1.0/max((lPN*(unitDelta & edgeNormal)), 0.05*lPN);
    }

    forAll(nonOrthDeltaCoeffs.boundaryField(), patchi)
    {
        vectorField delta(mesh().boundary()[patchi].delta());

        nonOrthDeltaCoeffs.boundaryFieldRef()[patchi] =
            1.0/max(mesh().boundary()[patchi].edgeNormals()
                    & delta, 0.05*mag(delta));
    }
}

void edgeInterpolation::makeNonOrthCorrectionVectors() const
{
    if (debug)
    {
        Info<< "edgeInterpolation::makeCorrectionVectors() : "
            << "Constructing non-orthogonal correction vectors"
            << endl;
    }

    nonOrthCorrectionVectors_ = new edgeVectorField
    (
        IOobject
        (
            "nonOrthCorrectionVectors",
            mesh()().pointsInstance(),
            mesh()(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh(),
        dimless
    );
    edgeVectorField& corrVecs = *nonOrthCorrectionVectors_;

    // Set local references to mesh data
    const areaVectorField& faceCentres = mesh().areaCentres();

    const labelUList& owner = mesh().owner();
    const labelUList& neighbour = mesh().neighbour();

    const edgeVectorField& lengths = mesh().Le();

    const edgeList& edges = mesh().edges();
    const pointField& points = mesh().points();

    const edgeScalarField& nonOrthoDeltaCoeffs = nonOrthDeltaCoeffs();

    forAll(owner, edgeI)
    {
        // Edge normal - area normal
        vector edgeNormal = lengths[edgeI] ^ edges[edgeI].vec(points);

        scalar magEN =  mag(edgeNormal);
        if (magEN != 0)
        {
            edgeNormal /= magEN;
        }
        else
        {
            edgeNormal = Zero;
        }

        // Unit delta vector
        vector unitDelta =
            faceCentres[neighbour[edgeI]]
          - faceCentres[owner[edgeI]];

        unitDelta -= edgeNormal*(edgeNormal & unitDelta);

//        unitDelta /= mag(unitDelta);

        // Edge normal - area tangent
        scalar magL = mag(lengths[edgeI]);
        if (magL != 0)
        {
            edgeNormal = lengths[edgeI]/magL;
        }
        else
        {
            edgeNormal = Zero;
        }

        // Edge correction vector
        corrVecs.primitiveFieldRef()[edgeI] =
            edgeNormal
          - nonOrthoDeltaCoeffs[edgeI]*unitDelta;
    }


    // Boundary correction vectors set to zero for boundary patches
    // and calculated consistently with internal corrections for
    // coupled patches

    forAll(corrVecs.boundaryField(), patchI)
    {
        faePatchVectorField& patchCorrVecs = corrVecs.boundaryFieldRef()[patchI];

        if (!patchCorrVecs.coupled())
        {
            patchCorrVecs = vector::zero;
        }
        else
        {
            const faePatchScalarField& patchNonOrthDeltaCoeffs
                = nonOrthoDeltaCoeffs.boundaryField()[patchI];

            const vectorField patchDeltas(mesh().boundary()[patchI].delta());

            forAll(patchCorrVecs, patchEdgeI)
            {
                // Unit delta vector
                const vector& unitDelta = patchDeltas[patchEdgeI];

                // Edge normal - area tangent
                vector edgeNormal = lengths.boundaryField()[patchI][patchEdgeI]
                    /mag(lengths.boundaryField()[patchI][patchEdgeI]);

                patchCorrVecs[patchEdgeI] =
                    edgeNormal - patchNonOrthDeltaCoeffs[patchEdgeI]*unitDelta;
            }

        }
    }

    if (debug)
    {
        Info<< "edgeInterpolation::makeCorrectionVectors() : "
            << "Finished constructing non-orthogonal correction vectors"
            << endl;
    }
}


void edgeInterpolation::makeSkewCorrectionVectors() const
{
    if (debug)
    {
        Info<< "edgeInterpolation::makeSkewCorrectionVectors() : "
            << "Constructing skew correction vectors"
            << endl;
    }

    skewCorrectionVectors_ = new edgeVectorField
    (
        IOobject
        (
            "skewCorrectionVectors",
            mesh()().pointsInstance(),
            mesh()(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh(),
        dimless
    );
    edgeVectorField& SkewCorrVecs = *skewCorrectionVectors_;

    // Set local references to mesh data
    const areaVectorField& C = mesh().areaCentres();
    const edgeVectorField& Ce = mesh().edgeCentres();

    const labelUList& owner = mesh().owner();
    const labelUList& neighbour = mesh().neighbour();

    const pointField& points = mesh().points();
    const edgeList& edges = mesh().edges();


    forAll(neighbour, edgeI)
    {
        vector P = C[owner[edgeI]];
        vector N = C[neighbour[edgeI]];
        vector S = points[edges[edgeI].start()];
        vector e = edges[edgeI].vec(points);

        if ((((N - P)^e)&((N - P)^e))<VSMALL)
        {
            SkewCorrVecs[edgeI] = vector::zero;
        }
        else
        {
            scalar alpha = - ( ( (N - P)^(S - P) )&( (N - P)^e ) )/
                    ( ( (N - P)^e )&( (N - P)^e ) );

            vector E = S + alpha*e;

            SkewCorrVecs[edgeI] = Ce[edgeI] - E;
        }
    }


    forAll(SkewCorrVecs.boundaryField(), patchI)
    {
        faePatchVectorField& patchSkewCorrVecs =
            SkewCorrVecs.boundaryFieldRef()[patchI];

        if (patchSkewCorrVecs.coupled())
        {
            const labelUList& edgeFaces =
                mesh().boundary()[patchI].edgeFaces();

            const edgeList::subList patchEdges =
                mesh().boundary()[patchI].patchSlice(edges);

            vectorField ngbC ( C.boundaryField()[patchI].patchNeighbourField() );

            forAll(patchSkewCorrVecs, edgeI)
            {
                vector P = C[edgeFaces[edgeI]];
                vector N = ngbC[edgeI];
                vector S = points[patchEdges[edgeI].start()];
                vector e = patchEdges[edgeI].vec(points);

                if ((((N - P)^e)&((N - P)^e))<VSMALL)
                   {
                    patchSkewCorrVecs[edgeI] = vector::zero;
                   }
                else
                {
                    scalar alpha = - ( ( (N - P)^(S - P) )&( (N - P)^e ) )/
                            ( ( (N - P)^e )&( (N - P)^e ) );

                    vector E = S + alpha*e;

                    patchSkewCorrVecs[edgeI] =
                            Ce.boundaryField()[patchI][edgeI] - E;
                }
            }
        }
        else
        {
            patchSkewCorrVecs = vector::zero;
        }
    }


    scalar skewCoeff = 0.0;

    // Calculating PN arc length
    scalarField lPN(owner.size());

    forAll(owner, edgeI)
    {
        lPN[edgeI] =
            mag
            (
                Ce[edgeI]
              - SkewCorrVecs[edgeI]
              - C[owner[edgeI]]
            )
          + mag
            (
                C[neighbour[edgeI]]
              - Ce[edgeI]
              + SkewCorrVecs[edgeI]
            );
    }

    if (lPN.size() > 0)
    {
        skewCoeff = max(mag(SkewCorrVecs.internalField())/mag(lPN));
    }

    if (debug)
    {
        Info<< "edgeInterpolation::makeSkewCorrectionVectors() : "
            << "skew coefficient = " << skewCoeff << endl;
    }

    if (skewCoeff < 0.1)
    {
        skew_ = false;
        deleteDemandDrivenData(skewCorrectionVectors_);
    }
    else
    {
        skew_ = true;
    }

    if (debug)
    {
        Info<< "edgeInterpolation::makeSkewCorrectionVectors() : "
            << "Finished constructing skew correction vectors"
            << endl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
