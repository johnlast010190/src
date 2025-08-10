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
    (c) 2011-2014 OpenFOAM Foundation
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "dynamicMultiCriterionRefineFvMesh/dynamicMultiCriterionRefineFvMesh.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "finiteVolume/fvc/fvcGrad.H"
#include "finiteVolume/fvc/fvcSnGrad.H"
#include "finiteVolume/fvc/fvcCurl.H"
#include "finiteVolume/fvc/fvcAverage.H"
#include "interpolation/volPointInterpolation/volPointInterpolation.H"
#include "sets/topoSetSource/topoSetSource.H"
#include "sets/topoSets/cellSet.H"
#include "fvMatrices/fvScalarMatrix/fvScalarMatrix.H"
#include "finiteVolume/fvm/fvm.H"
#include "fields/fvPatchFields/basic/zeroGradient/zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dynamicMultiCriterionRefineFvMesh, 0);
    addToRunTimeSelectionTable(dynamicFvMesh, dynamicMultiCriterionRefineFvMesh, IOobject);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::dynamicMultiCriterionRefineFvMesh::selectRefineCandidates
(
    const scalar lowerRefineLevel,
    const scalar upperRefineLevel,
    const scalarField& vFld,
    PackedBoolList& candidateCell
) const
{
    if (enableRefinementControl_)
    {
        forAll(vFld, cellI)
        {
            if ((vFld[cellI] >= lowerRefineLevel) && (vFld[cellI] <= upperRefineLevel))
            {
                candidateCell.set(cellI, 1);
            }
        }
    }
    else
    {
        dynamicRefineFvMesh::selectRefineCandidates
        (
            lowerRefineLevel,
            upperRefineLevel,
            vFld,
            candidateCell
        );
    }
}


// Get max of connected cell (overwrite minCellField with maxCellField)
Foam::scalarField
Foam::dynamicMultiCriterionRefineFvMesh::minCellField(const volScalarField& vFld) const
{
    scalarField pFld(nPoints(), -GREAT);

    if (enableRefinementControl_)
    {
        forAll(pointCells(), pointi)
        {
            const labelList& pCells = pointCells()[pointi];

            forAll(pCells, i)
            {
                pFld[pointi] = max(pFld[pointi], vFld[pCells[i]]);
            }
        }
    }
    else
    {
        dynamicRefineFvMesh::minCellField(vFld);
    }

    return pFld;
}


void Foam::dynamicMultiCriterionRefineFvMesh::updateRefinementField()
{
    Info<< "Calculating internal refinement field" << endl;

    volScalarField& intRefFld = *internalRefinementFieldPtr_;
    volScalarField& targetFld = *targetLevelPtr_;
    volScalarField& currFld = *isLevelPtr_;

    // Set the internal refinement field to zero to start with
    intRefFld = dimensionedScalar("zero",dimless,0.0);

    // Get the cell level field from dynamicRefineFvMesh
    const labelList& cellLevel = meshCutter().cellLevel();

    // Initialize list for target refinement level per cell
    labelList targetLevel (this->nCells(), 0);

    if (fixedTargetField_)
    {
        forAll(targetFld, cellI)
        {
            targetLevel[cellI] = targetFld[cellI];
        }
    }
    else
    {

        // First fields
        List<word> fieldNames = fields_.toc();
        Field<scalar> refFld(nCells(),0.0);

        forAll(fieldNames, i)
        {
            word dictName = fieldNames[i];
            const word fldName(fields_[dictName].lookup("field"));
            const scalar minValue = readScalar(fields_[dictName].lookup("minValue"));
            const scalar maxValue = readScalar(fields_[dictName].lookup("maxValue"));
            const label refineLevel = readLabel(fields_[dictName].lookup("refineLevel"));

            const volScalarField& fld = this->lookupObject<volScalarField>(fldName);

            // Limit the value of fld based on its max level
            forAll(fld, cellI)
            {
                if ((fld[cellI] >= minValue) && (fld[cellI] <= maxValue))
                {
                    // increase targetLevel up to refineLevel
                    // BUT: do not decrease if cell already marked for higher refinement level by previous criterion
                    targetLevel[cellI] = max(targetLevel[cellI], refineLevel);
                }
            }
        }


        // Then gradients
        List<word> gradFieldNames = gradFields_.toc();

        Field<scalar> cubeRtV = Foam::pow(this->V(),1.0/3.0);

        forAll(gradFieldNames, i)
        {
            word fldName = gradFieldNames[i];
            scalar minValue = gradFields_[fldName][0];
            scalar maxValue = gradFields_[fldName][1];
            label refineLevel = static_cast<label>(gradFields_[fldName][2]);

            const volScalarField& fld = this->lookupObject<volScalarField>(fldName);

            refFld = mag(fvc::grad(fld)) * cubeRtV;

            // Limit the value of refFld based on its max level
            forAll(refFld, cellI)
            {
                if ((refFld[cellI] >= minValue) && (refFld[cellI] <= maxValue))
                {
                    targetLevel[cellI] = max(targetLevel[cellI], refineLevel);
                }
            }
        }


        // Then curls
        List<word> curlFieldNames = curlFields_.toc();

        forAll(curlFieldNames, i)
        {
            word fldName = curlFieldNames[i];
            scalar minValue = curlFields_[fldName][0];
            scalar maxValue = curlFields_[fldName][1];
            label refineLevel = static_cast<label>(curlFields_[fldName][2]);

            const volVectorField& fld = this->lookupObject<volVectorField>(fldName);

            refFld = mag(fvc::curl(fld));

            // Limit the value of refFld based on its max level
            forAll(refFld, cellI)
            {
                if ((refFld[cellI] >= minValue) && (refFld[cellI] <= maxValue))
                {
                    targetLevel[cellI] = max(targetLevel[cellI], refineLevel);
                }
            }
        }


        // At the interface, assumed to be always the maximum refinement level
        List<word> interfaceRefineField = interface_.toc();

        forAll(interfaceRefineField, i)
        {
            word fldName = interfaceRefineField[i];

            // read region of maximum refinement levels inside and outside of interface indicator field
            // (does not need to be alpha, can also be a concentration field)
            scalar innerRefLayers = readScalar(interface_[fldName].lookup("innerRefLayers"));
            scalar outerRefLayers = readScalar(interface_[fldName].lookup("outerRefLayers"));

            Switch interfaceLayers(false);
            scalar nAddLayers(0);
            if (interface_[fldName].found("interfaceLayers"))
            {
                interfaceLayers = Switch(interface_[fldName].lookup("interfaceLayers"));
                nAddLayers = readScalar(interface_[fldName].lookup("nAddLayers"));
            }

            word addLayerMethod("constant");
            if (interface_[fldName].found("addLayerMethod"))
            {
                addLayerMethod = interface_[fldName]["addLayerMethod"][0].wordToken();
                Info<< " addLayerMethod : " << addLayerMethod << endl;
            }


            label refineLevel = maxRefinement_;

            if (interface_[fldName].found("maxRefineLevel"))
            {
                refineLevel = readScalar(interface_[fldName].lookup("maxRefineLevel"));

                // to avoid user input mistakes, limit the value with the maximum allowed
                refineLevel = min(refineLevel, maxRefinement_);
            }


            //TODO: add a min max value of the specified field to be marked for refinement
            // NOW: hard-coded method snGradAlpha: checks mag(alpha_N - alpha_P) > value ? 1:0

            const volScalarField& fld = this->lookupObject<volScalarField>(fldName);

            volScalarField isInterface(intRefFld * 0.0);
            isInterface = dimensionedScalar("isInterface",dimless,0.0);

            surfaceScalarField deltaAlpha ( mag(fvc::snGrad(fld) / this->deltaCoeffs()) );

            const labelUList& owner = this->owner();
            const labelUList& neighbour = this->neighbour();

            scalar threshold = interface_[fldName].lookupOrDefault<scalar>("threshold", 0.05);

            forAll(deltaAlpha, faceI)
            {
                if (deltaAlpha[faceI] > threshold) // currently a fixed prescribed value; should be read from díctionary
                {
                    label own = owner[faceI];
                    label nei = neighbour[faceI];

                    // set isInterface field to one
                    isInterface[own] = 1.0;
                    isInterface[nei] = 1.0;
                }
            }

            // assumed fld=0.5*(fldMax+fldMin) defines the interface
            dimensionedScalar fldInterfaceValue(0.5*(max(fld)+min(fld)));

            //-DD: old implementation based on face interpolation
            //     which results in slower transport in diagonal direction
            // add inner refinement layers
            //for(label i=0; i < innerRefLayers; i++)
            //{
            //    isInterface += neg(- fvc::average(fvc::interpolate(isInterface)) * pos0(fld - fldInterfaceValue));
            //    isInterface = neg(- isInterface);
            //}
            //
            // add outer refinement layers
            //for(label i=0; i < outerRefLayers; i++)
            //{
            //    isInterface += neg(- fvc::average(fvc::interpolate(isInterface)) * pos0(fldInterfaceValue - fld));
            //    isInterface = neg(- isInterface);
            //}
            //
            //forAll(isInterface, cellI)
            //{
            //    if (isInterface[cellI] > 0.5)
            //    {
            //        targetLevel[cellI] = max(targetLevel[cellI], refineLevel);
            //    }
            //}

            //-DD: new version using volPointInterpolation (direction independent buffer layer)
            const volPointInterpolation& pInterp = volPointInterpolation::New(*this);
            const fvMesh& mesh = fld.mesh();

            // add inner refinement layers
            for (label i=0; i < innerRefLayers; i++)
            {
                volScalarField markInner(isInterface*pos0(fld - fldInterfaceValue));
                pointScalarField markLayerP(pInterp.interpolate(markInner));

                forAll(mesh.C(), cellI)
                {
                    scalar sum = 0.;
                    label nPoints = 0;

                    forAll(mesh.cellPoints()[cellI], pointI)
                    {
                        sum += markLayerP[mesh.cellPoints()[cellI][pointI]];
                        nPoints++;
                    }
                    if (nPoints > 0)
                    {
                        sum /= nPoints;
                    }
                    isInterface[cellI] += sum;
                }
            }
            isInterface = pos0(isInterface - SMALL);

            // add outer refinement layers
            for (label i=0; i < outerRefLayers; i++)
            {
                volScalarField markOuter(isInterface*pos0(fldInterfaceValue - fld));
                pointScalarField markLayerP(pInterp.interpolate(markOuter));

                forAll(mesh.C(), cellI)
                {
                    scalar sum = 0.;
                    label nPoints = 0;

                    forAll(mesh.cellPoints()[cellI], pointI)
                    {
                        sum += markLayerP[mesh.cellPoints()[cellI][pointI]];
                        nPoints++;
                    }
                    if (nPoints > 0)
                    {
                        sum /= nPoints;
                    }
                    isInterface[cellI] += sum;
                }
            }
            isInterface = pos0(isInterface - SMALL);

            forAll(isInterface, cellI)
            {
                if (isInterface[cellI] > 0.5)
                {
                    targetLevel[cellI] = max(targetLevel[cellI], refineLevel);
                }
            }

            //-DD: old implementation based on face interpolation
            //     which results in slower transport in diagonal direction
            // expand additional layers if specified:
            if (nAddLayers > 0)
            {
                // loop over additional layers
                for (label i=1; i < refineLevel; i++)
                {
                    // #nAddLayers cells per additional level
                    for (label j=0; j < nAddLayers; j++)
                    {
                        isInterface += neg(- fvc::average(fvc::interpolate(isInterface)));
                        isInterface = neg(- isInterface);
                    }

                    forAll(isInterface, cellI)
                    {
                        if (isInterface[cellI] == 1.0)
                        {
                            targetLevel[cellI] = max(targetLevel[cellI], (refineLevel - i));
                        }
                    }
                }
            }
        }

        // Set refinement in physical regions (force the mesh to stay refined
        // near key features)
        forAll(refinedRegions_, regionI)
        {
            const entry& region = refinedRegions_[regionI];

            autoPtr<topoSetSource> source =
                topoSetSource::New(region.keyword(), *this, region.dict());

            cellSet selectedCellSet
            (
                *this,
                "cellSet",
                nCells()/10+1 //estimate
            );

            source->applyToSet
            (
                topoSetSource::NEW,
                selectedCellSet
            );

            const labelList cells = selectedCellSet.toc();

            label minLevel = readLabel(region.dict().lookup("minLevel"));

            forAll(cells, i)
            {
                const label& cellI = cells[i];
                targetLevel[cellI] = max(targetLevel[cellI], minLevel);
            }
        }

        //-DD: add buffer layer based on targetLevel field to prevent 2-to-1 refinement
        if (nBufferLayers_ > 1)
        {
            volScalarField blockedLevel ( targetFld * 0. );

            for (label currLayer=maxRefinement_; currLayer>=1; currLayer--)
            {
                forAll(targetLevel, cellI)
                {
                    if (targetLevel[cellI] >= currLayer)
                    {
                        blockedLevel[cellI] = targetLevel[cellI];
                    }
                }

                for (label i=1; i<nBufferLayers_; i++)
                {
                    blockedLevel = max(blockedLevel, pos0(fvc::average(fvc::interpolate(blockedLevel)) - SMALL)*(currLayer-1));
                }

                labelList blockLev(blockedLevel.internalField().size(), 0);
                forAll(blockLev, cellI)
                {
                    blockLev[cellI] = blockedLevel[cellI];
                }
                targetLevel = max(targetLevel, blockLev);
            }
        }

        // smooth target refinement level field
        if (targetLevelSmoothing_)
        {
            // create a field with non-default boundary conditions in order to solve relaxation
            volScalarField targetLevelField
            (
                IOobject
                (
                    "targetLevelField",
                    time().timeName(),
                    *this,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                *this,
                dimensionedScalar
                (
                    "zero",
                    targetFld.dimensions(),
                    scalar(0)
                ),
                zeroGradientFvPatchScalarField::typeName
            );
            forAll(targetLevelField.internalField(), cellI)
            {
                targetLevelField.primitiveFieldRef()[cellI] = targetLevel[cellI];
            }
            targetLevelField.correctBoundaryConditions();

            // call smoothing algorithm
            smoothTargetLevel(targetLevelField, targetLevel);
        }
    }

    // mark cells to be refined/unrefined based on above criteria:
    // simply check, if targetLevel lower or higher cellLevel
    forAll(intRefFld.internalField(), cellI)
    {
        intRefFld.primitiveFieldRef()[cellI] = targetLevel[cellI] - cellLevel[cellI];
        targetFld.primitiveFieldRef()[cellI] = targetLevel[cellI];
        currFld.primitiveFieldRef()[cellI] = cellLevel[cellI];
    }

    intRefFld.correctBoundaryConditions();

    //Info<<"Min,max refinement field = " << Foam::min(intRefFld).value() << ", "
    //    << Foam::max(intRefFld).value() << endl;
}


// smooth the target level field by elliptical relaxation
void Foam::dynamicMultiCriterionRefineFvMesh::smoothTargetLevel
(
    const volScalarField& targetLevelField,
    labelList& targetLevel
)
{
    volScalarField targetLevelSmooth("targetLevelSmooth", targetLevelField);
    volScalarField targetLevelDiscrete("targetLevelDiscrete", targetLevelField);

    // needed to discretize implicit source term
    dimensionedScalar unity("unity", dimless, 1);

  // define diffusion coefficient based on local mesh size
    surfaceScalarField diffCoeff ( pow(diffCoeff_/targetLevelField.mesh().deltaCoeffs(),2) );

  // set diffusion coefficient to zero away from refinement jumps
    // compute an indicator field around jumps in between refinement levels
    surfaceScalarField hasRefLevelJump(neg(-mag(fvc::snGrad(targetLevelField))));

    // expand the indicator field by N layers
    for (label count=0; count < nLayerSmooth_; count++)
    {
        hasRefLevelJump = neg(-mag(fvc::interpolate(fvc::average(hasRefLevelJump))));
    }
    // penalize diffusion
    diffCoeff *= hasRefLevelJump;

  // scale target refinement level field
    // unique sort indices (get number of different target levels)
    labelList targetLevelUnique(targetLevel);
    uniqueOrder(targetLevel, targetLevelUnique);
    label nTargetLevels(targetLevelUnique.size());

    // get unique sorted list of different levels
    labelList sortedLevels(nTargetLevels, 0);
    forAll(sortedLevels, idx)
    {
        sortedLevels[idx] = targetLevel[targetLevelUnique[idx]];
    }
    labelList mappedLevels(identity(nTargetLevels)+1);

    Info<< " sorted levels: " << sortedLevels << endl;
    Info<< " mapped levels: " << mappedLevels << endl;

    // scale target refinement level field
    forAll(targetLevelDiscrete, cellI)
    {
        forAll(sortedLevels, idx)
        {
            if (targetLevelDiscrete[cellI] == sortedLevels[idx])
            {
                targetLevelDiscrete[cellI] = mappedLevels[idx];
            }
        }
    }
    targetLevelDiscrete.correctBoundaryConditions();
    targetLevelSmooth = targetLevelDiscrete;

  // elliptical relaxation loop
    for (label count=0; count < 2; count++)
    {
        Info<<"solving elliptical relaxation loop " << count << endl;

        fvScalarMatrix smoothRefinementEqn
        (
          - fvm::laplacian(diffCoeff,targetLevelSmooth)
          + fvm::Sp(unity, targetLevelSmooth)
          ==
            targetLevelDiscrete
        );
        smoothRefinementEqn.solve();

        // re-set target source term based on smooth field solution
        forAll(targetLevelDiscrete, cellI)
        {
            label value = round(targetLevelSmooth[cellI]);
            targetLevelDiscrete[cellI] = value;
        }
        targetLevelDiscrete.correctBoundaryConditions();
    }

 // scale target refinement level field back!!
    forAll(targetLevelDiscrete, cellI)
    {
        for (label idx = nTargetLevels-1; idx >= 0; idx--)
        {
            if (targetLevelDiscrete[cellI] == mappedLevels[idx])
            {
                targetLevelDiscrete[cellI] = sortedLevels[idx];
            }
        }
    }

    if (allowOnlyLevelGrowth_)
    {
        // allow only layer growth!
        forAll(targetLevelDiscrete, cellI)
        {
            targetLevelDiscrete[cellI] = max(targetLevelDiscrete[cellI], targetLevel[cellI]);
        }
    }

    // overwrite targetLevel
    forAll(targetLevelDiscrete, cellI)
    {
        targetLevel[cellI] = targetLevelDiscrete[cellI];
    }

    Info<<"finished smoothing target refinement level field! " << endl;
}


void Foam::dynamicMultiCriterionRefineFvMesh::readDict()
{
    IOdictionary dynamicMeshDict
    (
        IOobject
        (
            "dynamicMeshDict",
            time().constant(),
            *this,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE,
            false
        )
    );
    dictionary& dict = dynamicMeshDict.subDict("dynamicRefineFvMeshCoeffs");

    nBufferLayers_ = dict.lookupOrDefault<label>("nBufferLayers",0);
    refineInterval_ = readLabel(dict.lookup("refineInterval"));
    maxRefinement_ = readLabel(dict.lookup("maxRefinement"));

    if (dynamicMeshDict.isDict("refinementControls"))
    {
        dictionary& refineControlDict =
            dynamicMeshDict.subDict("refinementControls");

        enableRefinementControl_ =
            refineControlDict.lookupOrDefault("enableRefinementControl", true);

        fixedTargetField_ =
            refineControlDict.lookupOrDefault<Switch>("fixedTargetField", false);

        targetLevelSmoothing_ =
            refineControlDict.lookupOrDefault<Switch>("targetLevelSmoothing", false);

        allowOnlyLevelGrowth_ =
            refineControlDict.lookupOrDefault<Switch>("allowOnlyLevelGrowth", false);

        diffCoeff_ =
            refineControlDict.lookupOrDefault<scalar>("diffCoeff", 2.5);

        nLayerSmooth_ =
            refineControlDict.lookupOrDefault<label>("nLayerSmooth", 3);

        if (enableRefinementControl_)
        {
            const word fieldName(dict.lookup("field"));

            // Overwrite field name entry in dynamicRefineFvMeshCoeffs if necessary
            if (fieldName != "internalRefinementField")
            {
                Info<<"Overwriting dictionary entries: \n"
                    <<"'field', 'lowerRefineLevel' and 'upperRefineLevel' \n"
                    <<" in dynamicRefineFvMeshCoeffs!" << endl;
                dict.set
                (
                  "field",
                  word("internalRefinementField")
                );

                dict.set
                (
                  "lowerRefineLevel",
                  scalar(0.5)
                );

                dict.set
                (
                  "upperRefineLevel",
                  scalar(maxRefinement_ + 0.5)
                );

                dict.set
                (
                  "unrefineLevel",
                  scalar(-0.5)
                );

                // write dictionary to file
                dynamicMeshDict.regIOobject::write();
            }

            // Read HashTable of field-refinement scalars
            if (refineControlDict.found("fields"))
            {
                fields_ = HashTable< dictionary >
                (
                    refineControlDict.lookup("fields")
                );
            }

            // Read HashTable of gradient-refinement scalars
            if (refineControlDict.found("gradients"))
            {
                gradFields_ = HashTable< List<scalar>>
                (
                    refineControlDict.lookup("gradients")
                );
            }

            // Read HashTable of curl-refinement vectors
            if (refineControlDict.found("curls"))
            {
                curlFields_ = HashTable< List<scalar>>
                (
                    refineControlDict.lookup("curls")
                );
            }

            // Read interface refinement data
            if (refineControlDict.found("interface"))
            {
                interface_ = HashTable< dictionary >
                (
                    refineControlDict.lookup("interface")
                );
            }

            // Read refinement regions
            if (refineControlDict.found("regions"))
            {
                refinedRegions_ = PtrList<entry>
                (
                    refineControlDict.lookup("regions")
                );
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicMultiCriterionRefineFvMesh::dynamicMultiCriterionRefineFvMesh
(
    const IOobject& io
)
:
    dynamicRefineFvMesh(io),
    internalRefinementFieldPtr_(nullptr),
    targetLevelPtr_(nullptr),
    isLevelPtr_(nullptr),
    fields_(),
    gradFields_(),
    curlFields_(),
    interface_(),
    refinedRegions_(),
    enableRefinementControl_(false),
    fixedTargetField_(false),
    targetLevelSmoothing_(false),
    allowOnlyLevelGrowth_(false),
    diffCoeff_(0),
    nLayerSmooth_(0),
    nBufferLayers_(0),
    refineInterval_(1),
    maxRefinement_(0)
{
    readDict();

    if (enableRefinementControl_)
    {
        internalRefinementFieldPtr_ = new volScalarField
        (
            IOobject
            (
                "internalRefinementField",
                this->time().timeName(),
                *this,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            *this,
            dimensionedScalar("zero", dimless, 0.0)
        );
        targetLevelPtr_ = new volScalarField
        (
            IOobject
            (
                "targetLevel",
                this->time().timeName(),
                *this,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            *this,
            dimensionedScalar("zero", dimless, 0.0)
        );
        isLevelPtr_ = new volScalarField
        (
            IOobject
            (
                "isLevel",
                this->time().timeName(),
                *this,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            *this,
            dimensionedScalar("zero", dimless, 0.0)
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dynamicMultiCriterionRefineFvMesh::~dynamicMultiCriterionRefineFvMesh()
{
    deleteDemandDrivenData(internalRefinementFieldPtr_);
    deleteDemandDrivenData(targetLevelPtr_);
    deleteDemandDrivenData(isLevelPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::dynamicMultiCriterionRefineFvMesh::update()
{
    //Part 1 - Update internally calculated refinement field
    readDict();

    if (time().value() >= refineStartTime_ && time().value() <= refineStopTime_)
    {
        if (enableRefinementControl_ && time().timeIndex() > 0 && time().timeIndex() % refineInterval_ == 0)
        {
            updateRefinementField();
        }
    }

    //Part 2 - Call normal update from dynamicRefineFvMesh
    bool hasChanged = dynamicRefineFvMesh::update();

    return hasChanged;
}

// ************************************************************************* //
