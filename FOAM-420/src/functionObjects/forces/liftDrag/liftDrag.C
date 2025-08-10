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
    (c) 2010-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "liftDrag/liftDrag.H"
#include "fields/volFields/volFields.H"
#include "db/dictionary/dictionary.H"
#include "db/Time/Time.H"
#include "primitives/ints/label/label.H"
#include "fields/fvsPatchFields/fvsPatchField/fvsPatchField.H"
#include "fields/fvPatchFields/basic/zeroGradient/zeroGradientFvPatchFields.H"
#include "surfaceMesh/surfaceMesh.H"
#include "turbulentTransportModels/turbulentTransportModel.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"
#include "cfdTools/general/porosityModel/porosityModel/porosityModel.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "meshes/meshShapes/cell/cell.H"
#include "meshes/polyMesh/polyMesh.H"
#include "centreOfMass/centreOfMass.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(liftDrag, 0);
    addToRunTimeSelectionTable(functionObject, liftDrag, dictionary);
}

template<>
const char* Foam::NamedEnum
<
    Foam::functionObjects::liftDrag::axisType,
    3
>::names[] =
{
    "global",
    "boundBox",
    "normalized"
};
}

const Foam::NamedEnum<Foam::functionObjects::liftDrag::axisType, 3>
    Foam::functionObjects::liftDrag::axisTypeNames_;


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::label Foam::functionObjects::liftDrag::binIndex
(
    const binData& ldBin,
    const point& p
)
{
    scalar bv = (binAxis(ldBin) & p);
    bv -= ldBin.binMin_;

    bv /= ldBin.binSize_;

    if (bv < 0 || bv > ldBin.nBins_)
    {
        FatalErrorInFunction
            << "Point (X = " << (p & binAxis(ldBin))
            << ") out of bin range: "
            <<  ldBin.binMin_ << "-" << ldBin.binMax_ << endl
            << "Bin index: " << bv << ", binSize: " << ldBin.binSize_
            << ", nBins: " << ldBin.nBins_
            << exit(FatalError);
    }

    return label(bv);
}


void Foam::functionObjects::liftDrag::calculateBins()
{
    forAll(liftDragBinPtr_, bI)
    {
        binData& ldBin = liftDragBinPtr_[bI];

        // Find min and max bin values
        ldBin.binMin_ = GREAT;
        ldBin.binMax_ = -GREAT;

        // Surface
        forAllConstIter(labelHashSet, patchSet(), iter)
        {
            label patchi = iter.key();

            scalarField axisDotCf
            (
                binAxis(ldBin) & mesh_.boundary()[patchi].Cf()
            );

            forAll(axisDotCf, i)
            {
                ldBin.binMin_ = min(ldBin.binMin_, axisDotCf[i]);
                ldBin.binMax_ = max(ldBin.binMax_, axisDotCf[i]);
            }
        }

        // Porous media
        if (porosity_)
        {
            const HashTable<const porosityModel*> models =
                obr_.lookupClass<porosityModel>();

            forAll(porousZoneNames_, i)
            {
                const porosityModel& pm =
                    const_cast<porosityModel&>(*models[porousZoneNames_[i]]);

                forAll(pm.cellZoneIDs(), zoneI)
                {
                    label cellZoneI = pm.cellZoneIDs()[zoneI];
                    const cellZone& cZone = mesh_.cellZones()[cellZoneI];
                    bool calculateForce = checkZone(cZone);

                    if (calculateForce)
                    {
                        const labelList& cells = mesh_.cellZones()[cellZoneI];

                        forAll(cells, czI)
                        {
                            label cellI = cells[czI];
                            scalar axisDotC =
                                binAxis(ldBin) & mesh_.C()[cellI];
                            ldBin.binMin_ = min(ldBin.binMin_, axisDotC);
                            ldBin.binMax_ = max(ldBin.binMax_, axisDotC);
                        }
                    }
                }
            }
        }
        reduce(ldBin.binMin_, minOp<scalar>());
        reduce(ldBin.binMax_, maxOp<scalar>());

        ldBin.binSize_ = (ldBin.binMax_ - ldBin.binMin_)/scalar(ldBin.nBins_);
        ldBin.binMax_ += 0.01*ldBin.binSize_;
        ldBin.binSize_ = (ldBin.binMax_ - ldBin.binMin_)/scalar(ldBin.nBins_);

    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::liftDrag::liftDrag
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    forceCoeffs(name, runTime, dict),
    CofP_refPoint_(false),
    cOfPressure_(false),
    porousZoneNames_(wordList(0)),
    refCofP_(Zero),
    wheelBase_(0),
    maxCp_(GREAT),
    minCp_(-GREAT),
    nAveragingSteps_(1),
    totalLift_(0),
    frontLift_(0),
    rearLift_(0),
    drag_(0),
    side_(0),
    moment_(Zero),
    CofP_(Zero),
    CofM_(Zero),
    averagingIndex_(0),
    liftDragRegions_(false),
    liftDragBinPtr_(0),
    liftDragFilePtr_(nullptr)
{
    liftDrag::read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::liftDrag::~liftDrag()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::vector Foam::functionObjects::liftDrag::binAxis
(
    const binData& data
) const
{
    if (coorFramePtr_ && definedInFrame_)
    {
        return coorFramePtr_->coorSys().globalVector(data.binAxis_);
    }
    return data.binAxis_;
}


bool Foam::functionObjects::liftDrag::execute()
{
    Log << type() << " " << name() <<  " execute:" << nl;

    // Calculate
    averagingIndex_++;

    calculate();

    if (averagingIndex_ == nAveragingSteps_)
    {
        averagingIndex_ = 0;
        // Write
        if (Pstream::master() || !Pstream::parRun())
        {
            liftDragFilePtr_() << mesh_.time().timeName() << tab
                               << totalLift_ << tab
                               << frontLift_ << tab
                               << rearLift_ << tab
                               << drag_ << tab
                               << side_ << tab
                               << moment_.x() << tab
                               << moment_.y() << tab
                               << moment_.z() << endl;

            forAll(liftDragBinPtr_, bI)
            {
                binData& ldBin = liftDragBinPtr_[bI];

                forAll(ldBin.filePtr_, fileI)
                {
                    ldBin.filePtr_[fileI]()
                        << mesh_.time().timeName() << tab << "0.0" << tab;

                    scalar cumForce = 0;
                    vector cumMom = Zero;

                    if (fileI == 0)
                    {
                        forAll(ldBin.forceCoeffs_, i)
                        {
                            cumForce += ldBin.forceCoeffs_[i] & liftDir();
                            ldBin.filePtr_[fileI]() << cumForce << tab;
                        }
                    }
                    else if (fileI == 1)
                    {
                        forAll(ldBin.forceCoeffs_, i)
                        {
                            cumForce += ldBin.forceCoeffs_[i] & dragDir();
                            ldBin.filePtr_[fileI]() << cumForce << tab;
                        }
                    }
                    else if (fileI == 2)
                    {
                        forAll(ldBin.forceCoeffs_, i)
                        {
                            cumForce += ldBin.forceCoeffs_[i] & pitchAxis();
                            ldBin.filePtr_[fileI]() << cumForce << tab;
                        }
                    }
                    else
                    {
                        forAll(ldBin.forceCoeffs_, i)
                        {
                            cumMom += ldBin.momentCoeffs_[i];
                            ldBin.filePtr_[fileI]() << cumMom << tab;
                        }
                    }

                    ldBin.filePtr_[fileI]() << endl;

                }
            }

            if (liftDragRegions_)
            {
                forAll(liftDragRegions_.filePtr_, fileI)
                {
                    liftDragRegions_.filePtr_[fileI]()
                        << mesh_.time().timeName() << tab;

                    if (fileI == 0)
                    {
                        forAll(liftDragRegions_.forceCoeffs_, i)
                        {
                            liftDragRegions_.filePtr_[fileI]()
                                << (liftDir()
                                    & liftDragRegions_.forceCoeffs_[i])
                                << tab;
                        }
                    }
                    else if (fileI == 1)
                    {
                        forAll(liftDragRegions_.forceCoeffs_, i)
                        {
                            liftDragRegions_.filePtr_[fileI]()
                                << (dragDir()
                                    & liftDragRegions_.forceCoeffs_[i])
                                << tab;
                        }
                    }
                    else if (fileI == 2)
                    {
                        forAll(liftDragRegions_.forceCoeffs_, i)
                        {
                            liftDragRegions_.filePtr_[fileI]()
                                << (pitchAxis()
                                    & liftDragRegions_.forceCoeffs_[i])
                                << tab;
                        }
                    }
                    else
                    {
                        forAll(liftDragRegions_.momentCoeffs_, i)
                        {
                            liftDragRegions_.filePtr_[fileI]()
                                << liftDragRegions_.momentCoeffs_[i]
                                << tab;
                        }
                    }

                    liftDragRegions_.filePtr_[fileI]() << endl;

                }
            }
        }
    }

    Info<< endl;

    return true;
}


void Foam::functionObjects::liftDrag::calculate()
{
    resetFields();

    // Check whether UName and pName exists, if not deactivate lift/Drag
    if
    (
        !foundObject<volVectorField>(UName())
     || !foundObject<volScalarField>(pName())
     || (
            rhoName_ != "rhoInf"
         && !foundObject<volScalarField>(rhoName_)
        )
    )
    {
        WarningInFunction
            <<"Not calculating "<< type() << " " << name()
            <<" Failure of lookup of one of following fields : "
            << UName() <<" "<< pName() <<" "<< rhoName_ << endl;
        return;
    }

    // Get U
    const volVectorField& U = lookupObject<volVectorField>(UName());

    // Get p (because p might be p/rho or simply p)
    tmp<volScalarField> p = getp();

    // Reset bins if moving case
    if (mesh_.moving() || coorFramePtr_)
    {
        calculateBins();
    }

    // Non-dimensionalizing coefficients
    scalar qRef = 0.5*rhoRef()*sqr(magUInf());
    scalar Fref = qRef*Aref();

    // Averaging coefficients
    scalar alpha = scalar(averagingIndex_ - 1.0)/scalar(averagingIndex_);
    scalar beta = 1.0/scalar(averagingIndex_);

    // Limit pressure contribution
    scalar maxP = maxCp_*qRef;
    scalar minP = minCp_*qRef;

    vector pressureForce = Zero;
    vector viscousForce = Zero;
    vector porousForce = Zero;
    vector pressureMoment = Zero;
    vector viscousMoment = Zero;
    vector porousMoment = Zero;

    // Calculate viscous and pressure forces and moments
    label regionI = 0;

    // Scale bin data for averaging purposes
    forAll(liftDragBinPtr_, bI)
    {
        binData& ldBin = liftDragBinPtr_[bI];
        ldBin.forceCoeffs_ *= alpha;
        ldBin.momentCoeffs_ *= alpha;
    }

    tmp<volSymmTensorField> tdevRhoReff = forces::devRhoReff();
    const volSymmTensorField::Boundary& devRhoReffb =
        tdevRhoReff().boundaryField();

    const surfaceVectorField::Boundary& Sfb =
        mesh_.Sf().boundaryField();

    forAllConstIter(labelHashSet, patchSet(), iter)
    {
        const label patchi = iter.key();

        vectorField Md
        (
            mesh_.C().boundaryField()[patchi] - csys().origin()
        );

        vectorField pf
        (
            max(min(p->boundaryField()[patchi], maxP), minP)
           *mesh_.Sf().boundaryField()[patchi]
        );

        const vector pressurePatchForce = gSum(pf);
        const vector pressurePatchMoment = gSum(Md ^ pf);

        vectorField vf(Sfb[patchi] & devRhoReffb[patchi]);

        // Add to force/moment fields
        vectorField fP(Md.size(), Zero);
        addToFields(patchi, Md, pf, vf, fP);

        const vector viscousPatchForce = gSum(vf);
        const vector viscousPatchMoment = gSum(Md ^ vf);

        // Calc regional lift/drag/moment coefficients
        if (liftDragRegions_)
        {
            liftDragRegions_.forceCoeffs_[regionI] *= alpha;
            liftDragRegions_.forceCoeffs_[regionI] +=
                beta*(viscousPatchForce + pressurePatchForce)/Fref;

            liftDragRegions_.momentCoeffs_[regionI] *= alpha;
            liftDragRegions_.momentCoeffs_[regionI] +=
                beta*(viscousPatchMoment + pressurePatchMoment)
                /Fref/lRef();
        }

        forAll(liftDragBinPtr_, bI)
        {
            binData& ldBin = liftDragBinPtr_[bI];

            const vectorField& Cf = mesh_.C().boundaryField()[patchi];
            forAll(Cf, fI)
            {

                label binI = binIndex(ldBin, Cf[fI]);

                ldBin.forceCoeffs_[binI] += beta * (pf[fI] + vf[fI])/Fref;

                ldBin.momentCoeffs_[binI] +=
                    beta*(((Cf[fI] - csys().origin())^pf[fI])
                  + ((Cf[fI] - csys().origin())^vf[fI]))/Fref/lRef();

            }
        }

        pressureForce += pressurePatchForce;
        pressureMoment += pressurePatchMoment;

        viscousForce += viscousPatchForce;
        viscousMoment += viscousPatchMoment;

        regionI++;
    }

    // Add porous forces
    if (porosity_)
    {
        const volScalarField rho(this->rho());
        const volScalarField mu(this->mu());

        const HashTable<const porosityModel*> models =
            obr_.lookupClass<porosityModel>();

        forAll(porousZoneNames_, i)
        {
            // Non-const access required if mesh is changing
            porosityModel& pm =
                const_cast<porosityModel&>(*models[porousZoneNames_[i]]);

            vectorField fPTot(pm.force(U, rho, mu));

            const labelList& cellZoneIDs = pm.cellZoneIDs();

            vector porousZoneForce = Zero;
            vector porousZoneMoment = Zero;

            forAll(cellZoneIDs, i)
            {
                label zonei = cellZoneIDs[i];
                const cellZone& cZone = mesh_.cellZones()[zonei];

                bool calculateForce = checkZone(cZone);

                if (calculateForce)
                {
                    const vectorField d(mesh_.C(), cZone);
                    const vectorField fP(fPTot, cZone);
                    const vectorField Md(d - csys().origin());

                    // Add porous force
                    const vectorField fDummy(Md.size(), Zero);
                    addToFields(cZone, Md, fDummy, fDummy, fP);

                    forAll(liftDragBinPtr_, bI)
                    {
                        binData& ldBin = liftDragBinPtr_[bI];

                        forAll(cZone, czI)
                        {
                            label cellI = cZone[czI];

                            const point& Cc = mesh_.C()[cellI];
                            label binI = binIndex(ldBin, Cc);

                            ldBin.forceCoeffs_[binI] +=
                                beta * fPTot[cellI]/Fref;

                            ldBin.momentCoeffs_[binI] +=
                                beta*((Cc - csys().origin())^fPTot[cellI])
                               /Fref/lRef();
                        }
                    }
                    porousZoneForce += sum(fP);
                    porousZoneMoment += sum(Md^fP);
                }
            }

            reduce(porousZoneForce, sumOp<vector>());
            reduce(porousZoneMoment, sumOp<vector>());

            porousForce += porousZoneForce;
            porousMoment += porousZoneMoment;

            if (liftDragRegions_)
            {
                liftDragRegions_.forceCoeffs_[regionI] *= alpha;
                liftDragRegions_.forceCoeffs_[regionI] +=
                    beta*porousZoneForce/Fref;

                liftDragRegions_.momentCoeffs_[regionI] *= alpha;
                liftDragRegions_.momentCoeffs_[regionI] +=
                    beta*(porousZoneMoment)/Fref/lRef();
            }
            regionI++;
        }
    }

    // Synchronise parallel bin data
    forAll(liftDragBinPtr_, bI)
    {
        binData& ldBin = liftDragBinPtr_[bI];

        combineReduce(ldBin.forceCoeffs_, plusEqOp<vectorField>());
        combineReduce(ldBin.momentCoeffs_, plusEqOp<vectorField>());
    }

    const vector pressureForceCoeff = pressureForce/Fref;
    const vector pressureMomentCoeff = pressureMoment/Fref/lRef();

    const vector viscousForceCoeff = viscousForce/Fref;
    const vector viscousMomentCoeff = viscousMoment/Fref/lRef();

    // If posoity is added later
    const vector porousForceCoeff = porousForce/Fref;
    const vector porousMomentCoeff = porousMoment/Fref/lRef();

    Info<< "Viscous drag: " << (viscousForceCoeff & dragDir())
         << ", pressure drag: " << (pressureForceCoeff & dragDir())
         << ", porous drag: "<< (porousForceCoeff & dragDir())
         << endl;

    Info<< "Viscous lift: " << (viscousForceCoeff & liftDir())
         << ", pressure lift: " << (pressureForceCoeff & liftDir())
         << ", porous lift: "<< (porousForceCoeff & liftDir())
         << endl;

    const vector totalForceCoeff =
        viscousForceCoeff + pressureForceCoeff + porousForceCoeff;

    const scalar drag = totalForceCoeff & dragDir();
    const scalar totalLift = totalForceCoeff & liftDir();
    const scalar side = totalForceCoeff & pitchAxis();

    const vector moment =
        (pressureMomentCoeff + viscousMomentCoeff + porousMomentCoeff);

    const scalar frontLift =
        totalLift/2.0 + (moment & pitchAxis())/(wheelBase_/lRef());
    const scalar rearLift =
        totalLift/2.0 - (moment & pitchAxis())/(wheelBase_/lRef());

    // Write to screen
    Info<< "Total lift: " << totalLift
         << ", Front lift: " << frontLift
         << ", Rear lift: " << rearLift
         << ", Drag: " << drag << endl;

    drag_ = alpha*drag_ + beta*drag;
    totalLift_ = alpha*totalLift_ + beta*totalLift;
    moment_ = alpha*moment_ + beta*moment;
    side_ = alpha*side_ + beta*side;
    frontLift_ = alpha*frontLift_ + beta*frontLift;
    rearLift_ = alpha*rearLift_ + beta*rearLift;

    // Write state/results information
    {
        setResult("Cm", moment_);
        setResult("Cd", drag_);
        setResult("Cl", totalLift_);
        setResult("Cl(f)", frontLift_);
        setResult("Cl(r)", rearLift_);
    }

    if (cOfPressure_)
    {
        // Calculation of Center of Pressure vector
        const vector totalForce = viscousForce + pressureForce + porousForce;

        if(mag(totalForce) > VSMALL)
        {
            const vector totalMoment =
                viscousMoment + pressureMoment + porousMoment;
            point clP(Zero);

            CofM_ = centreOfMass::centre(mesh_, patchSet());
            Log << "Location of Center of Mass: " << CofM_ <<endl;

            // If explicity specified in dict, use referencePointCP
            // otherwise, use the center of mass by default
            if (CofP_refPoint_)
            {
                clP = refCofP_;
                if (coorFramePtr_ && definedInFrame_)
                {
                    clP = coorFramePtr_->coorSys().globalVector(refCofP_);
                }
            }
            else
            {
                clP = CofM_;
            }

            // Find the closest point to the reference point for CofP
            const vector tCofP =
                (totalForce ^ totalMoment)/magSqr(totalForce) + csys().origin();
            const vector totalForce_dir = totalForce/mag(totalForce);
            const vector vDiff = tCofP - clP;
            const vector vPerpendicular =
                vDiff - (totalForce_dir & vDiff)*totalForce_dir;

            CofP_ = clP + vPerpendicular;

            Log << "Location of Center Of Pressure: " << CofP_ <<endl;
            Log << endl;
        }
    }
}


bool Foam::functionObjects::liftDrag::write()
{
    return forces::write();
}


bool Foam::functionObjects::liftDrag::read(const dictionary& dict)
{
    forceCoeffs::read(dict);

    Log << type() << " " << name() <<  " read:" << nl;

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    // Deprecated exactNamed/partialNamed support
    if (dict.found("liftDragPatches"))
    {
        Warning << "Use of 'liftDragPatches' keyword in liftDrag "
                << "functionObject is deprecated." << nl << "New usage: "
                << "patches(<words/regExs/Groups>);" << endl;

        const dictionary& patchDict = dict.subDict("liftDragPatches");

        wordList namedPatches(patchDict.lookup("exactNamed"));
        const wordList allPatchNames(mesh_.boundaryMesh().names());

        forAll(namedPatches, nI)
        {
            label patchi = mesh_.boundaryMesh().findPatchID(namedPatches[nI]);

            if (patchi != -1 && !mesh_.boundaryMesh()[patchi].coupled())
            {
                if (!patchSet().found(patchi))
                {
                    patchSet().insert(patchi);
                }
            }
            else
            {
                WarningInFunction
                    << "Could not find Exact patch name " << namedPatches[nI]
                    << " excluding from liftDrag calculation"
                    << endl;
            }
        }

        wordList partialNamePatches(patchDict.lookup("partialNamed"));

        forAll(partialNamePatches, nI)
        {
            word substring = partialNamePatches[nI];

            forAll(pbm, pI)
            {
                word name = pbm[pI].name();

                if (!pbm[pI].coupled())
                {
                    if (name.find(substring, 0) != string::npos)
                    {
                        if (!patchSet().found(pI))
                        {
                            patchSet().insert(pI);
                        }
                    }
                }
            }
        }
    }

    // Porous switch
    if (porosity_)
    {
        const HashTable<const porosityModel*> models =
            obr_.lookupClass<porosityModel>();
        if (models.size() != 0)
        {
            porousZoneNames_ = models.sortedToc();
        }
    }

    maxCp_ = GREAT;
    if (dict.found("maxCp"))
    {
        maxCp_ = readScalar(dict.lookup("maxCp"));
    }

    minCp_ = -GREAT;
    if (dict.found("minCp"))
    {
        minCp_ = readScalar(dict.lookup("minCp"));
    }

    averagingIndex_ = 0;
    nAveragingSteps_ = 1;
    if (dict.found("nAveragingSteps"))
    {
        nAveragingSteps_ = readLabel(dict.lookup("nAveragingSteps"));
    }

    if (dict.found("referencePointCP"))
    {
        CofP_refPoint_ = true;
        refCofP_ = dict.lookup("referencePointCP");
    }

    if (dict.found("cOfPressure"))
    {
        cOfPressure_ = readBool(dict.lookup("cOfPressure"));
    }

    wheelBase_ = readScalar(dict.lookup("wheelbase"));
    if (wheelBase_ <= 0.0)
    {
        FatalErrorInFunction
            << "wheelbase must be greater than zero."
            << exit(FatalError);
    }


    Info<< "    Creating lift/drag file." << nl << endl;

    // File update
    if (Pstream::master() || !Pstream::parRun())
    {
        fileName liftDragFileName("liftDrag");

        // Open new file at startup
        liftDragFilePtr_ = createFile(liftDragFileName);

        // Add headers
        writeCommented(liftDragFilePtr_(), "Time");
        writeDelimited(liftDragFilePtr_(), "totalLift");
        writeDelimited(liftDragFilePtr_(), "frontLift");
        writeDelimited(liftDragFilePtr_(), "rearLift");
        writeDelimited(liftDragFilePtr_(), "drag");
        writeDelimited(liftDragFilePtr_(), "side");
        writeDelimited(liftDragFilePtr_(), "Xmoment");
        writeDelimited(liftDragFilePtr_(), "Ymoment");
        writeDelimited(liftDragFilePtr_(), "Zmoment");

        liftDragFilePtr_() << endl;
    }

    // Create liftDragRegions_ struct entries
    if (dict.found("outputRegionData"))
    {
        dynamic_cast<Switch&>(liftDragRegions_) =
            Switch(dict.lookup("outputRegionData"));

        if (liftDragRegions_)
        {
            label nRegions = patchSet().size();

            if (porosity_)
            {
                nRegions += porousZoneNames_.size();
            }

            liftDragRegions_.forceCoeffs_ = vectorField(nRegions, Zero);
            liftDragRegions_.momentCoeffs_ = vectorField(nRegions, Zero);

            // File update
            if (Pstream::master() || !Pstream::parRun())
            {
                // File for each of lift, drag and side force and for moment
                // components
                liftDragRegions_.filePtr_.setSize(4);
                liftDragRegions_.filePtr_[0] = createFile("regionalLift");
                liftDragRegions_.filePtr_[1] = createFile("regionalDrag");
                liftDragRegions_.filePtr_[2] = createFile("regionalSide");
                liftDragRegions_.filePtr_[3] = createFile("regionalMoment");

                // Add headers
                forAll(liftDragRegions_.filePtr_, fileI)
                {
                    writeCommented
                    (
                        liftDragRegions_.filePtr_[fileI](),
                        "Time"
                    );

                    forAllConstIter(labelHashSet, patchSet(), iter)
                    {
                        label patchi = iter.key();
                        word patchName = pbm[patchi].name();

                        writeDelimited
                        (
                            liftDragRegions_.filePtr_[fileI](),
                            patchName
                        );
                    }

                    if (porosity_)
                    {
                        const HashTable<const porosityModel*> models =
                            obr_.lookupClass<porosityModel>();
                        forAll(porousZoneNames_, zoneI)
                        {
                            const porosityModel& pm =
                                const_cast<porosityModel&>
                                (*models[porousZoneNames_[zoneI]]);

                            writeDelimited
                            (
                                liftDragRegions_.filePtr_[fileI](),
                                pm.name()
                            );
                        }
                    }

                    liftDragRegions_.filePtr_[fileI]()<< endl;
                }
            }
        }
    }

    // Create liftDragBinPtr_ struct entries (disabled for dynamic runs)
    if (dict.found("binData"))
    {
        const dictionary& binDict = dict.subDict("binData");

        // Check if old format is used and write warning
        const wordList binNames(binDict.toc());
        bool oldFormat(true);
        liftDragBinPtr_.resize(binNames.size());

        vector normalDir = Zero;
        const twoDPointCorrector& twoDCorrector =
            twoDPointCorrector::New(mesh_);
        if (twoDCorrector.required())
        {
            normalDir = twoDCorrector.planeNormal();
        }

        label nonEmptyIndex = 0;

        forAll(binNames, bI)
        {
            if (binDict.isDict(binNames[bI]))
            {
                const dictionary& binDictAxis = binDict.subDict(binNames[bI]);
                const vector binAxis = binDictAxis.lookup("axis");
                vector checkBinAxis = binAxis;
                if (coorFramePtr_ && definedInFrame_)
                {
                    checkBinAxis =
                        coorFramePtr_->coorSys().globalVector(binAxis);
                }

                if (mag(normalDir & normalised(checkBinAxis)) < SMALL)
                {
                    liftDragBinPtr_.set(nonEmptyIndex, new binData(true));

                    binData& ldBin = liftDragBinPtr_[nonEmptyIndex];

                    ldBin.name_ = binNames[bI] ;

                    ldBin.binAxis_ = binAxis;
                    ldBin.nBins_ = readLabel(binDictAxis.lookup("nBins"));

                    if (binDictAxis.found("axisFormat"))
                    {
                        ldBin.axisFormat_ =
                            axisTypeNames_.read
                            (
                                binDictAxis.lookup("axisFormat")
                            );
                    }

                    // Set data list sizes
                    ldBin.forceCoeffs_ = vectorField(ldBin.nBins_, Zero);
                    ldBin.momentCoeffs_ = vectorField(ldBin.nBins_, Zero);

                    oldFormat = false;
                    nonEmptyIndex++;
                }
                else
                {
                    Info<< "liftDrag: Omitting " << binNames[bI]
                        << " as its axis direction is aligned "
                        << "with the empty direction!" << endl;
                }
            }
            else
            {
                break;
            }
        }

        // Remove empty directions
        liftDragBinPtr_.resize(nonEmptyIndex);

        if (oldFormat)
        {
            WarningInFunction
                << "Old format used in binData." << endl;

            liftDragBinPtr_.resize(1);

            liftDragBinPtr_.set(0, new binData(true));

            liftDragBinPtr_[0].name_ = word("bin");

            liftDragBinPtr_[0].binAxis_ = binDict.lookup("axis");
            liftDragBinPtr_[0].nBins_ = readLabel(binDict.lookup("nBins"));

            // Set data list sizes
            liftDragBinPtr_[0].forceCoeffs_ =
                vectorField(liftDragBinPtr_[0].nBins_, Zero);
            liftDragBinPtr_[0].momentCoeffs_ =
                vectorField(liftDragBinPtr_[0].nBins_, Zero);
        }

        calculateBins();

        // File ptr creation
        if (Pstream::master() || !Pstream::parRun())
        {
            forAll(liftDragBinPtr_, bI)
            {
                binData& ldBin = liftDragBinPtr_[bI];

                // File for each of lift, drag and side force and for moment
                // components
                ldBin.filePtr_.setSize(4);
                ldBin.filePtr_[0] = createFile(ldBin.name_ + "BinLift");
                ldBin.filePtr_[1] = createFile(ldBin.name_ + "BinDrag");
                ldBin.filePtr_[2] = createFile(ldBin.name_ + "BinSide");
                ldBin.filePtr_[3] = createFile(ldBin.name_ + "BinMoment");

                // Add headers
                forAll(ldBin.filePtr_, fileI)
                {
                    writeCommented(ldBin.filePtr_[fileI](), "Time");
                    writeDelimited(ldBin.filePtr_[fileI](), "0.0");

                    switch (ldBin.axisFormat_)
                    {
                        case global :
                        {
                            for (label binI = 0; binI < ldBin.nBins_; binI++)
                            {
                                scalar binLen
                                (
                                    (1.0 + scalar(binI))*(ldBin.binSize_)
                                  + ldBin.binMin_
                                );
                                ldBin.filePtr_[fileI]() << tab << binLen;
                            }
                            break;
                        }

                        case boundBox :
                        {
                            for (label binI = 0; binI < ldBin.nBins_; binI++)
                            {
                                scalar binLen
                                (
                                    (1.0 + scalar(binI))*(ldBin.binSize_)
                                );
                                ldBin.filePtr_[fileI]() << tab << binLen;
                            }
                            break;
                        }

                        case normalized :
                        {
                            for (label binI = 0; binI < ldBin.nBins_; binI++)
                            {
                                scalar binLen
                                (
                                    (1.0 + scalar(binI))
                                   *(ldBin.binSize_)
                                   /(ldBin.binMax_ - ldBin.binMin_)
                                );
                                ldBin.filePtr_[fileI]() << tab << binLen;
                            }
                        }
                    }

                    ldBin.filePtr_[fileI]() << endl;
                }
            }
        }
    }

    return true;
}

// ************************************************************************* //
