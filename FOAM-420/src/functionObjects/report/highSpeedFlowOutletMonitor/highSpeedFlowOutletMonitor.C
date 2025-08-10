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
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "highSpeedFlowOutletMonitor/highSpeedFlowOutletMonitor.H"
#include "meshes/polyMesh/syncTools/syncTools.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "basicThermo/basicThermo.H"
#include "fluidThermo/fluidThermo.H"
#include "materialModels/baseModels/materialModels.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(highSpeedFlowOutletMonitor, 0);
    addToRunTimeSelectionTable(functionObject, highSpeedFlowOutletMonitor, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::labelList functionObjects::highSpeedFlowOutletMonitor::patchIDs
(
    const wordList& patchNames
)
{
    labelList pIDs(patchNames.size());
    forAll(pIDs, i)
    {
        pIDs[i] = mesh_.boundaryMesh().findPatchID(patchNames[i]);
        if (pIDs[i] == -1)
        {
            FatalErrorInFunction
                << "Function object: " << this->typeName
                << " didn't find boundary named:" << patchNames[i] << nl
                << exit(FatalError);
        }
    }
    return pIDs;
}


Foam::labelList functionObjects::highSpeedFlowOutletMonitor::zoneIDs
(
    const wordList& zoneNames
)
{
    labelList zIDs(zoneNames.size());
    forAll(zIDs, i)
    {
        zIDs[i] = mesh_.faceZones().findZoneID(zoneNames[i]);
        if (zIDs[i] == -1)
        {
            FatalErrorInFunction
                << "Function object: " << this->typeName
                << " didn't find boundary named:" << zoneNames[i] << nl
                << exit(FatalError);
        }
    }
    return zIDs;
}


Foam::scalar
Foam::functionObjects::highSpeedFlowOutletMonitor::sumPatchAreas
(
    const labelList& patchIDs
)
{
    scalar patchAreaSum = 0;
    forAll(patchIDs, i)
    {
        patchAreaSum += gSum(magSf_.boundaryField()[patchIDs[i]]);
    }
    return patchAreaSum;
}


Foam::scalar
Foam::functionObjects::highSpeedFlowOutletMonitor::sumZoneAreas
(
    const labelList& zoneIDs
)
{
    scalar zoneAreaSum = 0;
    forAll(zoneIDs, i)
    {
        const faceZone& fZone = mesh_.faceZones()[zoneIDs[i]];
        zoneAreaSum += gSum(fZone().magFaceAreas());

    }
    return zoneAreaSum;
}


void Foam::functionObjects::highSpeedFlowOutletMonitor::faceInterpolate()
{
    if (inFaceZoneIDs_.size() > 0 || outFaceZoneIDs_.size() > 0)
    {
        pFaces_ = fvc::interpolate(p_);
        TFaces_ = fvc::interpolate(T_);
        UFaces_ = fvc::interpolate(U_);
        rhoFaces_ = fvc::interpolate(rho_);
    }
}


void Foam::functionObjects::highSpeedFlowOutletMonitor::clearTmpFields()
{
    pFaces_.clear();
    TFaces_.clear();
    UFaces_.clear();
    rhoFaces_.clear();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

functionObjects::highSpeedFlowOutletMonitor::highSpeedFlowOutletMonitor
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name, typeName, dict),
    p_(obr_.lookupObject<volScalarField>("p")),
    T_(obr_.lookupObject<volScalarField>("T")),
    U_(obr_.lookupObject<volVectorField>("U")),
    phi_(obr_.lookupObject<surfaceScalarField>("phi")),
    rho_(obr_.lookupObject<volScalarField>("rho")),
    magSf_(mesh_.magSf()),
    mat_
    (
        obr_.subRegistry("materialModels").lookupObject<materialTables>
        (
            "materialTables"
        )
    )
{
    read(dict);
    writeFileHeader(file());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool functionObjects::highSpeedFlowOutletMonitor::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    // Function object works only for compressible cases
    // Identification might be trickier for mixtures/mixed psi/rho models
    // It requires that Cp isn't equal to Cv
    const scalar gammaMinusOne =
        gSum(mat_(gammaModel::typeName).primitiveField()*mesh_.V())/gSum(mesh_.V()) - 1.0;

    if (mag(gammaMinusOne) < SMALL)
    {
        FatalErrorInFunction
            << "\"" << this->name()
            << "\" requires compressible case. "
            << "Ratio of specific heats (Cp/Cv) can't be 1."
            << exit(FatalError);
    }

    // Inlet patches
    inletPatches_ =
        dict.lookupOrDefault<wordList>("inletPatches", wordList::null());
    inPatchIDs_ = patchIDs(inletPatches_);
    inletFaceZones_ =
        dict.lookupOrDefault<wordList>("inletZones", wordList::null());
    inFaceZoneIDs_ = zoneIDs(inletFaceZones_);

    // Outlet patches
    outletPatches_ =
        dict.lookupOrDefault<wordList>("outletPatches", wordList::null());
    outPatchIDs_ = patchIDs(outletPatches_);
    outletFaceZones_ =
        dict.lookupOrDefault<wordList>("outletZones", wordList::null());
    outFaceZoneIDs_ = zoneIDs(outletFaceZones_);

    if
    (
        (inPatchIDs_.size() + inFaceZoneIDs_.size()) < 1
     || (outPatchIDs_.size() + outFaceZoneIDs_.size()) < 1
    )
    {
        FatalErrorInFunction
            << "\"" << this->name()
            << "\" requires at least one inlet and one outlet (patch or zone)."
            << exit(FatalError);
    }

    // Calculate areas of inlet patches
    sumInAreas_ = sumPatchAreas(inPatchIDs_) + sumZoneAreas(inFaceZoneIDs_);

    // Calculate areas of outlet patches
    sumOutAreas_ = sumPatchAreas(outPatchIDs_) + sumZoneAreas(outFaceZoneIDs_);

    return true;
}


bool Foam::functionObjects::highSpeedFlowOutletMonitor::execute()
{
    Log << type() << " " << name() <<  " execute:" << nl;
    calculate();
    return true;
}


void functionObjects::highSpeedFlowOutletMonitor::calculate()
{
    // If the face zones are available create interpolated fields
    faceInterpolate();

    // Update areas
    if (mesh_.changing())
    {
        sumInAreas_ =
            sumPatchAreas(inPatchIDs_) + sumZoneAreas(inPatchIDs_);
        sumOutAreas_ =
            sumPatchAreas(outPatchIDs_) + sumZoneAreas(outPatchIDs_);
    }

    // Lookup reference pressure
    const scalar pRef =
            obr_.lookupObject<refScalarField>("pRef").offset().value();

    // Calculate surface area averages
    const scalar inP(areaAverage(p_, pFaces_) + pRef);
    const scalar outP(areaAverage(p_, pFaces_, false) + pRef);
    const scalar inU(areaAverage(U_, UFaces_));
    const scalar outU(areaAverage(U_, UFaces_, false));
    const scalar inletT(areaAverage(T_, TFaces_));
    const scalar inRho(areaAverage(rho_, rhoFaces_));
    const scalar outRho(areaAverage(rho_, rhoFaces_, false));

    // Calculate inlet mass flow rate
    const scalar phi(calcPhi(phi_, inPatchIDs_, inFaceZoneIDs_));

    // Specific gas constant [m^2*s^-2*K^-1]
    const scalar R(mat_(RModel::typeName)[0]);

    // Volume averaged gamma
    const scalar gamma(calcGamma());

    // Upstream stagnation temperature
    const scalar T0(calcT0());

    // Upstream stagnation density
    const scalar rho0 = inRho*pow(T0/inletT, 1.0/(gamma - 1.0));

    // Absolute pressure ratio
    const scalar Pr(outP/(inP + 0.5*inRho*sqr(inU)));

    // Trashold pressure ratio
    const scalar PrCrit
    (
        pow(2.0/(1.0 + gamma), gamma/(gamma - 1.0))
    );

    // Isentropic velocity/density at the throat
    scalar Uis, rhois;
    if (Pr <= PrCrit)
    {
        rhois = rho0*pow(2.0/(gamma + 1.0), 1.0/(gamma - 1.0));
        Uis = sqrt(2.0*gamma*R*T0/(gamma + 1.0));
    }
    else
    {
        rhois = rho0*pow(Pr, 1.0/gamma);

        // The max value has to be usually used under wrong setup
        // It is there just to not cause division by zero.
        // However, results under wrong setup will not be correct.
        Uis =
            sqrt
            (
                R*T0*2*gamma/(gamma - 1.0)
               *max((1.0 - pow(Pr, (gamma - 1.0)/gamma)), SMALL)
            );
    }

    // Effective flow area
    const scalar Aeff(phi/(rhois*Uis));

    const scalar backFlowProc = (calcBackFlowProc());

    Log << tab << "Aeff " << Aeff << " [m^2], "
        << " area in backflow " << backFlowProc << " [%]"<< endl;

    // Log to the file
    file() << obr_.time().timeName() << tab
           << Aeff << tab
           << backFlowProc << tab
           << phi << tab
           << inU << tab
           << outU << tab
           << inP << tab
           << outP << tab
           << inletT << tab
           << inRho << tab
           << outRho << tab
           << T0 << tab
           << rho0 << tab
           << Pr << tab
           << PrCrit << tab
           << rhois << tab
           << Uis << tab
           << endl;

    // Clear out tmp fields
    clearTmpFields();
}


bool functionObjects::highSpeedFlowOutletMonitor::write()
{
    return true;
}


Foam::scalar functionObjects::highSpeedFlowOutletMonitor::areaAverage
(
    const volScalarField& field,
    const tmp<surfaceScalarField>& surfaceField,
    const bool isInlet
)
{
    // Decide if inlet/outlet
    scalar area = sumInAreas_;
    labelList patchIds = inPatchIDs_;
    labelList zoneIds = inFaceZoneIDs_;
    if (!isInlet)
    {
        area = sumOutAreas_;
        patchIds = outPatchIDs_;
        zoneIds = outFaceZoneIDs_;
    }

    scalar average = 0;
    forAll(patchIds, i)
    {
        const scalarField& pField = field.boundaryField()[patchIds[i]];
        const scalarField& areas = magSf_.boundaryField()[patchIds[i]];
        average += gSum(areas*pField);
    }

    forAll(zoneIds, i)
    {
        const labelList& faces = mesh_.faceZones()[zoneIds[i]];
        scalarField zField(surfaceField.ref(), faces);
        scalarField areas(magSf_, faces);
        average += gSum(areas*zField);
    }
    return average/area;
}


Foam::scalar functionObjects::highSpeedFlowOutletMonitor::areaAverage
(
    const volVectorField& field,
    const tmp<surfaceVectorField>& surfaceField,
    const bool isInlet
)
{
    // Decide if inlet/outlet
    scalar area = sumInAreas_;
    labelList patchIds = inPatchIDs_;
    labelList zoneIds = inFaceZoneIDs_;
    if (!isInlet)
    {
        area = sumOutAreas_;
        patchIds = outPatchIDs_;
        zoneIds = outFaceZoneIDs_;
    }

    scalar average = 0;
    forAll(patchIds, i)
    {
        scalarField pField(mag(field.boundaryField()[patchIds[i]]));
        const scalarField& areas = magSf_.boundaryField()[patchIds[i]];
        average += gSum(areas*pField);
    }
    forAll(zoneIds, i)
    {
        const labelList& faces = mesh_.faceZones()[zoneIds[i]];
        vectorField zVectorField(surfaceField.ref(), faces);
        scalarField zField(mag(zVectorField));
        scalarField areas(magSf_, faces);
        average += gSum(areas*zField);
    }
    return average/area;
}


Foam::scalar functionObjects::highSpeedFlowOutletMonitor::calcPhi
(
    const surfaceScalarField& field,
    const labelList& patchIds,
    const labelList& zoneIds
)
{
    scalar phi = 0;
    forAll(patchIds, i)
    {
        phi -= gSum(field.boundaryField()[patchIds[i]]);
    }
    // Need the flip map as well
    forAll(zoneIds, i)
    {
        scalarField fieldPhi(field, mesh_.faceZones()[zoneIds[i]]);
        const boolList& fMap = mesh_.faceZones()[zoneIds[i]].flipMap();
        forAll(fMap, j)
        {
            if (fMap[j])
            {
                fieldPhi[j] *= -1;
            }
        }
        phi -= gSum(fieldPhi);
    }
    return phi;
}


Foam::scalar functionObjects::highSpeedFlowOutletMonitor::calcGamma()
{
    // Simplification internal field averaged gamma
    // To avoid complicated calculation and possible mistake in implementation
    // For analytical calculations gamma is just one number
    return
        gSum(mat_(gammaModel::typeName).primitiveField()*scalarField(mesh_.V()))
       /gSum(scalarField(mesh_.V()));
}


Foam::scalar functionObjects::highSpeedFlowOutletMonitor::calcT0()
{
    // Compressible case is: h0 = h + mag(U)^2/2
    // First h0 is calculated than T0 by inversion from energy field
    const fluidThermo& thermo =
        obr_.lookupObject<fluidThermo>(basicThermo::dictName);

    volScalarField& he = const_cast<volScalarField&>(thermo.he());

    scalar T0 = 0.0;

    forAll(inPatchIDs_, i)
    {
        const label patchi(inPatchIDs_[i]);
        scalarField h0
        (
            mat_(HsModel::typeName).boundaryField()[patchi]
          + 0.5*sqr(mag(U_.boundaryField()[patchi]))
        );
        scalarField heOld(he.boundaryField()[patchi]);
        he.boundaryFieldRef()[patchi].forceAssign(h0);

        T0 +=
            gSum
            (
                mat_(THsModel::typeName).boundaryField()[patchi]
               *magSf_.boundaryField()[patchi]
            );

        he.boundaryFieldRef()[patchi].forceAssign(heOld);
    }

    // HACK treatment of surface field through the volume field
    // This has to be done for T and p as well
    tmp<surfaceScalarField> surfaceHe;
    tmp<scalarField> heOld;
    tmp<scalarField> pOld;
    tmp<scalarField> TOld;
    volScalarField& p = const_cast<volScalarField&>(p_);
    volScalarField& T = const_cast<volScalarField&>(T_);
    if (inFaceZoneIDs_.size() > 0)
    {
        heOld = new scalarField(he.primitiveField());
        pOld = new scalarField(p_.primitiveField());
        TOld = new scalarField(T_.primitiveField());
        surfaceHe = fvc::interpolate(he);
    }

    forAll(inFaceZoneIDs_, i)
    {
        const label zonei(inFaceZoneIDs_[i]);
        const labelList& faces = mesh_.faceZones()[zonei];
        if (faces.size() > mesh_.nCells())
        {
            FatalErrorInFunction
                << "The number of faces in the faceZone: "
                << mesh_.faceZones()[zonei].name()
                << " is larger than number of cells in the mesh."
                << " This isn't courrently supported."
                << exit(FatalError);
        }

        // Initialise list 0,1,2 etc.
        labelList cells(faces.size());
        forAll(cells, celli)
        {
            cells[celli] = celli;
        }

        // Prepare face values in begining positions
        // of the cell values
        forAll(faces, j)
        {
            he[j] = surfaceHe.ref()[faces[j]];
            p[j] = pFaces_.ref()[faces[j]];
            T[j] = TFaces_.ref()[faces[j]];
        }
        scalarField h0
        (
            scalarField(mat_(HsModel::typeName).primitiveField(), cells)
          + 0.5*sqr(mag(vectorField(UFaces_, faces)))
        );

        forAll(cells, j)
        {
            he.primitiveFieldRef()[j] = h0[j];
        }
        scalarField THs(mat_(THsModel::typeName).primitiveField(), cells);
        forAll(faces, j)
        {
            THs[j] *= magSf_[faces[j]];
        }
        T0 += gSum(THs);
    }

    if (inFaceZoneIDs_.size() > 0)
    {
        he.primitiveFieldRef() = heOld.ref();
        p.primitiveFieldRef() = pOld.ref();
        T.primitiveFieldRef() = TOld.ref();
    }
    return T0/sumInAreas_;
}


Foam::scalar functionObjects::highSpeedFlowOutletMonitor::calcBackFlowProc()
{
    scalar backflowArea = 0;

    // Loop through outlet patches
    forAll(outPatchIDs_, i)
    {
        const label patchi = outPatchIDs_[i];
        const scalarField& outPatchPhi = phi_.boundaryField()[patchi];
        forAll(outPatchPhi, faceI)
        {
            if (outPatchPhi[faceI] < 0)
            {
                backflowArea += magSf_.boundaryField()[patchi][faceI];
            }
        }
    }

    // Loop through outlet face zones
    forAll(inFaceZoneIDs_, i)
    {
        const label zonei = outFaceZoneIDs_[i];
        const labelList& faces = mesh_.faceZones()[zonei];
        scalarField outZonePhi(phi_, mesh_.faceZones()[zonei]);
        const boolList& fMap = mesh_.faceZones()[zonei].flipMap();
        forAll(outZonePhi, faceI)
        {
            if (fMap[faceI])
            {
                outZonePhi[faceI] *= -1;
            }
            if (outZonePhi[faceI] < 0)
            {
                backflowArea += magSf_[faces[faceI]];
            }
        }
    }

    if (Pstream::parRun())
    {
        reduce(backflowArea, sumOp<scalar>());
    }

    scalar backFlowProcentage = 0;
    if (sumOutAreas_ > SMALL)
    {
        backFlowProcentage = (backflowArea/sumOutAreas_)*100.0;
    }
    return backFlowProcentage;
}


void Foam::functionObjects::highSpeedFlowOutletMonitor::writeFileHeader
(
    Ostream& os
) const
{
    Log << "    Logging surface statistics to file: " << fileName_
        << endl;

    writeCommented(os, "Time");
    writeDelimited(os, "Aeff");
    writeDelimited(os, "BackflowProcentage");
    writeDelimited(os, "InletFlowRate");
    writeDelimited(os, "UInlet");
    writeDelimited(os, "UOutlet");
    writeDelimited(os, "pInlet");
    writeDelimited(os, "pOutlet");
    writeDelimited(os, "TInlet");
    writeDelimited(os, "rhoInlet");
    writeDelimited(os, "rhoOutlet");
    writeDelimited(os, "T0");
    writeDelimited(os, "rho0");
    writeDelimited(os, "Pr");
    writeDelimited(os, "PrCritical");
    writeDelimited(os, "rhois");
    writeDelimited(os, "Uis");

    os << endl;
}

// ************************************************************************* //
