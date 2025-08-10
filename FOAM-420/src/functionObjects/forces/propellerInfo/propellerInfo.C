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
    (c) 2021 OpenCFD Ltd.
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/fvOptions/fvOptions.H"
#include "propellerInfo/propellerInfo.H"
#include "fvMesh/fvMesh.H"
#include "coordinate/systems/cylindricalCS.H"
#include "interpolation/interpolation/interpolation/interpolation.H"
#include "algorithms/indexedOctree/treeDataCell.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(propellerInfo, 0);
    addToRunTimeSelectionTable(functionObject, propellerInfo, dictionary);
}
}

const Foam::Enum<Foam::functionObjects::propellerInfo::rotationMode>
Foam::functionObjects::propellerInfo::rotationModeNames_
({
    { rotationMode::SPECIFIED, "specified" }
});


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::propellerInfo::setCoordinateSystem
(
    const dictionary& dict
)
{
    vector origin(Zero);
    vector axis(Zero);

    switch (rotationMode_)
    {
        case rotationMode::SPECIFIED:
        {
            origin = dict.lookup("CofR");
            axis = dict.lookup("axis");
            axis.normalise();

            n_ = readScalar(dict.lookup("n"));
            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unhandled enumeration "
                << rotationModeNames_[rotationMode_]
                << abort(FatalError);
        }
    }

    vector alphaAxis;
    if (!dict.readIfPresent("alphaAxis", alphaAxis))
    {
        // Value has not been set - find vector orthogonal to axis

        vector cand(Zero);
        scalar minDot = GREAT;
        for (direction d = 0; d < 3; ++d)
        {
            vector test(Zero);
            test[d] = 1;
            scalar dotp = mag(test & axis);
            if (dotp < minDot)
            {
                minDot = dotp;
                cand = test;
            }
        }

        alphaAxis = axis ^ cand;
    }

    alphaAxis.normalise();

    coordSys_ = coordSystem::cylindrical(origin, axis, alphaAxis);
}


void Foam::functionObjects::propellerInfo::setRotationalSpeed()
{
    switch (rotationMode_)
    {
        case rotationMode::SPECIFIED:
        {
            // Set on dictionary re-read
            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unhandled enumeration "
                << rotationModeNames_[rotationMode_]
                << abort(FatalError);
        }
    }
}


void Foam::functionObjects::propellerInfo::createFiles()
{
    if (!writeToFile())
    {
        return;
    }

    if (writePropellerPerformance_ && !propellerPerformanceFilePtr_)
    {
        propellerPerformanceFilePtr_ = createFile("propellerPerformance");
        auto& os = propellerPerformanceFilePtr_();

        writeHeader(os, "Propeller performance");
        writeHeaderValue(os, "CofR", csys().origin());
        writeHeaderValue(os, "Radius", radius_);
        writeHeaderValue(os, "Axis", csys().e3());

        writeHeader(os, "");

        writeCommented(os, "Time");
        writeDelimited(os, "n");
        writeDelimited(os, "URef");
        writeDelimited(os, "J");
        writeDelimited(os, "KT");
        writeDelimited(os, "10*KQ");
        writeDelimited(os, "eta0");
        os  << nl;
    }

    if (writeWakeFields_)
    {
        if (!wakeFilePtr_) wakeFilePtr_ = createFile("wake");
        if (!axialWakeFilePtr_) axialWakeFilePtr_ = createFile("axialWake");
    }
}


const Foam::volVectorField& Foam::functionObjects::propellerInfo::U() const
{
    const auto* UPtr = mesh_.cfindObject<volVectorField>(UName_);

    if (!UPtr)
    {
        FatalErrorInFunction
            << "Unable to find velocity field " << UName_
            << " . Available vector fields are: "
            << mesh_.names<volVectorField>()
            << exit(FatalError);
    }

    return *UPtr;
}


Foam::scalar Foam::functionObjects::propellerInfo::getURefValue() const
{
    if (!URefPtr_.valid())
    {
        const IOdictionary& state = mesh_.time().functionObjects().stateDict();

        bool hasResultObject =
            (
                state.found("results")
            ? state.subDict("results").found(name())
            : false
            );

        if (!hasResultObject)
        {
            if (hasDefaultValue_)
            {
                DebugInformation
                    << "    Function object " << name()
                    << " not found; using default value " << defaultValue_
                    << endl;

                return defaultValue_;
            }

            wordList objectResultNames;
            if (state.found("results"))
            {
                objectResultNames = state.subDict("results").sortedToc();
            }

            FatalErrorInFunction
                << "Function object " << name()
                << " results not found. Valid objects with results include: "
                << objectResultNames
                << exit(FatalError);
        }

        const dictionary& objectDict =
            state.subDict("results").subDict(name());

        DynamicList<word> resultEntries;
        bool hasResultObjectEntry = false;

        for (const entry& dEntry : objectDict)
        {
            const dictionary& dict = dEntry.dict();

            resultEntries.append(dict.toc());

            if (dict.found(foResultName_))
            {
                hasResultObjectEntry = true;
            }
        }

        if (!hasResultObjectEntry)
        {
            if (hasDefaultValue_)
            {
                DebugInformation
                    << "    Function object " << name()
                    << " result " << foResultName_
                    << " not found; using default value " << defaultValue_
                    << endl;

                return defaultValue_;
            }

            FatalErrorInFunction
                << "Function object " << name()
                << " does not have a result field " << foResultName_ << nl
                << " Available result fields include: "
                << resultEntries
                << exit(FatalError);
        }

        scalar value = scalar(Zero);
        const word& dictTypeName = pTraits<scalar>::typeName;

        if (objectDict.found(dictTypeName))
        {
            const dictionary& resultTypeDict =
                objectDict.subDict(dictTypeName);

            (void)resultTypeDict.readIfPresent<scalar>(foResultName_, value);
        }

        DebugInformation
            << "    Using " << name() << " function object value: " << value
            << endl;

        return value;
    }
    else
    {
        return URefPtr_->value(time_.timeOutputValue());
    }
}


void Foam::functionObjects::propellerInfo::setSampleDiskGeometry
(
    const coordinateSystem& coordSys,
    const scalar r1,
    const scalar r2,
    const scalar nTheta,
    const label nRadius,
    faceList& faces,
    pointField& points
) const
{
    label nPoint = nRadius*nTheta;
    if (r1 < SMALL)
    {
        nPoint += 1;    // 1 for origin
    }
    else
    {
        nPoint += nTheta;
    }
    const label nFace = nRadius*nTheta;

    points.setSize(nPoint);
    faces.setSize(nFace);

    const point& origin = coordSys.origin();
    const scalar zCoord = 0;
    label pointi = 0;

    for (int radiusi = 0; radiusi <= nRadius; ++radiusi)
    {
        if (r1 < SMALL && radiusi == 0)
        {
            points[pointi++] = origin;
        }
        else
        {
            const scalar r = r1 + radiusi*(r2 - r1)/nRadius;

            for (label i = 0; i < nTheta; ++i)
            {
                point p
                (
                    r,
                    (i/scalar(nTheta))*constant::mathematical::twoPi,
                    zCoord
                );

                points[pointi++] = coordSys.globalPosition(p);
            }
        }
    }

    const List<label> ptIDs(identity(nTheta));

    // Faces
    label facei = 0;
    label pointOffset0 = 0;
    label radiusOffset = 0;
    DynamicList<label> facePts(4);

    for (int radiusi = 0; radiusi < nRadius; ++radiusi)
    {
        if (r1 < SMALL && radiusi == 0)
        {
            radiusOffset = 1;
            pointOffset0 = 1;

            // Adding faces as a fan about the origin
            for (label thetai = 1; thetai <= nTheta; ++thetai)
            {
                facePts.clear();

                // Append triangle about origin
                facePts.append(0);
                facePts.append(thetai);
                facePts.append(1 + ptIDs.fcIndex(thetai - 1));

                faces[facei++] = face(facePts);
            }
        }
        else
        {
            for (label thetai = 0; thetai < nTheta; ++thetai)
            {
                facePts.clear();

                label offset = pointOffset0 + (radiusi-radiusOffset)*nTheta;

                // Inner
                facePts.append(offset + ptIDs.fcIndex(thetai - 1));
                facePts.append(offset + ptIDs.fcIndex(thetai));

                // Outer
                facePts.append(offset + nTheta + ptIDs.fcIndex(thetai));
                facePts.append(offset + nTheta + ptIDs.fcIndex(thetai - 1));

                faces[facei++] = face(facePts);
            }
        }
    }
}


void Foam::functionObjects::propellerInfo::setSampleDiskSurface
(
    const dictionary& dict
)
{
    const dictionary& sampleDiskDict(dict.subDict("sampleDisk"));

    const scalar r1 = readScalar(sampleDiskDict.lookup("r1"));
    const scalar r2 = readScalar(sampleDiskDict.lookup("r2"));

    nTheta_ = readLabel(sampleDiskDict.lookup("nTheta"));
    nRadial_ = readLabel(sampleDiskDict.lookup("nRadial"));

    setSampleDiskGeometry
    (
        csys(),
        r1,
        r2,
        nTheta_,
        nRadial_,
        faces_,
        points_
    );

    errorOnPointNotFound_ =
        sampleDiskDict.lookupOrDefault("errorOnPointNotFound", false);

    updateSampleDiskCells();
}


void Foam::functionObjects::propellerInfo::updateSampleDiskCells()
{
    if (!writeWakeFields_)
    {
        return;
    }

    treeBoundBox bb(points_);
    bb.inflate(0.05);
    DynamicList<label> treeCellIDs(10*points_.size());

    const auto& meshCells = mesh_.cells();
    const auto& meshFaces = mesh_.faces();
    const auto& meshPoints = mesh_.points();

    forAll(meshCells, celli)
    {
        bool found = false;

        for (const label facei : meshCells[celli])
        {
            for (const label fpi : meshFaces[facei])
            {
                if (bb.contains(meshPoints[fpi]))
                {
                    found = true;
                    break;
                }
            }

            if (found)
            {
                treeCellIDs.append(celli);
                break;
            }
        }
    }

    indexedOctree<treeDataCell> tree
    (
        treeDataCell(true, mesh_, std::move(treeCellIDs), polyMesh::CELL_TETS),
        bb,
        10,
        100,
        10
    );

    cellIds_.setSize(points_.size(), -1);
    pointMask_.setSize(points_.size(), false);

    // Kick the tet base points calculation to avoid parallel comms later
    (void)mesh_.tetBasePtIs();

    const auto& cellLabels = tree.shapes().cellLabels();

    forAll(points_, pointi)
    {
        const vector& pos = points_[pointi];

//        label meshCelli = mesh_.findCell(pos);
        label treeCelli = tree.findInside(pos);

        label proci = treeCelli >= 0 ? Pstream::myProcNo() : -1;

        reduce(proci, maxOp<label>());

        pointMask_[pointi] = treeCelli != -1;

        if (proci >= 0)
        {
            cellIds_[pointi] =
                proci == Pstream::myProcNo()
              ? cellLabels[treeCelli]
              : -1;
        }
        else
        {
            if (errorOnPointNotFound_)
            {
                FatalErrorInFunction
                    << "Position " << pos << " not found in mesh"
                    << abort(FatalError);
            }
            else
            {
                DebugInformation
                    << "Position " << pos << " not found in mesh"
                    << endl;
            }
        }
    }

    Pstream::listCombineReduce(pointMask_, orOp<bool>());
}


Foam::scalar Foam::functionObjects::propellerInfo::meanSampleDiskField
(
    const scalarField& field
) const
{
    if (field.size() != points_.size())
    {
        FatalErrorInFunction
            << "Inconsistent field sizes: input:" << field.size()
            << " points:" << points_.size()
            << abort(FatalError);
    }

    scalar sumArea = 0;
    scalar sumFieldArea = 0;
    forAll(faces_, facei)
    {
        const face& f = faces_[facei];

        bool valid = true;
        scalar faceValue = 0;
        for (const label pti : f)
        {
            // Exclude contributions where sample cell for point was not found
            if (!pointMask_[pti])
            {
                valid = false;
                break;
            }
            faceValue += field[pti];
        }

        if (valid)
        {
            scalar area = f.mag(points_);
            sumArea += area;
            sumFieldArea += faceValue/f.size()*area;
        }
    }

    return sumFieldArea/(sumArea + ROOTVSMALL);
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::functionObjects::propellerInfo::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& psi,
    const Type& defaultValue
) const
{
    auto tfield = tmp<Field<Type>>::New(points_.size(), defaultValue);
    auto& field = tfield.ref();

    autoPtr<interpolation<Type>> interpolator
    (
        interpolation<Type>::New(interpolationScheme_, psi)
    );

    forAll(points_, pointi)
    {
        const label celli = cellIds_[pointi];

        if (cellIds_[pointi] != -1)
        {
            const point& position = points_[pointi];
            field[pointi] = interpolator().interpolate(position, celli);
        }
    }

    Pstream::listCombineReduce(field, maxOp<Type>());

    return tfield;
}


void Foam::functionObjects::propellerInfo::writePropellerPerformance()
{
    if (!writePropellerPerformance_)
    {
        return;
    }

    // Update n_
    setRotationalSpeed();

    const vector sumForce(sum(force_[0]) + sum(force_[1]) + sum(force_[2]));
    const vector sumMoment
    (
        sum(moment_[0]) + sum(moment_[1]) + sum(moment_[2])
    );

    const scalar diameter = 2*radius_;
    const scalar URef = getURefValue();
    const scalar j = -URef/mag(n_ + ROOTVSMALL)/diameter;
    const scalar denom = rhoRef_*sqr(n_)*pow4(diameter);
    const scalar kt = (sumForce & csys().e3())/denom;
    const scalar kq =
        -sign(n_)*(sumMoment & csys().e3())/(denom*diameter);
    const scalar etaO = kt*j/(kq*constant::mathematical::twoPi + ROOTVSMALL);

    if (writeToFile())
    {
        auto& os = propellerPerformanceFilePtr_();

        writeTime(os);
        os  << tab << n_
            << tab << URef
            << tab << j
            << tab << kt
            << tab << 10*kq
            << tab << etaO
            << nl;

        os.flush();
    }

    Log << type() << " " << name() <<  " output:" << nl
        << "    Revolutions per second, n : " << n_ << nl
        << "    Reference velocity, URef  : " << URef << nl
        << "    Advance coefficient, J    : " << j << nl
        << "    Thrust coefficient, Kt    : " << kt << nl
        << "    Torque coefficient, 10*Kq : " << 10*kq << nl
        << "    Efficiency, etaO          : " << etaO << nl
        << nl;

    // Write state/results information
    setResult("n", n_);
    setResult("URef", URef);
    setResult("Kt", kt);
    setResult("Kq", kq);
    setResult("J", j);
    setResult("etaO", etaO);
}


void Foam::functionObjects::propellerInfo::writeWake
(
    const vectorField& U,
    const scalar URef
)
{
    if (!Pstream::master()) return;

    // Velocity
    auto& os = wakeFilePtr_();

    const pointField propPoints(csys().localPosition(points_));
    const label offset =
        mag(propPoints[1][0] - propPoints[0][0]) < SMALL ? 0 : 1;
    const scalar rMax = propPoints.last()[0];

    const scalar UzMean = meanSampleDiskField(U.component(2));

    writeHeaderValue(os, "Time", time_.timeOutputValue());
    writeHeaderValue(os, "Reference velocity", URef);
    writeHeaderValue(os, "Direction", csys().e3());
    writeHeaderValue(os, "Wake", 1 - UzMean/URef);
    writeHeader(os, "");
    writeCommented(os, "r/R");
    writeDelimited(os, "alpha");
    writeDelimited(os, "(x y z)");
    writeDelimited(os, "(Ur Utheta Uz)");
    os << nl;

    for (label thetai = 0; thetai < nTheta_; ++thetai)
    {
        const scalar deg = 360*thetai/scalar(nTheta_);

        for (label radiusi = 0; radiusi <= nRadial_; ++radiusi)
        {
            label pointi = radiusi*nTheta_ + thetai + offset;

            if (radiusi == 0 && offset == 1)
            {
                // Only a single point at the centre - repeat for all thetai
                pointi = 0;
            }

            if (pointMask_[pointi])
            {
                const scalar rR = propPoints[radiusi*nTheta_][0]/rMax;

                os  << rR << tab << deg << tab
                    << points_[pointi] << tab << U[pointi] << nl;
            }
        }
    }

    writeHeader(os, "===");

    os  << endl;
}


void Foam::functionObjects::propellerInfo::writeAxialWake
(
    const vectorField& U,
    const scalar URef
)
{
    if (!Pstream::master()) return;

    // Alternative common format - axial wake component
    auto& os = axialWakeFilePtr_();

    const pointField propPoints(csys().localPosition(points_));
    const label offset =
        mag(propPoints[1][0] - propPoints[0][0]) < SMALL ? 0 : 1;
    const scalar rMax = propPoints.last()[0];

    writeHeaderValue(os, "Time", time_.timeOutputValue());

    os  << "# angle";
    for (label radiusi = 0; radiusi <= nRadial_; ++radiusi)
    {
        label pointi = radiusi*nTheta_;
        scalar r = propPoints[pointi][0];
        os  << tab << "r/R=" << r/rMax;
    }
    os  << nl;

    for (label thetai = 0; thetai < nTheta_; ++thetai)
    {
        os  << 360*thetai/scalar(nTheta_);

        for (label radiusi = 0; radiusi <= nRadial_; ++radiusi)
        {
            label pointi = radiusi*nTheta_ + thetai + offset;

            if (radiusi == 0 && offset == 1)
            {
                // Only a single point at the centre - repeat for all thetai
                pointi = 0;
            }

            if (pointMask_[pointi])
            {
                os << tab << 1 - U[pointi][2]/URef;
            }
            else
            {
                os << tab << "undefined";
            }
        }

        os  << nl;
    }

    writeHeader(os, "===");

    os  << endl;
}


void Foam::functionObjects::propellerInfo::writeWakeFields(const scalar URef)
{
    if (!writeWakeFields_)
    {
        return;
    }

    // Normalised velocity
    const vectorField UDisk(interpolate(U(), vector::uniform(nanValue_))());
    const vectorField UrDisk(csys().localVector(UDisk));

    // Write wake text files
    writeWake(UrDisk, URef);
    writeAxialWake(UrDisk, URef);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::propellerInfo::propellerInfo
(
    const word& name,
    const Time& runTime,
    const dictionary& dict,
    bool readFields
)
:
    forces(name, runTime, dict, true),
    radius_(0),
    URefPtr_(nullptr),
    foType_(dict.lookup<word>("URef")),
    foResultName_(dict.lookup<word>("functionObjectResult")),
    defaultValue_(Zero),
    hasDefaultValue_(dict.readIfPresent("defaultValue", defaultValue_)),
    rotationMode_(rotationMode::SPECIFIED),
    n_(0),
    writePropellerPerformance_(true),
    propellerPerformanceFilePtr_(nullptr),
    writeWakeFields_(true),
    nTheta_(0),
    nRadial_(0),
    points_(),
    errorOnPointNotFound_(false),
    faces_(),
    cellIds_(),
    pointMask_(),
    interpolationScheme_("cell"),
    wakeFilePtr_(nullptr),
    axialWakeFilePtr_(nullptr),
    nanValue_(pTraits<scalar>::min)
{
    if (readFields)
    {
        read(dict);
        Log << endl;
    }
}


Foam::functionObjects::propellerInfo::propellerInfo
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    bool readFields
)
:
    propellerInfo(name, obr.time(), dict, false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::propellerInfo::read(const dictionary& dict)
{
    if (forces::read(dict))
    {
        radius_ = readScalar(dict.lookup("radius"));
        rotationMode_ = rotationModeNames_.lookup("rotationMode", dict);

        // Set a scalar Function1, if not of functionObjectValue type
        if (foType_ != "functionObjectValue")
        {
            URefPtr_.reset(Function1<scalar>::New("URef", dict));
        }

        // Must be set before setting the surface
        setCoordinateSystem(dict);

        writePropellerPerformance_ = dict.lookup("writePropellerPerformance");

        writeWakeFields_ = dict.lookup("writeWakeFields");
        if (writeWakeFields_)
        {
            setSampleDiskSurface(dict);

            dict.readIfPresent("interpolationScheme", interpolationScheme_);
            dict.readIfPresent("nanValue", nanValue_);
        }

        return true;
    }

    return false;
}


bool Foam::functionObjects::propellerInfo::execute()
{
    calcForcesMoment();

    createFiles();

    if (writeWakeFields_)
    {
        // Only setting mean axial velocity result during execute:
        // wake fields are 'heavy' and controlled separately using the
        // writeControl setting
        const vectorField
            UDisk
            (
                csys().localVector
                (
                    interpolate
                    (
                        U(),
                        vector::uniform(nanValue_)
                    )()
                )
            );
        const scalar UzMean = meanSampleDiskField(UDisk.component(2));

        setResult("UzMean", UzMean);
    }

    writePropellerPerformance();

    return true;
}


bool Foam::functionObjects::propellerInfo::write()
{
    const scalar URef = getURefValue();
    writeWakeFields(URef);

    return true;
}


void Foam::functionObjects::propellerInfo::UpdateMesh(const mapPolyMesh& mpm)
{
    updateSampleDiskCells();
}


void Foam::functionObjects::propellerInfo::movePoints(const polyMesh& mesh)
{
    updateSampleDiskCells();
}


// ************************************************************************* //
