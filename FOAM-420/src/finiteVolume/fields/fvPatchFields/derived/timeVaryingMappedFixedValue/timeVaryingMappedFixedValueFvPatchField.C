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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2015 OpenCFD Ltd.
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/timeVaryingMappedFixedValue/timeVaryingMappedFixedValueFvPatchField.H"
#include "db/Time/Time.H"
#include "fields/fvPatchFields/derived/timeVaryingMappedFixedValue/AverageField.H"
#include "db/IOstreams/Fstreams/IFstream.H"
#include "f4st/read/f4stReader.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/fvPatchFields/fvPatchField/fieldMappers/gibFvPatchFieldMapper.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::timeVaryingMappedFixedValueFvPatchField<Type>::
timeVaryingMappedFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    fieldTableName_(iF.name()),
    setAverage_(false),
    cylindricalCoords_(false),
    perturb_(0),
    delabella_(true),
    isPtr_(nullptr),
    f4stFormat_(false),
    mapperPtr_(nullptr),
    samplePoints_(nullptr),
    sampleTimes_(0),
    startSampleTime_(-1),
    startSampledValues_(0),
    startAverage_(Zero),
    endSampleTime_(-1),
    endSampledValues_(0),
    endAverage_(Zero),
    offset_(),
    timeOffset_(0)
{}


template<class Type>
Foam::timeVaryingMappedFixedValueFvPatchField<Type>::
timeVaryingMappedFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF, dict, false),
    fieldTableName_(iF.name()),
    setAverage_(dict.lookupOrDefault("setAverage", false)),
    cylindricalCoords_(dict.lookupOrDefault("cylindricalCoords", false)),
    perturb_(dict.lookupOrDefault("perturb", 1e-6)),
    delabella_(dict.lookupOrDefault("delabella", true)),
    mapMethod_
    (
        dict.lookupOrDefault<word>
        (
            "mapMethod",
            "planarInterpolation"
        )
    ),
    isPtr_(nullptr),
    f4stFormat_(dict.lookupOrDefault("f4st", false)),
    mapperPtr_(nullptr),
    samplePoints_(nullptr),
    sampleTimes_(0),
    startSampleTime_(-1),
    startSampledValues_(0),
    startAverage_(Zero),
    endSampleTime_(-1),
    endSampledValues_(0),
    endAverage_(Zero),
    offset_(),
    timeOffset_(dict.lookupOrDefault<scalar>("timeOffset", 0.0))
{
    if (dict.found("offset"))
    {
        offset_ = Function1<Type>::New("offset", dict);
    }

    if
    (
        mapMethod_ != "planarInterpolation"
     && mapMethod_ != "nearest"
    )
    {
        FatalIOErrorInFunction(dict)
            << "mapMethod should be one of 'planarInterpolation'"
            << ", 'nearest'" << exit(FatalIOError);
    }

    dict.readIfPresent("fieldTable", fieldTableName_);

    if (dict.found("value"))
    {
        fvPatchField<Type>::forceAssign(Field<Type>("value", dict, p.size()));
    }
    else
    {
        // Note: we use evaluate() here to trigger updateCoeffs followed
        //       by re-setting of fvatchfield::updated_ flag. This is
        //       so if first use is in the next time step it retriggers
        //       a new update.
        this->evaluate(Pstream::commsTypes::blocking);
    }

    if (f4stFormat_)
    {
        word name_ = "surfaceData";
        char dataFile[80];
        sprintf(dataFile, "f4st/%s.f4", name_.c_str());

        isPtr_.reset(new IFstream(dataFile, IOstream::BINARY));
    }
}


template<class Type>
Foam::timeVaryingMappedFixedValueFvPatchField<Type>::
timeVaryingMappedFixedValueFvPatchField
(
    const timeVaryingMappedFixedValueFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    fieldTableName_(ptf.fieldTableName_),
    setAverage_(ptf.setAverage_),
    cylindricalCoords_(ptf.cylindricalCoords_),
    perturb_(ptf.perturb_),
    delabella_(ptf.delabella_),
    mapMethod_(ptf.mapMethod_),
    isPtr_(nullptr),
    f4stFormat_(ptf.f4stFormat_),
    mapperPtr_(nullptr),
    samplePoints_(nullptr),
    sampleTimes_(0),
    startSampleTime_(-1),
    startSampledValues_(0),
    startAverage_(Zero),
    endSampleTime_(-1),
    endSampledValues_(0),
    endAverage_(Zero),
    offset_(ptf.offset_, false),
    timeOffset_(ptf.timeOffset_)
{
    if (ptf.isPtr_.valid())
    {
        isPtr_.reset(new IFstream(ptf.isPtr_->name(), IOstream::BINARY));
    }
}


template<class Type>
Foam::timeVaryingMappedFixedValueFvPatchField<Type>::
timeVaryingMappedFixedValueFvPatchField
(
    const timeVaryingMappedFixedValueFvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf),
    fieldTableName_(ptf.fieldTableName_),
    setAverage_(ptf.setAverage_),
    cylindricalCoords_(ptf.cylindricalCoords_),
    perturb_(ptf.perturb_),
    delabella_(ptf.delabella_),
    mapMethod_(ptf.mapMethod_),
    isPtr_(nullptr),
    f4stFormat_(ptf.f4stFormat_),
    mapperPtr_(nullptr),
    samplePoints_(nullptr),
    sampleTimes_(ptf.sampleTimes_),
    startSampleTime_(ptf.startSampleTime_),
    startSampledValues_(ptf.startSampledValues_),
    startAverage_(ptf.startAverage_),
    endSampleTime_(ptf.endSampleTime_),
    endSampledValues_(ptf.endSampledValues_),
    endAverage_(ptf.endAverage_),
    offset_(ptf.offset_, false),
    timeOffset_(ptf.timeOffset_)
{
    if (ptf.isPtr_.valid())
    {
        isPtr_.reset(new IFstream(ptf.isPtr_->name(), IOstream::BINARY));
    }
}


template<class Type>
Foam::timeVaryingMappedFixedValueFvPatchField<Type>::
timeVaryingMappedFixedValueFvPatchField
(
    const timeVaryingMappedFixedValueFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    fieldTableName_(ptf.fieldTableName_),
    setAverage_(ptf.setAverage_),
    cylindricalCoords_(ptf.cylindricalCoords_),
    perturb_(ptf.perturb_),
    delabella_(ptf.delabella_),
    mapMethod_(ptf.mapMethod_),
    isPtr_(nullptr),
    f4stFormat_(ptf.f4stFormat_),
    mapperPtr_(nullptr),
    samplePoints_(nullptr),
    sampleTimes_(ptf.sampleTimes_),
    startSampleTime_(ptf.startSampleTime_),
    startSampledValues_(ptf.startSampledValues_),
    startAverage_(ptf.startAverage_),
    endSampleTime_(ptf.endSampleTime_),
    endSampledValues_(ptf.endSampledValues_),
    endAverage_(ptf.endAverage_),
    offset_(ptf.offset_, false),
    timeOffset_(ptf.timeOffset_)
{
    if (ptf.isPtr_.valid())
    {
        isPtr_.reset(new IFstream(ptf.isPtr_->name(), IOstream::BINARY));
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::timeVaryingMappedFixedValueFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchField<Type>::autoMap(m);
    if (startSampledValues_.size())
    {
        m(startSampledValues_, startSampledValues_);
        m(endSampledValues_, endSampledValues_);
    }
    // Clear interpolator
    mapperPtr_.clear();
    startSampleTime_ = -1;
    endSampleTime_ = -1;
}


template<class Type>
void Foam::timeVaryingMappedFixedValueFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchField<Type>::rmap(ptf, addr);

    const timeVaryingMappedFixedValueFvPatchField<Type>& tiptf =
        refCast<const timeVaryingMappedFixedValueFvPatchField<Type>>(ptf);

    startSampledValues_.rmap(tiptf.startSampledValues_, addr);
    endSampledValues_.rmap(tiptf.endSampledValues_, addr);

    // Clear interpolator
    mapperPtr_.clear();
    startSampleTime_ = -1;
    endSampleTime_ = -1;
}


template<class Type>
void Foam::timeVaryingMappedFixedValueFvPatchField<Type>::autoMapGIB
(
    const gibFvPatchFieldMapper& mapper
)
{
    fixedValueFvPatchField<Type>::autoMapGIB(mapper);

    if (startSampledValues_.size())
    {
        mapper.map(startSampledValues_, pTraits<Type>::zero);
        mapper.map(endSampledValues_, pTraits<Type>::zero);
    }
    // Clear interpolator
    mapperPtr_.clear();
    startSampleTime_ = -1;
    endSampleTime_ = -1;
}


template<class Type>
Foam::tmp<Foam::pointField>
Foam::timeVaryingMappedFixedValueFvPatchField<Type>::makePointsGlobal
(
    pointField& samplePoints
)
{
    tmp<pointField> globalPoints(samplePoints);
    if (this->coorFramePtr() && this->isDefinedInFrame())
    {
        if (cylindricalCoords_)
        {
            FatalErrorInFunction
                << "Cylindrical coordinates not allowed "
                << "for use with reference frame. "
                << exit(FatalError);
        }
        globalPoints =
            this->coorFramePtr()->coorSys().globalPosition(samplePoints);
    }

    return globalPoints;
}


template<class Type>
void Foam::timeVaryingMappedFixedValueFvPatchField<Type>::checkTable()
{
    // Initialise
    if (mapperPtr_.empty())
    {
        // Reread values and interpolate
        fileName samplePointsFile
        (
            this->db().time().path()
           /this->db().time().caseConstant()
           /"boundaryData"
           /this->patch().name()
           /"points"
        );

        samplePoints_.reset(new pointField(IFstream(samplePointsFile)()));

        DebugInformation
            << "timeVaryingMappedFixedValueFvPatchField :"
            << " Read " << samplePoints_().size() << " sample points from "
            << samplePointsFile << endl;


        if (cylindricalCoords_)
        {
            // x is r(m), y is theta(radians) and z is z(m)
            forAll(samplePoints_(), pI)
            {
                point& pntI = samplePoints_()[pI];

                scalar x = pntI.x() * Foam::cos(pntI.y());
                scalar y = pntI.x() * Foam::sin(pntI.y());

                pntI.x() = x;
                pntI.y() = y;
            }

            DebugInformation
                << "timeVaryingMappedFixedValueFvPatchField :"
                << " Converting cylindrical to cartesian coordinates "
                << endl;
        }


        // tbd: run-time selection
        bool nearestOnly =
        (
           !mapMethod_.empty()
         && mapMethod_ != "planarInterpolation"
        );

        // Allocate the interpolator
        mapperPtr_.reset
        (
            new pointToPointPlanarInterpolation
            (
                makePointsGlobal(samplePoints_()),
                this->patch().patch().faceCentres(),
                perturb_,
                delabella_,
                nearestOnly
            )
        );

        // Read the times for which data is available
        const fileName samplePointsDir = samplePointsFile.path();
        sampleTimes_ = Time::findTimes(samplePointsDir);

        DebugInformation
            << "timeVaryingMappedFixedValueFvPatchField : In directory "
            << samplePointsDir << " found times "
            << pointToPointPlanarInterpolation::timeNames(sampleTimes_)
            << endl;
    }

    if
    (
        this->coorFramePtr()
     && this->isDefinedInFrame()
     && this->coorFramePtr()->anyDynamic()
    )
    {
        mapperPtr_.reset
        (
            new pointToPointPlanarInterpolation
            (
                makePointsGlobal(samplePoints_()),
                this->patch().patch().faceCentres(),
                perturb_,
                delabella_,
                mapMethod_ != "planarInterpolation"
            )
        );
    }


    // Find current time in sampleTimes
    label lo = -1;
    label hi = -1;

    bool foundTime = mapperPtr_().findTime
    (
        sampleTimes_,
        startSampleTime_,
        this->db().time().value() + timeOffset_,
        lo,
        hi
    );

    if (!foundTime)
    {
        FatalErrorInFunction
            << "Cannot find starting sampling values for current time "
            << this->db().time().value() + timeOffset_ << nl
            << "Have sampling values for times "
            << pointToPointPlanarInterpolation::timeNames(sampleTimes_) << nl
            << "In directory "
            <<  this->db().time().constant()/"boundaryData"/this->patch().name()
            << "\n    on patch " << this->patch().name()
            << " of field " << fieldTableName_
            << exit(FatalError);
    }


    // Update sampled data fields.

    if (lo != startSampleTime_)
    {
        startSampleTime_ = lo;

        if (startSampleTime_ == endSampleTime_)
        {
            // No need to reread since are end values
            if (debug)
            {
                Pout<< "checkTable : Setting startValues to (already read) "
                    <<   "boundaryData"
                        /this->patch().name()
                        /sampleTimes_[startSampleTime_].name()
                    << endl;
            }
            startSampledValues_ = endSampledValues_;
            startAverage_ = endAverage_;
        }
        else
        {
            if (debug)
            {
                Pout<< "checkTable : Reading startValues from "
                    <<   "boundaryData"
                        /this->patch().name()
                        /sampleTimes_[lo].name()
                    << endl;
            }


            // Reread values and interpolate
            fileName valsFile
            (
                this->db().time().path()
               /this->db().time().caseConstant()
               /"boundaryData"
               /this->patch().name()
               /sampleTimes_[startSampleTime_].name()
               /fieldTableName_
            );

            Field<Type> vals;

            if (setAverage_)
            {
                AverageField<Type> avals((IFstream(valsFile)()));
                vals = avals;
                startAverage_ = avals.average();
            }
            else
            {
                IFstream(valsFile).operator()() >> vals;
            }

            if (vals.size() != mapperPtr_().sourceSize())
            {
                FatalErrorInFunction
                    << "Number of values (" << vals.size()
                    << ") differs from the number of points ("
                    <<  mapperPtr_().sourceSize()
                    << ") in file " << valsFile << exit(FatalError);
            }

            startSampledValues_ = mapperPtr_().interpolate(vals);
        }
    }

    if (hi != endSampleTime_)
    {
        endSampleTime_ = hi;

        if (endSampleTime_ == -1)
        {
            // endTime no longer valid. Might as well clear endValues.
            if (debug)
            {
                Pout<< "checkTable : Clearing endValues" << endl;
            }
            endSampledValues_.clear();
        }
        else
        {
            if (debug)
            {
                Pout<< "checkTable : Reading endValues from "
                    <<   "boundaryData"
                        /this->patch().name()
                        /sampleTimes_[endSampleTime_].name()
                    << endl;
            }

            // Reread values and interpolate
            fileName valsFile
            (
                this->db().time().path()
               /this->db().time().caseConstant()
               /"boundaryData"
               /this->patch().name()
               /sampleTimes_[endSampleTime_].name()
               /fieldTableName_
            );

            Field<Type> vals;

            if (setAverage_)
            {
                AverageField<Type> avals((IFstream(valsFile)()));
                vals = avals;
                endAverage_ = avals.average();
            }
            else
            {
                IFstream(valsFile).operator()() >> vals;
            }

            if (vals.size() != mapperPtr_().sourceSize())
            {
                FatalErrorInFunction
                    << "Number of values (" << vals.size()
                    << ") differs from the number of points ("
                    <<  mapperPtr_().sourceSize()
                    << ") in file " << valsFile << exit(FatalError);
            }

            endSampledValues_ = mapperPtr_().interpolate(vals);
        }
    }
}

template<class Type>
void Foam::timeVaryingMappedFixedValueFvPatchField<Type>::checkTableF4st()
{
    // Initialise
    word name_ = "surfaceData";
    char dataFile[80];
    sprintf(dataFile, "f4st/%s.f4", name_.c_str());

    if (mapperPtr_.empty())
    {
        // Reread values and interpolate
        samplePoints_.reset(new pointField(f4stReader::getPoints(isPtr_())));

        DebugInformation
            << "timeVaryingMappedFixedValueFvPatchField :"
            << " Read " << samplePoints_().size() << " sample points from "
            << dataFile << endl;


        if (cylindricalCoords_)
        {
            // x is r(m), y is theta(radians) and z is z(m)
            forAll(samplePoints_(), pI)
            {
                point& pntI = samplePoints_()[pI];

                scalar x = pntI.x() * Foam::cos(pntI.y());
                scalar y = pntI.x() * Foam::sin(pntI.y());

                pntI.x() = x;
                pntI.y() = y;
            }

            DebugInformation
                << "timeVaryingMappedFixedValueFvPatchField :"
                << " Converting cylindrical to cartesian coordinates "
                << endl;
        }


        // tbd: run-time selection
        bool nearestOnly =
            (
                !mapMethod_.empty()
                    && mapMethod_!="planarInterpolation"
            );

        // Allocate the interpolator
        mapperPtr_.reset
        (
            new pointToPointPlanarInterpolation
            (
                makePointsGlobal(samplePoints_()),
                this->patch().patch().faceCentres(),
                perturb_,
                delabella_,
                nearestOnly
            )
        );

        // Read the times for which data is available
        const fileName samplePointsDir = "file";
        sampleTimes_ = f4stReader::getSampleTimes(dataFile);
        DebugInformation
            << "timeVaryingMappedFixedValueFvPatchField : In directory "
            << samplePointsDir << " found times "
            << pointToPointPlanarInterpolation::timeNames(sampleTimes_)
            << endl;
    }

    if
    (
        this->coorFramePtr()
     && this->isDefinedInFrame()
     && this->coorFramePtr()->anyDynamic()
    )
    {
        mapperPtr_.reset
        (
            new pointToPointPlanarInterpolation
            (
                makePointsGlobal(samplePoints_()),
                this->patch().patch().faceCentres(),
                perturb_,
                delabella_,
                mapMethod_ != "planarInterpolation"
            )
        );
    }

    // Find current time in sampleTimes
    label lo = -1;
    label hi = -1;

    bool foundTime = mapperPtr_().findTime
    (
        sampleTimes_,
        startSampleTime_,
        this->db().time().value() + timeOffset_,
        lo,
        hi
    );

    if (!foundTime)
    {
        FatalErrorInFunction
            << "Cannot find starting sampling values for current time "
            << this->db().time().value() + timeOffset_ << nl
            << "Have sampling values for times "
            << pointToPointPlanarInterpolation::timeNames(sampleTimes_) << nl
            << "In directory "
            <<  this->db().time().constant()/"boundaryData"/this->patch().name()
            << "\n    on patch " << this->patch().name()
            << " of field " << fieldTableName_
            << exit(FatalError);
    }


    // Update sampled data fields.

    if (lo != startSampleTime_)
    {
        startSampleTime_ = lo;

        if (startSampleTime_ == endSampleTime_)
        {
            // No need to reread since are end values
            if (debug)
            {
                Pout<< "checkTable : Setting startValues to (already read) "
                    <<   "boundaryData"
                        /this->patch().name()
                        /sampleTimes_[startSampleTime_].name()
                    << endl;
            }
            startSampledValues_ = endSampledValues_;
            startAverage_ = endAverage_;
        }
        else
        {
            if (debug)
            {
                Pout<< "checkTable : Reading startValues from "
                    <<   "boundaryData"
                        /this->patch().name()
                        /sampleTimes_[lo].name()
                    << endl;
            }
            word fieldName =
                this->internalField().name() + "_"
                  + sampleTimes_[startSampleTime_].name();

            Field<Type> vals
            (
                f4stReader::getField<Type>
                (
                    fieldName, isPtr_(), mapperPtr_().sourceSize()
                )
            );

            if (setAverage_)
            {
                Type average = gAverage(vals);
                AverageField<Type> avals(vals, average);
                vals = avals;
                startAverage_ = avals.average();
            }

            if (vals.size() != mapperPtr_().sourceSize())
            {
                FatalErrorInFunction
                    << "Number of values (" << vals.size()
                    << ") differs from the number of points ("
                    <<  mapperPtr_().sourceSize()
                    << ") in file " << dataFile << exit(FatalError);
            }

            startSampledValues_ = mapperPtr_().interpolate(vals);
        }
    }

    if (hi != endSampleTime_)
    {
        endSampleTime_ = hi;

        if (endSampleTime_ == -1)
        {
            // endTime no longer valid. Might as well clear endValues.
            if (debug)
            {
                Pout<< "checkTable : Clearing endValues" << endl;
            }
            endSampledValues_.clear();
        }
        else
        {
            if (debug)
            {
                Pout<< "checkTable : Reading endValues from "
                    <<   "boundaryData"
                        /this->patch().name()
                        /sampleTimes_[endSampleTime_].name()
                    << endl;
            }

            // Reread values and interpolate
            word fieldName =
                this->internalField().name() + "_"
                    + sampleTimes_[endSampleTime_].name();
            Field<Type> vals
            (
                f4stReader::getField<Type>
                    (
                        fieldName, isPtr_(), mapperPtr_().sourceSize()
                    )
            );

            if (setAverage_)
            {
                Type average = gAverage(vals);
                AverageField<Type> avals(vals, average);
                vals = avals;
                endAverage_ = avals.average();
            }

            endSampledValues_ = mapperPtr_().interpolate(vals);
        }
    }
}

template<class Type>
void Foam::timeVaryingMappedFixedValueFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (f4stFormat_)
    {
        checkTableF4st();
    }
    else
    {
        checkTable();
    }

    // Interpolate between the sampled data

    Type wantedAverage;

    if (endSampleTime_ == -1)
    {
        // Only start value
        if (debug)
        {
            Pout<< "updateCoeffs : Sampled, non-interpolated values"
                << " from start time:"
                << sampleTimes_[startSampleTime_].name() << nl;
        }

        this->forceAssign(startSampledValues_);
        wantedAverage = startAverage_;
    }
    else
    {
        const scalar start = sampleTimes_[startSampleTime_].value();
        const scalar end = sampleTimes_[endSampleTime_].value();

        const scalar s =
            (this->db().time().value() + timeOffset_ - start)/(end - start);

        if (debug)
        {
            Pout<< "updateCoeffs : Sampled, interpolated values"
                << " between start time:"
                << sampleTimes_[startSampleTime_].name()
                << " and end time:" << sampleTimes_[endSampleTime_].name()
                << " with weight:" << s << endl;
        }

        this->forceAssign((1 - s)*startSampledValues_ + s*endSampledValues_);
        wantedAverage = (1 - s)*startAverage_ + s*endAverage_;
    }

    // Enforce average. Either by scaling (if scaling factor > 0.5) or by
    // offsetting.
    if (setAverage_)
    {
        const Field<Type>& fld = *this;

        const Type averagePsi =
            gSum(this->patch().magSf()*fld)
           /gSum(this->patch().magSf());

        if (debug)
        {
            Pout<< "updateCoeffs :"
                << " actual average:" << averagePsi
                << " wanted average:" << wantedAverage
                << endl;
        }

        if (mag(averagePsi) < VSMALL)
        {
            // Field too small to scale. Offset instead.
            const Type offset = wantedAverage - averagePsi;
            if (debug)
            {
                Pout<< "updateCoeffs :"
                    << " offsetting with:" << offset << endl;
            }
            this->forceAssign(fld + offset);
        }
        else
        {
            const scalar scale = mag(wantedAverage)/mag(averagePsi);

            if (debug)
            {
                Pout<< "updateCoeffs :"
                    << " scaling with:" << scale << endl;
            }
            this->forceAssign(scale*fld);
        }
    }

    // Apply offset to mapped values
    if (offset_.valid())
    {
        const scalar t = this->db().time().timeOutputValue();
        this->forceAssign(*this + offset_->value(t));
    }
    this->makeVectorGlobal(*this);
    this->forceAssign(*this + getFrameVelocity());

    if (debug)
    {
        Pout<< "updateCoeffs : set fixedValue to min:" << gMin(*this)
            << " max:" << gMax(*this)
            << " avg:" << gAverage(*this) << endl;
    }

    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::timeVaryingMappedFixedValueFvPatchField<Type>::write
(
    Ostream& os
) const
{
    fixedValueFvPatchField<Type>::write(os);
    this->writeEntryIfDifferent(os, "setAverage", Switch(false), setAverage_);
    this->writeEntryIfDifferent(os, "perturb", scalar(1e-6), perturb_);
    this->writeEntryIfDifferent(os, "delabella", Switch(true), delabella_);
    this->writeEntryIfDifferent(os, "f4st", Switch(false), f4stFormat_);
    this->writeEntryIfDifferent
    (
        os,
        "fieldTable",
        this->internalField().name(),
        fieldTableName_
    );
    this->writeEntryIfDifferent
    (
        os,
        "mapMethod",
        word("planarInterpolation"),
        mapMethod_
    );

    if (offset_.valid())
    {
        offset_->writeData(os);
    }

    this->writeEntryIfDifferent(os, "timeOffset", scalar(0.0), timeOffset_);
}


// ************************************************************************* //
