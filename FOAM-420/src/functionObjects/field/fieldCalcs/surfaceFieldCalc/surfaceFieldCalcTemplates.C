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
    (c) 2019 Esi Ltd.
    (c) 2015-2017 OpenCFD Ltd.
    (c) 2011-2017 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "fieldCalcs/surfaceFieldCalc/surfaceFieldCalc.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "surfFields/surfFields/surfFields.H"
#include "fields/volFields/volFields.H"
#include "sampledSurface/sampledSurface/sampledSurface.H"
#include "sampledSurface/writers/surfaceWriter.H"
#include "interpolation/interpolation/interpolationCellPoint/interpolationCellPoint.H"
#include "global/unitConversion/unitConversion.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //
namespace Foam
{

template<class T>
class isNotEqOp
{
public:

    void operator()(T& x, const T& y) const
    {
        const T unsetVal(-VGREAT*pTraits<T>::one);

        if (x != unsetVal)
        {
            // Keep x.

            // Note:chould check for y != unsetVal but multiple sample cells
            // already handled in read().
        }
        else
        {
            // x is not set. y might be.
            x = y;
        }
    }
};

}

template<class Type>
bool Foam::functionObjects::fieldCalcs::surfaceFieldCalc::validField
(
    const word& fieldName
) const
{
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> sf;
    typedef GeometricField<Type, fvPatchField, volMesh> vf;
    typedef DimensionedField<Type, surfGeoMesh> smt;

    return
    (
        foundObject<smt>(fieldName)
     || foundObject<vf>(fieldName)
     || (regionType_ != stSampledSurface && foundObject<sf>(fieldName))
    );
}

void Foam::functionObjects::fieldCalcs::surfaceFieldCalc::field
(
    const dictionary& dict,
    const vectorField& Sf,
    const meshedSurf& surfToWrite,
    writeFile& writer,
    autoPtr<writeFile>& probesWriter
)
{
    word fieldName(dict.lookup("fieldName"));

    Switch report(dict.lookupOrDefault<Switch>("writeReport", true));
    Switch write(dict.lookupOrDefault<Switch>("writeField", true));

    if (!report && !write)
    {
        WarningInFunction
            <<"writeReport and writeField not set deactivating : "
            << dict.lookup("operation")<<endl;
        return;
    }

    if (validField<scalar>(fieldName))
    {
        Field<scalar> values(getFieldCalcs<scalar>(fieldName, true));

        word operation(dict.lookup("operation"));
        word outputName(operation + "_" + fieldName);

        if (report)
        {
            writeReport(values, outputName, Sf, writer, probesWriter);
        }

        if (write)
        {
            writeField(values, outputName, surfToWrite, writer);
        }

    }
    else if (validField<vector>(fieldName))
    {
        Field<vector> values(getFieldCalcs<vector>(fieldName, true));

        word operation(dict.lookup("operation"));
        word outputName(operation + "_" + fieldName);

        if (report)
        {
            writeReport(values, outputName, Sf, writer, probesWriter);
        }

        if (write)
        {
            writeField(values, outputName, surfToWrite, writer);
        }
    }
}


void Foam::functionObjects::fieldCalcs::surfaceFieldCalc::axialDeviation
(
    const dictionary& dict,
    const vectorField& Sf,
    const meshedSurf& surfToWrite,
    writeFile& writer,
    autoPtr<writeFile>& probesWriter
)
{
    word fieldName(dict.lookup("fieldName"));

    Switch report(dict.lookupOrDefault<Switch>("writeReport", true));
    Switch write(dict.lookupOrDefault<Switch>("writeField", true));

    if (!report && !write)
    {
        WarningInFunction
            <<"writeReport and writeField not set deactivating : "
            << dict.lookup("operation")<<endl;
        return;
    }

    if (validField<scalar>(fieldName))
    {
        Field<scalar> values(getFieldCalcs<scalar>(fieldName, true));

        scalar meanVal = 0;
        if (dict.found("meanValue"))
        {
            meanVal = readScalar(dict.lookup("meanValue"));
        }
        else
        {
            scalar areaSum = 0;
            forAll(filterFaces_, i)
            {
                label facei = filterFaces_[i];
                scalar fa = mag(Sf[facei]);
                meanVal += (values[facei]*fa);
                areaSum += fa;
            }

            meanVal = returnReduce(meanVal, sumOp<scalar>());
            areaSum = returnReduce(areaSum, sumOp<scalar>());

            if (areaSum > SMALL)
            {
                meanVal /= areaSum;
            }
        }

        Field<scalar> deviation(values.size(), -1);

        if (mag(meanVal) > SMALL)
        {
            forAll(filterFaces_, i)
            {
                label facei = filterFaces_[i];
                deviation[facei] = 100*(values[facei]-meanVal)/meanVal;
            }
        }
        else
        {
            WarningInFunction
                <<"Not calculating axial deviation as zero mean value : "<<endl;
        }

        word operation(dict.lookup("operation"));
        word outputName(operation + "_" + fieldName);

        if (report)
        {
            writeReport(deviation, outputName, Sf, writer, probesWriter);
        }

        if (write)
        {
            writeField(deviation, outputName, surfToWrite, writer);
        }

    }
    else if (validField<vector>(fieldName))
    {
        Field<vector> values(getFieldCalcs<vector>(fieldName, true));

        //Calculate average normal
        scalar areaSum = 0;
        vector aveNorm = vector::zero;

        forAll(filterFaces_, i)
        {
            label facei = filterFaces_[i];
            scalar fa = mag(Sf[facei]);
            aveNorm += Sf[facei];
            areaSum += fa;
        }
        aveNorm = returnReduce(aveNorm, sumOp<vector>());
        areaSum = returnReduce(areaSum, sumOp<scalar>());
        if (areaSum > SMALL)
        {
            aveNorm /= areaSum;
        }

        scalar meanVal = 0;
        if (dict.found("meanValue"))
        {
            meanVal = readScalar(dict.lookup("meanValue"));
        }
        else
        {
            forAll(filterFaces_, i)
            {
                label facei = filterFaces_[i];
                scalar fa = mag(Sf[facei]);
                meanVal += fa*(values[facei]&aveNorm);
            }
            meanVal = returnReduce(meanVal, sumOp<scalar>());

            if (areaSum > SMALL)
            {
                meanVal /= areaSum;
            }
        }

        Field<scalar> deviation(values.size(), -1);

        if (mag(meanVal) > SMALL)
        {
            forAll(filterFaces_, i)
            {
                label facei = filterFaces_[i];
                deviation[facei] =
                    100*((values[facei]&aveNorm)-meanVal)/meanVal;
            }
        }
        else
        {
            WarningInFunction
                <<"Not calculating deviation as zero mean value : "<<endl;
        }

        word operation(dict.lookup("operation"));
        word outputName(operation + "_" + fieldName);

        if (report)
        {
            writeReport(deviation, outputName, Sf, writer, probesWriter);
        }

        if (write)
        {
            writeField(deviation, outputName, surfToWrite, writer);
        }
    }
}


void Foam::functionObjects::fieldCalcs::surfaceFieldCalc::swirlAngle
(
    const dictionary& dict,
    const vectorField& Cf,
    const vectorField& Sf,
    const meshedSurf& surfToWrite,
    writeFile& writer,
    autoPtr<writeFile>& probesWriter
)
{
    word uName(dict.lookupOrDefault<word>("U", "U"));
    point swirlOrigin(dict.lookup("origin"));

    word fieldName(dict.lookupOrDefault<word>("fieldName", "swirlAngle"));
    Switch report(dict.lookupOrDefault<Switch>("writeReport", true));
    Switch write(dict.lookupOrDefault<Switch>("writeField", true));

    if (!report && !write)
    {
        WarningInFunction
            <<"writeReport and writeField not set deactivating : "
            << dict.lookup("operation")<<endl;
        return;
    }

    vector dir = vector::zero;
    vector averagePt = vector::zero;
    scalar areaSum = 0;
    forAll(filterFaces_, i)
    {
        label facei = filterFaces_[i];
        scalar fa = mag(Sf[facei]);
        dir += Sf[facei];
        areaSum += fa;
        averagePt += fa*Cf[facei];
    }

    areaSum = returnReduce(areaSum, sumOp<scalar>());
    dir = returnReduce(dir, sumOp<vector>());
    averagePt = returnReduce(averagePt, sumOp<vector>());

    if (areaSum > SMALL)
    {
        dir /= areaSum;
        dir /= mag(dir) + SMALL;
        averagePt /= areaSum;
    }

    plane pl(averagePt, dir);
    swirlOrigin = pl.nearestPoint(swirlOrigin);

    Field<vector> uValues(getFieldCalcs<vector>(uName, true));
    Field<scalar> swirl(uValues.size(), -1);

    forAll(filterFaces_, i)
    {
        label facei = filterFaces_[i];
        point fc = Cf[facei];

        vector dR = fc - swirlOrigin;
        scalar radius = mag(dR);

        vector axialVel = (uValues[facei] & dir) * dir;

        if (radius > SMALL)
        {
            vector radialDir(dR/radius);
            vector tanDir(dir ^ radialDir);
            tanDir /= (mag(tanDir) + SMALL);
            scalar tangentialVel(uValues[facei] & tanDir);
            swirl[facei] = tangentialVel;
        }
        else
        {
            swirl[facei] = mag(axialVel-uValues[facei]);
        }
    }

    scalar axialVelSum = 0;
    if (dict.found("meanAxialVelocity"))
    {
        axialVelSum = readScalar(dict.lookup("meanAxialVelocity"));
    }
    else
    {
        forAll(filterFaces_, i)
        {
            label facei = filterFaces_[i];
            scalar fa = mag(Sf[facei]);
            axialVelSum += (uValues[facei] & dir)*fa;
        }

        axialVelSum = returnReduce(axialVelSum, sumOp<scalar>());
        if (areaSum > SMALL)
        {
            axialVelSum /= areaSum;
        }
    }

    if (mag(axialVelSum) > SMALL)
    {
        forAll(filterFaces_, i)
        {
            label facei = filterFaces_[i];
            swirl[facei] = radToDeg
            (
                Foam::atan
                (
                    swirl[facei]
                    /mag(axialVelSum)
                )
            );
        }
    }
    else
    {
        forAll(filterFaces_, i)
        {
            label facei = filterFaces_[i];
            swirl[facei] = scalar(-1);
        }
        WarningInFunction
            << " Swirl axial compont is zero."
            << " Disabling swirl angle calculation." <<endl;
    }

    word operation(dict.lookup("operation"));
    word outputName(operation + "_" + fieldName);

    if (report)
    {
        writeReport(swirl, outputName, Sf, writer, probesWriter);
    }

    if (write)
    {
        writeField(swirl, outputName, surfToWrite, writer);
    }
}


void Foam::functionObjects::fieldCalcs::surfaceFieldCalc::bulkSwirlAngle
(
    const dictionary& dict,
    const vectorField& Cf,
    const vectorField& Sf,
    const meshedSurf& surfToWrite,
    writeFile& writer,
    autoPtr<writeFile>& probesWriter
)
{
    word uName(dict.lookupOrDefault<word>("U", "U"));
    word rhoName(dict.lookupOrDefault<word>("rho", "rho"));
    point swirlOrigin(dict.lookup("origin"));
    scalar swirlRadius(readScalar(dict.lookup("radius")));

    word fieldName(dict.lookupOrDefault<word>("fieldName", "swirlAngle"));
    Switch report(dict.lookupOrDefault<Switch>("writeReport", true));
    Switch write(dict.lookupOrDefault<Switch>("writeField", true));

    if (!report && !write)
    {
        WarningInFunction
            <<"writeReport and writeField not set deactivating : "
            << dict.lookup("operation") <<endl;
        return;
    }

    vector dir = vector::zero;
    vector averagePt = vector::zero;
    scalar areaSum = 0;
    forAll(filterFaces_, i)
    {
        label facei = filterFaces_[i];
        scalar fa = mag(Sf[facei]);
        dir += Sf[facei];
        averagePt += fa*Cf[facei];
        areaSum += fa;
    }

    areaSum = returnReduce(areaSum, sumOp<scalar>());
    dir = returnReduce(dir, sumOp<vector>());
    averagePt = returnReduce(averagePt, sumOp<vector>());

    if (areaSum > SMALL)
    {
        dir /= areaSum;
        dir /= mag(dir) + SMALL;
        averagePt /= areaSum;
    }

    plane pl(averagePt, dir);
    //Make sure swirl origin is on the plane
    swirlOrigin = pl.nearestPoint(swirlOrigin);

    Field<vector> uValues(getFieldCalcs<vector>(uName, true));
    Field<scalar> swirl(uValues.size(), -1);
    Field<scalar> rho(getFieldCalcs<scalar>(rhoName, false));

    scalar angMomSum = 0;
    scalar inertialSum = 0;

    forAll(filterFaces_, i)
    {
        label facei = filterFaces_[i];
        scalar fa = mag(Sf[facei]);
        point fc = Cf[facei];

        vector dR = fc - swirlOrigin;
        scalar radius = mag(dR);

        vector axialVel = (uValues[facei] & dir) * dir;

        if (radius > SMALL)
        {
            vector radialDir(dR/radius);
            vector radialVel = (uValues[facei] & radialDir) *radialDir;

            vector tangentialVel = uValues[facei] - radialVel
                - axialVel;

            scalar dM = mag(axialVel)*fa;
            if (rho.size() > 0)
            {
                dM *= rho[facei];
            }

            angMomSum += ((dR ^ tangentialVel) & dir) * dM;
            inertialSum += pow(radius,2) * dM;
        }
    }

    scalar axialVelSum = 0;
    if (dict.found("meanAxialVelocity"))
    {
        axialVelSum = readScalar(dict.lookup("meanAxialVelocity"));
    }
    else
    {
        forAll(filterFaces_, i)
        {
            label facei = filterFaces_[i];
            scalar fa = mag(Sf[facei]);
            axialVelSum += (uValues[facei] & dir)*fa;
        }

        axialVelSum = returnReduce(axialVelSum, sumOp<scalar>());
        if (areaSum > SMALL)
        {
            axialVelSum /= areaSum;
        }
    }

    angMomSum = returnReduce(angMomSum, sumOp<scalar>());
    inertialSum = returnReduce(inertialSum, sumOp<scalar>());

    if (inertialSum > SMALL && mag(axialVelSum) > SMALL)
    {
        scalar angVel = mag(angMomSum/inertialSum);
        scalar swirlAngle = radToDeg
        (
            Foam::atan
            (
                angVel*swirlRadius
                /mag(axialVelSum)
            )
        );

        forAll(filterFaces_, i)
        {
            label facei = filterFaces_[i];
            swirl[facei] = swirlAngle;
        }
    }
    else
    {
        if (inertialSum < SMALL)
        {
            WarningInFunction
                << " Bulk swirl Inertial component is zero."
                << " Disabling bulk swirl calculation." <<endl;
        }
        else if (mag(axialVelSum) < SMALL)
        {
            WarningInFunction
                << " Bulk swirl axial compont is zero."
                << " Disabling bulk swirl calculation." <<endl;
        }
    }

    word operation(dict.lookup("operation"));
    word outputName(operation + "_" + fieldName);

    if (report)
    {
        writeReport(swirl, outputName, Sf, writer, probesWriter);
    }

    if (write)
    {
        writeField(swirl, outputName, surfToWrite, writer);
    }
}


Foam::scalar Foam::functionObjects::fieldCalcs::surfaceFieldCalc::integrate
(
    const Field<scalar> values,
    const vectorField& Sf
)
{
    scalar intergralVal = Zero;

    forAll(filterFaces_, i)
    {
        label facei = filterFaces_[i];
        intergralVal += (values[facei]*mag(Sf[facei]));
    }

    reduce(intergralVal, sumOp<scalar>());

    return intergralVal;
}


Foam::scalar Foam::functionObjects::fieldCalcs::surfaceFieldCalc::integrate
(
    const Field<vector> values,
    const vectorField& Sf
)
{
    scalar intergralVal = Zero;

    forAll(filterFaces_, i)
    {
        label facei = filterFaces_[i];
        intergralVal += (values[facei]&Sf[facei]);
    }

    reduce(intergralVal, sumOp<scalar>());

    return intergralVal;
}


template<class Type>
void Foam::functionObjects::fieldCalcs::surfaceFieldCalc::writeReport
(
    const Field<Type> values,
    const word& fieldName,
    const vectorField& Sf,
    writeFile& writer,
    autoPtr<writeFile>& probesWriter
)
{
    Type aveVal = Zero;
    scalar aveMag = 0;

    scalar maxVal = -GREAT;
    scalar minVal = GREAT;
    scalar stdDev = 0;
    scalar area = 0;

    bool nonScalar = true;
    if (string(pTraits<Type>::typeName) == string("scalar"))
    {
        nonScalar = false;
    }

    forAll(filterFaces_, i)
    {
        label facei = filterFaces_[i];
        Type fv = values[facei];
        scalar fa = mag(Sf[facei]);

        aveVal += (fv*fa);
        aveMag += (mag(fv)*fa);

        area += fa;
        if (nonScalar)
        {
            maxVal = max(mag(fv),maxVal);
            minVal = min(mag(fv),minVal);
        }
        else
        {
            maxVal = max(cmptAv(fv),maxVal);
            minVal = min(cmptAv(fv),minVal);
        }
    }

    reduce(area, sumOp<scalar>());
    reduce(aveMag, sumOp<scalar>());
    reduce(aveVal, sumOp<Type>());
    reduce(maxVal, maxOp<scalar>());
    reduce(minVal, minOp<scalar>());

    if (area)
    {
        aveVal /= area;
        aveMag /= area;
    }
    else
    {
        aveVal = 0.;
        aveMag = 0.;
    }

    forAll(filterFaces_, i)
    {
        label facei = filterFaces_[i];
        Type fv = values[facei];
        scalar fa = mag(Sf[facei]);

        stdDev += (sqr(mag(fv- aveVal))*fa);
    }

    reduce(stdDev, sumOp<scalar>());

    if (area > SMALL)
    {
        stdDev /= area;
        stdDev = sqrt(stdDev);
    }
    else
    {
        stdDev = 0.;
    }

    //Calculate value integral
    scalar integral = integrate(values,Sf);

    writer.file()<< minVal << tab << maxVal << tab << aveMag << tab
                 << aveVal << tab << stdDev << tab << integral <<endl;

    // Write state/results information
    word prefix, suffix;
    {
        prefix += fieldName;
        prefix += '(';
        suffix += ')';
    }

    Log << "    " << prefix << regionName_ << suffix
        << " of " << fieldName
        << " Min " << minVal << " Max "<< maxVal
        << " Ave(mag) " << aveMag <<" Ave "<<aveVal
        << " StdDev "<< stdDev << " Integral "<< integral
        << endl;

     // Write state/results information
    word resultMin = prefix + regionName_ + ',' + fieldName + suffix + "Min";
    this->setResult(resultMin, minVal);

    word resultMax = prefix + regionName_ + ',' + fieldName + suffix + "Max";
    this->setResult(resultMax, maxVal);

    word resultAve = prefix + regionName_ + ',' + fieldName + suffix + "Ave";
    this->setResult(resultAve, aveVal);

    word resultAveMag = prefix + regionName_ + ',' + fieldName + suffix
        + "Ave(Mag)";
    this->setResult(resultAve, mag(aveVal));

    word resultStdDev = prefix + regionName_ + ','
        + fieldName + suffix + "StdDev";
    this->setResult(resultStdDev, stdDev);

    word resultIntegral = prefix + regionName_ + ','
        + fieldName + suffix + "Integral";
    this->setResult(resultIntegral, integral);

    //write optional probes date
    if (probesWriter.valid())
    {
        const Type unsetVal(-VGREAT*pTraits<Type>::one);

        tmp<Field<Type>> tValues
        (
            new Field<Type>(probeIndex_.size(), unsetVal)
        );

        Field<Type>& probeValues = tValues.ref();

        forAll(probeIndex_, probei)
        {
            if (probeIndex_[probei] >= 0)
            {
                probeValues[probei] = values[probeIndex_[probei]];
            }
        }
        Pstream::listCombineGather(probeValues, isNotEqOp<Type>());
        Pstream::listCombineScatter(probeValues);

        probesWriter().file() << obr_.time().timeName();
        forAll(probeIndex_, probei)
        {
            probesWriter().file()<< tab << probeValues[probei];
        }
        probesWriter().file()<<endl;
    }
}


template<class Type>
bool Foam::functionObjects::fieldCalcs::surfaceFieldCalc::writeField
(
    const Field<Type> values,
    const word& fieldName,
    const meshedSurf& surfToWrite,
    const writeFile& writer
)
{
    // Write raw values on surface if specified
    if (surfaceWriterPtr_.valid())
    {
        Field<Type> allValues(values);
        combineFields(allValues);

        fileName outputDir =
            writer.baseFileDir()/name()/"surface"/time_.timeName();

        if (Pstream::master())
        {
            surfaceWriterPtr_->write
            (
                outputDir,
                regionTypeNames_[regionType_] + ("_" + regionName_),
                surfToWrite,
                fieldName,
                allValues,
                false
            );
        }

        return true;
    }

    return false;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::functionObjects::fieldCalcs::surfaceFieldCalc::getFieldCalcs
(
    const word& fieldName,
    const bool mustGet
) const
{
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> sf;
    typedef GeometricField<Type, fvPatchField, volMesh> vf;
    typedef DimensionedField<Type, surfGeoMesh> smt;

    if (foundObject<smt>(fieldName))
    {
        return lookupObject<smt>(fieldName);
    }
    else if (regionType_ != stSampledSurface && foundObject<sf>(fieldName))
    {
        return filterField(lookupObject<sf>(fieldName));
    }
    else if (foundObject<vf>(fieldName))
    {
        const vf& fld = lookupObject<vf>(fieldName);

        if (surfacePtr_.valid())
        {
            if (surfacePtr_().interpolate())
            {
                const interpolationCellPoint<Type> interp(fld);
                tmp<Field<Type>> tintFld(surfacePtr_().interpolate(interp));
                const Field<Type>& intFld = tintFld();

                // Average
                const faceList& faces = surfacePtr_().faces();
                tmp<Field<Type>> tavg
                (
                    new Field<Type>(faces.size(), Zero)
                );
                Field<Type>& avg = tavg.ref();

                forAll(faces, facei)
                {
                    const face& f = faces[facei];
                    forAll(f, fp)
                    {
                        avg[facei] += intFld[f[fp]];
                    }
                    avg[facei] /= f.size();
                }

                return tavg;
            }
            else
            {
                return surfacePtr_().sample(fld);
            }
        }
        else
        {
            return filterField(fld);
        }
    }

    if (mustGet)
    {
        FatalErrorInFunction
            << "Field " << fieldName << " not found in database"
            << abort(FatalError);
    }

    return tmp<Field<Type>>(new Field<Type>(0));
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::functionObjects::fieldCalcs::surfaceFieldCalc::filterField
(
    const GeometricField<Type, fvPatchField, volMesh>& field
) const
{
    tmp<Field<Type>> tvalues(new Field<Type>(faceId_.size()));
    Field<Type>& values = tvalues.ref();

    forAll(values, i)
    {
        label facei = faceId_[i];
        label patchi = facePatchId_[i];
        if (patchi >= 0)
        {
            values[i] = field.boundaryField()[patchi][facei];
        }
        else
        {
            FatalErrorInFunction
                << type() << " " << name() << ": "
                << regionTypeNames_[regionType_] << "(" << regionName_ << "):"
                << nl
                << "    Unable to process internal faces for volume field "
                << field.name() << nl << abort(FatalError);
        }
    }

    // No need to flip values - all boundary faces point outwards

    return tvalues;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::functionObjects::fieldCalcs::surfaceFieldCalc::filterField
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& field
) const
{
    tmp<Field<Type>> tvalues(new Field<Type>(faceId_.size()));
    Field<Type>& values = tvalues.ref();

    forAll(values, i)
    {
        label facei = faceId_[i];
        label patchi = facePatchId_[i];
        if (patchi >= 0)
        {
            values[i] = field.boundaryField()[patchi][facei];
        }
        else
        {
            values[i] = field[facei];
        }
    }

    if (debug)
    {
        Pout<< "field " << field.name() << " oriented: "
            << field.oriented()() << endl;
    }

    if (field.oriented()())
    {
        forAll(values, i)
        {
            if (faceFlip_[i])
            {
                values[i] *= -1;
            }
        }
    }

    return tvalues;
}


// ************************************************************************* //
