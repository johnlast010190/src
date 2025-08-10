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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2015-2016 OpenCFD Ltd.
    (c) 2017 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "sampledSurface/sampledSurfaces/sampledSurfaces.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "containers/Lists/ListListOps/ListListOps.H"
#include "primitives/strings/lists/stringListOps.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class GeoField, class Type>
void Foam::sampledSurfaces::returnField
(
    List<Field<Type>>& field
)
{
    const IOobjectList objects(obr(), obr().time().timeName());

    sampleAndSet<GeoField, Type>(objects, field);

    forAll(field, fldI)
    {
        if (Pstream::parRun())
        {
            // Collect values from all processors
            List<Field<Type>> gatheredValues(Pstream::nProcs());
            gatheredValues[Pstream::myProcNo()] = field[fldI];
            Pstream::gatherList(gatheredValues);

            if (Pstream::master())
            {
                // Combine values into single field
                Field<Type> allValues
                (
                    ListListOps::combine<Field<vector>>
                    (
                        gatheredValues,
                        accessOp<Field<Type>>()
                    )
                );

                field[fldI] = allValues;
            }

            Pstream::scatter(field[fldI]);
        }
    }
}


template<class GeoField, class Type>
Foam::List<Type> Foam::sampledSurfaces::returnAverage
(
    List<Field<Type>>& field
)
{
    const IOobjectList objects(obr(), obr().time().timeName());

    sampleAndSet<GeoField, Type>(objects, field);

    List<Type> avg(this->size());

    forAll(*this, surfI)
    {
        const sampledSurface& s = operator[](surfI);

        avg[surfI] = s.average(field[surfI]);
    }

    return avg;
}


template<class Type>
void Foam::sampledSurfaces::writeSurface
(
    const Field<Type>& values,
    const label surfI,
    const word& fieldName,
    const fileName& outputDir
)
{
    const sampledSurface& s = operator[](surfI);

    if (Pstream::parRun())
    {
        // Collect values from all processors
        List<Field<Type>> gatheredValues(Pstream::nProcs());
        gatheredValues[Pstream::myProcNo()] = values;
        Pstream::gatherList(gatheredValues);

        fileName sampleFile;
        if (Pstream::master())
        {
            // Combine values into single field
            Field<Type> allValues
            (
                ListListOps::combine<Field<Type>>
                (
                    gatheredValues,
                    accessOp<Field<Type>>()
                )
            );

            // Renumber (point data) to correspond to merged points
            if (mergedList_[surfI].pointsMap().size() == allValues.size())
            {
                inplaceReorder(mergedList_[surfI].pointsMap(), allValues);
                allValues.setSize(mergedList_[surfI].points().size());
            }

            // Write to time directory under outputPath_
            // skip surface without faces (eg, a failed cut-plane)
            if (mergedList_[surfI].size())
            {
                fileName inputName = outputDir;
                if (formatter_->sepFile()) inputName = mesh_.time().timeName();

                sampleFile = formatter_->write
                (
                    inputName,
                    s.name(),
                    mergedList_[surfI],
                    fieldName,
                    allValues,
                    s.interpolate()
                );
            }
        }

        Pstream::scatter(sampleFile);
        if (sampleFile.size())
        {
            dictionary propsDict;
            propsDict.add("file", sampleFile);
            setProperty(fieldName, propsDict);
        }
    }
    else
    {
        // Write to time directory under outputPath_
        // skip surface without faces (eg, a failed cut-plane)
        if (s.faces().size())
        {
            fileName inputName = outputDir;
            if (formatter_->sepFile()) inputName = mesh_.time().timeName();

            fileName fName = formatter_->write
            (
                inputName,
                s.name(),
                s,
                fieldName,
                values,
                s.interpolate()
            );

            dictionary propsDict;
            propsDict.add("file", fName);
            setProperty(fieldName, propsDict);
        }
    }
}


template<class Type>
void Foam::sampledSurfaces::sampleAndWrite
(
    const GeometricField<Type, fvPatchField, volMesh>& vField
)
{
    // interpolator for this field
    autoPtr<interpolation<Type>> interpolatorPtr;

    const word& fieldName = vField.name();

    fileName outputDir = outputPath_/vField.time().timeName();
    if (!writeTimeName_)
    {
        outputDir = outputPath_;
    }

    forAll(*this, surfI)
    {
        const sampledSurface& s = operator[](surfI);

        Field<Type> values;

        if (s.interpolate())
        {
            if (interpolatorPtr.empty())
            {
                interpolatorPtr = interpolation<Type>::New
                (
                    interpolationScheme_,
                    vField
                );
            }

            values = s.interpolate(interpolatorPtr());
        }
        else
        {
            values = s.sample(vField);
        }


        if (writeStats())
        {
            if (s.interpolate())
            {
                s.statistics(outputDir, fieldName, Field<Type>(s.sample(vField)));
            }
            else
            {
                s.statistics(outputDir, fieldName, values);
            }
        }

        if (writeFields())
        {
            writeSurface<Type>(values, surfI, fieldName, outputDir);
        }
    }
}


template<class Type>
void Foam::sampledSurfaces::sampleAndWrite
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& sField
)
{
    const word& fieldName = sField.name();
    fileName outputDir = outputPath_/sField.time().timeName();
    if (!writeTimeName_)
    {
        outputDir = outputPath_;
    }

    forAll(*this, surfI)
    {
        const sampledSurface& s = operator[](surfI);
        Field<Type> values(s.sample(sField));
        writeSurface<Type>(values, surfI, fieldName, outputDir);
    }
}


template<class GeoField>
void Foam::sampledSurfaces::sampleAndWrite(const IOobjectList& objects)
{
    wordList fieldNames;
    if (loadFromFiles_)
    {
        fieldNames = objects.sortedNames(GeoField::typeName, fieldSelection_);
    }
    else
    {
        fieldNames = obr().sortedNames<GeoField>(fieldSelection_);

        writeOriginalIds();
    }

    forAll(fieldNames, fieldi)
    {
        const word& fieldName = fieldNames[fieldi];

        if ((Pstream::master()) && verbose_)
        {
            Pout<< "sampleAndWrite: " << fieldName << endl;
        }

        if (loadFromFiles_)
        {
            const GeoField fld
            (
                IOobject
                (
                    fieldName,
                    time_.timeName(),
                    obr(),
                    IOobject::MUST_READ
                ),
                mesh_
            );

            sampleAndWrite(fld);
        }
        else
        {
            sampleAndWrite
            (
                obr().lookupObject<GeoField>(fieldName)
            );
        }
    }
}


template<class Type>
void Foam::sampledSurfaces::sampleAndSet
(
    const GeometricField<Type, fvPatchField, volMesh>& vField,
    List<Field<Type>>& field
)
{
    // interpolator for this field
    autoPtr<interpolation<Type>> interpolatorPtr;

    field.setSize(this->size());

    forAll(*this, surfI)
    {
        const sampledSurface& s = operator[](surfI);

        Field<Type> values;

        if (s.interpolate())
        {
            if (interpolatorPtr.empty())
            {
                interpolatorPtr = interpolation<Type>::New
                (
                    interpolationScheme_,
                    vField
                );
            }

            values = s.interpolate(interpolatorPtr());
        }
        else
        {
            values = s.sample(vField);
        }

        field[surfI] = values;
    }
}


template<class Type>
void Foam::sampledSurfaces::sampleAndSet
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& sField,
    List<Field<Type>>& field
)
{
    field.setSize(this->size());

    forAll(*this, surfI)
    {
        const sampledSurface& s = operator[](surfI);
        Field<Type> values(s.sample(sField));
        field[surfI] = values;
    }
}


template<class GeoField, class Type>
void Foam::sampledSurfaces::sampleAndSet
(
    const IOobjectList& objects,
    List<Field<Type>>& field
)
{
    wordList fieldNames;
    if (loadFromFiles_)
    {
        fieldNames = objects.sortedNames(GeoField::typeName, fieldSelection_);
    }
    else
    {
        fieldNames = obr().sortedNames<GeoField>(fieldSelection_);
    }

    forAll(fieldNames, fieldi)
    {
        const word& fieldName = fieldNames[fieldi];

        if ((Pstream::master()) && verbose_)
        {
            Pout<< "sampleAndSet: " << fieldName << endl;
        }

        if (loadFromFiles_)
        {
            const GeoField fld
            (
                IOobject
                (
                    fieldName,
                    time_.timeName(),
                    obr(),
                    IOobject::MUST_READ
                ),
                mesh_
            );

            sampleAndSet(fld, field);
        }
        else
        {
            sampleAndSet
            (
                obr().lookupObject<GeoField>(fieldName),
                field
            );
        }
    }
}


// ************************************************************************* //
