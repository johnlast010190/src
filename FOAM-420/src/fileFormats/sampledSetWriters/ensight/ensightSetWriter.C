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

\*---------------------------------------------------------------------------*/

#include "sampledSetWriters/ensight/ensightSetWriter.H"
#include "coordSet/coordSet.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "db/IOstreams/IOstreams/IOmanip.H"
#include "global/foamVersion.H"
#include "ensight/type/ensightPTraits.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::ensightSetWriter<Type>::ensightSetWriter()
:
    writer<Type>()
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::ensightSetWriter<Type>::~ensightSetWriter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const char* Foam::ensightSetWriter<Type>::ext() const {
    return ".case";
}

template<class Type>
void Foam::ensightSetWriter<Type>::write
(
    const coordSet& points,
    const wordList& valueSetNames,
    const List<const Field<Type>*>& valueSets,
    Ostream& os
) const
{
    const fileName base(os.name().lessExt());
    const fileName meshFile(base + ".mesh");

    // Write .case file
    os  << "FORMAT" << nl
        << "type: ensight gold" << nl
        << nl
        << "GEOMETRY" << nl
        << "model:        1     " << meshFile.name().c_str() << nl
        << nl
        << "VARIABLE"
        << nl;
    forAll(valueSetNames, setI)
    {
        fileName dataFile(base + ".***." + valueSetNames[setI]);

        os.setf(ios_base::left);
        os  << ensightPTraits<Type>::typeName
            << " per node:            1       "
            << setw(15) << valueSetNames[setI]
            << " " << dataFile.name().c_str()
            << nl;
    }
    os  << nl
        << "TIME" << nl
        << "time set:                      1" << nl
        << "number of steps:               1" << nl
        << "filename start number:         0" << nl
        << "filename increment:            1" << nl
        << "time values:" << nl
        << "0.00000e+00" << nl;

    // Write .mesh file
    {
        string desc = string("written by OpenFOAM-") + Foam::FOAMversion;
        OFstream os(meshFile);
        os.setf(ios_base::scientific, ios_base::floatfield);
        os.precision(5);

        os  << "EnSight Geometry File" << nl
            << desc.c_str() << nl
            << "node id assign" << nl
            << "element id assign" << nl
            << "part" << nl
            << setw(10) << 1 << nl
            << "internalMesh" << nl
            << "coordinates" << nl
            << setw(10) << points.size() << nl;

        for (direction cmpt = 0; cmpt < vector::nComponents; cmpt++)
        {
            forAll(points, pointi)
            {
                const scalar comp = points[pointi][cmpt];
                if (mag(comp) >= scalar(floatScalarVSMALL))
                {
                    os  << setw(12) << comp << nl;
                }
                else
                {
                    os  << setw(12) << scalar(0) << nl;
                }
            }
        }
        os  << "point" << nl
            << setw(10) << points.size() << nl;
        forAll(points, pointi)
        {
            os  << setw(10) << pointi+1 << nl;
        }
    }

    // Write data files
    forAll(valueSetNames, setI)
    {
        fileName dataFile(base + ".000." + valueSetNames[setI]);
        OFstream os(dataFile);
        os.setf(ios_base::scientific, ios_base::floatfield);
        os.precision(5);
        {
            os  << ensightPTraits<Type>::typeName << nl
                << "part" << nl
                << setw(10) << 1 << nl
                << "coordinates" << nl;

            for (direction i=0; i < pTraits<Type>::nComponents; ++i)
            {
                label cmpt = ensightPTraits<Type>::componentOrder[i];

                const scalarField fld(valueSets[setI]->component(cmpt));
                forAll(fld, i)
                {
                    if (mag(fld[i]) >= scalar(floatScalarVSMALL))
                    {
                        os  << setw(12) << fld[i] << nl;
                    }
                    else
                    {
                        os  << setw(12) << scalar(0) << nl;
                    }
                }
            }
        }
    }
}


template<class Type>
void Foam::ensightSetWriter<Type>::write
(
    const bool writeTracks,
    const PtrList<coordSet>& tracks,
    const wordList& valueSetNames,
    const List<List<Field<Type>>>& valueSets,
    Ostream& os
) const
{
    const fileName base(os.name().lessExt());
    const fileName meshFile(base + ".mesh");

    // Write .case file
    os  << "FORMAT" << nl
        << "type: ensight gold" << nl
        << nl
        << "GEOMETRY" << nl
        << "model:        1     " << meshFile.name().c_str() << nl
        << nl
        << "VARIABLE"
        << nl;
    forAll(valueSetNames, setI)
    {
        fileName dataFile(base + ".***." + valueSetNames[setI]);

        os.setf(ios_base::left);
        os  << ensightPTraits<Type>::typeName
            << " per node:            1       "
            << setw(15) << valueSetNames[setI]
            << " " << dataFile.name().c_str()
            << nl;
    }
    os  << nl
        << "TIME" << nl
        << "time set:                      1" << nl
        << "number of steps:               1" << nl
        << "filename start number:         0" << nl
        << "filename increment:            1" << nl
        << "time values:" << nl
        << "0.00000e+00" << nl;

    // Write .mesh file
    {
        string desc = string("written by OpenFOAM-") + Foam::FOAMversion;
        OFstream os(meshFile);
        os.setf(ios_base::scientific, ios_base::floatfield);
        os.precision(5);
        os  << "EnSight Geometry File" << nl
            << desc.c_str() << nl
            << "node id assign" << nl
            << "element id assign" << nl;

        forAll(tracks, trackI)
        {
            const coordSet& points = tracks[trackI];

            os  << "part" << nl
                << setw(10) << trackI+1 << nl
                << "internalMesh" << nl
                << "coordinates" << nl
                << setw(10) << points.size() << nl;

            for (direction cmpt = 0; cmpt < vector::nComponents; cmpt++)
            {
                forAll(points, pointi)
                {
                    const scalar comp = points[pointi][cmpt];
                    if (mag(comp) >= scalar(floatScalarVSMALL))
                    {
                        os  << setw(12) << comp << nl;
                    }
                    else
                    {
                        os  << setw(12) << scalar(0) << nl;
                    }
                }
            }

            if (writeTracks)
            {
                os  << "bar2" << nl
                    << setw(10) << points.size()-1 << nl;
                for (label i = 0; i < points.size()-1; i++)
                {
                    os  << setw(10) << i+1
                        << setw(10) << i+2
                        << nl;
                }
            }
        }
    }


    // Write data files
    forAll(valueSetNames, setI)
    {
        fileName dataFile(base + ".000." + valueSetNames[setI]);
        OFstream os(dataFile);
        os.setf(ios_base::scientific, ios_base::floatfield);
        os.precision(5);
        {
            os  << ensightPTraits<Type>::typeName << nl;

            const List<Field<Type>>& fieldVals = valueSets[setI];
            forAll(fieldVals, trackI)
            {
                os  << "part" << nl
                    << setw(10) << trackI+1 << nl
                    << "coordinates" << nl;

                for (direction i=0; i < pTraits<Type>::nComponents; ++i)
                {
                    label cmpt = ensightPTraits<Type>::componentOrder[i];

                    const scalarField fld(fieldVals[trackI].component(cmpt));
                    forAll(fld, i)
                    {
                        if (mag(fld[i]) >= scalar(floatScalarVSMALL))
                        {
                            os  << setw(12) << fld[i] << nl;
                        }
                        else
                        {
                            os  << setw(12) << scalar(0) << nl;
                        }
                    }
                }
            }
        }
    }
}


// ************************************************************************* //
