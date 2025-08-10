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
    (c) 2010-2012 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "volumeReport/volumeReport.H"
#include "sets/topoSetSource/topoSetSource.H"
#include "sets/topoSets/cellSet.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "dynamicFvMesh/dynamicFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(volumeReport, 0);
    addToRunTimeSelectionTable(functionObject, volumeReport, dictionary);
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::volumeReport::writeFileHeader() const
{
    fileWriter_.writeCommented("Time");

    forAll(fields_, fI)
    {
        // For backwards-compatability use the old header names
        string min_location_header = "at location";
        string max_location_header = "at location";
        if (!fileWriter_.isDatFile())
        {
            // But for other formats, use more descriptive (and non-duplicated!)
            // colum names
            min_location_header = fields_[fI] + "min location";
            max_location_header = fields_[fI] + "max location";
        }
        fileWriter_.writeDelimited
        (
            fields_[fI]+"min",
            min_location_header,
            fields_[fI]+"max",
            max_location_header,
            fields_[fI]+"mean",
            fields_[fI]+"stDev"
        );
    }

    if (isA<dynamicFvMesh>(mesh_) || mesh_.hasChangers())
    {
        fileWriter_.writeDelimited("totalVolume");
    }

    fileWriter_.endLine();
}

template<class GeoField>
void Foam::functionObjects::volumeReport::calcNonScalarData
(
    const scalarField& V,
    const volVectorField& C,
    const GeoField& f,
    label fI
)
{
    localData_[fI][0] = fieldInstance(GREAT);
    localData_[fI][1] = fieldInstance(-GREAT);
    globalData_[fI][0] = 0;
    globalData_[fI][1] = 0;

    forAll(reportCells_, rcI)
    {
        if (reportCells_.get(rcI) == 1)
        {
            const scalar magF = mag(f[rcI]);
            // look for minimum and record location
            if (magF < localData_[fI][0].value())
            {
                localData_[fI][0].setValue(magF);
                localData_[fI][0].setPosition(C[rcI]);
            }

            // look for maximum and record location
            if (magF > localData_[fI][1].value())
            {
                localData_[fI][1].setValue(magF);
                localData_[fI][1].setPosition(C[rcI]);
            }


            globalData_[fI][0] += V[rcI]*magF;
        }
    }
    reduce(
        std::tie(
            localData_[fI][0],
            localData_[fI][1],
            globalData_[fI][0]
        ),
        ParallelOp<minFIOp, maxFIOp, sumOp<scalar>>{}
    );
    globalData_[fI][0] /= totalVolume_;

    forAll(reportCells_, rcI)
    {
        if (reportCells_.get(rcI) == 1)
        {
            globalData_[fI][1] += V[rcI]*sqr(mag(f[rcI])-globalData_[fI][0]);
        }
    }

    reduce(globalData_[fI][1], sumOp<scalar>());
    globalData_[fI][1] /= totalVolume_;
    globalData_[fI][1] = sqrt(globalData_[fI][1]);
}

void Foam::functionObjects::volumeReport::calcData
(
    const scalarField& V,
    const volVectorField& C,
    const volScalarField& f,
    label fI
)
{
    localData_[fI][0] = fieldInstance(GREAT);
    localData_[fI][1] = fieldInstance(-GREAT);
    globalData_[fI][0] = 0;
    globalData_[fI][1] = 0;

    forAll(reportCells_, rcI)
    {
        if (reportCells_.get(rcI) == 1)
        {
            // look for minimum and record location
            if (f[rcI] < localData_[fI][0].value())
            {
                localData_[fI][0].setValue(f[rcI]);
                localData_[fI][0].setPosition(C[rcI]);
            }

            // look for maximum and record location
            if (f[rcI] > localData_[fI][1].value())
            {
                localData_[fI][1].setValue(f[rcI]);
                localData_[fI][1].setPosition(C[rcI]);
            }

            globalData_[fI][0] += V[rcI]*f[rcI];
        }
    }

    reduce(
        std::tie(
            localData_[fI][0],
            localData_[fI][1],
            globalData_[fI][0]
        ),
        ParallelOp<minFIOp, maxFIOp, sumOp<scalar>>{}
    );
    globalData_[fI][0] /= totalVolume_;

    forAll(reportCells_, rcI)
    {
        if (reportCells_.get(rcI) == 1)
        {
            globalData_[fI][1] += V[rcI]*sqr(f[rcI]-globalData_[fI][0]);
        }
    }

    reduce(globalData_[fI][1], sumOp<scalar>());
    globalData_[fI][1] /= totalVolume_;
    globalData_[fI][1] = sqrt(globalData_[fI][1]);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::volumeReport::volumeReport
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    fields_(0),
    dict_(dict),
    reportCells_(mesh_.nCells(), 0),
    totalVolume_(GREAT),
    globalData_(0),
    localData_(0),
    convergenceChecks_(),
    fileWriter_
    (
        mesh_.time(),
        fileName
        (
            outputFileDir() + "/" + name + "/" +
            Time::timeName
            (
                mesh_.time().timeToUserTime(mesh_.time().startTime().value())
            )
        ),
        typeName,
        dict
    )
{
    read(dict);
    writeFileHeader();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::volumeReport::~volumeReport()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::volumeReport::execute()
{
    Log << type() << " " << name() <<  " execute:" << nl;

    if (mesh_.changing())
    {
        updateReportCells();
        computeTotalVolume();
    }

    //calculate
    calculate();

    //Write to screen

    Ostream& rpInfo(Info);

    rpInfo.precision(5);

    forAll(fields_, fI)
    {
        rpInfo<< fields_[fI]
               << ": min " << localData_[fI][0].value()
               << " pos: " << localData_[fI][0].position()
               << ", max " << localData_[fI][1].value()
               << " pos: " << localData_[fI][1].position()
               << ", mean " << globalData_[fI][0]
               << ", stdDev " << globalData_[fI][1] << endl;
    }

    if (isA<dynamicFvMesh>(mesh_) || mesh_.hasChangers())
    {
        rpInfo<< "totalVolume: " << totalVolume_ << endl;
    }

    forAll(fields_, fI)
    {
        word fieldName = fields_[fI];
        forAll(convergenceChecks_, checkI)
        {
            convergenceTermination& check = convergenceChecks_[checkI];

            if (check.fieldName() == fieldName)
            {
                if (check.calculate(globalData_[fI][0]))
                {
                    Info<<"Convergence obtained for field : "<< fieldName
                        <<" in function object : "<< name() <<endl;
                }
            }
        }
    }

    rpInfo.precision(IOstream::defaultPrecision());

    fileWriter_.writeTime();

    forAll(fields_, fI)
    {
        scalar minVal = localData_[fI][0].value();
        vector minLoc = localData_[fI][0].position();
        scalar maxVal = localData_[fI][1].value();
        vector maxLoc = localData_[fI][1].position();
        scalar meanVal = globalData_[fI][0];
        scalar stdDev = globalData_[fI][1];

        if (fileWriter_.isDatFile())
        {
            // The legacy volReport format is a little weird - the body is
            // space-separated (which is unhelpful if it contains points...),
            // but the header is tab-separated.
            fileWriter_.file() << minVal << " " << minLoc << " " << maxVal << " "
               << maxLoc << " " << meanVal << " " << stdDev << " ";
        }
        else
        {
            fileWriter_.writeDelimited
            (
                minVal,
                minLoc,
                maxVal,
                maxLoc,
                meanVal,
                stdDev
            );
        }


        word nameStr('(' + fields_[fI] + ')');
        this->setResult("min" + nameStr, minVal);
        this->setResult("min" + nameStr + "_position", minLoc);
        this->setResult("max" + nameStr, maxVal);
        this->setResult("max" + nameStr + "_position", maxLoc);
        this->setResult("mean" + nameStr, meanVal);
        this->setResult("stdDev" + nameStr, stdDev);
    }

    if (isA<dynamicFvMesh>(mesh_) || mesh_.hasChangers())
    {
        if (fileWriter_.isDatFile())
        {
            // The legacy volReport format is a little weird - the body is
            // space-separated (which is unhelpful if it contains points...),
            // but the header is tab-separated.
            fileWriter_.file() << totalVolume_ << " ";
        }
        else
        {
            fileWriter_.writeDelimited(totalVolume_);
        }
        this->setResult("totalVolume", totalVolume_);
    }
    fileWriter_.file() << endl;
    Info<< endl;

    return true;
}


void Foam::functionObjects::volumeReport::calculate()
{
    const scalarField& V = mesh_.V();
    const volVectorField& C = mesh_.C();

    forAll(fields_, fI)
    {
        if (foundObject<volScalarField>(fields_[fI]))
        {
            calcData(V, C, lookupObject<volScalarField>(fields_[fI]), fI);
        }
        else if (foundObject<volVectorField>(fields_[fI]))
        {
            calcNonScalarData
                (V, C, lookupObject<volVectorField>(fields_[fI]), fI);
        }
        else if (foundObject<volTensorField>(fields_[fI]))
        {
            calcNonScalarData
                (V, C, lookupObject<volTensorField>(fields_[fI]), fI);
        }
        else if (foundObject<volSymmTensorField>(fields_[fI]))
        {
            calcNonScalarData
                (V, C, lookupObject<volSymmTensorField>(fields_[fI]), fI);
        }
        else if (foundObject<volSphericalTensorField>(fields_[fI]))
        {
            calcNonScalarData
                (V, C, lookupObject<volSphericalTensorField>(fields_[fI]), fI);
        }
        else
        {
            WarningInFunction
                    << "Field " << fields_[fI]
                    << " could not be found in the object registry."
                    << " No volume report compiled." << endl;
        }
    }
}

bool Foam::functionObjects::volumeReport::write()
{
    return true;
}


bool Foam::functionObjects::volumeReport::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    Info<< type() << " " << name() << ":" << nl;

    // This is one way to make the volumeReport format look like the legacy
    // format, but it doesn't work - in the legacy format, the header is
    // tab-separated but the body is space-separated.
    //
    // if (formats::dat == format_)
    //  {
    //     // The legacy volReport format is a little weird
    //     delimiter_ = token::SPACE;
    // }

    //field names
    {
        wordList tfns = dict.lookup("fields");
        fields_.setSize(tfns.size());
        fields_ = tfns;
        localData_.setSize(tfns.size());
        globalData_.setSize(tfns.size());

        forAll(localData_, dI)
        {
            //data_[dI].setSize(4);
            localData_[dI] = -1;
        }
        forAll(globalData_, dI)
        {
            //data_[dI].setSize(4);
            globalData_[dI] = -1;
        }
    }

    const dictionary* convSubDictPtr = dict.subDictPtr
    (
        "convergence"
    );

    if (convSubDictPtr)
    {
        label nConvChecks = 0;
        forAll(fields_, fieldI)
        {
            word fieldName = fields_[fieldI];
            if (convSubDictPtr->found(fieldName))
            {
                nConvChecks++;
            }
        }
        convergenceChecks_.setSize(nConvChecks);

        nConvChecks = 0;
        forAll(fields_, fieldI)
        {
            word fieldName = fields_[fieldI];
            if (convSubDictPtr->found(fieldName))
            {
                const dictionary convDict
                (
                    convSubDictPtr->subDict(fieldName)
                );

                convergenceChecks_.set
                (
                    nConvChecks,
                    new convergenceTermination
                    (
                        obr_,
                        convDict
                    )
                 );
                nConvChecks++;
            }
        }
    }

    updateReportCells();
    computeTotalVolume();

    return true;
}


void Foam::functionObjects::volumeReport::computeTotalVolume()
{
    totalVolume_= 0;
    forAll(mesh_.V(), cI)
    {
        if (reportCells_.get(cI) == 1)
        {
            totalVolume_ += mesh_.V()[cI];
        }
    }
    reduce(totalVolume_, sumOp<scalar>());
}


void Foam::functionObjects::volumeReport::updateReportCells()
{
    //set list size and initialise with 0
    reportCells_ = PackedList<1>(mesh_.nCells(), 0);

    //select cells
    //if sets entry is not present, select all cells
    if (dict_.found("sets"))
    {
        PtrList<entry> regions(dict_.lookup("sets"));

        forAll(regions, regioni)
        {
            const entry& region = regions[regioni];

            autoPtr<topoSetSource> cellSelector =
                topoSetSource::New(region.keyword(), mesh_, region.dict());

            cellSet selectedCellSet
            (
                mesh_,
                "cellSet",
                mesh_.nCells()/10+1  // Reasonable size estimate.
            );

            cellSelector->applyToSet
            (
                topoSetSource::NEW,
                selectedCellSet
            );

            forAllConstIter(labelHashSet, selectedCellSet, iter)
            {
                reportCells_.set(iter.key(), 1);
            }
        }
    }
    else
    {
        forAll(mesh_.C(), cI)
        {
            reportCells_.set(cI, 1);
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
