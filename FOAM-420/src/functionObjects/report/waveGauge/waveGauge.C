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
    (c) 2011-2013 OpenFOAM Foundation
    (c) 2016, 2020 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "waveGauge/waveGauge.H"
#include "db/dictionary/dictionary.H"
#include "db/Time/Time.H"
#include "fields/volFields/volFields.H"
#include "containers/Lists/ListListOps/ListListOps.H"
#include "containers/Lists/SortableList/SortableList.H"
#include "interpolation/volPointInterpolation/volPointInterpolation.H"
#include "db/dictionary/dictionaryEntry/dictionaryEntry.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(waveGauge, 0);
    addToRunTimeSelectionTable(functionObject, waveGauge, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::waveGauge::waveGauge
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    sampledSets(name, runTime, dict),
    writeFile(mesh(), name, typeName, dict),
    vofField_("alpha1"),
    contourLevel_(0.5)
{
    read(dict);
    writeFileHeader(file());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::waveGauge::~waveGauge()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::functionObjects::waveGauge::writeFileHeader(Ostream& os) const
{
    //make sure file directories exist and create OFstream
    writeCommented(os,"Time");

    forAll(gaugeNames_, i)
    {
        writeDelimited(os,gaugeNames_[i]);
    }

    os << endl;
}


bool Foam::functionObjects::waveGauge::read(const dictionary& dict)
{
    writeFile::read(dict);

    Log << type() << " " << name() <<  " read:" << nl;

    // Assemble dictionary for base class sample object
    dictionary sampleDict;
    dict.lookup("vofField") >> vofField_;
    wordReList fields(1);
    fields[0] = vofField_;
    sampleDict.add("fields", fields);
    sampleDict.add("interpolationScheme", "cell");
    sampleDict.add("setFormat", "raw");

    PtrList<entry> gaugesList(dict.lookup("gauges"));
    gaugeNames_.resize(gaugesList.size());
    PtrList<dictionaryEntry> setsDicts(gaugesList.size());

    Log << gaugeNames_ << endl;

    forAll(gaugesList, i)
    {
        dictionary setsDict;
        setsDict.add("type", "midPointAndFace"); //To get the boundary faces
        setsDict.add("axis", "distance");
        setsDict.add("start", gaugesList[i].dict().lookup("start"));
        setsDict.add("end", gaugesList[i].dict().lookup("end"));
        setsDicts.set
        (
            i,
            new dictionaryEntry(gaugesList[i].keyword(), sampleDict, setsDict)
        );
        gaugeNames_[i] = gaugesList[i].keyword();
    }
    sampleDict.add("sets", setsDicts);
    dict.lookup("contourLevel") >> contourLevel_;

    // Call base class
    sampledSets::read(sampleDict);

    return true;
}


bool Foam::functionObjects::waveGauge::execute()
{
    Log << type() << " " << name() <<  " execute:" << nl;

    // Storage for interpolated values
    PtrList<volFieldSampler<scalar>> sampledFields(1);

    sampledFields.set
    (
        0,
        new volFieldSampler<scalar>
        (
            lookupObject<volScalarField>(vofField_),
            *this
        )
    );

    // Combine sampled fields from processors.
    // Note: only master results are valid

    PtrList<volFieldSampler<scalar>> masterFields(1);
    combineSampledValues(sampledFields, indexSets(), masterFields);

    scalarList heights(masterFields[0].size());

    if (Pstream::master())
    {
        forAll(masterFields[0], gaugeI)
        {
            // Find highest intersection
            const List<scalar>& h = masterSampledSets()[gaugeI].curveDist();

            const scalarField& alpha(masterFields[0][gaugeI]); //already sorted
            heights[gaugeI] = h[0];

            for (label i = 1; i < alpha.size(); i++)
            {
                if
                (
                    (alpha[i-1] >= contourLevel_ && alpha[i] < contourLevel_)
                    || (alpha[i-1] < contourLevel_ && alpha[i] >= contourLevel_)
                )
                {
                    scalar r = (alpha[i]-contourLevel_)
                        /stabilise(alpha[i]-alpha[i-1], SMALL);
                    heights[gaugeI] = r*h[i-1]+(1-r)*h[i];
                }
            }
        }
    }

    //Write to screen
    Log << "Wave height at:"<<nl;
    forAll(heights, i)
    {
        Log << tab << gaugeNames_[i] << ": " << heights[i] << endl;
    }
    Log << endl;

    file() << obr_.time().timeName() << " ";

    forAll(heights, i)
    {
        file() << heights[i] << " ";
    }
    file() << endl;

    return true;
}

bool Foam::functionObjects::waveGauge::write()
{
    return true;
}



// ************************************************************************* //
