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
    (c) 2020 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "pedestrianComfort/pedestrianComfort.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "db/Time/Time.H"
#include "db/solutionInstanceRegistry/solutionInstanceRegistry.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace functionObjects
    {
        defineTypeNameAndDebug(pedestrianComfort, 0);
        addToRunTimeSelectionTable(functionObject, pedestrianComfort, dictionary);
    }
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::pedestrianComfort::calculateTotalProbabilityExceedance()
{
    // loop over thresholds given by the criteria
    forAll(criteriaThresholds_, thresholdI)
    {
        scalar Ulim = criteriaThresholds_[thresholdI];

        volScalarField& tpe = *totProbExceedPtr_[thresholdI];

        if (windDirections_.size())
        {
            forAll(windDirections_, wdirI)
            {
                word instanceRegionName = instanceNames_[wdirI];
                word windDirRegion = windDirections_[wdirI];

                const solutionInstanceRegistry& instanceRegistry =
                    time_.lookupObject<solutionInstanceRegistry>
                    (
                        instanceRegionName
                    );
                if (instanceRegistry.isActive())
                {
                    const objectRegistry& regionRegistry = mesh_.lookupObject<objectRegistry>
                    (
                        windDirRegion
                    );
                    const volVectorField& U =
                        regionRegistry.lookupObject<volVectorField>(UName_);

                    // get scale factor for provided reference height
                    scalar sf = scaleFactor_->value(hRef_[wdirI]);
                    // get Weibull distribution parameters
                    scalar pw = WeibullParameters_[wdirI](0);
                    scalar cw = WeibullParameters_[wdirI](1);
                    scalar kw = WeibullParameters_[wdirI](2);

                    volScalarField Uscaled
                    (
                        IOobject
                        (
                            "Uscaled",
                            time_.timeName(),
                            obr_,
                            IOobject::NO_READ,
                            IOobject::NO_WRITE,
                            false
                        ),
                        mesh_,
                        dimensionedScalar("totProbExceed", dimVelocity, 0.0)
                    );
                    if (gustCorr_)
                    {
                        const volScalarField& tke =
                            regionRegistry.lookupObject<volScalarField>(kName_);

                        if (meanVelCorr_)
                        {
                            volScalarField UMeanCorr2(magSqr(U) + 2*tke);
                            volScalarField Ugust(sqrt(UMeanCorr2)+gustG_*sqrt(2/3*tke));
                            Uscaled = max(sqrt(UMeanCorr2), Ugust/gustT_)/URef_[wdirI]*sf*cw;
                        }
                        else
                        {
                            volScalarField Ugust(mag(U)+gustG_*sqrt(2/3*tke));
                            Uscaled = max(mag(U), Ugust/gustT_)/URef_[wdirI]*sf*cw;
                        }

                    }
                    else
                    {
                        if (meanVelCorr_)
                        {
                            const volScalarField& tke =
                                regionRegistry.lookupObject<volScalarField>("k");
                            volScalarField UMeanCorr2(magSqr(U) + 2*tke);
                            Uscaled = sqrt(UMeanCorr2)/URef_[wdirI]*sf*cw;
                        }
                        else
                        {
                            Uscaled = mag(U)/URef_[wdirI]*sf*cw;
                        }
                    }


                    forAll(tpe, cellI)
                    {
                        tpe[cellI] +=
                            pw*Foam::exp(-Foam::pow(Ulim/(Uscaled[cellI]+VSMALL), kw));
                    }
                }
            }
        }

        totProbExceedPtr_[thresholdI]->write();
    }

    assessComfort();
}

void Foam::functionObjects::pedestrianComfort::assessComfort()
{
    volScalarField& PC = pComfort_();

    forAll(PC, cellI)
    {
        // marker corresponding to fist category
        // i.e. lowest velocity
        label categoryMarker = 1;
        forAll(criteriaThresholds_, thresholdI)
        {
            const volScalarField& tpe = *totProbExceedPtr_[thresholdI];
            // exceedance given by criteria
            scalar crEx = criteriaExceedance_[thresholdI];
            // local total exceedance
            scalar prEx = tpe[cellI];

            if (prEx <= crEx)
            {
                PC[cellI] = categoryMarker;
                break; // assign lowest category marker
            }
            categoryMarker++;
        }
    }

    PC.write();
}

// * * * * * * * * * * * * *      Constructors         * * * * * * * * * * * //

Foam::functionObjects::pedestrianComfort::pedestrianComfort
(
    const Foam::word &name,
    const Foam::Time &runTime,
    const Foam::dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    UName_(),
    kName_("k"),
    meanVelCorr_(),
    gustCorr_(),
    gustG_(),
    gustT_(),
    windDirections_(),
    WeibullParameters_(),
    hRef_(),
    URef_(),
    scaleFactor_(),
    criteriaThresholds_(),
    criteriaExceedance_(),
    totProbExceedPtr_(0),
    pComfort_(nullptr)
{
    read(dict);
    totProbExceedPtr_.setSize(criteriaThresholds_.size());
    forAll(totProbExceedPtr_, thresholdI)
    {
        word fieldName = "tPEx_"+Foam::word(this->name())+"_"+Foam::word(Foam::name(criteriaThresholds_[thresholdI]));
        if (mesh_.foundObject<volScalarField>(fieldName))
        {
            totProbExceedPtr_[thresholdI] = mesh_.lookupObjectRefPtr<volScalarField>(fieldName);
        }
        else
        {
            totProbExceedPtr_[thresholdI] =
                new volScalarField
                (
                    IOobject
                    (
                        fieldName,
                        time_.timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh_,
                    dimensionedScalar
                    (
                        "totProbExceed",
                        dimless,
                        0.0
                    )
                );
            totProbExceedPtr_[thresholdI]->store();
        }
    }

    // init with marker from comfort criteria last category
    // i.e. total excedance > last threshold
    scalar  nThresh = totProbExceedPtr_.size() + 1.0;

    pComfort_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "PC-"+Foam::word(this->name()),
                time_.timeName(),
                obr_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("PC", dimless, nThresh)
         )
    );
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::pedestrianComfort::~pedestrianComfort()
= default;

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
bool Foam::functionObjects::pedestrianComfort::execute()
{
    return true;
}

bool Foam::functionObjects::pedestrianComfort::write()
{
    calculateTotalProbabilityExceedance();
    return true;
}

bool Foam::functionObjects::pedestrianComfort::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    Info<< type() << " " << name() << ":" << nl;

    dict.lookup("UName") >> UName_;
    meanVelCorr_ = dict.lookupOrDefault<bool>("meanVelCorr", false);
    gustCorr_ = dict.lookupOrDefault<bool>("gust", false);
    if (gustCorr_)
    {
        kName_ = dict.lookupOrDefault<word>("kName", "k");
    }
    gustG_ = dict.lookupOrDefault<scalar>("peakFactor", 3.5);
    gustT_ = dict.lookupOrDefault<scalar>("meanFactor", 1.85);

    dict.lookup("windDirections") >> windDirections_;
    numberOfWindDirections_ = windDirections_.size();

    dict.lookup("WeibullParameters") >> WeibullParameters_;
    dict.lookup("instances") >> instanceNames_;

    if (instanceNames_.size() != numberOfWindDirections_)
    {
        FatalErrorInFunction
            << "instanceNames.size() "<<instanceNames_.size()<<
            "  not equal to windDirections_.size() "<<numberOfWindDirections_
            <<exit(FatalError);
    }
    forAll(instanceNames_, instI)
    {
        word instanceRegionName = instanceNames_[instI];
        const solutionInstanceRegistry& instanceRegistry =
            time_.lookupObject<solutionInstanceRegistry>
            (
                 instanceRegionName
            );
        if (!instanceRegistry.isSelected())
        {
            WarningInFunction<<"Instances are defined but are not selected "<<endl;
            break;
        }
    }

    if (WeibullParameters_.size() != numberOfWindDirections_)
    {
        FatalErrorInFunction
            << "WeibullParameters_.size() "<<WeibullParameters_.size()<<
            "  not equal to windDirections_.size() "<<numberOfWindDirections_
            <<exit(FatalError);
    }

    dict.lookup("referenceHeight") >> hRef_;

    if (hRef_.size() != numberOfWindDirections_)
    {
        FatalErrorInFunction
            << "hRef_.size() "<<hRef_.size()<<
            "  not equal to windDirections_.size() "<<numberOfWindDirections_
            <<exit(FatalError);
    }

    dict.lookup("referenceVelocity") >> URef_;

    if (URef_.size() != numberOfWindDirections_)
    {
        FatalErrorInFunction
            << "URef_.size() "<<URef_.size()<<
            "  not equal to windDirections_.size() "<<numberOfWindDirections_
            <<exit(FatalError);
    }

    scaleFactor_ = Function1<scalar>::New("scaleFactor", dict);

    dict.lookup("criteriaExceedance") >> criteriaExceedance_;
    dict.lookup("criteriaThresholds") >> criteriaThresholds_;

    if (criteriaExceedance_.size() != criteriaThresholds_.size())
    {
        FatalErrorInFunction
            << "criteriaExceedance_.size() "<<criteriaExceedance_.size()<<
            "  not equal to criteriaThresholds_.size() "<<criteriaThresholds_.size()
            <<exit(FatalError);
    }

    return true;
}
