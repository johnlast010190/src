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
     (c) 2020-2020 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "motionSmoother/polyMeshGeometry/meshStatistics.H"
#include "fields/fvPatchFields/basic/zeroGradient/zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(meshStatistics, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::meshStatistics::meshStatistics
(
    const fvMesh& mesh,
    const word& name,
    const dictionary& dict
)
:
    name_(name),
    dict_(dict),
    quality_
    (
        IOobject
        (
            name_,
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    frequency_(0),
    nSamples_(0),
    minRange_(-GREAT),
    binSize_(0),
    minVal_(-GREAT),
    maxVal_(GREAT),
    aveVal_(0),
    stdDev_(0),
    active_(false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::meshStatistics::calcStats
(
    const List<scalar>& metrics
)
{
    nSamples_ = metrics.size();
    reduce(nSamples_, sumOp<label>());

    active_ = dict_.lookupOrDefault<bool>("active", true);

    if (nSamples_ == 0 || !active_)
    {
        active_ = false;
        return;
    }

    minRange_ = dict_.lookupOrDefault<scalar>("minRange", GREAT);
    scalar maxRange(dict_.lookupOrDefault<scalar>("maxRange", -GREAT));
    label nIntervals(dict_.lookupOrDefault<label>("nIntervals", 20));

    if (!(dict_.found("minRange") && dict_.found("maxRange")))
    {
        forAll(metrics, i)
        {
            scalar mv = metrics[i];
            minRange_ = min(minRange_, mv);
            maxRange = max(maxRange, mv);
        }

        reduce(
            std::tie(minRange_, maxRange),
            ParallelOp<minOp<scalar>, maxOp<scalar>>{}
        );
    }

    if ((maxRange - minRange_) < 1e-8 || nIntervals < 0)
    {
        nIntervals = 1;
    }

    frequency_.setSize(nIntervals, 0);
    binSize_ = (maxRange-minRange_)/ nIntervals;

    maxVal_ = -GREAT;
    minVal_ = GREAT;

    aveVal_ =  scalar(0);
    stdDev_ =  scalar(0);

    forAll(metrics, i)
    {
        scalar mv = metrics[i];

        maxVal_ = max(mv,maxVal_);
        minVal_ = min(mv,minVal_);
        aveVal_ += mv;

        label bini = 0;

        if (binSize_ > SMALL)
        {
            mv -= minRange_;
            mv /= binSize_;

            if (mv < 0)
            {
                bini = 0;
            }
            else if (mv >= scalar(nIntervals))
            {
                bini = nIntervals-1;
            }
            else
            {
                bini = label(mv);
            }
        }
        frequency_[bini]++;
    }

    reduce(
        std::tie(aveVal_, maxVal_, minVal_),
        ParallelOp<sumOp<scalar>, maxOp<scalar>, minOp<scalar>>{}
    );
    aveVal_ /= nSamples_;
    Pstream::listCombineGather(frequency_, plusEqOp<label>());

    forAll(metrics, i)
    {
        scalar mv = metrics[i];
        stdDev_ += sqr(mv - aveVal_);
    }
    reduce(stdDev_, sumOp<scalar>());

    stdDev_ /= nSamples_;
    stdDev_ =sqrt(stdDev_);
}

void Foam::meshStatistics::writeStats()
{
    if (active_ && Pstream::master())
    {
        Info<< endl;

        Info<<"Calculating " << name_ <<" mesh statistics : "<<endl;
        Info<< " Min : " << minVal_ << " Max : " << maxVal_
             << " Ave : " << aveVal_ << " Std Dev : " << stdDev_
             << nl << endl;

        Info<< "Interval, frequency, percentage" <<endl;
        Info<< "-------------------------------" <<endl;

        scalar prevVal = minRange_;
        forAll(frequency_, i)
        {
            scalar percent = 100.*frequency_[i]/nSamples_;
            scalar nextVal = minRange_ + (i+1)*binSize_;

            Info<<prevVal<<" - "<<nextVal
                <<" ( "<<frequency_[i]<<" / "<<nSamples_<<" ) "
                << percent << " % " << endl;
            prevVal = nextVal;
        }
        Info<< endl;
    }
}



// ************************************************************************* //
