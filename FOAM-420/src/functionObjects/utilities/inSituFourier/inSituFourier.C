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
    (c) 2015 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "inSituFourier/inSituFourier.H"
#include "fields/Fields/vector2DField/vector2DField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(inSituFourier, 0);
    addToRunTimeSelectionTable(functionObject, inSituFourier, dictionary);
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::inSituFourier::createFields()
{
    const Time& runTime = mesh_.time();

    const volScalarField& Phi = lookupObject<volScalarField>(fieldName_);

    //storage for each band
    label numBands = fm_.size();
    avgPower_.resize(numBands);
    integralPhiSinWinDt_.resize(numBands);
    integralPhiCosWinDt_.resize(numBands);
    integralSinWinDt_.resize(numBands);
    integralCosWinDt_.resize(numBands);
    dB_.resize(numBands);
    avgPhi_.resize(numBands);

    forAll(avgPower_,m)
    {
        OStringStream bandNumStr;
        bandNumStr << m;

        if (!avgPower_.set(m))
        {
            avgPower_.set
            (
                m,
                new volScalarField
                (
                    IOobject
                    (
                        name()+"_avgPower_"+bandNumStr.str(),
                        runTime.timeName(),
                        mesh_,
                        IOobject::READ_IF_PRESENT,
                        IOobject::NO_WRITE
                    ),
                    mesh_,
                    dimensionedScalar
                    (
                        "0",
                        sqr(Phi.dimensions()*PhiRescaling_->dimensions()),
                        0.0
                    )
                )
            );
            deactivateFields(avgPower_[m]);
        }

        if (!dB_.set(m))
        {
            dB_.set
            (
                m,
                new volScalarField
                (
                    IOobject
                    (
                        name()+"_dB_"+bandNumStr.str(),
                        runTime.timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh_,
                    dimensionedScalar("0", dimless, 0.0)
                )
            );
            deactivateFields(dB_[m]);
        }

        if (!integralPhiSinWinDt_.set(m))
        {
            integralPhiSinWinDt_.set
            (
                m,
                new volScalarField
                (
                    IOobject
                    (
                        "integralPhiSinWinDt_"+bandNumStr.str(),
                        runTime.timeName(),
                        name()+"Restart",
                        mesh_,
                        IOobject::READ_IF_PRESENT,
                        IOobject::NO_WRITE
                    ),
                    mesh_,
                    dimensionedScalar("0", Phi.dimensions()*dimTime, 0.0)
                )
            );
            deactivateFields(integralPhiSinWinDt_[m]);
        }

        if (!integralPhiCosWinDt_.set(m))
        {
            integralPhiCosWinDt_.set
            (
                m,
                new volScalarField
                (
                    IOobject
                    (
                        "integralPhiCosWinDt_"+bandNumStr.str(),
                        runTime.timeName(),
                        name()+"Restart",
                        mesh_,
                        IOobject::READ_IF_PRESENT,
                        IOobject::NO_WRITE
                    ),
                    mesh_,
                    dimensionedScalar("0", Phi.dimensions()*dimTime, 0.0)
                )
            );
            deactivateFields(integralPhiCosWinDt_[m]);
        }

        if (!integralSinWinDt_.set(m))
        {
            integralSinWinDt_.set
            (
                m,
                new uniformDimensionedScalarField
                (
                    IOobject
                    (
                        "integralSinWinDt_"+bandNumStr.str(),
                        runTime.timeName(),
                        name()+"Restart",
                        mesh_,
                        IOobject::READ_IF_PRESENT,
                        IOobject::NO_WRITE
                    ),
                    dimensionedScalar("0", dimTime, 0.0)
                )
            );
        }
        if (!integralCosWinDt_.set(m))
        {
            integralCosWinDt_.set
            (
                m,
                new uniformDimensionedScalarField
                (
                    IOobject
                    (
                        "integralCosWinDt_"+bandNumStr.str(),
                        runTime.timeName(),
                        name()+"Restart",
                        mesh_,
                        IOobject::READ_IF_PRESENT,
                        IOobject::NO_WRITE
                    ),
                    dimensionedScalar("0", dimTime, 0.0)
                )
            );
        }

        if (!avgPhi_.set(m))
        {
            avgPhi_.set
            (
                m,
                new volScalarField
                (
                    IOobject
                    (
                        "avgPhi_"+bandNumStr.str(),
                        runTime.timeName(),
                        name()+"Restart",
                        mesh_,
                        IOobject::READ_IF_PRESENT,
                        IOobject::NO_WRITE
                    ),
                    mesh_,
                    dimensionedScalar("0", Phi.dimensions(), 0.0)
                )
            );
            deactivateFields(avgPhi_[m]);
        }

    }

}

void Foam::functionObjects::inSituFourier::deactivateFields(volScalarField& f)
{
    // Resize the unused internal field and/or patches to
    // 1 to get a uniform write.
    // A bit of a hack, but makes I/O for restarts much less clunky
    if (!doVolume_)
    {
        f.primitiveFieldRef().resize(1);
    }
    forAll(f.boundaryField(),patchi)
    {
        if (!patchSet_.found(patchi))
        {
            f.boundaryFieldRef()[patchi].resize(1);
        }
    }
}


void Foam::functionObjects::inSituFourier::doIntegration
(
    const label bandNo,
    const scalarField& Phi,
    scalarField& integralPhiSinWinDt,
    scalarField& integralPhiCosWinDt,
    scalar& integralSinWinDt,
    scalar& integralCosWinDt,
    scalarField& avgPhi,
    scalarField& avgPower,
    scalarField& dB
)
{
    const label m = bandNo;
    scalar t = obr_.time().value();
    scalar dt = obr_.time().deltaTValue();


    scalar T = 1.0/(fhigh_[m]-flow_[m]);
    int n = int((t-t0_)/T);
    int nold = int((t-dt-t0_)/T);
    scalar tseg = t-t0_-n*T;

    if (n != nold) //end of window
    {
        scalarField ca( (integralPhiSinWinDt-avgPhi*integralSinWinDt)*PhiRescaling_->value() );
        scalarField cb( (integralPhiCosWinDt-avgPhi*integralCosWinDt)*PhiRescaling_->value() );
        scalarField power( 2*(sqr(ca)+sqr(cb))/sqr(T) );
        avgPower = (power+n*avgPower)/(n+1);
        avgPower = max(avgPower,SMALL);
        dB = 10*log10(avgPower/sqr(dBRefVal_().value()));

        integralSinWinDt = 0;
        integralCosWinDt = 0;
        integralPhiSinWinDt = 0;
        integralPhiCosWinDt = 0;
        avgPhi = Phi;
    }

    //integration for this timestep

    const scalar pi = constant::mathematical::pi;

    //window function (Hanning)
    scalar win = ecf_*0.5*(1-cos(2*pi*tseg/T));

    scalar sinWinDt = sin(-2*pi*fm_[m]*tseg)*win*dt;
    scalar cosWinDt = cos(-2*pi*fm_[m]*tseg)*win*dt;

    integralSinWinDt += sinWinDt;
    integralCosWinDt += cosWinDt;

    integralPhiSinWinDt += Phi*sinWinDt;
    integralPhiCosWinDt += Phi*cosWinDt;

    avgPhi = (tseg*avgPhi + dt*Phi)/(tseg+dt);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::inSituFourier::inSituFourier
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict)
{
    read(dict);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::inSituFourier::~inSituFourier()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::inSituFourier::execute()
{
    Log << type() << " " << name() <<  " execute:" << nl;

    const volScalarField& Phi = lookupObject<volScalarField>(fieldName_);

    if (!PhiRescaling_.valid())
    {
        //autodetect an incompressible 'pressure' that needs rescaling
        PhiRescaling_.reset(new dimensionedScalar("PhiRescaling", dimless, 1.0));
        if
        (
            Phi.dimensions() == dimPressure/dimDensity
        )
        {
            IOdictionary transportProperties
            (
                IOobject
                (
                    "transportProperties" ,
                    mesh_.time().constant( ) ,
                    mesh_,
                    IOobject::MUST_READ ,
                    IOobject::NO_WRITE
                )
            );
            PhiRescaling_.reset
            (
                new dimensionedScalar(transportProperties.lookup("rho"))
            );
        }
    }

    //create fields if not already done
    createFields();

    //check dimensions
    if
    (
        dimensionSet::debug &&
        Phi.dimensions()*PhiRescaling_->dimensions()
        != dBRefVal_->dimensions()
    )
    {
        FatalErrorInFunction
            << "Dimensions of " << Phi.name() << " and "
            << dBRefVal_->name() << nl
            << "are not consistent."
            << Foam::abort(FatalError);
    }

    for (label m = 0; m < fm_.size(); m++)
    {
        scalar integralSinWinDt = integralSinWinDt_[m].value();
        scalar integralCosWinDt = integralCosWinDt_[m].value();
        if (doVolume_)
        {
            doIntegration
            (
                m,
                Phi.primitiveField(),
                integralPhiSinWinDt_[m].primitiveFieldRef(),
                integralPhiCosWinDt_[m].primitiveFieldRef(),
                integralSinWinDt,
                integralCosWinDt,
                avgPhi_[m].primitiveFieldRef(),
                avgPower_[m].primitiveFieldRef(),
                dB_[m].primitiveFieldRef()
             );
        }

        forAllConstIter(labelHashSet, patchSet_, iter)
        {
            label patchi = iter.key();
            integralSinWinDt = integralSinWinDt_[m].value();
            integralCosWinDt = integralCosWinDt_[m].value();
            doIntegration
            (
                m,
                Phi.boundaryField()[patchi],
                integralPhiSinWinDt_[m].boundaryFieldRef()[patchi],
                integralPhiCosWinDt_[m].boundaryFieldRef()[patchi],
                integralSinWinDt,
                integralCosWinDt,
                avgPhi_[m].boundaryFieldRef()[patchi],
                avgPower_[m].boundaryFieldRef()[patchi],
                dB_[m].boundaryFieldRef()[patchi]
            );
        }
        integralSinWinDt_[m].value() = integralSinWinDt;
        integralCosWinDt_[m].value() = integralCosWinDt;
    }

    return true;
}


bool Foam::functionObjects::inSituFourier::write()
{
    Log << type() << " " << name() <<  " write:" << nl;

    forAll(avgPower_,m)
    {
        if (avgPower_.set(m))
        {
            avgPower_[m].write();
        }

        if (dB_.set(m))
        {
            dB_[m].write();
        }

        if (integralPhiSinWinDt_.set(m))
        {
            integralPhiSinWinDt_[m].write();
        }

        if (integralPhiCosWinDt_.set(m))
        {
            integralPhiCosWinDt_[m].write();
        }

        if (integralSinWinDt_.set(m))
        {
            integralSinWinDt_[m].write();
        }

        if (integralCosWinDt_.set(m))
        {
            integralCosWinDt_[m].write();
        }

        if (avgPhi_.set(m))
        {
            avgPhi_[m].write();
        }
    }

    return true;
}

bool Foam::functionObjects::inSituFourier::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    Log << type() << " " << name() <<  " read:" << nl;

    fieldName_ = dict.lookupOrDefault<word>("field", "p");

    doVolume_ = dict.lookupOrDefault<Switch>("doVolume", true);

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    patchSet_.clear();
    if (dict.found("patches"))
    {
        patchSet_ =
            pbm.patchSet(wordReList(dict.lookup("patches")), false, true);
    }
    else
    {
        forAll(pbm, patchi)
        {
            if (!pbm[patchi].coupled())
            {
                patchSet_.insert(patchi);
            }
        }
    }

    //we can, otherwise use user-specified value
    t0_ = dict.lookupOrDefault("timeStart",mesh_.time().startTime().value());

    dBRefVal_.set
    (
        new dimensionedScalar
        (
            dict.lookupOrDefault
            (
                "dBRefVal",
                dimensionedScalar("dBRefVal", dimPressure, 20e-6)
            )
         )
    );
    ecf_ = dict.lookupOrDefault("convFactor", sqrt(8.0/3.0));

    fm_.clear();
    flow_.clear();
    fhigh_.clear();
    if (readBool(dict.lookup("generateFrequencyBands")))
    {
        scalar fl = readScalar(dict.lookup("fl"));
        if (fl <= 0)
        {
            FatalErrorInFunction
                << "fl entry must be positive."
                << endl << exit(FatalError);
        }
        scalar fU = readScalar(dict.lookup("fU"));

        scalar bandsPerOctave = dict.lookupOrDefault("bandsPerOctave", 3.0);

        scalar cbrt2 = Foam::pow(2.0, 1.0/bandsPerOctave);
        scalar sixthrt2 = Foam::sqrt(cbrt2);

        label pp = 0;
        do
        {
            fm_.append(fl*pow(cbrt2,pp));
            flow_.append(1.0/sixthrt2*fm_.last());
            fhigh_.append(flow_.last()*cbrt2);
            pp++;
        } while (fm_.last() < fU);

    }
    else
    {
        const vector2DField bands(dict.lookup("frequencyBands"));
        forAll(bands,m)
        {
            if (bands[m].y() <= bands[m].x())
            {
                FatalErrorInFunction
                    << "band entry incorrectly ordered."
                    << endl << exit(FatalError);
            }
            flow_.append(bands[m].x());
            fhigh_.append(bands[m].y());
            fm_.append(sqrt(flow_.last()*fhigh_.last()));
        }
        if (dict.found("centreFrequencies"))
        {
            const scalarField centreFrequencies
            (
                dict.lookup("centreFrequencies")
            );
            if (centreFrequencies.size() != bands.size())
            {
                FatalErrorInFunction
                    << "frequencyBands and centreFrequencies entries"
                    <<" must be same length."
                    << endl << exit(FatalError);
            }
            forAll(centreFrequencies,m)
            {
                fm_[m] = centreFrequencies[m];
            }
        }

    }

    Info<<"Band number"<<tab<<"f"<<tab<<"flow"<<tab<<"fhigh"<<endl;
    forAll(fm_, i)
    {
        Info<<tab<<i<<tab<<fm_[i]<<tab<<flow_[i]
            <<tab<<fhigh_[i]<<endl;
    }

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
