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
    (c) 2011-2014 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "acousticAnalogy/FfowcsWilliamsHawkings.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(FfowcsWilliamsHawkings, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        FfowcsWilliamsHawkings,
        dictionary
    );

}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::FfowcsWilliamsHawkings::writeFileHeader
(
    Ostream& os
)
const
{
    writeHeader(os, "FfowcsWilliamsHawkings Acoustic Analogy");
    writeCommented(os, "Time");
    forAll(observerSets_, obsSetI)
    {
        if (observerSets_[obsSetI].write())
        {
            writeDelimited(os,observerSets_[obsSetI].name());
            if (writeComponent_)
            {
				if (modelName_=="FWH")
				{
					writeDelimited(os,"(thickness loading) ");
				}
				else
				{
					writeDelimited(os,"(thickness loading nonlinear) ");
				}
			}
            writeDelimited(os,"Size ");
            os << observerSets_[obsSetI].size();
        }
    }
    os << endl;
}

void
Foam::functionObjects::FfowcsWilliamsHawkings::writeFfowcsWilliamsHawkings()
{
    Log << "Sum of instantaneous acoustic pressure: " << endl;
    forAll(observerSets_, obsSetI)
    {
        observer& obs = observerSets_[obsSetI];
        scalar pPrimeSum = gSum(obs.pPrime());
        Log << "   -" << obs.name() << ": " << pPrimeSum << " [Pa]" << endl;
    }

    // File output
	label currentTimeStep = floor
    (
        (time_.time().value() - startSamplingTime_) / dtAco_
    );
    scalar stime=currentTimeStep*dtAco_;

   // file() << obr_.time().timeName() << tab << setw(1) << "   ";
    file() <<stime<< tab << setw(1) << "   ";
    forAll(observerSets_, obsSetI)
    {
        observer& obs = observerSets_[obsSetI];
        if (obs.write())
        {
            file() << obs.name() << "    ";
            const scalarField& pPrime = obs.pPrime();
			const scalarField& pthick = obs.pThickness();
			scalarField pload(pPrime);
			if (writeComponent_)
			{
				pload-=pthick;
			}
            forAll(obs.positions(), obsI)
            {
				if (!writeComponent_)
				{
					file() << pPrime[obsI] << "  ";
				}
				else
				{
                    file() << pPrime[obsI] <<" ("<< pthick[obsI]
					<<"  "<<pload[obsI]<<")"
					<< "  ";
				}
            }
        }
    }
    file() << endl;
}


void Foam::functionObjects::FfowcsWilliamsHawkings::initialise()
{
    if (initialised_)
    {
        return;
    }

    if
    (
        !foundObject<volVectorField>(UName_)
     || !foundObject<volScalarField>(pName_)
    )
    {
        FatalErrorInFunction
            << "Could not find " << UName_ << ", " << pName_
            << exit(FatalError);
    }

    scalar totalNumberOfTimeSteps = endSamplingTime_-startSamplingTime_;
    totalNumberOfTimeSteps = min
    (
       floor(totalNumberOfTimeSteps/dtAco_) + 1,

        scalar(1E7)
    );

    if (totalNumberOfTimeSteps == 1E7)
    {
        WarningInFunction
            <<"Calculated total number of time steps exceeds 1E7 limit."<<nl
            <<"This maybe because of default controlDict time settings."<<nl
            <<"Try resetting the functionObject startSamplingTime or "
            <<"endSamplingTime keywords or increase the physical timestep."<<nl
            <<"The simulation will continue assuming number of timesteps will "
            <<"not exceed 1E7 limit"<<nl
            <<endl;
    }

	label sizeTime=floor(totalNumberOfTimeSteps);


	label	sizeThk=sizeTime;

    forAll(observerSets_, obsI)
    {
        pPrime_[obsI].setSize(
                                observerSets_[obsI].size(),
                                scalarField(sizeTime, 0.)
                             );

		pThickness_[obsI].setSize(
                                observerSets_[obsI].size(),
                                scalarField(sizeThk, 0.)
                             );
    }


    initialised_ = true;
}


Foam::tmp<Foam::volScalarField>
Foam::functionObjects::FfowcsWilliamsHawkings::p() const
{
    const volScalarField& p = lookupObject<volScalarField>(pName_);
    if (p.dimensions() == dimPressure)
    {
        return p;
    }
    else
    {
        return
        (
            rhoRef_ * p
        );
    }
}


Foam::tmp<Foam::volScalarField>
Foam::functionObjects::FfowcsWilliamsHawkings::dpdt() const
{
    const volScalarField& p = lookupObject<volScalarField>(pName_);

    if (p.dimensions() == dimPressure)
    {
        return
        (
            Foam::fvc::ddt(p)
        );
    }
    else
    {
        return
        (
            rhoRef_ * Foam::fvc::ddt(p)
        );
    }
}


Foam::tmp<Foam::volVectorField>
Foam::functionObjects::FfowcsWilliamsHawkings::U() const
{
    return
    (
        lookupObject<volVectorField>(UName_)
    );
}


Foam::tmp<Foam::volVectorField>
Foam::functionObjects::FfowcsWilliamsHawkings::dUdt() const
{
    return
    (
        Foam::fvc::ddt(lookupObject<volVectorField>(UName_))
    );
}


Foam::tmp<Foam::surfaceVectorField>
Foam::functionObjects::FfowcsWilliamsHawkings::n() const
{
    return
    (
        -mesh_.Sf() / mesh_.magSf()
    );
}

Foam::tmp<Foam::scalarField>
Foam::functionObjects::FfowcsWilliamsHawkings::getSurfaceNormalVelocityGradient
(
    const Foam::label &patchI
) const
{
    vectorField np( -mesh_.boundary()[patchI].nf() );
    vectorField np0( -mesh_.Sf().oldTime().boundaryField()[patchI] / mesh_.magSf().oldTime().boundaryField()[patchI] );
    vectorField np00( - mesh_.Sf().oldTime().oldTime().boundaryField()[patchI] / mesh_.magSf().oldTime().oldTime().boundaryField()[patchI] );

    volVectorField U( this->U() );

    // Get time stepping size
    scalar deltaT = mesh_.time().deltaTValue();
    scalar deltaT0;
    if (mesh_.time().timeIndex() < 2)
    {
        // Revert to Euler
        deltaT0 = GREAT;
    }
    else
    {
        deltaT0 = mesh_.time().deltaT0Value();
    }

    const vectorField& Up = U.boundaryField()[patchI];
    const vectorField& Up0 = U.oldTime().boundaryField()[patchI];
    const vectorField& Up00 = U.oldTime().oldTime().boundaryField()[patchI];

    // Calculate coefficients
    scalar coefft   = 1.0 + deltaT/(deltaT + deltaT0);
    scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    scalar coefft0  = coefft + coefft00;

    // Second-order backward-differencing time derivative
    // assuming constant time stepping
    // see: backwardDdtScheme.C

    return
    (
        (
            (
                coefft*(Up&np)
              - coefft0*(Up0&np0)
              + coefft00*(Up00&np00)
            )
            / mesh_.time().deltaTValue()
        )
    );
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::functionObjects::FfowcsWilliamsHawkings::FfowcsWilliamsHawkings
(
    const word& name,
    const Time& runTime,
    const dictionary& dict,
    const word &modelName
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name, typeName, dict),
    initialised_(false),
    patches_(0),
    pName_(word::null),
    UName_(word::null),
    rhoRef_(-1),
    cRef_(-1),
    observerSets_(0),
    pPrime_(0),

    startSamplingTime_
    (
        dict.lookupOrDefault<scalar>
        (
            "startSamplingTime",
            runTime.startTime().value()
        )
    ),
    endSamplingTime_
    (
        dict.lookupOrDefault<scalar>
        (
            "endSamplingTime",
            runTime.endTime().value()
        )
    )
{
	tAcc_=2;
	hasOfield_=false;
    hasOOfield_=false;
    dt_= mesh_.time().deltaTValue();
	dto_=dt_;
	modelName_=modelName;
    writeComponent_=false;
    dtAco_=1.0e-4;
    modelType_="simplified";

    read(dict);
    pPrime_.setSize(observerSets_.size());
	pThickness_.setSize(observerSets_.size());
	setFields();
    writeFileHeader(file());
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


Foam::functionObjects::FfowcsWilliamsHawkings::~FfowcsWilliamsHawkings()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool Foam::functionObjects::FfowcsWilliamsHawkings::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    Log << type() << " " << name() <<  " write:" << nl;

	if (modelName_!="KFWH")
	{
		patches_ =
			mesh_.boundaryMesh().patchSet(wordReList(dict.lookup("patches")));
	}

    pName_ = dict.lookupOrDefault<word>("pName", "p");

    UName_ = dict.lookupOrDefault<word>("UName", "U");

    rhoRef_ = readScalar(dict.lookup("rhoRef"));

    cRef_ = readScalar(dict.lookup("cRef"));

	scalar dt0= mesh_.time().deltaTValue();
	dtAco_=dict.lookupOrDefault<scalar>("dt", dt0);
	modelType_=dict.lookupOrDefault<word>("modelType", "simplified");
	writeComponent_=dict.lookupOrDefault<bool>("writeComponent", false);//write noise component thickness, loading

	if (modelType_=="full")
	{
		writeComponent_=true;
	}

    // Read observers
    {
        const dictionary& obsDict = dict.subDict("observers");
        wordList obsNames = obsDict.toc();

        observerSets_.setSize(obsNames.size());
        label index(0);

        forAll(obsNames, obsI)
        {
            word oName = obsNames[obsI];

            observerSets_.set
            (
                index++,
                observer::New
                (
                    oName,
                    mesh_,
                    dict.subDict("observers").subDict(oName)
                )
            );
        }
    }

    return true;
}


bool Foam::functionObjects::FfowcsWilliamsHawkings::execute()
{
    if
    (
        time_.time().value() >= startSamplingTime_
        && time_.time().value() < endSamplingTime_
    )
    {
        Log << type() << " " << name() <<  " execute:" << nl;

        calculate();
        writeFfowcsWilliamsHawkings();
    }
    return true;
}

bool Foam::functionObjects::FfowcsWilliamsHawkings::write()
{
    return true;
}

void Foam::functionObjects::FfowcsWilliamsHawkings::setFields()
{
	volScalarField p( this->p() );
	volVectorField U( this->U() );
	forAllConstIter(labelHashSet, patches_, iter)
	{
		// Get patch ID
		label patchI = iter.key();
		pp_[patchI]=p.boundaryField()[patchI];
		pp0_[patchI]=pp_[patchI];
		pp00_[patchI]=pp0_[patchI];
		Up_[patchI]=U.boundaryField()[patchI];
		Up0_[patchI]=Up_[patchI];
		Up00_[patchI]=Up0_[patchI];
	}
}

void Foam::functionObjects::FfowcsWilliamsHawkings::updateFields()
{
	volScalarField p( this->p() );
	volVectorField U( this->U() );
	forAllConstIter(labelHashSet, patches_, iter)
	{
		// Get patch ID
		label patchI = iter.key();
		pp00_[patchI]=pp0_[patchI];
		Up00_[patchI]=Up0_[patchI];

		pp0_[patchI]=pp_[patchI];
		Up0_[patchI]=Up_[patchI];
		pp_[patchI]=p.boundaryField()[patchI];
		Up_[patchI]=U.boundaryField()[patchI];
	}
}

void Foam::functionObjects::FfowcsWilliamsHawkings::calUdotPdot
(
	std::map<label,vectorField>& Udot,
	std::map<label,scalarField>& pdot
)
{
	dt_=mesh_.time().deltaTValue();
	std::map<label,vectorField>::iterator it;
	if (tAcc_==1)
	{
		tmp<volScalarField> tdpdt(this->dpdt());
		tmp<volVectorField> tdUdt(this->dUdt());

		const volScalarField& dpdt = tdpdt;
		const volVectorField& dUdt = tdUdt;

		for (it=Up_.begin(); it!=Up_.end();it++)
		{
			label patchI=it->first;
			const scalarField& dotp = dpdt.boundaryField()[patchI];
			pdot[patchI]=dotp;
			const vectorField &dUdtp = dUdt.boundaryField()[patchI];
			Udot[patchI]=dUdtp;
		}
		return;
	}


	scalar coefft   = 1.0/dt_ + 1.0/(dt_ + dto_);
    scalar coefft00 = dt_/(dto_*(dt_ + dto_));
    scalar coefft0  = coefft + coefft00;


    for (it=Up_.begin(); it!=Up_.end();it++)
    {
		label patchI=it->first;
		const vectorField &Up=Up_[patchI];
		const vectorField &Up0=Up0_[patchI];
		const vectorField &Up00=Up00_[patchI];
		vectorField dotUp( coefft*Up- coefft0*Up0+ coefft00*Up00 );
		Udot[patchI]=dotUp;

		const scalarField &pp=pp_[patchI];
		const scalarField &pp0=pp0_[patchI];
		const scalarField &pp00=pp00_[patchI];
		scalarField dotp( coefft*pp- coefft0*pp0+ coefft00*pp00 );
		pdot[patchI]=dotp;
	}
}

void Foam::functionObjects::FfowcsWilliamsHawkings::updateState()
{
	dto_=dt_;
	if (hasOfield_)
	{
		hasOOfield_=true;
	}
	hasOfield_=true;
}

void Foam::functionObjects::FfowcsWilliamsHawkings::calculate()
{
    initialise();
	updateFields();

	//second-order time derivatives
	std::map<label,vectorField> Udot;
	std::map<label,scalarField> pdot;
	calUdotPdot( Udot,pdot);

    // Pressure field and its time derivative
    volScalarField p( this->p() );
    volScalarField dpdt( this->dpdt() );

    // Velocity field and its time derivative
    volVectorField U( this->U() );
    volVectorField dUdt( this->dUdt() );

    // Surface normal vector and its time derivative
    surfaceVectorField n( this->n() );

    // Calculate constant
    scalar coeff = 1.0 / ( 4.0 * Foam::constant::mathematical::pi );

    // Loop over all observers
    forAll(observerSets_, obsSetI)
    {
        observer& obsSet = observerSets_[obsSetI];
        const vectorField& pos = obsSet.positions();
        List<scalarList>&  pPrimeSet = pPrime_[obsSetI];
		List<scalarList>&  pThicknessSet = pThickness_[obsSetI];

        scalarField& pPrimeObsCur = obsSet.pPrime();
		scalarField& pThicknessObsCur = obsSet.pThickness();

        forAll(pos, obsI)
        {
            scalar pPrime = 0.0;
            scalarList& pPrimeObs = pPrimeSet[obsI];
			scalarList& pThicknessObs = pThicknessSet[obsI];
            // Surface integral - loop over all patches
            forAllConstIter(labelHashSet, patches_, iter)
            {
                // Get patch ID
                label patchI = iter.key();

                // Surface area vector and face center at patch
                scalarField magSf = mesh_.magSf().boundaryField()[patchI];
                vectorField Cf = mesh_.Cf().boundaryField()[patchI];

                // Normal vector pointing towards fluid
                vectorField np = n.boundaryField()[patchI];

                // Pressure field and time derivative at patch
                const scalarField& pp = p.boundaryField()[patchI];

				const scalarField& dpdtp= pdot[patchI];
                // Velocity field ant time derivative at patch
                vectorField Up = U.boundaryField()[patchI];
               // vectorField dUdtp = dUdt.boundaryField()[patchI];
				const vectorField dUdtp =Udot[patchI];

                // Distance surface-observer
                vectorField rHat( pos[obsI] - Cf );
                scalarField r( mag(rHat) );
                rHat /= r;

                // Surface normal velocity
                scalarField Un( Up & np );
                scalarField UnDot( getSurfaceNormalVelocityGradient(patchI) );

                // Mach number
                vectorField M( Up / cRef_ );
                scalarField Mr( M & rHat );
                scalarField Mn( M & np );
                scalarField MrDot( rHat & (dUdtp / cRef_) );

                // Intermediate variables
                vectorField Li( pp*np );

				vectorField dLidt( dpdtp*np );

                forAll(p.boundaryField()[patchI], fI)
                {
                    scalar distance = r[fI];
                    const scalar& currentTime = time_.time().value();
                    scalar timeDelay = distance/cRef_;
                    scalar advancedTime = currentTime + timeDelay;
                    if (advancedTime < endSamplingTime_)
                    {
                        // Thickness noise

                        scalar pthk1=UnDot[fI] * magSf[fI]
							/ (r[fI]*sqr(1.0-Mr[fI]))*rhoRef_*coeff;

                        pPrime = pthk1;

						scalar pthk2=Un[fI]*magSf[fI]
                         *(r[fI]*MrDot[fI] + cRef_*(Mr[fI]-magSqr(M[fI])))
                            / (sqr(r[fI])*pow3(1.0-Mr[fI])) * rhoRef_ *coeff;

                        pPrime += pthk2;

                        // Loading noise
                        pPrime += magSf[fI]*coeff*(dLidt[fI]&rHat[fI])
                            / (r[fI]*sqr(1.0-Mr[fI])) / cRef_;

                        pPrime += magSf[fI]*coeff*((Li[fI]&rHat[fI])
                                                   - (Li[fI]&M[fI]))
                          / (sqr(r[fI])*sqr(1.0 - Mr[fI]));

                        pPrime += magSf[fI]*coeff
                          *((Li[fI] & rHat[fI])
                            *(r[fI]*MrDot[fI] + cRef_*(Mr[fI]-magSqr(M[fI]))))
                           / (sqr(r[fI])*pow3(1.0 - Mr[fI])) / cRef_;

                        label previousTimeStep =
                            floor
                            (
                                (advancedTime-startSamplingTime_)
                              / dtAco_
                            );
                        scalar previousTime = startSamplingTime_
                            + previousTimeStep*dtAco_;

                        scalar weight = 1. -
                            (advancedTime - previousTime)/dtAco_;
                        pPrimeObs[previousTimeStep] += weight*pPrime;

                        if (writeComponent_)
                        {
							pThicknessObs[previousTimeStep] += weight*(pthk1+pthk2);
						}
                        if (previousTimeStep + 1 < pPrimeObs.size())
                        {
                            pPrimeObs[previousTimeStep+1] += (1-weight)*pPrime;
                            if (writeComponent_)
							{
								pThicknessObs[previousTimeStep+1] += (1-weight)*(pthk1+pthk2);
							}
                        }
                    }
                }
            }

            label currentTimeStep = floor
            (
                (time_.time().value() - startSamplingTime_) / dtAco_
            );

            pPrime = pPrimeObs[currentTimeStep];

            reduce(pPrime, sumOp<scalar>());
            pPrimeObsCur[obsI] = pPrime;

            if (writeComponent_)
            {
				scalar pThk= pThicknessObs[currentTimeStep];
				reduce(pThk, sumOp<scalar>());
				pThicknessObsCur[obsI] = pThk;
			}
        }
    }
    updateState();
}


// ************************************************************************* //
