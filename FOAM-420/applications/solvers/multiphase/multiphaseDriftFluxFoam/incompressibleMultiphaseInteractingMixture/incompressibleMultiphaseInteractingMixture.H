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
    (c) 2014-2015 OpenFOAM Foundation
    (c) 2010-2020 Esi Ltd.

Class
    Foam::incompressibleMultiphaseInteractingMixture

Description
    A two-phase incompressible transportModel for interacting phases
    requiring the direct evaluation of the mixture viscosity,
    e.g. activated sludge or slurry.

SourceFiles
    incompressibleMultiphaseInteractingMixture.C

\*---------------------------------------------------------------------------*/

#ifndef incompressibleMultiphaseInteractingMixture_H
#define incompressibleMultiphaseInteractingMixture_H

#include "compressibleTransportModel/compressibleTransportModel.H"
#include "incompressible/viscosityModels/viscosityModel/viscosityModel.H"
#include "mixtureViscosityModel/mixtureViscosityModel.H"
#include "db/IOobjects/IOdictionary/IOdictionary.H"
#include "containers/Dictionaries/PtrDictionary/PtrDictionary.H"
#include "phase/phase.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
         Class incompressibleMultiphaseInteractingMixture Declaration
\*---------------------------------------------------------------------------*/

class incompressibleMultiphaseInteractingMixture
:
    public IOdictionary,
    public compressibleTransportModel
{
protected:

    // Protected data

        // disperse (solids) phase fraction
        volScalarField alphad_;

        //- Dictionary of phases
        PtrDictionary<phase> phases_;

        // continuous phase name
        word continuousPhase_;

        // gaseous phase name
        word gasPhase_;

        // mixture viscosity model
        autoPtr<mixtureViscosityModel> muModel_;

        //- Optional maximum dispersed phase-fraction (e.g. packing limit)
        scalar alphaMax_;

        const volVectorField& U_;
        const surfaceScalarField& phi_;

        // mixture viscosity
        volScalarField mu_;

        // drift flux ddt argument (store to get access to oldTime)
        volVectorField Hdm_;


public:

    TypeName("incompressibleMultiphaseInteractingMixture");


    // Constructors

        //- Construct from components
        incompressibleMultiphaseInteractingMixture
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~incompressibleMultiphaseInteractingMixture()
    {}


    // Member Functions

        //- Optional maximum phase-fraction (e.g. packing limit)
        //  Defaults to 1
        scalar alphaMax() const
        {
            return alphaMax_;
        }

        //- Return access to the mixture VOF fields
        PtrDictionary<phase>& phases()
        {
            return phases_;
        }

        //- Return access to continuous phase
        volScalarField& alphac()
        {
            for
            (
                PtrDictionary<phase>::iterator iter = phases_.begin();
                iter != phases_.end(); ++iter
            )
            {
                if (iter().name() == continuousPhase_)
                {
                    return iter();
                }
            }

            FatalErrorInFunction
                << "Phase " << continuousPhase_ << " does not exist."
                << abort(FatalError);

            return alphad_;
        }

        //- Return access to continuous phase
        phase& phasec()
        {
            for
            (
                PtrDictionary<phase>::iterator iter = phases_.begin();
                iter != phases_.end(); ++iter
            )
            {
                if (iter().name() == continuousPhase_)
                {
                    return iter();
                }
            }

            FatalErrorInFunction
                << "Phase " << continuousPhase_ << " does not exist."
                << abort(FatalError);

            return phases_[0];
        }

        //- Return const-access to the mixture phase
        const volScalarField& alphad() const
        {
            return alphad_;
        }

        //- Return continuous phase name
        const word& continuousPhase() const
        {
            return continuousPhase_;
        }

        //- Return gas phase name
        const word& gasPhase() const
        {
            return gasPhase_;
        }

        //- Return const-access to the mixture velocity
        const volVectorField& U() const
        {
            return U_;
        }

        //- Return the dynamic mixture viscosity
        tmp<volScalarField> mu() const
        {
            return mu_;
        }

        //- Return the dynamic mixture viscosity for patch
        virtual tmp<scalarField> mu(const label patchi) const
        {
            return mu_.boundaryField()[patchi];
        }

        //- Return the mixture density
        virtual tmp<volScalarField> rho() const;

        //- Return the mixture density for patch
        virtual tmp<scalarField> rho(const label patchi) const;

        //- Return the mixture viscosity
        virtual tmp<volScalarField> nu() const
        {
            return mu_/rho();
        }

        //- Return the mixture viscosity for patch
        virtual tmp<scalarField> nu(const label patchi) const
        {
            return mu_.boundaryField()[patchi]/rho(patchi);
        }

        //- Return the mixture heat capacity
        virtual tmp<volScalarField> Cp() const
        {
            NotImplemented;
        }

        //- Return the mixture conductivity
        virtual tmp<volScalarField> lambda() const
        {
            NotImplemented;
        }

        //- Return the mixture Prandtl number
        virtual tmp<volScalarField> Prt() const
        {
            NotImplemented;
        }

        //- Correct the laminar viscosity
        virtual void correct();

        //- Correct the drift velocities
        virtual void correctUdm()
        {
            forAllIter(PtrDictionary<phase>, phases_, iter)
            {
                iter().correctUdm();
            }
            correctHdm();
        }

        //- Return drift velocity addition
        tmp<volVectorField> Hvol() const;

        //- Return drift velocity addition
        tmp<volVectorField> Hmass() const;

        //- Return drift velocity ddt term
        const volVectorField& Hdm() const
        {
            return Hdm_;
        }

        //- compute drift velocity ddt term
        void correctHdm()
        {
            tmp<volVectorField> Hm = Hmass();
            tmp<volVectorField> H = Hvol();
            Hdm_ = rho()*H-Hm;
        }

        // add drift flux momentum sources
        tmp<volSymmTensorField> tauDrift() const;

        //- Read base transportProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
