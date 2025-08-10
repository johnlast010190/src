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
    (c) held by original author

Class
    Foam::curvature

SourceFiles
    curvature.C

Authors
    Daniel Rettenmaier < rettenmaier@gsc.tu-darmstadt.de>
    Daniel Deising     < deising@mma.tu-darmstadt.de>
    All rights reserved.

Description
    Base class for different curvature models for implicit or explicit
    interface representations
    Provides smoothing functionality for curvature fields to reduce
    parasitic currents

    You may refer to this software as :
    //- full bibliographic data to be provided

    This code has been developed by :
        Daniel Rettenmaier <rettenmaier@gsc.tu-darmstadt.de> (main developer).

    Method Development and Intellectual Property :
        Daniel Rettenmaier <rettenmaier@gsc.tu-darmstadt.de>
        Daniel Deising <deising@mma.tu-darmstadt.de>
        Holger Marschall <marschall@csi.tu-darmstadt.de>
        Dieter Bothe <bothe@csi.tu-darmstadt.de>
        Cameron Tropea <ctropea@sla.tu-darmstadt.de>

        Mathematical Modeling and Analysis
        Institute for Fluid Mechanics and Aerodynamics
        Center of Smart Interfaces
        Technische Universitaet Darmstadt

    If you use this software for your scientific work or your publications,
    please don't forget to acknowledge explicitly the use of it.

\*---------------------------------------------------------------------------*/


#ifndef curvature_H
#define curvature_H

#include "db/typeInfo/typeInfo.H"
#include "db/runTimeSelection/construction/runTimeSelectionTables.H"
#include "memory/tmp/tmp.H"
#include "fields/volFields/volFieldsFwd.H"
#include "fields/volFields/volFields.H"
#include "memory/autoPtr/autoPtr.H"
#include "db/IOobjects/IOdictionary/IOdictionary.H"

#include "finiteVolume/fvc/fvc.H"

#include "meshes/polyMesh/syncTools/syncTools.H"
#include "interface/distributeField/distributeField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    class curvature;

    Ostream& operator<<(Ostream&, curvature&);

/*---------------------------------------------------------------------------*\
                         Class curvature Declaration
\*---------------------------------------------------------------------------*/

class curvature
{

protected:

    // Protected data

        volScalarField kappa_;

        word name_;
        const fvMesh& mesh_;
        const volScalarField& alpha1_;
        const volVectorField& nHatv_;
        const surfaceVectorField& nHatfv_;
        const volScalarField& interfaceDensity_;
        const List<bool>& isWallPatch_;
        const volScalarField& isInterface_;

        //- shall cells at symmetryPlanes be smoothened?
        Switch smoothSym_;

        //- special smoothing treatment for symmetry boundaries
        List<bool> attachedToSym_;

        //- shall cells at wallPlanes be smoothened?
        Switch smoothWall_;

        //- special smoothing treatment for wall boundaries
        List<bool> attachedToWall_;

        //- number of smoothing cycles on interfacDesity > 0 cells
        scalar preSmoothCycles_;

        //- number of smoothing cycles on isInterface stencil
        scalar postSmoothCycles_;

        const distributeField distributeField_;

        //- switch for numerical divergence correction
        Switch complexDivCorrection_;

    // Protected Member Functions

        //- set smoothing marker fields
        void initialize();

        void smooth(const boolList& set, volScalarField& kappa, int nLoops);

        //- -fvc::div(nHatf_) + (nHat & fvc::grad(nHatfv) & nHat)*complexDivCorrection_;
        void calcK();

        //- Disallow default bitwise copy construct
        curvature(const curvature&);

        //- Disallow default bitwise assignment
        void operator=(const curvature&);

public:

    // Runtime type name
        TypeName("curvature");

    // Declare run-time selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            curvature,
            dictionary,
            (
                const word& name,
                const volScalarField& alpha,
                const volVectorField& nHatv,
                const surfaceVectorField& nHatfv,
                const volScalarField& interfaceDensity,
                const List<bool>& isWallPatch,
                const dictionary& transpProp,
                const volScalarField& isInterface
            ),
            (
                name,
                alpha,
                nHatv,
                nHatfv,
                interfaceDensity,
                isWallPatch,
                transpProp,
                isInterface
            )
        );

    // Selectors

        //- Return a reference to the selected curvature model
        static autoPtr<curvature> New
        (
            const word& name,
            const volScalarField& alpha,
            const volVectorField& nHatv,
            const surfaceVectorField& nHatfv,
            const volScalarField& interfaceDensity,
            const List<bool>& isWallPatch,
            const dictionary& transpProp,
            const volScalarField& isInterface
        );

    // Constructors

        //- Construct from components
        curvature
        (
            const word& name,
            const volScalarField& alpha,
            const volVectorField& nHatv,
            const surfaceVectorField& nHatfv,
            const volScalarField& interfaceDensity,
            const List<bool>& isWallPatch,
            const dictionary& transpProp,
            const volScalarField& isInterface
        );

    //- Destructor
    virtual ~curvature(){}

    // Member Functions

        //- calculates the curvature
        virtual void calculateK(){}

        //- calculates the force due to the surfacetension
        virtual tmp<surfaceScalarField> surfaceTensionForce(dimensionedScalar sigma) const
        {
            return fvc::interpolate(sigma * kappa_) * fvc::snGrad(alpha1_);
        };

        //- reduce kappa on real interface, smooth, propagate on isInterface
        //  and smooth again
        void prePostSmooth(volScalarField& kappa);

        // Access
        const volScalarField& kappa() const
        {
            return kappa_;
        }

        const surfaceScalarField kf()
        {
            return fvc::interpolate(kappa_)();
        }

        word name() const
        {
            return name_;
        }

        friend Ostream& operator<<(Ostream& os, curvature&);

        void write(Ostream& os) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
