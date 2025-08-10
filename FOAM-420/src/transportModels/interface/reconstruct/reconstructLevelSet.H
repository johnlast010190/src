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
    Foam::reconstructLevelSet

SourceFiles
    reconstruct.C

Authors
    Daniel Rettenmaier < rettenmaier@gsc.tu-darmstadt.de>
    Vignesh Thammanna Gurumurthy < vignesh@sla.tu-darmstadt.de >

    All rights reserved.


S-SCLSVOF - Density Balanced - Method Developement:

    Albadawi A, Donoghue DB, Robinson AJ, Murray DB, Delauré YMC.
    Influence of surface tension implementation in volume of fluid and coupled volume of fluid with level set methods for bubble growth and detachment. International
    Journal of Multiphase Flow 2013; 53:11–28.

    T. YAMAMOTO, Y. OKANO AND S. DOST
    Validation of the S-CLSVOF method with the density-scaled balanced continuum surface force model in multiphase systems coupled with thermocapillary flows
    Int. J. Numer. Meth. Fluids 2017; 83:223-244.

   Yokoi K.
    A density-scaled continuum surface force model within a balanced force formulation. Journal of
    Computational Physics 2014; 278:221–228.

Description
    Reconstructs the interface with an 0.5 alpha iso surface by interpolating
    cell centered alpha values to the cell points. At edges of the cell a
    0.5-alpha point is interpolated linearely if the alpha value from the two
    edge defining points is higher and lower then 0.5.
    Calculates the a distance field to the interface in a band of cells near
    the interface (if selected)

    You may refer to this software as :
    //- full bibliographic data to be provided

    This code has been developed by :
        Rettenmaier adapting the code bases of Kunkelmann (main developer).

    Method Development and Intellectual Property :
        Rettenmaier<rettenmaier@gsc.tu-darmstadt.de>
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


#ifndef reconstructLevelSet_H
#define reconstructLevelSet_H

#include "reconstruct/reconstruct.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class reconstructLevelSet Declaration
\*---------------------------------------------------------------------------*/

class reconstructLevelSet
:   public reconstruct
{

private:

    //- const reference to the mesh
    const fvMesh& mesh_;

    volScalarField psi_;

    surfaceVectorField gradPsif_;

    volScalarField psi0_;
    volScalarField delta_;
    volScalarField H_;
    volScalarField Hscaled_;

    //- Parameter for S-CLSVOF method (deltaX,gamma,epsilon,deltaTau)

    //- Grid spacing
    volScalarField      deltaX_;

    //- artificial time
    volScalarField      deltaTau_;

    //- non-dim number depended on grid spacing
    volScalarField      gamma_;

    //- interface thickness coefficient
    //  optimal value lies between 2 and 3
    scalar              epsC_;

    //- value of the interface-thickness depending on deltaX
    volScalarField      epsilon_;

    //- scale delta and H to reduce the parasitic current
    //  dependency on the density ratio
    Switch scaled_;

    const scalar pi = Foam::constant::mathematical::pi;

    //- Field on where to apply the surface tension force
    void calculateDelta();
    void calculateH();

    //- Initializes constants, deltaX, gamma, epsilon, deltaTau
    void initConstants();

    //- Recalculates the LevelSet psi
    void solveReinitializationEqn();


public:

    // Runtime type name
        TypeName("reconstructLevelSet");

    // Constructors

        //- Construct from components
        reconstructLevelSet
        (
            const word& name,
            const volScalarField& alpha,
            const dictionary& transpProp,
            const List<bool>& isWallPatch,
            const volScalarField& isInterface
        );

    //- Destructor
    virtual ~reconstructLevelSet(){}

    // Member Functions

        //- calculates interface normal nHatv and nHatfv
        void reconstructInterface();

    //- returns nHatf projected on the face normals
    //  OVERWRIDES nHatf() in reconstruct.H
    //  The via LevelSet reconstructed nHatf may not be fit for compression!
    const surfaceScalarField& nHatf()
    {

/*
        const volVectorField gradAlpha(fvc::grad(alpha1_, "nHat"));

        //- Interpolated face-gradient of alpha
        surfaceVectorField gradAlphaf(fvc::interpolate(gradAlpha));

        nHatf_  = (gradAlphaf/(mag(gradAlphaf) + deltaN_)) & mesh_.Sf();

        nHatf_.boundary FieldRef() = nHatfv_.boundaryField() & mesh_.Sf().boundaryField();
*/

        nHatf_  = nHatfv_ & mesh_.Sf();

        return nHatf_;
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
