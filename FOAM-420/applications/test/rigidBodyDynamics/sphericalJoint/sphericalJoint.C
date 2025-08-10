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
    (c) 2016 OpenFOAM Foundation

Application
    sphericalJoint

Description
    Simple spherical-joint pendulum.

\*---------------------------------------------------------------------------*/

#include "rigidBodyMotion/rigidBodyMotion.H"
#include "bodies/masslessBody/masslessBody.H"
#include "bodies/sphere/sphere.H"
#include "joints.H"
#include "restraints/restraint/rigidBodyRestraint.H"
#include "rigidBodyModelState/rigidBodyModelState.H"
#include "db/IOstreams/Fstreams/IFstream.H"
#include "db/IOstreams/Fstreams/OFstream.H"

using namespace Foam;
using namespace RBD;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    dictionary sphericalJointDict(IFstream("sphericalJoint")());

    // Create the sphericalJoint model from dictionary
    rigidBodyMotion sphericalJoint(sphericalJointDict);

    label nIter(readLabel(sphericalJointDict.lookup("nIter")));

    Info<< sphericalJoint << endl;

    // Create the joint-space force field
    scalarField tau(sphericalJoint.nDoF(), Zero);

    // Create the external body force field
    Field<spatialVector> fx(sphericalJoint.nBodies(), Zero);

    // Set the angle of the pendulum to 0.3rad
    sphericalJoint.joints()[1].unitQuaternion
    (
        quaternion(quaternion::ZYX, vector(0.3, 0, 0)),
        sphericalJoint.state().q()
    );

    // Set the gravitational acceleration
    sphericalJoint.g() = vector(0, -9.81, 0);

    OFstream omegaFile("omegaVsTime");

    // Integrate the motion of the sphericalJoint for 4.1s
    scalar deltaT = 0.01;
    for (scalar t=0; t<4.1; t+=deltaT)
    {
        sphericalJoint.newTime();

        for (label i=0; i<nIter; i++)
        {
            sphericalJoint.solve(deltaT, tau, fx);
        }

        // Write the results for graph generation
        // using 'gnuplot sphericalJoint.gnuplot'
        omegaFile
            << t << " "
            << sphericalJoint.joints()[1].unitQuaternion
               (
                   sphericalJoint.state().q()
               ).eulerAngles(quaternion::ZYX).x()
            << endl;
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
