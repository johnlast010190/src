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
    (c) 2017-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "dynamicGIBFvMesh/geometryTransformations/deformingTyreTransformation/deformingTyreTransformation.H"
#include "global/constants/mathematical/mathematicalConstants.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"


namespace Foam
{
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

    defineTypeNameAndDebug(deformingTyreTransformation, 0);
    addToRunTimeSelectionTable
    (
        geometryTransformation,
        deformingTyreTransformation,
        dictionary
    );

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void deformingTyreTransformation::initializeCBox()
{
    // get box
    controlBox& box = vNurbs->getBox();

    // for every depth zone of the box
    for (label k = 0; k < 3; k++)
    {
        // define the box dimensions

        scalar Length = radius_*1.1;
        scalar Depth = tWidth_*1.1;

        scalar aInit = 3.*constant::mathematical::pi/2.-alphaContact_/2.;
        scalar aEnd  = 3.*constant::mathematical::pi/2.+alphaContact_/2.;

        scalar zeta = -Depth/2. + scalar(k)*Depth/2.;

        // define the 27 points on the tyre patch to be deformed

        for (label i = 1; i < 4; i++)
        {
            scalar angle = aInit + scalar(i-1)*(aEnd-aInit)/2.;

            for (label j = 1; j < 4; j++)
            {
                scalar r = radius_ -scalar(j-1)*hTread_/2.;

                vector dv = vX_*r*cos(angle)
                          + gNormal_*r*sin(angle)
                          + tNormal_*zeta
                          + center_;

                box.setDisplacement(i, j, k, dv);
            }
        }

        // define the rest 16 points that create the buffer zones

        // i = 0
        for (label j = 0; j < 5; j++)
        {
            vector v = box.getDisplacement(1, j, k);
            scalar y = (v-center_)&gNormal_;
            vector dv = center_
                      + gNormal_*y//(-Length + scalar(j)*Length/2.)
                      - Length*vX_
                      + tNormal_*zeta;



            box.setDisplacement(0, j, k, dv);
        }

        // i = 4
        for (label j = 0; j < 5; j++)
        {
            vector v = box.getDisplacement(3, j, k);
            scalar y = (v-center_)&gNormal_;
            vector dv = center_
                      + gNormal_*y//(-Length + scalar(j)*Length/2.)
                      + Length*vX_
                      + tNormal_*zeta;

            box.setDisplacement(4, j, k, dv);
        }

        // j = 0
        for (label i = 0; i < 5; i++)
        {
            vector v = box.getDisplacement(i, 1, k);
            scalar x = (v-center_)&vX_;
            vector dv = center_
                      + vX_*x//(-Length + scalar(i)*Length/2.)
                      - Length*gNormal_
                      + tNormal_*zeta;

            box.setDisplacement(i, 0, k, dv);
        }

        // j = 4
        for (label i = 0; i < 5; i++)
        {
            vector v = box.getDisplacement(i, 3, k);
            scalar x = (v-center_)&vX_;
            vector dv = center_
                      + vX_*x//(-Length + scalar(i)*Length/2.)
                      + Length*gNormal_
                      + tNormal_*zeta;

            box.setDisplacement(i, 4, k, dv);
        }
    }

    box.updateBoxPosition();

    return;
}


void deformingTyreTransformation::createDisplacementFields()
{
    //the radius reduction at the contact point has to be less
    //than the tread height
    scalar dy = radius_-radius_*cos(alphaContact_/2.);

    if (dy > hTread_)
    {
        FatalErrorInFunction
            << "Bad geometry settings in tyre definition result in a "
            << "deformation bigger than the tread."
            << exit(FatalError);
    }

    controlBox& box = vNurbs->getBox();

    // Deformation along the ground normal
    for (label k = 0; k < 3; k++)
    {
        vector dv1 = dy*gNormal_;
        vector dv2 = dy/2.*gNormal_;
        box.setDisplacement(2, 1, k, dv1);
        box.setDisplacement(2, 2, k, dv2);
    }

    for (label i = 0; i < 5; i++)
    {
        for (label k = 0; k < 3; k++)
        {
            vector dv = dy/2.*gNormal_;

            box.setDisplacement(i, 0, k, dv);
        }
    }

    // Deformation along the axis of the tyre
    vector b1 = dy/2.*tNormal_;

    box.setDisplacement(2, 2, 2, b1);

    b1 *= -1;

    box.setDisplacement(2, 2, 0, b1);

    return;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

deformingTyreTransformation::deformingTyreTransformation
(
    const dictionary& dict
)
:
    geometryTransformation(dict),
    center_(vector::zero),
    tNormal_(vector::zero),
    gNormal_(vector::zero),
    vX_(vector::zero),
    radius_(0.0),
    hTread_(0.0),
    tWidth_(0.0),
    alphaContact_(0.0),
    distToGround_(0.0),
    vNurbs(nullptr)
{
    //- After zero Initialization attempt to generate the morpher
    bool done = this->read(dict);

    if (!done)
    {
        WarningInFunction
            << "possible warning: dictionary not read "
            << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool deformingTyreTransformation::read
(
    const dictionary& dict
)
{
    // Read the Parametes and allocate the pointer to the morpher
    center_ = dict.lookup("center");
    tNormal_ = dict.lookup("tyreNormal");
    tNormal_ /= mag(tNormal_);
    gNormal_ = dict.lookup("groundNormal");
    gNormal_ /= mag(gNormal_);
    vX_ = gNormal_^tNormal_;
    vX_ /= mag(vX_);
    radius_ = readScalar(dict.lookup("radius"));
    hTread_ = readScalar(dict.lookup("treadHeight"));
    tWidth_ = readScalar(dict.lookup("tyreWidth"));
    distToGround_ = readScalar(dict.lookup("distToGround"));
    alphaContact_ = readScalar(dict.lookup("contactAngle"));

    alphaContact_ *= constant::mathematical::pi/180.;

    if (vNurbs!=nullptr)
    {
        delete vNurbs;
    }

    vNurbs = new volumetricNurbs
    (
        point(0,0,0),
        point(0,0,0),
        5, 2,
        5, 2,
        3, 1
    );

    // Create the initialBox based on tyre dimensions

    this->initializeCBox();

    // Create the user required displacements to the control points

    this->createDisplacementFields();

    return true;
}


tmp<pointField> deformingTyreTransformation::transformPoints
(
    const pointField& stlPoints
) const
{
    // Get a copy of the morpher
    volumetricNurbs *vNurbsCopy = new volumetricNurbs(*vNurbs);

    // Training of the morpher

    scalarField U(stlPoints.size(), 0);
    scalarField V(stlPoints.size(), 0);
    scalarField W(stlPoints.size(), 0);

//    label mark(-1);
//    scalar vMax = 2;

    forAll(stlPoints, pI)
    {
        const point& p = stlPoints[pI];

        point params = vNurbsCopy->invert(p);
        U[pI] = params.x();
        V[pI] = params.y();
        W[pI] = params.z();
/*
        if (V[pI] < vMax)
        {
            mark = pI;
            vMax = V[pI];
        }
        */
    }

//    scalar distToGround = vNurbsCopy->getDisplacement(U[mark], V[mark], W[mark])&gNormal_;

    // After training is complete update the cBox

    controlBox& box = vNurbsCopy->getBox();

    box.updateBoxPosition();

    // Translate the Assembly so that the tyre touches the ground


    for (label i = 0; i < 5; i++)
    {
        for (label j = 0; j < 5; j++)
        {
            for (label k = 0; k < 3; k++)
            {
                vector dv = -distToGround_*gNormal_;
                box.setDisplacement(i, j, k, dv);
            }
        }
    }
    box.updateBoxPosition();


    // return the result

    tmp<pointField> tresult(new pointField(stlPoints.size(), vector::zero));

    pointField& result = tresult();

    forAll(result, pI)
    {
        result[pI] = vNurbsCopy->getPoint(U[pI], V[pI], W[pI]);
    }

    delete vNurbsCopy;

    vNurbsCopy = nullptr;

    return tresult;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// ************************************************************************* //
