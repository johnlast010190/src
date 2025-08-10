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
    (c) 2010-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "ventilationDrag/ventilationDrag.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/derived/wallVelocity/wallVelocityFvPatchVectorField.H"
#include "referenceFrames/motionCoordinateFrame/motionCoordinateFrame.H"
#include "solidBodyMotionFunctions/rotatingWheelMotion/rotatingWheelMotion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(ventilationDrag, 0);
    addToRunTimeSelectionTable(functionObject, ventilationDrag, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::ventilationDrag::ventilationDrag
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    forces(name, runTime, dict, true),
    Uname_("U"),
    ventilationDragFilePtr_(nullptr)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::ventilationDrag::~ventilationDrag()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::ventilationDrag::read(const dictionary& dict)
{
    // The inherited forces function object requires a reference point /
    // center of rotation (CofR) to be defined in the dictionary, but this
    // reference point / center of rotation (CofR) is not a member variable of
    // forces. The following meets the requirement of a referencePoint or CofR
    // entry in the dictionary.
    dictionary newDict(dict.parent(), dict);
    newDict.add("CofR", vector::zero);

    forces::read(newDict);


    // Get access to the velocity field
    Uname_ = dict.lookupOrDefault<word>("Uname", "U");

    const volVectorField& U = lookupObject<volVectorField>(UName_);

    // Initialise prevCoorFramePtr
    const coordinateFrame* prevCoorFramePtr = nullptr;

    // check that patch(es) is/are indeed specified for the ventilation drag 
    // calculation
    if (patchSet_.size() == 0)
    {
        FatalErrorInFunction
                << "No patches are specified for the Ventilation Drag "
                << "calculation."
                << exit(FatalError);
    }


    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        label patchi = iter.key();

        const fvPatchField<vector>& pf = U.boundaryField()[patchi];
        const wallVelocityFvPatchVectorField& wvpvf =
            refCast<const wallVelocityFvPatchVectorField>(pf);

        coorFramePtr_ = wvpvf.coorFramePtr();

        // check that the patch(es) do indeed belong to a reference frame
        if (coorFramePtr_ == nullptr)
        {
            FatalErrorInFunction
                << "Patches used with the Ventilation Drag function object "
                << "all have to belong to a (the same) rotational reference "
                << "frame."
                << exit(FatalError);
        }

        // check that all patches belong to the same reference frame
        if (prevCoorFramePtr == nullptr || prevCoorFramePtr == coorFramePtr_)
        {
            prevCoorFramePtr = coorFramePtr_;
        }
        else if (prevCoorFramePtr != coorFramePtr_)
        {
            FatalErrorInFunction
                << "All the patches specified for use with the "
                << "ventilationDrag function object must belong to the same "
                << "reference frame"
                << exit(FatalError);
        }
    };


    // Free stream velocity magnitude (from input dictionary)
    vector Uinf = dict.lookup("Uinf");
    if (mag(Uinf) == 0.0)
    {
        FatalErrorInFunction
            << "mag(Uinf) must be greater than zero."
            << exit(FatalError);
    }
    magUInf_ = mag(Uinf);


    // Reference area (from input dictionary)
    dict.lookup("referenceArea") >> Aref_;
    if (Aref_ <= 0.0)
    {
        FatalErrorInFunction
            << "referenceArea must be greater than zero."
            << exit(FatalError);
    }



    Info<< "    Creating ventilation drag file." << nl << endl;

    //file update
    if (Pstream::master() || !Pstream::parRun())
    {
        fileName ventilationDragFileName("ventilationDrag");

        // Open new file at startup
        ventilationDragFilePtr_ = createFile(ventilationDragFileName);

        //add headers
        writeCommented(ventilationDragFilePtr_(),"Time");
        writeDelimited(ventilationDragFilePtr_(),"ventilationDragPower");
        writeDelimited(ventilationDragFilePtr_(),"ventilationDragCoefficient");

        ventilationDragFilePtr_() <<endl;
    }

    return true;
}


bool Foam::functionObjects::ventilationDrag::execute()
{
    // Moment computation

    // Ventilation drag is computated using the moment about the center of
    // rotation i.e. axis of the rotating reference frame of the patches,
    //
    // M_ventDrag = sumOf ( r_CofR x F ).
    //
    // The definition of the position vector is r_CofR = c_f - c_CofR , where
    // c_f is the cell centers on the patches and c_CofR is the center of
    // rotation of the patches. Substituting this into the equation above and
    // applying the distributive property of the cross product gives
    //
    // M_ventDrag = sumOf [ (c_f - c_CofR) x F ]
    // M_ventDrag = sumOf (c_f x F) - sumOf (c_CofR x F) .
    //
    // Since all the patches in a ventilationDrag function object belong to the
    // same rotating reference frame with the same center of rotation, c_CofR is
    // constant and can be moved out of the summation
    //
    // M_ventDrag = sumOf [ (c_f x F) ] - [ c_CofR x sumOf (F) ].
    //
    // The variables in the second term, c_CofR and sumOf (F), are
    // available from the reference frame and the forces function object
    // respectively. The first term can be rewritten as
    //
    // sumOf [ (c_f - vector::zero) x F ]
    //
    // which is equal to the moment around the origin computed by the forces
    // function object, since the forces function object computes the moment
    // about the (0, 0, 0) point that was added as the CofR in the dictionary
    // in the read() function.


    // Use the available member function of the forces function object
    // to compute the moment about the CofR read in from the dictionary in the
    // read() function as (0,0,0)
    forces::calcForcesMoment();
    vector momentAboutZero = forces::momentEff();

    // CofR of the reference frame
    vector centerOfRotation = (*coorFramePtr_).CofR();

    // Angular velocity of the reference frame multiplied with a unit vector
    // in the direction of the axis
    vector wheelAngularVelocity = (*coorFramePtr_).Omega();

    // Use the available member function of the forces function object to
    // compute the sum of the forces on the patches
    vector forceSumPatches = forces::forceEff();

    // Compute the second term in the ventilation drag moment by taking the
    // cross product of the reference frame CofR and sum of the forces on the
    // patches
    vector momentSecondTerm = centerOfRotation ^ forceSumPatches;

    // Ventilation drag moment
    //     = M_ventDrag = (c_f - c_CofR) x F_sumPatches
    //     = M_ventDrag = (c_f x F_sumPatches) - (c_CofR x F_sumPatches)
    //     = M_ventDrag = ((c_f - vector::zero) x F_sumPatches)
    //                      - (c_CofR x F_sumPatches)
    vector momentVentDrag = momentAboutZero - momentSecondTerm;

    // Dynamic pressure
    scalar pDyn = 0.5*rhoRef()*magUInf_*magUInf_;

    // Ventilation drag power
    //      = moment perpendicular to the wheel axis x rotational speed
    //      = moment dot (angular velocity x unit vector in the axis direction)
    scalar ventilationDragPower = momentVentDrag & -wheelAngularVelocity;

    // Ventilation drag coefficient
    // = ventilation drag power / (0.5 x rhoInf x Uinf^3 x referenceArea)
    // = ventilation drag power / (pDyn x Uinf x Aref)
    scalar ventilationDragCoeff = ventilationDragPower/(pDyn*magUInf_*Aref_);

    // Output results to the terminal
    Log << type() << " " << name() << " execute:" << nl;

    Log << " Ventilation drag power         : " << ventilationDragPower << nl
        << " Ventilation drag coefficient   : " << ventilationDragCoeff << nl
        << endl;


    // Write results to file
    if (writeToFile() && Pstream::master())
    {
        writeTime(ventilationDragFilePtr_());
        ventilationDragFilePtr_()
        << tab << ventilationDragPower
        << tab << ventilationDragCoeff
        << endl;
    }


    // Write state/results information
    {
        setResult("Ventilation drag power", ventilationDragPower);
        setResult("Ventilation drag coefficient", ventilationDragCoeff);
    }

    return true;
}