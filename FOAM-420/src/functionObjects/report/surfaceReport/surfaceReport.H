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
    (c) 2010-2012 Esi Ltd.

Class
    Foam::functionObjects::surfaceReport

Group
    grpReportFunctionObjects

Description
    List of patches and constants to use for lift/drag calculations
    Call surfaceReport.write() to calculate&write files.

See also
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    surfaceReport.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_surfaceReport_H
#define functionObjects_surfaceReport_H

#include "functionObjects/fvMeshFunctionObject/fvMeshFunctionObject.H"
#include "db/Time/Time.H"
#include "cfdTools/general/include/fvCFD.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "containers/Lists/PackedList/PackedList.H"
#include "convergenceTermination/convergenceTermination.H"
#include "global/columnatedFileWriter/columnatedFileWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                           Class surfaceReport Declaration
\*---------------------------------------------------------------------------*/

class surfaceReport
:
    public fvMeshFunctionObject
{
protected:

    // Protected data

        // Read from dictonary

            dictionary dict_;

            Switch homogeneity_;

            Switch fluxWeighting_;
            Switch directionAware_;
            Switch addDirectionMeanOutput_;
            word fluxName_;

            wordList fields_;

            //- Swirl setup
            Switch swirl_;
            vector swirlUnitNorm_;
            vector swirlCentre_;

            //- Tumble setup
            Switch tumble_;
            // z-Direction
            vector tumbleDir_;
            // Origin/centre of tumble plane
            vector tumbleOrig_;
            // rho name
            word rhoName_;
            // rho reference
            scalar rhoRef_;
            // Engine stroke length
            scalar Ls_;
            // Average over mass or area
            bool massAveraged_;

            autoPtr<vector> surfaceDirectionPtr_;

        //- Result values
            PackedList<2> reportFaces_;
            Switch reportContainsInternalFaces_;
            Switch compressible_;

            Switch backFlowReport_;
            scalar backFlowArea_;
            scalar backFlowPer_;

            scalar totalSurface_;
            scalar totalFlux_;
            scalar relativeFlux_;
            List< FixedList<scalar, 5>> data_;

            //Check for convergence to abort calculation
            PtrList<convergenceTermination> convergenceChecks_;

            scalar swirlData_;
            scalar tumbleValue_;

        // Whether we need to update for changes in mesh
        bool haveIndirectPatches_;
        bool meshChanged_;
        bool meshMoved_;

        columnatedFileWriter fileWriter_;


    // Private Member Functions
        virtual void writeFileHeader() const;

        void reduceMaxVectorMag(vector& v);

        tmp<surfaceScalarField> correctedFlux();
        void totalFlux();
        void relativeFlux();
        void backFlowArea();

        template<class GeoField>
        void calcNonScalarData
            (const GeoField& f, FixedList<scalar, 5>&);

        void calcScalarData
            (const volScalarField& f, FixedList<scalar, 5>&);

        void calcSwirlData
            (const volVectorField& f);

        void calcTumble
            (const volVectorField& U);

        //- Return rho if rhoName is specified otherwise rhoRef
        tmp<volScalarField> getRho() const;

        void markFaceZone(const dictionary&);
        void markCellZonePair(const dictionary&);
        void markPatchFaces(const dictionary&);
        void markPlaneFaces(const dictionary&);

        //Update faces to report on
        void updateFaces(const dictionary&);

        //Mark faces to report on
        void markFaces(const dictionary&);

        //Unmark duplicate processor and constrained faces
        void unmarkFaces(const dictionary&);

        void setDirection();

        void calcWeighting();

        //- Disallow default bitwise copy construct and assignment
        surfaceReport(const surfaceReport&);
        void operator=(const surfaceReport&);


public:

    //- Runtime type information
    TypeName("surfaceReport");


    // Constructors

        //- Construct for given objectRegistry and dictionary
        surfaceReport
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    // Destructor

        virtual ~surfaceReport();


    // Member Functions

        //- Update solver
        void calculate();

        //- Execute, currently does nothing
        virtual bool execute();

        //- Sample and write
        virtual bool write();

        //- Read the parameters
        virtual bool read(const dictionary&);

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh& mpm);

        //- Update for changes of mesh
        virtual void movePoints(const polyMesh& mesh);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
