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
    Foam::functionObjects::liftDrag

Group
    grpForcesFunctionObjects

Description
    List of patches and constants to use for lift/drag calculations
    Call liftDrag.write() to calculate&write files.

Usage
    Example of function object specification:
    \verbatim
    liftDrag1
    {
        type            liftDrag;
        functionObjectLibs ( "libforces.so" );
        ...
        patches         ( walls  );
        writeFields     false;
        porosity        false;
        nAveragingSteps 1;
        binData
        {
            xAxis
            {
                nBins           10;
                axis            ( 1 0 0 );
                axisFormat      boundBox; //global; //normalized;
            }
            zAxis
            {
                nBins           10;
                axis            ( 0 0 1 );
                axisFormat      normalized;
            }
        }
        maxCp           1e+15;
        minCp           -1e+15;
        Uinf            ( 50 0 0 );
        rhoInf          1.0;
        liftDirection   ( 0 0 1 );
        dragDirection   ( 1 0 0 );
        referenceArea   2.2;
        referenceLength 3.5;
        referencePoint  ( 0. 0. 0. );
        wheelbase       1;
        Pref            0;
    }
    \endverbatim


SourceFiles
    liftDrag.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_liftDrag_H
#define functionObjects_liftDrag_H

#include "forceCoeffs/forceCoeffs.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "fields/volFields/volFieldsFwd.H"
#include "containers/HashTables/HashPtrTable/HashPtrTable.H"
#include "fields/Fields/vectorField/vectorField.H"
#include "primitives/bools/lists/boolList.H"
#include "meshes/meshShapes/cell/cell.H"
#include "meshes/meshShapes/cell/pyramidPointFaceRef.H"
#include "meshes/meshShapes/face/faceList.H"
#include "meshes/polyMesh/zones/cellZone/cellZone.H"
#include "twoDPointCorrector/twoDPointCorrector.H"
#include "primitives/enums/NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    class forceCoeffs;

/*---------------------------------------------------------------------------*\
                           Class liftDrag Declaration
\*---------------------------------------------------------------------------*/

class liftDrag
:
    public forceCoeffs
{
    public:
        enum axisType
        {
            global,
            boundBox,
            normalized
        };

    private:
        static const NamedEnum<axisType, 3> axisTypeNames_;

    protected:

    // Protected data

        // Read from dictonary

            Switch CofP_refPoint_;
            Switch cOfPressure_;

            //- Sorted list of porous zones
            List<word> porousZoneNames_;

            //- Constants
            vector refCofP_;
            scalar wheelBase_;

            scalar maxCp_;
            scalar minCp_;

            label nAveragingSteps_;

        //- Result values

            scalar totalLift_;
            scalar frontLift_;
            scalar rearLift_;
            scalar drag_;
            scalar side_;
            vector moment_;
            vector CofP_;
            vector CofM_;

            label averagingIndex_;


        //- All region output related information
            class regionData
            :
                public Switch
            {
                public:
                vectorField forceCoeffs_;
                vectorField momentCoeffs_;
                List<autoPtr<OFstream>> filePtr_;

                regionData(bool b):Switch(b){};

            } liftDragRegions_;


        //- All bin output related information
            class binData
            :
                public Switch
            {
                public:
                word   name_;
                vector binAxis_;
                label nBins_;
                scalar binMin_;
                scalar binMax_;
                scalar binSize_;
                vectorField forceCoeffs_;
                vectorField momentCoeffs_;

                List<autoPtr<OFstream>> filePtr_;
                axisType axisFormat_;

                binData(bool b):Switch(b), axisFormat_(boundBox){};
            };

            PtrList<binData> liftDragBinPtr_;

            //- liftDrag file ptr
            autoPtr<OFstream> liftDragFilePtr_;

    // Protected Member Functions

        //- hash bin index
        label binIndex(const binData&, const vector&);

        //- Calculate bins
        void calculateBins();

        //- Disallow default bitwise copy construct and assignment
        liftDrag(const liftDrag&);
        void operator=(const liftDrag&);


public:

    //- Runtime type information
    TypeName("liftDrag");

    // Constructors

        //- Construct for given objectRegistry and dictionary
        liftDrag
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

    // Destructor

        virtual ~liftDrag();


    // Member Functions

        //- Frame bin direction
        vector binAxis(const binData& data) const;

        //- Update solver
        void calculate();

        //- Execute, currently does nothing
        virtual bool execute() override;

        //- Sample and write
        virtual bool write() override;

        //- Read the parameters
        virtual bool read(const dictionary&) override;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
