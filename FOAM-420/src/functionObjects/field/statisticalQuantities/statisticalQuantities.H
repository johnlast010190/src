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
    (c) 2016 OpenCFD Ltd.
    (c) 2015 OpenFOAM Foundation
    (c) 2022 Esi Ltd

Class
    Foam::functionObjects::statisticalQuantities

Group
    grpFieldFunctionObjects

Description
    This function object calculates the average value from the output of
    function objects that generate singular values.

Usage
    Example of function object specification:
    \verbatim
    statisticalQuant
    {
        type        statisticalQuantities;
        libs        ("libfieldFunctionObjects.so");
        ...
        writeToFile yes;
        log         yes;
        functionObject forceCoeffs1;
        fields      (Cm Cd Cl);
        sampleSize  10;
        CoVaV true;
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description             | Required    | Default value
        type         | type name: statisticalQuantities | yes         |
        writeToFile  | write average data to file |  no      | yes
        log          | write average data to standard output | no | yes
        fields       | list of fields to process | yes       |
        functionObject | Name of function object to retrieve data from | yes |
        sampleSize  | Moving Average window | no | 1e+6;

        movAverage | Moving Average  | no | no
        movStdDeviation | StdDev (sigma) | no | no
        movAverage2 | Moving Average switch using x2 data points | no | no
        movStdDevOfMovAverage | StdDev(movAverage) | no | no
        CoV | movStdDeviation/movAverage | no | no
        CoDaV | abs(movAverage2-movAverage)/movAverage | no | no
        CoVaV | movStdDevOfMovAverage/movAverage | no | no
    \endtable

    Output data is written to the file \<timeDir\>/statisticalQuantities.dat

See also
    Foam::functionObject
    Foam::functionObjects::stateFunctionObject
    Foam::functionObjects::writeFile

SourceFiles
    statisticalQuantities.C
    statisticalQuantitiesTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_statisticalQuantities_H
#define functionObjects_statisticalQuantities_H

#include "db/functionObjects/regionFunctionObject/regionFunctionObject.H"
#include "db/functionObjects/writeFile/writeFile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                        Class statisticalQuantities Declaration
\*---------------------------------------------------------------------------*/

class statisticalQuantities
:
    public regionFunctionObject,
    public writeFile
{
    // Private data

        //- Reference to uniform/functionObjects/functionObjectProperties
        //  Needed for restarts
        dictionary& dict_;

        //- Name of function object to retrieve data from
        word functionObjectName_;

        //- List of fields on which to operate
        wordList fieldNames_;

        //- Modes
        //  Here some quantities are needed to compute other quantities
        //  For example if sCoVav is true, then the sAver sDevAver become true
        //  Check at the constructor

        // Moving Average switch
        Switch sAver;
        // StdDev switch
        Switch sDev;
        // Moving Average switch using x2 data points
        Switch sAver2;
        // StdDev of the moving Average switch
        Switch sDevAver;
        // Coefficient of variation
        Switch sCov;
        // Coefficient of delta average
        Switch sCoDav;
        // Coefficient of variation average
        Switch sCoVav;

        //- sample size
        label sampleSize_;

        //- deltaTs sample. It might be needed on the future to add weighted
        //  averaged
        DLList<scalar> runTimeLists_;

        //- Data points stored. Needed for running average
        labelList scalarMap_;
        List<DLList<scalar>> runFieldScalarLists_;
        List<DLList<scalar>> runFieldMeanScalarLists_;
        /*
        labelList vectorMap_;
        labelList sphTensorMap_;
        labelList symmTensorMap_;
        labelList tensorMap_;
        List<DLList<vector>> runFieldVectorLists_;
        List<DLList<sphericalTensor>> runFieldsphTensorLists_;
        List<DLList<symmTensor>> runFieldsymmTensorLists_;
        List<DLList<tensor>> runFieldTensorLists_;
        */

    // Private Member Functions

        //- Read data. Needed for restarts
        dictionary& setDict();

        //- Templated function to resize the lists
        void populateLists();


        //- Templated function to decide if pop is needed
        //  This will depend on the sampleSize and the sAver2
        //  if sAver2 is true we need 2xsampleSize data stored
        template<class Type>
        bool pop(const DLList<Type>&);

        //- Templated function to calculate the statistical quantities
        template<class Type>
        void calc
        (
            List<DLList<Type>>&,
            const labelList&
        );

        //- Templated function to resize the lists
        template<class Type>
        void populate
        (
            List<DLList<Type>>&,
            labelList&
        );

        //- Output file header information
        virtual void writeFileHeader(Ostream& os) const;

        //- Disallow default bitwise copy construct
        statisticalQuantities(const statisticalQuantities&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const statisticalQuantities&) = delete;


public:

    //- Runtime type information
    TypeName("statisticalQuantities");

    //- Constructor
    statisticalQuantities
    (
        const word& name,
        const Time& runTime,
        const dictionary& dict
    );

    //- Destructor
    virtual ~statisticalQuantities();


    // Public Member Functions

        //- Read the field value average data
        virtual bool read(const dictionary&);

        //- Execute
        virtual bool execute();

        //- Write the value average
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "statisticalQuantities/statisticalQuantitiesTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
