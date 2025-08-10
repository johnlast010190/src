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
    (c) 2014 Hannes Kroeger <hannes@kroegeronline.net>
 *
 */

#include "twoPointCorrelation/twoPointCorrelation.H"
#include "db/dictionary/dictionary.H"
#include "db/dynamicLibrary/dlLibraryTable/dlLibraryTable.H"
#include "db/Time/Time.H"
#include "interpolation/interpolation/interpolation/interpolation.H"

#include "containers/Lists/SortableList/SortableList.H"
#include "containers/Lists/List/List.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(twoPointCorrelation, 0);
    addToRunTimeSelectionTable(functionObject, twoPointCorrelation, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::twoPointCorrelation::twoPointCorrelation
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name, typeName, dict),
    p0_(point::zero),
    directionSpan_(vector::zero),
    np_(0),
    homogeneousTranslationUnit_(vector::zero),
    nph_(0),
    totalTime_(0.0)
{
    searchEngine_.reset(new meshSearch(mesh_));

    read(dict);
    writeFileHeader(file());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::twoPointCorrelation::~twoPointCorrelation()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::functionObjects::twoPointCorrelation::writeFileHeader(Ostream& os) const
{
    writeCommented(os, "Time");
    writeDelimited(os, "correlation values");
}

void Foam::functionObjects::twoPointCorrelation::resetAveraging()
{
  totalTime_=0.0;
  filePtr_.reset();
  correlationCoeffs_.reset(new tensorField(np_, tensor::zero));
}


template<class T>
void Foam::functionObjects::twoPointCorrelation::combineSampledValues
(
    const volFieldSampler<T>& sampledField,
    const labelListList& indexSets,
    autoPtr<volFieldSampler<T>>& masterField
)
{
    List<Field<T>> masterValues(indexSets.size());

    forAll(indexSets, setI)
    {
        // Collect data from all processors
        List<Field<T>> gatheredData(Pstream::nProcs());
        gatheredData[Pstream::myProcNo()] = sampledField[setI];
        Pstream::gatherList(gatheredData);

        if (Pstream::master())
        {
            Field<T> allData
            (
                ListListOps::combine<Field<T>>
                (
                    gatheredData,
                    Foam::accessOp<Field<T>>()
                )
            );

            masterValues[setI] = UIndirectList<T>
                                 (
                                     allData,
                                     indexSets[setI]
                                 )();
        }
    }

    masterField.reset
    (
        new volFieldSampler<T>
        (
            masterValues,
            sampledField.name()
        )
    );
}

template<class Type>
Foam::autoPtr<Foam::volFieldSampler<Type>>
Foam::functionObjects::twoPointCorrelation::sample
(
    PtrList<sampledSet>& sets,
    const GeometricField<Type, fvPatchField, volMesh>& field,
    const labelListList& indexSets
)
{
    // Storage for interpolated values
    volFieldSampler<Type> sampledField
    (
        "cellPointFace",
        field,
        sets
    );
    autoPtr<volFieldSampler<Type>> masterField;
    combineSampledValues(sampledField, indexSets, masterField);

    return masterField;
}

void Foam::functionObjects::twoPointCorrelation::combineSampledSets
(
    PtrList<coordSet>& masterSampledSets,
    labelListList& indexSets
)
{
    // Combine sampleSets from processors. Sort by curveDist. Return
    // ordering in indexSets.
    // Note: only master results are valid

    masterSampledSets_.clear();
    masterSampledSets_.setSize(lines_.size());
    indexSets_.setSize(lines_.size());

    const PtrList<sampledSet>& sampledSets = lines_;

    forAll(sampledSets, setI)
    {
        const sampledSet& samplePts = sampledSets[setI];

        // Collect data from all processors
        List<List<point>> gatheredPts(Pstream::nProcs());
        gatheredPts[Pstream::myProcNo()] = samplePts;
        Pstream::gatherList(gatheredPts);

        List<labelList> gatheredSegments(Pstream::nProcs());
        gatheredSegments[Pstream::myProcNo()] = samplePts.segments();
        Pstream::gatherList(gatheredSegments);

        List<scalarList> gatheredDist(Pstream::nProcs());
        gatheredDist[Pstream::myProcNo()] = samplePts.curveDist();
        Pstream::gatherList(gatheredDist);


        // Combine processor lists into one big list.
        List<point> allPts
        (
            ListListOps::combine<List<point>>
            (
                gatheredPts, accessOp<List<point>>()
            )
        );
        labelList allSegments
        (
            ListListOps::combine<labelList>
            (
                gatheredSegments, accessOp<labelList>()
            )
        );
        scalarList allCurveDist
        (
            ListListOps::combine<scalarList>
            (
                gatheredDist, accessOp<scalarList>()
            )
        );


        if (Pstream::master() && allCurveDist.size() == 0)
        {
            WarningInFunction
                << "Sample set " << samplePts.name()
                << " has zero points." << endl;
        }

        // Sort curveDist and use to fill masterSamplePts
        SortableList<scalar> sortedDist(allCurveDist);
        indexSets[setI] = sortedDist.indices();

        masterSampledSets.set
        (
            setI,
            new coordSet
            (
                samplePts.name(),
                samplePts.axis(),
                List<point>(UIndirectList<point>(allPts, indexSets[setI])),
                sortedDist
            )
        );
    }
}


bool Foam::functionObjects::twoPointCorrelation::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);
    Log << type() << " " << name() <<  " read:" << nl;

    p0_ = point(dict.lookup("p0"));
    directionSpan_=vector(dict.lookup("directionSpan"));
    np_=readLabel(dict.lookup("np"));

    homogeneousTranslationUnit_=
        vector(dict.lookup("homogeneousTranslationUnit"));
    nph_=readLabel(dict.lookup("nph"));

    dictionary csysDict(dict.subDict("csys"));
    csys_=coordinateSystem::New
    (
        obr_,
        csysDict
    );

    Log <<"Definition of twoPointCorrelation "<<name()<<":"<<nl
        <<"    from point "<<p0_<<" on "<<np_<<" points along "
        <<directionSpan_<<nl
        <<"    averaged over "<<nph_<<" copies, translated by "
        <<homogeneousTranslationUnit_<<endl;

    createInterpolators();


    IOobject propsDictHeader
    (
        "twoPointCorrelationProperties",
        obr_.time().timeName(obr_.time().startTime().value()),
        "uniform",
        obr_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );

    if (propsDictHeader.typeHeaderOk<IOdictionary>())
    {
        IOdictionary propsDict(propsDictHeader);

        if (Pstream::master())
        {
            totalTime_= obr_.time().deltaTValue();

            if (propsDict.found(name()))
            {
                Info<< "    Restarting averaging for twoPointCorrelation site "
                    << name() << nl;
                totalTime_ =
                    readScalar(propsDict.subDict(name()).lookup("totalTime"));

                for (label i=0; i<pTraits<tensor>::nComponents; i++)
                {
                    Istream& is=
                        propsDict.subDict(name()).lookup("correlationCoeffs");
                    label num=readLabel(is); // overread list size
                    if (num!=np_)
                        FatalErrorInFunction
                            << "number of sampling points does not match"
                            << Foam::abort(FatalError);
                    correlationCoeffs_.reset
                        (new tensorField(readList<tensor>(is)));
                }
            }
        }
    }

    return true;
}


bool Foam::functionObjects::twoPointCorrelation::execute()
{
    Log << type() << " " << name() <<  " execute:" << nl;

    if (!foundObject<volVectorField>("UMean"))
    {
        //resetAveraging();
        WarningInFunction
            << "No mean velocity field in registry, delaying until "
            << "next timestep"<<endl;
    }
    else
    {
        const volVectorField& U = lookupObject<volVectorField>("U");
        const volVectorField& Umean = lookupObject<volVectorField>("UMean");
        volVectorField uPrime( U-Umean );

        autoPtr<OFstream> dbgFile;
        if (debug && Pstream::master())
        {
            dbgFile.reset(new OFstream("twoPointCorrelation_"+name()+".csv"));
            dbgFile() << "X,Y,Z,Vx,Vy,Vz,Vr,Vtheta,Vz" <<nl;
        }

        combineSampledSets(masterSampledSets_, indexSets_);
        autoPtr<volFieldSampler<vector>> vfs =
            sample(lines_, uPrime, indexSets_);

        if (Pstream::master())
        {
            tensorField cCoeffs(correlationCoeffs_().size(), tensor::zero);
            forAll(vfs(), i)
            {
                //const cloudSet& samples = lines_[i];
                const vectorField& values=vfs()[i];

                //vectorField values(np_, vector::zero);
                // Fixed size according to input params!
                for (label j=0; j<np_; j++)
                {
                    if (dbgFile.valid())
                    {
                        const point& pt = masterSampledSets_[i][j];
                        const vector& v= values[j]; // in local CS
                        // in local CS
                        const vector& lv= csys_().localVector(values[j]);
                        Info<<j<<" "<<pt<<" "<<v<<endl;
                        dbgFile() << pt.x()<<","<<pt.y()<<","<<pt.z()<<","
                                  <<v.x()<<","<<v.y()<<","<<v.z()<<","
                                  <<lv.x()<<","<<lv.y()<<","<<lv.z()<<nl;
                    }
                    cCoeffs[j] += csys_().localVector
                        (values[0]) * csys_().localVector(values[j]);
                }
            }

            if (dbgFile.valid())
            {
                dbgFile.reset();
            }

            // averaging over homogeneous directions
            scalar dt = obr_.time().deltaTValue();
            totalTime_ += dt;
            scalar Dt = totalTime_;
            scalar alpha = (Dt - dt)/Dt;
            scalar beta = dt/Dt;

            correlationCoeffs_() = alpha * correlationCoeffs_()
                + beta*(cCoeffs/scalar(lines_.size()));

        }

        Pstream::scatter(totalTime_);
        Pstream::scatter(correlationCoeffs_());


        if (obr_.time().outputTime())
        {
            IOdictionary propsDict
            (
                IOobject
                (
                    "twoPointCorrelationProperties",
                    obr_.time().timeName(),
                    "uniform",
                    obr_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE,
                    false
                 )
            );

            {
                propsDict.add(name(), dictionary());

                propsDict.subDict(name()).add("totalTime", totalTime_);

                propsDict.subDict(name()).add
                (
                    "correlationCoeffs",
                    correlationCoeffs_()
                );
                propsDict.regIOobject::write();
            }
        }
    }

    if (correlationCoeffs_.valid())
    {
        file() << obr_.time().timeName() <<token::TAB;
        for (label k=0; k < pTraits<tensor>::nComponents; k++)
        {
            for (label i=0; i<correlationCoeffs_().size(); i++)
            {
                file()<<correlationCoeffs_()[i][k]<<token::SPACE;
            }
            file()<<token::TAB;
        }

        tensor L=tensor::zero;
        for (label l=0; l<correlationCoeffs_().size()-1; l++)
        {
            L += 0.5*(correlationCoeffs_()[l]+correlationCoeffs_()[l+1])
                * (x_()[l+1]-x_()[l]);
        }
        if (mag(correlationCoeffs_()[0]) > SMALL)
        {
            L=cmptDivide( L, correlationCoeffs_()[0]);
        }
        for (label k=0; k < pTraits<tensor>::nComponents; k++)
        {
            file()<<L[k]<<token::SPACE;
        }

        file()<<endl;
    }

    return true;
}


bool Foam::functionObjects::twoPointCorrelation::write()
{
    return true;
}


void Foam::functionObjects::twoPointCorrelation::createInterpolators()
{
    Info<< "Building interpolators for twoPointCorrelation "<<name() << endl;

    lines_.resize(nph_);
    x_.reset(new scalarField(np_, 0.0));

    for (label i=0; i<nph_; i++)
    {
        pointField pts(np_);

        for (label j=0; j<np_; j++)
        {
            pts[j]=csys_().globalPosition
            (
                p0_
                + scalar(i) * homogeneousTranslationUnit_
                + scalar(j)/scalar(np_-1)  * directionSpan_
            );
            if ((i==0)&&(j>0))
                x_()[j] = x_()[j-1] + mag(pts[j]-pts[j-1]);
        }

        lines_.set
        (
            i,
            new cloudSet
            (
                name(),
                mesh_,
                searchEngine_(),
                "distance",
                pts
            )
        );
    }

    combineSampledSets(masterSampledSets_, indexSets_);

    bool reset=false;
    if (!correlationCoeffs_.valid())
    {
        reset=true;
    }
    else if (correlationCoeffs_().size()!=np_)
    {
        Info<< "Reset averaging because parameters became incompatible "
             << "to previous averaging."<<endl;
        reset=true;
    }

    if (reset)
    {
        resetAveraging();
    }

}


// ************************************************************************* //
