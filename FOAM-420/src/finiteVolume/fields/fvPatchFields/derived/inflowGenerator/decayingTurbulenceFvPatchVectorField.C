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
    (c) 2014 Hannes Kroeger
    (c) 2017-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/inflowGenerator/decayingTurbulenceFvPatchVectorField.H"
#include "fields/volFields/volFields.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/fvPatchFields/fvPatchField/fieldMappers/gibFvPatchFieldMapper.H"
#include "containers/Lists/ListListOps/ListListOps.H"
#include "db/IOstreams/Pstreams/PstreamReduceOps.H"

#include <ctime>

#if !defined( WIN32 ) && !defined( WIN64 )
Foam::random11 ranGen(Foam::label(time(0)));
#else
Foam::Random ranGen(Foam::label(time(0)));
#endif

const Foam::scalar OVERLAP = 0.3;

const Foam::label NUM_VORT = 10;

inline Foam::scalar lspot(Foam::scalar l)
{
    return 3*l;
}

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::decayingTurbulenceFvPatchVectorField::initialiseFromDictionary
(
    const dictionary& dict
)
{
    if (dict.found("vortons"))
    {
        vortons_ = SLList<decayingVorton>(dict.lookup("vortons"));
    }

    // safety check on Reynolds stress input
    forAll(R_, I)
    {
        if
        (
            (R_[I].xx() < 0)
         || (R_[I].yy() < 0)
         || (R_[I].zz() < 0)
        )
        {
            FatalErrorInFunction
                << "Some RMS are negative"<<abort(FatalError);
        }
    }

    // generate Lund tensor based on mean Reynolds stress
    {
        Lund_.replace
        (
            tensor::XX,
            sqrt(R_.component(symmTensor::XX))
        );
        Lund_.replace
        (
            tensor::YX,
            R_.component(symmTensor::XY)/Lund_.component(tensor::XX)
        );
        Lund_.replace
        (
            tensor::ZX,
            R_.component(symmTensor::XZ)/Lund_.component(tensor::XX)
        );
        Lund_.replace
        (
            tensor::YY,
            sqrt(R_.component(symmTensor::YY)-sqr(Lund_.component(tensor::YX)))
        );
        Lund_.replace
        (
            tensor::ZY,
            (R_.component(symmTensor::YZ) - Lund_.component(tensor::YX)
            *Lund_.component(tensor::ZX) )/Lund_.component(tensor::YY)
        );
        Lund_.replace
        (
            tensor::ZZ,
            sqrt(R_.component(symmTensor::ZZ)
            - sqr(Lund_.component(tensor::ZX))-sqr(Lund_.component(tensor::ZY)))
        );
    }

    // Read or set instantaneous Reynolds stress tensor
    if (dict.found("R_inst"))
    {
        Rinst_ = Field<symmTensor>("R_inst", dict, this->size());
    }
    else
    {
        Rinst_ = R_;
    }

    // relaxation coefficient
    if (dict.found("ind"))
    {
        ind_ = readScalar(dict.lookup("ind"));
    }

    direction_ = dict.lookupOrDefault<label>("direction", 1);

    // vorton cut-off size
    if (dict.found("minVortonLength"))
    {
        minVortonLength_ = readScalar(dict.lookup("minVortonLength"));
    }

    const scalarField& sf = this->patch().magSf();

    int numsf = 0;

    forAll(sf, I)
    {
        if (sf[I] > L_[I]*L_[I])
        {
            numsf++;
        }
    }

    if (numsf > 0)
    {
        Pout<< "Warning: there are " << numsf <<" inlet cells (out of "
            << patch().size()
            << ") where the integral length is smaller than the cell size."
            << endl;
    }

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::decayingTurbulenceFvPatchVectorField::decayingTurbulenceFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    L_      (p.size(), 0.0),
    U_    (p.size(), vector::zero),
    R_      (p.size(), symmTensor::zero),
    curTimeIndex_(-1),
    vortons_     (),
    Lund_        (p.size(), tensor::zero),
    Rinst_           (p.size(), symmTensor::zero),
    ind_         (2),
    direction_   (1),
    minVortonLength_(0)
{
}

Foam::decayingTurbulenceFvPatchVectorField::decayingTurbulenceFvPatchVectorField
(
    const decayingTurbulenceFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    L_(mapper(ptf.L_)),
    U_(mapper(ptf.U_)),
    R_(mapper(ptf.R_)),
    curTimeIndex_(-1),
    vortons_(ptf.vortons_),
    Lund_(mapper(ptf.Lund_)),
    Rinst_(mapper(ptf.Rinst_)),
    ind_(ptf.ind_),
    direction_(ptf.direction_),
    minVortonLength_(ptf.minVortonLength_)
{}

Foam::decayingTurbulenceFvPatchVectorField::decayingTurbulenceFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF),
    L_      ("L", dict, p.size()),
    U_    ("U", dict, p.size()),
    R_      ("R", dict, p.size()),
    curTimeIndex_(-1),
    Lund_        (p.size(), pTraits<tensor>::zero),
    Rinst_           (R_),
    ind_         (2),
    direction_   (readScalar(dict.lookup("direction"))),
    minVortonLength_(dict.lookupOrDefault<scalar>("minVortonLength", 0.0))
{
    if (dict.found("vortons"))
        vortons_ = SLList<decayingVorton>(dict.lookup("vortons"));

    if (dict.found("value"))
        fixedValueFvPatchField<vector>::forceAssign(Field<vector>("value", dict, p.size()));
    else
        fixedValueFvPatchField<vector>::forceAssign(U_);

    forAll(R_, I)
    {
        if ((R_[I].xx() < 0) || (R_[I].yy() < 0) || (R_[I].zz() < 0))
        FatalErrorInFunction << "Some RMS are negative"<<abort(FatalError);
    }

    // generate Lund tensor based on mean Reynolds stress
    {
        Lund_.replace
        (
            tensor::XX,
            sqrt(R_.component(symmTensor::XX))
        );
        Lund_.replace
        (
            tensor::YX,
            R_.component(symmTensor::XY)/Lund_.component(tensor::XX)
        );
        Lund_.replace
        (
            tensor::ZX,
            R_.component(symmTensor::XZ)/Lund_.component(tensor::XX)
        );
        Lund_.replace
        (
            tensor::YY,
            sqrt(R_.component(symmTensor::YY)-sqr(Lund_.component(tensor::YX)))
        );
        Lund_.replace
        (
            tensor::ZY,
            (R_.component(symmTensor::YZ) - Lund_.component(tensor::YX)
            *Lund_.component(tensor::ZX) )/Lund_.component(tensor::YY)
        );
        Lund_.replace
        (
            tensor::ZZ,
            sqrt(R_.component(symmTensor::ZZ)
            - sqr(Lund_.component(tensor::ZX))-sqr(Lund_.component(tensor::ZY)))
        );
    }

    if (dict.found("R"))
    {
        Rinst_ = Field<symmTensor>("R", dict, p.size());
    }

    if (dict.found("ind"))
    {
        ind_ = readScalar(dict.lookup("ind"));
    }

    const scalarField& sf = patch().magSf();

    int numsf = 0;

    forAll(sf, I)
    {
        if (sf[I] > L_[I]*L_[I])
        numsf++;
    }

    if (numsf > 0)
    {
        Pout<<"Warning: there are "<<numsf<<" inlet cells (out of "
            <<patch().size()
            <<") where the integral length is smaller than the cell size."<<endl;
    }
}

Foam::decayingTurbulenceFvPatchVectorField::decayingTurbulenceFvPatchVectorField
(
    const decayingTurbulenceFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    L_      (ptf.L_),
    U_    (ptf.U_),
    R_      (ptf.R_),
    curTimeIndex_(-1),
    vortons_     (ptf.vortons_),
    Lund_        (ptf.Lund_),
    Rinst_           (ptf.Rinst_),
    ind_         (ptf.ind_),
    direction_   (ptf.direction_),
    minVortonLength_(ptf.minVortonLength_)
{
}

Foam::decayingTurbulenceFvPatchVectorField::decayingTurbulenceFvPatchVectorField
(
    const decayingTurbulenceFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    L_      (ptf.L_),
    U_    (ptf.U_),
    R_      (ptf.R_),
    curTimeIndex_(-1),
    vortons_     (ptf.vortons_),
    Lund_        (ptf.Lund_),
    Rinst_           (ptf.Rinst_),
    ind_         (ptf.ind_),
    direction_   (ptf.direction_),
    minVortonLength_(ptf.minVortonLength_)
{}

void Foam::decayingTurbulenceFvPatchVectorField::autoMap(const fvPatchFieldMapper& m)
{
    m(*this, *this);
    m(L_, L_);
    m(U_, U_);
    m(R_, R_);
    m(Rinst_, Rinst_);
}

void Foam::decayingTurbulenceFvPatchVectorField::rmap(const fvPatchField<vector>& ptf, const labelList& addr)
{
    fixedValueFvPatchField<vector>::rmap(ptf, addr);

    const decayingTurbulenceFvPatchVectorField& tiptf = refCast<const decayingTurbulenceFvPatchVectorField>(ptf);

    L_.rmap(tiptf.L_, addr);
    U_.rmap(tiptf.U_, addr);
    R_.rmap(tiptf.R_, addr);
    Rinst_.rmap(tiptf.Rinst_, addr);
}


void Foam::decayingTurbulenceFvPatchVectorField::autoMapGIB
(
    const gibFvPatchFieldMapper& mapper
)
{
    fixedValueFvPatchField<vector>::autoMapGIB(mapper);
    mapper.map(L_, pTraits<scalar>::zero);
    mapper.map(U_, pTraits<vector>::zero);
    mapper.map(R_, pTraits<symmTensor>::zero);
    mapper.map(Rinst_, pTraits<symmTensor>::zero);
}


void Foam::decayingTurbulenceFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
        return;

    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        doUpdate();
        curTimeIndex_ = this->db().time().timeIndex();
    }

    fixedValueFvPatchField<vector>::updateCoeffs();
}

void Foam::decayingTurbulenceFvPatchVectorField::doUpdate()
{
    Field<vector>& patchField = *this;

    Field<vector> tloc = patch().Cf();

    List<List<scalar>> l(Pstream::nProcs());
    List<List<vector>> cf(Pstream::nProcs());
    List<List<vector>> rf(Pstream::nProcs());

    l[Pstream::myProcNo()].setSize(L_.size());
    cf[Pstream::myProcNo()].setSize(tloc.size());
    rf[Pstream::myProcNo()].setSize(U_.size());

    label i = 0;
    forAll(tloc, I)
    {
        l[Pstream::myProcNo()][i]  = L_[I];
        cf[Pstream::myProcNo()][i] = tloc[I];
        rf[Pstream::myProcNo()][i] = U_[I];

        ++i;
    }

    Pstream::gatherList(l);
    Pstream::gatherList(cf);
    Pstream::gatherList(rf);

    List<scalar> L  = ListListOps::combine<List<scalar>>(l, accessOp<List<scalar>>());
    List<vector> CF = ListListOps::combine<List<vector>>(cf, accessOp<List<vector>>());
    List<vector> RF = ListListOps::combine<List<vector>>(rf, accessOp<List<vector>>());

    List<List<decayingVorton>> v(Pstream::nProcs());

    if (Pstream::master())
    {
        Pout<<"Starting inflgen, time = "<<this->db().time().elapsedClockTime()<<" s"<<endl;

        forAll(CF, I)
        {
            if (L[I] > minVortonLength_)
            {
                scalar x    = CF[I].x();
                scalar ls   = lspot(L[I]);
                scalar ymin = CF[I].y() - L[I];
                scalar zmin = CF[I].z() - L[I];

                for (label k = 0; k < NUM_VORT; k++)
                {
                    vector v((direction_ > 0) ? x-ls : x+ls, ymin+ranGen.scalar01()*2*L[I], zmin+ranGen.scalar01()*2*L[I]);
                    bool allowed = true;
                    for (SLList<decayingVorton>::const_iterator it = vortons_.begin(); it != vortons_.end(); ++it)
                    {
                        if (mag(v - it().location()) < OVERLAP*ls)
                        {
                            allowed = false;
                            break;
                        }
                    }
                    if (!allowed)
                    {
                        continue;
                    }
                    else
                    {
                        vortons_.insert(decayingVorton(L[I], v, RF[I], (direction_ > 0) ? x+ls : x-ls));
                    }
                }
            }
        }

        Pout<< "The number of vortons is " << vortons_.size() << endl;

        v[Pstream::myProcNo()].setSize(vortons_.size());
        i = 0;
        for (SLList<decayingVorton>::iterator it = vortons_.begin(); it != vortons_.end(); ++it)
            v[Pstream::myProcNo()][i++] = it();
    }

    Pstream::scatterList(v);

    List<decayingVorton> V
        = ListListOps::combine<List<decayingVorton>>(v, accessOp<List<decayingVorton>>());

    if (size() > 0)
    {
        Field<vector> turbulent(U_.size(), pTraits<vector>::zero);

        forAll(patchField, I)
        {
            const vector& pt = tloc[I];

            forAll(V, vI)
            {
                const decayingVorton& vorton = V[vI];

                if (magSqr(vorton.location() - pt) < sqr(lspot(vorton.length())))
                    turbulent[I] += vorton.velocityAt(pt);
            }
        }

        Rinst_ = ((ind_ - 1.0)/ind_)*Rinst_ + (1.0/ind_)*sqr(turbulent);

        Field<tensor> C_(Rinst_.size(), pTraits<tensor>::zero);
        forAll(C_, I)
        {
            C_[I].xx() = 1.0;
            C_[I].yy() = 1.0;
            C_[I].zz() = 1.0;
        }

        forAll(Rinst_, I)
        {
            scalar D1 = Rinst_[I].xx();
            if (D1 > 0)
                C_[I].xx() = 1.0/sqrt(D1);

            scalar D2 = Rinst_[I].xx()*Rinst_[I].yy() - Rinst_[I].xy()*Rinst_[I].xy();
            if (D1 > 0 && D2 > 0)
            {
                C_[I].yx() = -Rinst_[I].xy()/sqrt(D1*D2);
                C_[I].yy() = sqrt(D1/D2);
            }

            scalar D3 = det(Rinst_[I]);
            if (D2 > 0 && D3 > 0)
            {
                C_[I].zx() = (Rinst_[I].xy()*Rinst_[I].yz()-Rinst_[I].yy()*Rinst_[I].xz())/sqrt(D2*D3);
                C_[I].zy() = -(Rinst_[I].xx()*Rinst_[I].yz()-Rinst_[I].xz()*Rinst_[I].xy())/sqrt(D2*D3);
                C_[I].zz() = sqrt(D2/D3);
            }
        }

        ind_++;

        turbulent = C_ & turbulent;
        turbulent = Lund_ & turbulent;

        fixedValueFvPatchField<vector>::forceAssign(U_ + turbulent);
    }

    if (Pstream::master())
    {
        for
        (
            SLList<decayingVorton>::iterator it = vortons_.begin();
            it != vortons_.end();
            ++it
        )
        {
            it().move(this->db().time().deltaT().value());
        }

        for (SLList<decayingVorton>::iterator it = vortons_.begin(); it!=vortons_.end();)
        {
            if (direction_ > 0)
            {
                if (it().location().x() > it().xmax())
                {
                    SLList<decayingVorton>::iterator cl(it);
                    ++it;
                    vortons_.remove(cl);
                }
                else
                {
                    ++it;
                }
            }
            else
            {
                if (it().location().x() < it().xmax())
                {
                    SLList<decayingVorton>::iterator cl(it);
                    ++it;
                    vortons_.remove(cl);
                }
            }
        }

        Pout<<"Finishing inflow generation, time = "<<this->db().time().elapsedClockTime()<<" s"<<endl;
    }
}

void Foam::decayingTurbulenceFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    L_.writeEntry("L", os);
    U_.writeEntry("U", os);
    R_.writeEntry("R", os);
    os.writeEntry("direction", direction_);
    os.writeEntry("minVortonLength", minVortonLength_);
    this->writeEntry("value", os);
    Rinst_.writeEntry("Rinst", os);
    os.writeEntry("ind", ind_);

    if (Pstream::master())
    {
        os.writeEntry("vortons", vortons_);
    }
}


namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        decayingTurbulenceFvPatchVectorField
    );
}
