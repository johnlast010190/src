class pyramid;

class faceTriangle
{
    friend class ::pyramid;
    const meshGeometry& geom_;
    const fvMesh& mesh_;
    const pointField& p_;
    const labelList& fPoints_;
    label nPoints_;
    const vectorField& fCtrs_;
    const List<tensor>& faceCentreDerivatives_;
    const label& initialLabel_;
    const label& faceLabel_;
    label currentLabel_;
    label nextLabel_;
    label previousLabel_;
    vector vectorA_;
    vector vectorB_;
    tensor derVectA_;
    tensor derVectB_;

protected:
    tensor enhance(const tensor& A)
    {
        vector A_Enh_1(1.0-A.xx(), -A.xy(), -A.xz());
        vector A_Enh_2(-A.yx(), 1.0-A.yy(), -A.yz());
        vector A_Enh_3(-A.zx(), -A.zy(), 1.0-A.zz());
        tensor B(A_Enh_1,A_Enh_2,A_Enh_3);
        return B;
    }

    void initialize()
    {
        vectorA_ = p_[fPoints_[currentLabel_]] - fCtrs_[faceLabel_];
        vectorB_ = p_[fPoints_[nextLabel_]] - fCtrs_[faceLabel_];
        derVectB_ = faceCentreDerivatives_[faceLabel_];
        derVectA_ = faceCentreDerivatives_[faceLabel_];
    }

    vector squareCrossProductDerivative()
    {
        const vector& a = vectorA_;
        const vector& b = vectorB_;

        scalar coef = 0.;
        vector calc = vector::zero;

        scalar squareSumB = b.x()*b.x()+b.y()*b.y()+b.z()*b.z();
        scalar mixedSum = a.x()*b.x()+a.y()*b.y()+a.z()*b.z();

        coef = (a.x()*a.x()+a.y()*a.y()+a.z()*a.z())*
               squareSumB-
               mixedSum*mixedSum;

        if (coef>VSMALL)
        {
            coef = 0.5/Foam::sqrt(coef);
        }
        else
        {
            return vector::zero;
        }

        calc.x() = (squareSumB*a.x() - mixedSum*b.x());

        calc.y() = (squareSumB*a.y() - mixedSum*b.y());

        calc.z() = (squareSumB*a.z() - mixedSum*b.z());

        return coef*calc;
    }

public:

    faceTriangle();

    faceTriangle
    (
        const meshGeometry& geom,
        const label& faceLabel,
        const label& pointLabel
    )
    :
        geom_(geom),
        mesh_(geom_.getMesh()),
        p_(geom_.getPoints()),
        fPoints_(mesh_.faces()[faceLabel]),
        nPoints_(fPoints_.size()),
        fCtrs_(geom_.fCtrs()),
        faceCentreDerivatives_(geom_.faceCentreDerivatives()),
        initialLabel_(pointLabel),
        faceLabel_(faceLabel),
        currentLabel_(pointLabel)
    {
        nextLabel_ = (currentLabel_ + 1) % nPoints_;
        previousLabel_ = (nextLabel_ + nPoints_ - 2) % nPoints_;
    }

    ~faceTriangle()
    {}

    void next()
    {
        currentLabel_ = nextLabel_;
        nextLabel_ = (nextLabel_ + 1) % nPoints_;
    }

    vector contribution()
    {
        vector nxtTri = nextTriangleContribution();
        vector prvTri = prevTriangleContribution();
        return nxtTri + prvTri;
    }

    vector nextTriangleContribution()
    {
        vectorA_ = p_[fPoints_[currentLabel_]] - fCtrs_[faceLabel_];
        vectorB_ = p_[fPoints_[nextLabel_]] - fCtrs_[faceLabel_];
        vector deriv = squareCrossProductDerivative();
        return deriv;
    }

    vector prevTriangleContribution()
    {
        vectorA_ = p_[fPoints_[currentLabel_]] - fCtrs_[faceLabel_];
        vectorB_ = p_[fPoints_[previousLabel_]] - fCtrs_[faceLabel_];
        vector deriv = squareCrossProductDerivative();
        return deriv;
    }

    bool isANewPoint()
    {
        if (nextLabel_ == initialLabel_)
        {
            return false;
        }
        return true;
    }
};

class pyramid
{
    const meshGeometry& geom_;
    const cellMetric& metric_;
    const pointField& p_;
    const faceTriangle& tri;
    const label& cellT_;
    const label& faceLabel_;
    const labelList& fPoints_;
    const vector& faceCtr_;
    const labelList& own;
    const labelList& nei;
    label owner_;
    label neigh_;
    label bound_;
    vector cellCtr_;
    vector vectorA_;
    vector vectorB_;
    vector vectorC_;
    vector derVecA_;
    label cL_;

    void initialize()
    {
        const tensor& trans = metric_.tTensor()[cL_];
        if (cellT_ == owner_ || cellT_ == bound_)
        {
            vectorA_ = trans&(faceCtr_ - cellCtr_);
            const tensor& t = tri.faceCentreDerivatives_[faceLabel_];
            derVecA_.x() = t.xx();
            derVecA_.y() = t.yy();
            derVecA_.z() = t.zz();
        }
        else
        {
            vectorA_ = trans&(cellCtr_ - faceCtr_);
            const tensor& t = tri.faceCentreDerivatives_[faceLabel_];
            derVecA_.x() = -t.xx();
            derVecA_.y() = -t.yy();
            derVecA_.z() = -t.zz();
        }
    }

public:

    pyramid
    (
        const meshGeometry& geom,
        const faceTriangle& triangle,
        const label& cellType,
        const cellMetric& metric
    )
    :
        geom_(geom),
        metric_(metric),
        p_(geom_.getPoints()),
        tri(triangle),
        cellT_(cellType),
        faceLabel_(tri.faceLabel_),
        fPoints_(geom_.getMesh().faces()[faceLabel_]),
        faceCtr_(geom_.fCtrs()[faceLabel_]),
        own(geom_.getMesh().faceOwner()),
        nei(geom_.getMesh().faceNeighbour()),
        owner_(1),
        neigh_(0),
        bound_(-1),
        cellCtr_(vector::zero),
        cL_(-1)
    {
        if (cellT_ == owner_ || cellT_ == bound_)
        {
            const label& ownl = own[faceLabel_];
            cellCtr_ = geom_.cellCtrs()[ownl];
            cL_ = ownl;
        }
        else if (cellT_ == neigh_)
        {
            const label& neil = nei[faceLabel_];
            cellCtr_ = geom_.cellCtrs()[neil];
            cL_ = neil;
        }
        initialize();
    }

    ~pyramid(){}

    vector contribution()
    {
           const tensor& trans = metric_.tTensor()[cL_];
        if (tri.currentLabel_==tri.initialLabel_)
        {
            vectorB_ = trans&(p_[fPoints_[tri.currentLabel_]] - cellCtr_);
            vectorC_ = trans&(p_[fPoints_[tri.nextLabel_]] - cellCtr_);
            const vector& a = vectorA_;
            const vector& b = vectorB_;
            const vector& c = vectorC_;
            const vector& a_d = derVecA_;
            vector calc = vector::zero;

            calc.x() = a_d.x()*(b.y()*c.z() - b.z()*c.y())
                     + (a.z()*c.y() - a.y()*c.z());

            calc.y() = a_d.y()*(b.z()*c.x() - b.x()*c.z())
                     + (a.x()*c.z() - a.z()*c.x());

            calc.z() = (a.y()*c.x() - a.x()*c.y())
                     + a_d.z()*(b.x()*c.y() - b.y()*c.x());

            return calc;
        }
        else if (tri.nextLabel_==tri.initialLabel_)
        {
            vectorB_ = trans&(p_[fPoints_[tri.currentLabel_]] - cellCtr_);
            vectorC_ = trans&(p_[fPoints_[tri.nextLabel_]] - cellCtr_);
            const vector& a = vectorA_;
            const vector& b = vectorB_;
            const vector& c = vectorC_;
            const vector& a_d = derVecA_;

            vector calc = vector::zero;
            calc.x() = a_d.x()*(b.y()*c.z() - b.z()*c.y())
                     + (a.y()*b.z() - a.z()*b.y());

            calc.y() = a_d.y()*(b.z()*c.x() - b.x()*c.z())
                     + (a.z()*b.x() - a.x()*b.z());

            calc.z() = a_d.z()*(b.x()*c.y() - b.y()*c.x())
                     + (a.x()*b.y() - a.y()*b.x());

            return calc;
        }
        else
        {
            vectorB_ = trans&(p_[fPoints_[tri.currentLabel_]] - cellCtr_);
            vectorC_ = trans&(p_[fPoints_[tri.nextLabel_]] - cellCtr_);
            const vector& b = vectorB_;
            const vector& c = vectorC_;
            const vector& a_d = derVecA_;
            vector calc = vector::zero;

            calc.x() = a_d.x()*(b.y()*c.z() - b.z()*c.y());
            calc.y() = a_d.y()*(b.z()*c.x() - b.x()*c.z());
            calc.z() = a_d.z()*(b.x()*c.y() - b.y()*c.x());
            return calc;
        }
    }
};
