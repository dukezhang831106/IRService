#include <vector>
#include <set>
#include <iostream>
#include <algorithm>
#include <ql/quantlib.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

using namespace QuantLib;

struct ProbabilityTree{
    std::vector<boost::numeric::ublas::vector<double>> R;
    std::vector<boost::numeric::ublas::vector<double>> Fwd;
    std::vector<boost::numeric::ublas::matrix<int>> Connect;
    std::vector<boost::numeric::ublas::matrix<double>> Prob;
    std::vector<boost::numeric::ublas::vector<int>> Branches;
    Date valDate;
    std::vector<Date> matDates;
    std::vector<double> dT;
    boost::numeric::ublas::vector<double> dr;
    double sigma, alpha, theta;
    ProbabilityTree(int num) : R(std::vector<boost::numeric::ublas::vector<double>>(num, boost::numeric::ublas::vector<double>())), Fwd(std::vector<boost::numeric::ublas::vector<double>>(num, boost::numeric::ublas::vector<double>())), Connect(std::vector<boost::numeric::ublas::matrix<int>>(num-1, boost::numeric::ublas::matrix<int>())), Prob(std::vector<boost::numeric::ublas::matrix<double>>(num-1, boost::numeric::ublas::matrix<double>())) {};
    ProbabilityTree(Date& vDate, std::vector<Date>& mDates, double& vol, double& level, double& speed, std::vector<double>& dtSpan) : valDate(vDate), matDates(mDates), sigma(vol), alpha(speed), theta(level), dT(dtSpan) {
        int num = dtSpan.size() - 1;
        R = std::vector<boost::numeric::ublas::vector<double>>(num, boost::numeric::ublas::vector<double>());
        Fwd = std::vector<boost::numeric::ublas::vector<double>>(num, boost::numeric::ublas::vector<double>());
        Branches = std::vector<boost::numeric::ublas::vector<int>>(num, boost::numeric::ublas::vector<int>());
        Connect = std::vector<boost::numeric::ublas::matrix<int>>(num-1, boost::numeric::ublas::matrix<int>());
        Prob = std::vector<boost::numeric::ublas::matrix<double>>(num-1, boost::numeric::ublas::matrix<double>());
    }
    ProbabilityTree(Date& vDate, std::vector<Date>& mDates, double& vol, double& speed, std::vector<double>& dtSpan) : valDate(vDate), matDates(mDates), sigma(vol), alpha(speed), theta(0.0), dT(dtSpan) {
        int num = dtSpan.size() - 1;
        R = std::vector<boost::numeric::ublas::vector<double>>(num, boost::numeric::ublas::vector<double>());
        Fwd = std::vector<boost::numeric::ublas::vector<double>>(num, boost::numeric::ublas::vector<double>());
        Branches = std::vector<boost::numeric::ublas::vector<int>>(num, boost::numeric::ublas::vector<int>());
        Connect = std::vector<boost::numeric::ublas::matrix<int>>(num-1, boost::numeric::ublas::matrix<int>());
        Prob = std::vector<boost::numeric::ublas::matrix<double>>(num-1, boost::numeric::ublas::matrix<double>());
    }
};

double findb(double& dT, double& X){
    double sqdt = std::sqrt(dT), xsqrtdt = X/std::sqrt(1.5*dT), bcom = X/sqdt;
    double base = xsqrtdt > 0 ? std::floor(xsqrtdt): -std::floor(-xsqrtdt), base2 = xsqrtdt+1 > 0 ? std::floor(xsqrtdt+1): -std::floor(-xsqrtdt+1);
    double be = bcom/base, bc = bcom/base2;
    return std::abs(bc - std::sqrt(1.5)) < std::abs(be - std::sqrt(1.5)) ? bc : be;
}

std::vector<double> solvedT(double& dT, double& X){
    int cnt = 0, Kmax = 100;
    double tol = 1e-10, prev0 = dT, curr0 = dT + tol, fprev0 = findb(prev0, X), fcurr0 = findb(curr0, X), prev1 = prev0, curr1 = curr0, fprev1 = fprev0, fcurr1 = fcurr0;
    std::vector<double> retval(2, 0.0);
    while(cnt < Kmax){
        double error0 = fcurr0 - std::sqrt(1.5), error1 = fcurr1 - 1.0;
        if (std::abs(error0) < tol && std::abs(error1) < tol){
            retval[0] = curr0;
            retval[1] = curr1;
            return retval;
        }
        else{
            if (std::abs(error0) > tol){
                double temp = curr0 - (curr0 - prev0) * error0/(fcurr0 - fprev0);
                prev0 = curr0;
                fprev0 = fcurr0;
                curr0 = temp;
                fcurr0 = findb(curr0, X);
            }
            if (std::abs(error1) > tol){
                double temp = curr1 - (curr1 - prev1) * error1/(fcurr1 - fprev1);
                prev1 = curr1;
                fprev1 = fcurr1;
                curr1 = temp;
                fcurr1 = findb(curr1, X);
            }
        }
    }
    retval[0] = std::numeric_limits<double>::infinity();
    retval[1] = std::numeric_limits<double>::infinity();
    return retval;
}

// remove elements from a vector if value is approximately near (up to 9 decimals)
void remove_duplicates(std::vector<double>& vec){
    int rounder = 1e8;
    std::map<double, double> mapping;
    std::for_each(vec.begin(), vec.end(), [&](double& x) {
        double key = std::round(rounder*x)/rounder;
        if (mapping.find(key) == mapping.end()){
            mapping[key] = x;
        }
            
    });
    vec.resize(mapping.size());
    int cnt = 0;
    for(std::map<double, double>::iterator it = mapping.begin(); it != mapping.end(); it++){
        vec[cnt] = it->first; cnt++;
    }
    std::sort(vec.begin(), vec.end(), std::greater<double>());
}

ProbabilityTree getLattice(){
    //std::vector<Period> maturities = {1*Years, 2*Years, 3*Years, 5*Years, 7*Years, 10*Years, 12*Years, 15*Years, 20*Years};
    Calendar calendar = UnitedStates();
    DayCounter daycounter = Actual365NoLeap();
    Date valDate(28, May, 2019);
    std::vector<Date> volDates = {Date(28, May, 2020), Date(28, May, 2021), Date(31, May, 2022), Date(28, May, 2024), 
        Date(28, May, 2026), Date(29, May, 2029), Date(28, May, 2031), Date(30, May, 2034), Date(31, May, 2039)};
    //std::transform(maturities.begin(), maturities.end(), volDates.begin(), [&](Period& p){ return calendar.advance(valDate, p); });
    Date alphaDate = volDates[volDates.size() - 1];

    double alphaCurve = 0.00568106, thetaCurve = 0.190994, volCurve = 0.049579;
    std::vector<double> timeFracs(volDates.size());
    std::transform(volDates.begin(), volDates.end(), timeFracs.begin(), [&](Date& d) { return daycounter.yearFraction(valDate, d); });
    
    std::vector<double> zRates = {0.0212836, 0.0206989, 0.0205021, 0.0203964, 0.020933, 0.0219579, 0.022595, 0.0232348, 0.0238165};
    std::vector<double> discountFactors(zRates.size()), nominalBond(zRates.size());
    int timeSteps = -1;
    std::transform(zRates.begin(), zRates.end(), discountFactors.begin(), [&](double& r) { timeSteps++; return std::pow(1+r, -timeFracs[timeSteps]); });

    double dT = timeFracs[timeFracs.size() - 1] / timeFracs.size();
    
    // linear interpolation the first grid
    LinearInterpolation Interp (timeFracs.begin(),timeFracs.end(),zRates.begin());

    double refRate = Interp(dT);
    timeSteps = 0;
    // Find Prices of nominal bond
    std::for_each(timeFracs.begin(), timeFracs.end(), [&](double& t) { timeSteps++; double rate = Interp(timeSteps*dT); nominalBond[timeSteps-1] = std::pow(1 + rate, -timeSteps*dT); });
    
    int numLevels = volDates.size();
    std::vector<double> tSpan(numLevels + 1, 0.0), dtSpan(numLevels + 1, 0.0); tSpan[0] = 0.0;
    std::copy(timeFracs.begin(), timeFracs.end(), tSpan.begin()+1);
    std::adjacent_difference(tSpan.begin(), tSpan.end(), dtSpan.begin());
    
    ProbabilityTree ptree(valDate, volDates, volCurve, thetaCurve, alphaCurve, dtSpan);

    double Xprev = 2*sqrt(refRate)/ptree.sigma, sqdt = std::sqrt(dT), b = findb(dT, Xprev);
    std::vector<double> vecXprev(1, Xprev);

    if (b > std::sqrt(2.0) || b < 1.0){
        // find dT for which b will make sense
        std::vector<double> dt = solvedT(dT, Xprev);

    }

    double sigsq = ptree.sigma * ptree.sigma, sigsq4 = sigsq / 4.0;

    ptree.R[0] = boost::numeric::ublas::vector<double>(1, refRate);
    // jump if necessary
    double Xju = 2*std::sqrt(ptree.alpha*ptree.theta*dT)/ptree.sigma;


    for (int level = 1; level < numLevels; level++){
        std::vector<double> Xu(vecXprev.size()), Xm(vecXprev.size()), Xd(vecXprev.size());
        std::vector<int> nomask;
        //std::cout << "Level " << level << std::endl;        
        for (int index = 0; index < vecXprev.size(); index++){
            double muxt = (0.5 * ptree.alpha * (4.0 * ptree.theta / sigsq - vecXprev[index] * vecXprev[index]) - 0.5) / vecXprev[index];
            double J = std::floor(muxt * sqdt / b + 1.0 / b / b);
            if (vecXprev[index] > DBL_EPSILON){    
                Xu[index] = vecXprev[index] + b * (J + 1) * sqdt;
                Xm[index] = vecXprev[index] + b * J * sqdt;
                Xd[index] = std::max(vecXprev[index] + b * (J - 1) * sqdt, 0.0);
                Xd[index] = Xd[index] < 100000 * DBL_EPSILON ? 0.0 : Xd[index];
            }
            else{
                Xu[index] = 0.0;
                Xm[index] = 0.0;
                Xd[index] = 0.0;
                nomask.push_back(index);
            }
            ptree.Prob[level - 1].resize(3, vecXprev.size());
            ptree.Prob[level - 1](0, index) = 0.5/b/b - J/2 + sqdt*muxt/2/b;
            ptree.Prob[level - 1](1, index) = 1 - 1/b/b;
            ptree.Prob[level - 1](2, index) = 0.5/b/b + J/2 - sqdt*muxt/2/b;
        }
        remove_duplicates(Xu);
        remove_duplicates(Xm);
        remove_duplicates(Xd);

        std::vector<double> vecX(Xu);
        vecX.insert(vecX.end(), Xm.begin(), Xm.end());
        vecX.insert(vecX.end(), Xd.begin(), Xd.end());

        remove_duplicates(vecX);

        // update rates

        ptree.R[level].resize(vecX.size());
        std::transform(vecX.begin(), vecX.end(), ptree.R[level].begin(), [&](double& x) { return sigsq4*x*x; });
        int idxJump = -1;
        if (std::abs(ptree.R[level-1][ptree.R[level-1].size()-1]) < DBL_EPSILON){
            std::vector<double>::reverse_iterator low = std::lower_bound(Xu.rbegin(), Xu.rend(), Xju, std::less<double>());
            idxJump = std::distance(low, Xu.rend()) - 1;
        }


        // update Connect
        boost::numeric::ublas::matrix<int> conn(3, Xu.size());
        
        for(int i = 0; i < Xu.size(); i++){
            std::vector<double>::iterator it = std::find(vecX.begin(), vecX.end(), Xu[i]); 
            conn(0, i) = std::distance(vecX.begin(), it);
            it = std::find(vecX.begin(), vecX.end(), Xm[i]); 
            conn(1, i) = std::distance(vecX.begin(), it);
            it = std::find(vecX.begin(), vecX.end(), Xd[i]); 
            conn(2, i) = std::distance(vecX.begin(), it);
            if (idxJump >= 0){
                std::for_each(nomask.begin(), nomask.end(), [&](int& pos) { conn(0, pos) = idxJump;});
            }
            
        }
        ptree.Connect[level-1] = conn;

        if (idxJump >= 0){
            boost::numeric::ublas::vector<double> pu(conn.size2()), temp_u(conn.size2()), temp_d(conn.size2());
            boost::numeric::ublas::vector<int> puconn(boost::numeric::ublas::matrix_row<boost::numeric::ublas::matrix<int>>(conn, 0));
            boost::numeric::ublas::vector<int> pdconn(boost::numeric::ublas::matrix_row<boost::numeric::ublas::matrix<int>>(conn, 2));
            std::transform(ptree.R[level-1].begin(), ptree.R[level-1].end(), pu.begin(), [&](double& r) { return r + dT*ptree.alpha*(ptree.theta - r); });
            std::transform(puconn.begin(), puconn.end(), temp_u.begin(), [&](int& idx) { return ptree.R[level][idx]; });
            std::transform(pdconn.begin(), pdconn.end(), temp_d.begin(), [&](int& idx) { return ptree.R[level][idx]; });
            pu = pu - temp_d;
            temp_u = temp_u - temp_d;
            std::transform(pu.begin(), pu.end(), temp_u.begin(), pu.begin(), std::divides<double>());
            std::for_each(nomask.begin(), nomask.end(), [&](int& pos) { ptree.Prob[level-1](0, pos) = pu(pos); ptree.Prob[level-1](1, pos) = 0.0; ptree.Prob[level-1](2, pos) = 1.0 - pu(pos);});
        }
        
        // update probabilities
        vecXprev = vecX;
    }
    
    // Now calculate forward rates
    ptree.Fwd = ptree.R;
    //std::cout << "Forward: " << std::endl;
    for(int level = 0; level < ptree.R.size(); level++){
        std::transform(ptree.R[level].begin(), ptree.R[level].end(), ptree.Fwd[level].begin(), [&](double& r) { return std::exp(r*dT); });
    }

    // The trees we have are for the initial rate (first time Step). We have to now
    // adjust the tree to the rest of the term structure. We do that by
    // calculating the price of pseudo-bonds with face=1 maturing at each tree level

    std::vector<boost::numeric::ublas::vector<double>> Pt(ptree.Fwd.size() + 1);
    std::vector<boost::numeric::ublas::matrix<double>> ProbT(Pt.size());
    std::copy(ptree.Fwd.begin(), ptree.Fwd.end(), Pt.begin());
    std::copy(ptree.Fwd.end() - 1, ptree.Fwd.end(), Pt.end() - 1);
    
    std::transform(Pt.begin(), Pt.end(), ProbT.begin(), [&](boost::numeric::ublas::vector<double>& row) { return boost::numeric::ublas::matrix<double>(numLevels, row.size(), 0.0); });
    
    boost::numeric::ublas::matrix<double> Face(numLevels , numLevels + 1, 0.0);
    boost::numeric::ublas::subrange(Face, 0, numLevels, 1, numLevels + 1) = boost::numeric::ublas::identity_matrix<double>(numLevels);
   
    int NPLevels = Pt.size(), NBranches = Pt[1].size();
    for(int level = NPLevels - 1; level >= 0; level--){
        if (level < NPLevels - 1){
            for(int i = 0; i < ProbT[level].size1(); i++){
                for(int j = 0; j < ProbT[level].size2(); j++)
                    ProbT[level](i, j) = ProbT[level](i, j)/ptree.Fwd[level][j];
            }
        }
        
        ProbT[level] += boost::numeric::ublas::prod(boost::numeric::ublas::subrange(Face, 0, numLevels, level, level + 1), boost::numeric::ublas::matrix<double>(1, Pt[level].size(), 1.0)); 

        if (level == NPLevels - 1){
            ProbT[level - 1] = ProbT[level];
        }
        else if (level >= 1){
            for(int branch = 0; branch < NBranches; branch++){
                boost::numeric::ublas::vector<int> connectTo(boost::numeric::ublas::matrix_row<boost::numeric::ublas::matrix<int>>(ptree.Connect[level - 1], branch));
                for(int i = 0; i < ProbT[level - 1].size1(); i++){
                   for(int j = 0; j < ProbT[level - 1].size2(); j++){
                        ProbT[level - 1](i, j) += ptree.Prob[level - 1](branch, j) * ProbT[level](i, connectTo[j]);
                   }
                }
            }
        }
    }

    boost::numeric::ublas::vector<double> Pdt(discountFactors.size()), delta(Pdt.size());
    std::transform(ProbT[0].begin1(), ProbT[0].end1(), nominalBond.begin(), Pdt.begin(), std::divides<double>());
    std::transform(Pdt.begin(), Pdt.end(), Pdt.begin(), [&](double& z) { return std::log(z)/dT; });
    
    std::adjacent_difference(Pdt.begin(), Pdt.end(), delta.begin());

    for(int level = 0; level < numLevels; level++){
        std::transform(ptree.R[level].begin(), ptree.R[level].end(), ptree.R[level].begin(), [&](double& r) { return r + delta[level]; });
        std::transform(ptree.R[level].begin(), ptree.R[level].end(), ptree.Fwd[level].begin(), [&](double& r) { return std::exp(r*dT); });
        std::cout << "Level " << level << std::endl;
        std::cout << "forward: " << ptree.Fwd[level] << std::endl;
        if (level < numLevels - 1){
            std::cout << "connection: " << ptree.Connect[level] << std::endl;
            std::cout << "prob: " << ptree.Prob[level] << std::endl;
        }
    }

    
    return ptree;
}

int main(){
    ProbabilityTree ptree = getLattice();
    std::cout << "Done" << std::endl;
    return EXIT_SUCCESS;
}