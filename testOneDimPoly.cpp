#include <iostream>
#include <math.h>
#include "src/oneDimPoly.cpp"

int main(){
    
    // Test equation: f(x) = 0.5 + 1.5*x + 2.5*x**2 + 3.5*x**3
    std::vector<double> testCoeffsVect(4);
    testCoeffsVect[0] = 0.5;
    testCoeffsVect[1] = 1.5;
    testCoeffsVect[2] = 2.5;
    testCoeffsVect[3] = 3.5;
    OneDimPoly testPoly(testCoeffsVect);

    double testVal = testPoly.evalAt(2.0);
    std::cout << "Test case: f(x) = 0.5 + 1.5*x + 2.5*x**2 + 3.5*x**3" << std::endl; 
    std::cout << "Known: f(2) = " << 0.5+1.5*2.0+2.5*pow(2.0,2.0)+3.5*pow(2.0,3.0)  << std::endl;
    std::cout << "Calculated: f(2) = " << testVal << std::endl << std::endl;

    OneDimPoly testPolyDeriv = testPoly.getDerivative();
    double derivTestVal = testPolyDeriv.evalAt(2.0);
    std::cout << "Test case: f'(x) = 1.5 + 5.0*x + 10.5*x**2" << std::endl; 
    std::cout << "Known: f'(2) = " << 1.5+5.0*2.0+10.5*pow(2.0,2.0)  << std::endl;
    std::cout << "Calculated: f'(2) = " << derivTestVal << std::endl << std::endl;

    OneDimPoly testPolyMult = testPoly.multOneDimPoly(testPolyDeriv);
    double multTestVal = testPolyMult.evalAt(2.0);
    std::cout << "Test case: [f*f'](x) = 0.75 + 4.75*x + 16.5*x**2 + 33.5*x**3 + 43.75*x**4 + 36.75*x**5" << std::endl;
    std::cout << "Known: [f*f'](2.0) = " <<  0.75+4.75*2.0+16.5*pow(2.0,2.0)+33.5*pow(2.0,3.0)+43.75*pow(2.0,4.0)+36.75*pow(2.0,5.0) << std::endl;
    std::cout << "Calculated: [f*f'](2.0) = " << multTestVal << std::endl << std::endl; 

    OneDimPoly testPolyAdd = testPoly.addOneDimPoly(testPolyDeriv);
    double addTestVal = testPolyAdd.evalAt(2.0);
    std::cout << "Test case: f(x) + f'(x) = 2.0 + 6.5*x + 13.0*x**2 + 3.5*x**3" << std::endl;
    std::cout << "Known: f(2.0) + f'(2.0) = " << 2.0+6.5*2.0+13.0*pow(2.0,2.0)+3.5*pow(2.0,3.0) << std::endl;
    std::cout << "Calculated: f(2.0) + f'(2.0) = " << addTestVal << std::endl << std::endl;

    OneDimPoly testPolySub = testPoly.subOneDimPoly(testPolyDeriv);
    double subTestVal = testPolySub.evalAt(2.0);
    std::cout << "Test case: f(x) - f'(x) = -1.0 - 3.5*x - 8.0*x**2 + 3.5*x**3" << std::endl;
    std::cout << "Known: f(2.0) - f'(2.0) = " << -1.0 - 3.5*2.0-8.0*pow(2.0,2.0)+3.5*pow(2.0,3.0) << std::endl;
    std::cout << "Calculated: f(2.0) - f'(2.0) = " << subTestVal << std::endl << std::endl;

    OneDimPoly testPolyScalarMult = testPoly.multScalar(3.0);
    double scalMultTestVal = testPolyScalarMult.evalAt(2.0);
    std::cout << "Test case: 3*f(x) = 1.5 + 4.5*x + 7.5*x**2 + 10.5*x**3" << std::endl; 
    std::cout << "Known: 3*f(2) = " << 1.5+4.5*2.0+7.5*pow(2.0,2.0)+10.5*pow(2.0,3.0)  << std::endl;
    std::cout << "Calculated: 3*f(2) = " << scalMultTestVal << std::endl << std::endl;

    OneDimPoly testPolyScalarDiv = testPoly.divScalar(2.0);
    double scalDivTestVal = testPolyScalarDiv.evalAt(2.0);
    std::cout << "Test case: f(x)/2 = 0.25 + 0.75*x + 1.25*x**2 + 1.75*x**3" << std::endl; 
    std::cout << "Known: f(2)/2 = " << 0.25+0.75*2.0+1.25*pow(2.0,2.0)+1.75*pow(2.0,3.0)  << std::endl;
    std::cout << "Calculated: f(2)/2 = " << scalDivTestVal << std::endl << std::endl;

    OneDimPoly testPolyScalarAdd = testPoly.addScalar(3.0);
    double scalAddTestVal = testPolyScalarAdd.evalAt(2.0);
    std::cout << "Test case: f(x) + 3.0 = 3.5 + 1.5*x + 2.5*x**2 + 3.5*x**3" << std::endl; 
    std::cout << "Known: f(2) + 3.0 = " << 3.5+1.5*2.0+2.5*pow(2.0,2.0)+3.5*pow(2.0,3.0)  << std::endl;
    std::cout << "Calculated: f(2) + 3.0 = " << scalAddTestVal << std::endl << std::endl;

    OneDimPoly testPolyScalarSub = testPoly.subScalar(3.0);
    double scalSubTestVal = testPolyScalarSub.evalAt(2.0);
    std::cout << "Test case: f(x) - 3.0 = -2.5 + 1.5*x + 2.5*x**2 + 3.5*x**3" << std::endl; 
    std::cout << "Known: f(2) - 3.0 = " << -2.5+1.5*2.0+2.5*pow(2.0,2.0)+3.5*pow(2.0,3.0)  << std::endl;
    std::cout << "Calculated: f(2) - 3.0 = " << scalSubTestVal << std::endl << std::endl;

    double integralTestVal = testPoly.findIntegralOverRange(4.0, 5.0);
    std::cout << "Test case: integrate(f(x),4.0,5.0) = 0.5*(5.0-4.0) + 0.75*(5.0**2-4.0**2) + 0.8333*(5.0**3-4.0**3) + 0.875*(5.0**4-4.0**4)"  << std::endl;
    std::cout << "Known: " << 0.5*(5.0-4.0)+0.75*(pow(5.0,2.0)-pow(4.0,2.0))+2.5/3.0*(pow(5.0,3.0)-pow(4.0,3.0))+0.875*(pow(5.0,4.0)-pow(4.0,4.0)) << std::endl;
    std::cout << "Calculated: " << integralTestVal << std::endl << std::endl;
};
