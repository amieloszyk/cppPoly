#include <iostream>
#include <math.h>
#include "src/twoDimPoly.cpp"

int main() {

    // Test equation: f(x,y) = (1.0 + 2.0*x) + (-0.5 + 4.5*x**2)*y + (10.0*x + 2.5*x**2)*y**2
    std::vector<OneDimPoly> testOneDimPolyVect;
    std::vector<double> varTwoOrderZeroCoeffs(2);
    std::vector<double> varTwoOrderOneCoeffs(3);
    std::vector<double> varTwoOrderTwoCoeffs(3);

    varTwoOrderZeroCoeffs[0] = 1.0;
    varTwoOrderZeroCoeffs[1] = 2.0;

    varTwoOrderOneCoeffs[0] = -0.5;
    varTwoOrderOneCoeffs[1] = 0.0;
    varTwoOrderOneCoeffs[2] = 4.5;

    varTwoOrderTwoCoeffs[0] = 0.0;
    varTwoOrderTwoCoeffs[1] = 10.0;
    varTwoOrderTwoCoeffs[2] = 2.5;

    testOneDimPolyVect.push_back(OneDimPoly(varTwoOrderZeroCoeffs));
    testOneDimPolyVect.push_back(OneDimPoly(varTwoOrderOneCoeffs));
    testOneDimPolyVect.push_back(OneDimPoly(varTwoOrderTwoCoeffs));

    TwoDimPoly testPoly(testOneDimPolyVect);

    double testVal = testPoly.evalAt(2.0, 4.0);
    std::cout << "Test case: f(x,y) = (1.0 + 2.0*x) + (-0.5 + 4.5*x**2)*y + (10.0*x + 2.5*x**2)*y**2" << std::endl;
    std::cout << "Known: f(2, 4) = " << (1.0 + 2.0*2.0) + (-0.5 + 4.5*pow(2.0,2.0))*4.0 + (10.0*2.0 + 2.5*pow(2.0,2.0))*pow(4.0,2.0) << std::endl;
    std::cout << "Calculated: f(2, 4) = " << testVal << std::endl << std::endl; 
    
    
    // Second test equation: f(x,y) = (2.0 + 3.0*x**2) + (0.5*x)*y
    std::vector<OneDimPoly> secondTestOneDimPolyVect;
    std::vector<double> secondVarTwoOrderZero(3); 
    secondVarTwoOrderZero[0] = 2.0;
    secondVarTwoOrderZero[1] =  0.0;
    secondVarTwoOrderZero[2] =  3.0;
    std::vector<double> secondVarTwoOrderOne(2); 
    secondVarTwoOrderOne[0] = 0.0;
    secondVarTwoOrderOne[1] = 0.5;
    secondTestOneDimPolyVect.push_back(secondVarTwoOrderZero);
    secondTestOneDimPolyVect.push_back(secondVarTwoOrderOne);
    TwoDimPoly secondTestPoly(secondTestOneDimPolyVect);
    double secondTestVal = secondTestPoly.evalAt(2.0, 4.0);

    std::cout << "Test case: g(x,y) = (2.0 + 3.0*x**2) + (0.5*x)*y" << std::endl;
    std::cout << "Known: g(2, 4) = " << (2.0 + 3.0*pow(2.0,2.0)) + (0.5*2.0)*4.0 << std::endl;
    std::cout << "Calculated: g(2, 4) = " << secondTestVal << std::endl << std::endl; 

    TwoDimPoly testMultPoly = testPoly.multTwoDimPoly(secondTestPoly);
    double multTestVal = testMultPoly.evalAt(2.0, 4.0);
    std::cout << "Known: f(2,4)*g(2, 4) = " << ((1.0 + 2.0*2.0) + (-0.5 + 4.5*pow(2.0,2.0))*4.0 + (10.0*2.0 + 2.5*pow(2.0,2.0))*pow(4.0,2.0))*((2.0 + 3.0*pow(2.0,2.0)) + (0.5*2.0)*4.0) << std::endl;
    std::cout << "Calculated: f(2,4)*g(2, 4) = " << multTestVal << std::endl << std::endl;

    TwoDimPoly testAddPoly = testPoly.addTwoDimPoly(secondTestPoly);
    double addTestVal = testAddPoly.evalAt(2.0, 4.0);
    std::cout << "Known: f(2,4)+g(2, 4) = " << ((1.0 + 2.0*2.0) + (-0.5 + 4.5*pow(2.0,2.0))*4.0 + (10.0*2.0 + 2.5*pow(2.0,2.0))*pow(4.0,2.0))+((2.0 + 3.0*pow(2.0,2.0)) + (0.5*2.0)*4.0) << std::endl;
    std::cout << "Calculated: f(2,4)+g(2, 4) = " << addTestVal << std::endl << std::endl; 

    TwoDimPoly testVarOneDerivPoly = testPoly.getDerivVarOne();
    double varOneDerivTestVal = testVarOneDerivPoly.evalAt(2.0, 4.0);
    std::cout << "Test case: df/dx = (2.0) + (9.0*x)*y + (10.0 + 5.0*x)*y**2" << std::endl;
    std::cout << "Known: df/dx(2, 4) = " << (2.0) + (9.0*2.0)*4.0 + (10.0 + 5.0*2.0)*pow(4.0,2.0) << std::endl;
    std::cout << "Calculated: df/dx(2, 4) = " << varOneDerivTestVal << std::endl << std::endl; 

    TwoDimPoly testVarTwoDerivPoly = testVarOneDerivPoly.getDerivVarTwo();
    double varTwoDerivTestVal = testVarTwoDerivPoly.evalAt(2.0, 4.0);
    std::cout << "Test case: d2f/dxdy = (9.0*x) + (20.0 + 10.0*x)*y" << std::endl;
    std::cout << "Known: d2f/dxdy(2, 4) = " << (9.0*2.0) + (20.0 + 10.0*2.0)*4.0 << std::endl;
    std::cout << "Calculated: d2f/dxdy(2, 4) = " << varTwoDerivTestVal << std::endl << std::endl; 

    TwoDimPoly testVarTwoScalarMult = testPoly.multScalar(6.0);
    double testScalarMultVal = testVarTwoScalarMult.evalAt(2.0, 4.0);
    std::cout << "Test case: f(x,y) * 6.0 = (6.0 + 12.0*x) + (-3.0 + 27.0*x**2)*y + (60.0*x + 15.0*x**2)*y**2" << std::endl;
    std::cout << "Known: f(2, 4) * 6.0 = " << (6.0 + 12.0*2.0) + (-3.0 + 27.0*pow(2.0,2.0))*4.0 + (60.0*2.0 + 15.0*pow(2.0,2.0))*pow(4.0,2.0) << std::endl;
    std::cout << "Calculated: f(2, 4) * 6.0 = " << testScalarMultVal << std::endl << std::endl; 

    // Debuggin the last test -- make a standard printing option
    //std::vector<OneDimPoly> varTwoArray = testAddPoly.getVarOneArray();
    //int varTwoOrder = testAddPoly.getVarTwoOrder();
    //int varOneOrder = 0;
    //for (int twoIdx=0; twoIdx <= varTwoOrder; ++twoIdx) {
    //    varOneOrder = varTwoArray[twoIdx].getPolyOrder();
    //    std::vector<double> coeffArray = varTwoArray[twoIdx].getCoeffArray();
    //    for (int oneIdx=0; oneIdx <= varOneOrder; ++oneIdx) {
    //        std::cout << twoIdx << "," << oneIdx << "," << coeffArray[oneIdx] << std::endl;
    //    };
    //};
};
