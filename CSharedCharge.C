/////////////////////////////////////////////////////////////////////
//
// A root macro to distribute the charge measured on each fiber among the
// cubes that are sharing the fiber.  Takes two input trees.  The first is a
// tree of the 2D hits.  The second is the tree of the associated 3D hits.  It
// will write a new tree of 3D hits which contains the charge sharing.  The
// output tree is in a new file, and has the same structure as the input tree
// of 3D hits.  This debugged using ROOT 6.10.  It is run using:
//
// root -q -l -b CHit2D.cxx+ CHit3D.cxx+ CSharedCharge.C++
//
/////////////////////////////////////////////////////////////////////

#include <string>
#include <iostream>
#include <cmath>

#include <TTree.h>
#include <TFile.h>

#include "CHit2D.hxx"
#include "CHit3D.hxx"

TFile* gFile2D = nullptr;
TTree* gTree2D = nullptr;

// Vectors of the 2D hits read from the tree.
std::vector<CHit2D>* gHitsYZ = nullptr;
std::vector<CHit2D>* gHitsXZ = nullptr;
std::vector<CHit2D>* gHitsXY = nullptr;

// Set up the tree for the 2D hits.
void Attach2D(std::string file2D) {

    if (!TClass::GetDict("CHit2D")) {
        gROOT->ProcessLine(".L CHit2D.cxx+");
    }

    gFile2D = new TFile(file2D.c_str(),"READ");
    gTree2D = (TTree*)gFile2D->Get("treeWith2DHits");

    gTree2D->SetBranchAddress("2DHitsXY",&gHitsXY);
    gTree2D->SetBranchAddress("2DHitsXZ",&gHitsXZ);
    gTree2D->SetBranchAddress("2DHitsYZ",&gHitsYZ);

}

TFile* gFile3D = nullptr;
TTree* gTree3D = nullptr;

// A pointer to a vector of 3d hits read from the tree.
std::vector<CHit3D>* gHits = nullptr;
std::vector<CHit2D>* gUnusedHits = nullptr;

// Set up the tree for the input 3D hits.
void Attach3D(std::string file3D) {

    if (!TClass::GetDict("CHit3D")) {
        gROOT->ProcessLine(".L CHit3D.cxx+");
    }

    gFile3D = new TFile(file3D.c_str(),"READ");
    gTree3D = (TTree*) gFile3D->Get("treeWith3DHitt");

    gTree3D->SetBranchAddress("3DHits",&gHits);
    gTree3D->SetBranchAddress("Unused2DHits",&gUnusedHits);

}

TFile* gOutputFile = nullptr;
TTree *gOutputTree = nullptr;
std::vector<CHit3D> gOutput3D;
std::vector<CHit2D> gOutput2D;

// Set up the tree for the output hits.
void AttachOutput(std::string output3D) {

    gOutputFile = new TFile(output3D.c_str(),"RECREATE");
    gOutputTree = new TTree("treeWith3DHitt",
                            "tree with 3D hits (w/ shared deposits)");

    gOutputTree->Branch("3DHits","std::vector<CHit3D>",&gOutput3D,64000,0);
    gOutputTree->Branch("Unused2DHits","std::vector<CHit2D>",&gOutput2D,64000,0);
}

// Find the amount of light attenuated in the fiber.  The attenuation values
// are taken from Guang's ElecSim.C macro which is tuned for the super-fgd
// prototype at CERN.
double FiberAttenuation(double ell) {
    const double LongCompFrac_FGD = 0.816;
    const double LongAtt_FGD = 11926.; //*CLHEP::mm;
    const double ShortAtt_FGD = 312.; //*CLHEP::mm;
    double arg = LongCompFrac_FGD*ell/LongAtt_FGD +
        (1.0-LongCompFrac_FGD)*ell/ShortAtt_FGD;
    return std::exp(-arg);
}

// Forward declare the Augmented objects.
struct AugmentedCube;
struct AugmentedDeposit;
struct AugmentedFiber;
std::vector<struct AugmentedCube> gAugmentedCubes;
std::vector<struct AugmentedDeposit> gAugmentedDeposits;
std::map<int, struct AugmentedFiber> gAugmentedFibers;

// An object to keep track of which deposits are contributing to a particular
// cube.
struct AugmentedCube {
    // The index of the associated cube in the vector of CHit3D objects.
    int Index;
    // The index of the associated deposits in the vector of AugmentedDeposit
    // objects.
    std::vector<int> Deposits;
    // The derivative of the q1 deposit based on the other deposits in the
    // cube.  This is actually just a function, but is here since it is
    // conceptually associated with the cube.  If one of the fibers is
    // missing, then the associated charge should be set to zero.  The formula
    // was generated using wxmaxima, and is the derivative of
    //
    // Qavg : (q1+q2+q3)/3
    // X2: [(q1-Qavg)^2 + (q2-Qavg)^2 + (q3+Qavg)^2]/(Qavg + 1)
    // f90(diff(X2,q1));
    //
    // The code is then tweaked by hand to make sure that it fits C++.  The
    // function artificially imposes the constraint that all of the charges
    // are greater or equal to zero.
    //
    double CubeDerivative(double q1, double q2, double q3) const {
        if (q1 < 0.0) q1 = 0.0;
        if (q2 < 0.0) q2 = 0.0;
        if (q3 < 0.0) q3 = 0.0;
        double deriv = (((-2.0)*(q3-(q3+q2+q1)/3.0))/3.0+((-2.0)*(q2-(q3+q2+q1)/3.0))/3.0+(4.0*(q1-(q3+q2+q1)/3.0))/3.0)/((q3+q2+q1)/3.0+1)-((std::pow(q3-(q3+q2+q1)/3.0,2.0)+std::pow(q2-(q3+q2+q1)/3.0,2.0)+std::pow(q1-(q3+q2+q1)/3.0,2.0))/std::pow((q3+q2+q1)/3.0+1.0,2.0))/3.0;
        // The returned derivative is reduced by 2/3 because this is
        // overcounting the actual number of degrees of freedom.  That's so
        // that all of the deposits are treated the same way.
        return 2.0*deriv/3.0;
    }

    // Get the total deposit in this cube.  This is implemented after the
    // AugmentedDeposit class is declared.
    double GetDeposit() const;

};

// An object to keep track of the deposits that go with a particular CHit2D
// object.  Once this is created, the calculation will not be accessing the
// CHit2D because the only information that is used will be the measured
// charge in the CHit2D.  This object keeps track of which deposits
// contributed to the CHit2D.
struct AugmentedFiber {
    // The key of the fiber in the AugmentedFiber map.
    int Index;
    // The measured number of photo electrons in the CHit2D (this is a copy so
    // we don't need to track the CHit2D.
    double Measurement;
    // A vector of the indexes of the objects in the AugmentedDeposit vector
    // that are connected to cubes on this fiber.
    std::vector<int> Deposits;
    // Check to see if the fiber is shared by multiple cubes.  If the fiber
    // isn't shared, then there is only one deposit associated with this
    // fiber, and it's value should be fixed to this measurement.
    bool IsSharedFiber() const {return Deposits.size()>1;}
    // Get the measurement for this fiber.
    double GetMeasurement() const {return Measurement;}
    // Get the current estimate of the deposits from different cubes into this
    // fiber.  When the calculation has converged the sum of the deposits will
    // be equal to the measurement.
    double GetDeposits() const; // Implemented after AugmentedDeposits!
};

// An object to keep track of the contribution of a particular CHit2D object
// to the deposited energy in a particular CHit3D object.  Since a CHit2D may
// be measuring several CHit3D deposits, this is used to track the linkage
// between the 2D and 3D hits.
struct AugmentedDeposit {
    // The index of this object in the vector of AugmentedDeposit objects.
    int Index;
    // The number of photo electrons generated in the CHit3D cube.
    double Deposit;
    // The attenuation between the cube and the MPPC.
    double Attenuation;
    // The index of the cube in the vector of AugmentedCubes;
    int Cube;
    // The key of the fiber in the map of AugmentedFibers;
    int Fiber;
    // The contribution of this deposit to the cube.
    double GetDeposit() const {return Deposit;}
    // Change the deposit by a (usually) small amount.
    void ChangeDeposit(double change) {
        while (std::abs(change) > 0.5*Deposit) change = 0.5*change;
        Deposit += change;
    }
    // The contribution of the deposit to the MPPC.
    double GetMeasurement() const {return Attenuation*Deposit;}
    // Set the deposit using the amount of charge measured at the MPPC.
    void SetMeasurement(double measurement) {
        Deposit = measurement/Attenuation;
    }
    // A convenient way to check if the fiber is shared between deposits.
    bool HasSharedFiber() {
        const struct AugmentedFiber& theFiber = gAugmentedFibers[Fiber];
        return theFiber.IsSharedFiber();
    }
    // A convenient way to get the fiber charge.
    double GetFiberMeasurement() const {
        const struct AugmentedFiber& theFiber = gAugmentedFibers[Fiber];
        return theFiber.GetMeasurement();
    }
    // A convenient way to get the fiber charge.
    double GetFiberDeposits() const {
        const struct AugmentedFiber& theFiber = gAugmentedFibers[Fiber];
        return theFiber.GetDeposits();
    }
    // The number of shared cubes
    double GetFiberCubes() const {
        const struct AugmentedFiber& theFiber = gAugmentedFibers[Fiber];
        return theFiber.Deposits.size();
    }
    // Calculate the derivative for this deposit.  This is broken into two
    // components.  The component for the fiber, and the component for the
    // cube.  Those two components are calculated in separate methods.  The
    // alpha parameter is a Lagrange multiplier.  The sum of the deposits on a
    // fiber *must* add up to the measurement, but that constraint makes the
    // basic equations being minimized singular.  The true minimum of the
    // calculation will meet the constraint.  The multiplier will start small,
    // and be increaseduntil it's (approximately) infinite.
    double  DepositDerivative(double alpha) const {
        const struct AugmentedFiber& theFiber = gAugmentedFibers[Fiber];
        if (!theFiber.IsSharedFiber()) {
            std::cout << "Calculating derivative for a fiber with only one hit"
                      << std::endl;
            return 0.0;
        }
        double fiberDerivative = FiberDerivative();
        double cubeDerivative = CubeDerivative();
        double deriv = cubeDerivative + alpha*fiberDerivative;
        return deriv;
    }
    // The contribution to the derivative by the fiber...
    double  FiberDerivative() const {
        const struct AugmentedFiber& theFiber = gAugmentedFibers[Fiber];
        double deposits = theFiber.GetDeposits();
        double measurement = theFiber.GetMeasurement();
        double deriv
            = Attenuation*(theFiber.GetDeposits()-theFiber.GetMeasurement());
        return deriv;
    }
    // The contribution to the derivative by the cube...
    double CubeDerivative() const {
        const struct AugmentedCube& theCube = gAugmentedCubes[Cube];
        std::vector<double> q;
        q.push_back(GetDeposit());
        for (std::vector<int>::const_iterator d = theCube.Deposits.begin();
             d != theCube.Deposits.end(); ++d) {
            if (*d == Index) continue;
            q.push_back(gAugmentedDeposits[*d].GetDeposit());
        }
        while (q.size() < 3) q.push_back(0.0);
        double deriv = theCube.CubeDerivative(q[0],q[1],q[2]);
        return deriv;
    }
};

double AugmentedFiber::GetDeposits() const {
    double deposits = 0.0;
    for (std::vector<int>::const_iterator d = Deposits.begin();
         d != Deposits.end(); ++d) {
        deposits += gAugmentedDeposits[*d].GetMeasurement();
    }
    return deposits;
}

double AugmentedCube::GetDeposit() const {
    double deposit = 0.0;
    for (std::vector<int>::const_iterator d = Deposits.begin();
         d != Deposits.end(); ++d) {
        deposit += gAugmentedDeposits[*d].GetDeposit();
    }
    return deposit;
}

// Routines to translate between the index in the hitXY, hitXZ, or hitYZ
// vectors and the key to the fiber map.
int YZHit2DIndex(int hit) {
    return 1000000 + hit;
}
int YZHit2DIndex(const CHit2D& hit) {
    return YZHit2DIndex(hit.GetId());
}
int XZHit2DIndex(int hit) {
    return 2000000 + hit;
}
int XZHit2DIndex(const CHit2D& hit) {
    return XZHit2DIndex(hit.GetId());
}
int XYHit2DIndex(int hit) {
    return 3000000 + hit;
}
int XYHit2DIndex(const CHit2D& hit) {
    return XYHit2DIndex(hit.GetId());
}

// Fill all of the augmented cubes, fibers and deposits.  The augmented
// objects are used to track some extra information used by the charge sharing
// calculation.  They are temporary, and all of the important information is
// transfered to the CHit3D object at the end of the calculation.
void FillAugmented(const std::vector<CHit3D>& hit3D,
                   const std::vector<CHit2D>& hitYZ,
                   const std::vector<CHit2D>& HitXZ,
                   const std::vector<CHit2D>& hitXY) {

    gAugmentedCubes.clear();
    gAugmentedFibers.clear();
    gAugmentedDeposits.clear();

    // Create the augmented fiber hits for the YZ hits.
    for (std::size_t hit = 0; hit < hitYZ.size(); ++hit) {
        struct AugmentedFiber newFiber;
        newFiber.Index = YZHit2DIndex(hitYZ[hit]);
        newFiber.Measurement = hitYZ[hit].GetCharge();
        gAugmentedFibers[newFiber.Index] = newFiber;
    }

    // Create the augmented fiber hits for the XZ hits.
    for (std::size_t hit = 0; hit < HitXZ.size(); ++hit) {
        struct AugmentedFiber newFiber;
        newFiber.Index = XZHit2DIndex(HitXZ[hit]);
        newFiber.Measurement = HitXZ[hit].GetCharge();
        gAugmentedFibers[newFiber.Index] = newFiber;
    }

    // Create the augmented fiber hits for the XY hits.
    for (std::size_t hit = 0; hit < hitXY.size(); ++hit) {
        struct AugmentedFiber newFiber;
        newFiber.Index = XYHit2DIndex(hitXY[hit]);
        newFiber.Measurement = hitXY[hit].GetCharge();
        gAugmentedFibers[newFiber.Index] = newFiber;
    }

    // Fill the augmented cubes and augmented deposits for each of the 3D
    // hits.
    for (std::size_t cube = 0; cube < hit3D.size(); ++cube) {
        const CHit3D& theCube = hit3D[cube];
        struct AugmentedCube newCube;
        newCube.Index = cube;
        // Fill the YZ fiber deposit.
        if (theCube.GetFiberCharge(0) > 0) {
            int fiber = YZHit2DIndex(theCube.Get2DConstituents(0).front());
            struct AugmentedFiber& theFiber = gAugmentedFibers[fiber];
            struct AugmentedDeposit Deposit;
            Deposit.Index = gAugmentedDeposits.size();
            Deposit.Fiber = fiber;
            Deposit.Cube = cube;
            newCube.Deposits.push_back(Deposit.Index);
            theFiber.Deposits.push_back(Deposit.Index);
#define MPPC_POSITION (0.0)
            double dist = theCube.GetPosition().X() - MPPC_POSITION;
            Deposit.Attenuation = FiberAttenuation(dist);
            Deposit.SetMeasurement(theFiber.Measurement);
            gAugmentedDeposits.push_back(Deposit);
        }
        // Fill the XZ fiber deposit.
        if (theCube.GetFiberCharge(1) > 0) {
            int fiber = XZHit2DIndex(theCube.Get2DConstituents(1).front());
            struct AugmentedFiber& theFiber = gAugmentedFibers[fiber];
            struct AugmentedDeposit Deposit;
            Deposit.Index = gAugmentedDeposits.size();
            Deposit.Fiber = fiber;
            Deposit.Cube = cube;
            newCube.Deposits.push_back(Deposit.Index);
            theFiber.Deposits.push_back(Deposit.Index);
#define MPPC_POSITION (0.0)
            double dist = theCube.GetPosition().Y() - MPPC_POSITION;
            Deposit.Attenuation = FiberAttenuation(dist);
            Deposit.SetMeasurement(theFiber.Measurement);
            gAugmentedDeposits.push_back(Deposit);
        }
        // Fill the XY fiber deposit.
        if (theCube.GetFiberCharge(2) > 0) {
            int fiber = XYHit2DIndex(theCube.Get2DConstituents(2).front());
            struct AugmentedFiber& theFiber = gAugmentedFibers[fiber];
            struct AugmentedDeposit Deposit;
            Deposit.Index = gAugmentedDeposits.size();
            Deposit.Fiber = fiber;
            Deposit.Cube = cube;
            newCube.Deposits.push_back(Deposit.Index);
            theFiber.Deposits.push_back(Deposit.Index);
#define MPPC_POSITION (0.0)
            double dist = theCube.GetPosition().Z() - MPPC_POSITION;
            Deposit.Attenuation = FiberAttenuation(dist);
            Deposit.SetMeasurement(theFiber.Measurement);
            gAugmentedDeposits.push_back(Deposit);
        }
        gAugmentedCubes.push_back(newCube);
    }

    // Set the deposit measurements so that they sum to the total fiber charge.
    for (std::vector<struct AugmentedDeposit>::iterator d
             = gAugmentedDeposits.begin();
         d !=  gAugmentedDeposits.end(); ++d) {
        if (d->GetFiberCubes()<1) continue;
        d->SetMeasurement(d->GetFiberMeasurement()/d->GetFiberCubes());
    }

    std::cout << "Augmented Cubes    " << gAugmentedCubes.size() << std::endl;
    std::cout << "Augmented Deposits " << gAugmentedDeposits.size()<< std::endl;
    std::cout << "Augmented Fibers   " << gAugmentedFibers.size() << std::endl;
}

double EvolveDeposits(double step, double alpha) {
    double totalChange = 0.0;
    int deposit = 0;
    for (std::vector<struct AugmentedDeposit>::iterator d
             = gAugmentedDeposits.begin();
         d !=  gAugmentedDeposits.end(); ++d) {
        // If this is the only deposit for the fiber, then the deposit is
        // fixed to the measurement on the fiber.
        if (!(*d).HasSharedFiber()) {
            double r = (*d).GetMeasurement() - (*d).GetFiberMeasurement();
            if (std::abs(r) > 0.0001) {
            std::cout << "delta " << r
                      << " force " << (*d).GetMeasurement()
                      << " to " << (*d).GetFiberMeasurement()
                      << std::endl;
            }
            (*d).SetMeasurement((*d).GetFiberMeasurement());
            continue;
        }

        double deriv = (*d).DepositDerivative(alpha);
        deriv = step * deriv / (*d).GetFiberCubes();
        totalChange += std::abs(deriv);
#ifdef DEBUG_CHANGES
        std::cout << "Deposit " << (*d).GetDeposit()
                  << " contributes " << (*d).GetMeasurement()
                  << " to " << (*d).GetFiberMeasurement()
                  << " out of " << (*d).GetFiberDeposits()
                  << " from " << (*d).GetFiberCubes()
                  << " change " << - deriv
                  << std::endl;
#endif
        // Never make a step of more than one photoelectron.
        if (deriv > 1.0) deriv = 1.0;
        if (deriv < -1.0) deriv = -1.0;
        // The derivative points away from the minimum so "step backwards"
        (*d).ChangeDeposit(-deriv);

    }
    return totalChange;
}

void CSharedCharge(std::string file2D = "../FileWith2DHits.root",
                   std::string file3D = "../FileWith3DHits.root",
                   std::string fileOut = "output.root") {


    Attach2D(file2D);
    Attach3D(file3D);
    AttachOutput("FileWithShared.root");

    int entries3D = gTree3D->GetEntries();
    int entries2D = gTree2D->GetEntries();

    if (entries3D != entries2D) {
        std::cout << "OOPS!  Do the files go together?" << std::endl;
        std::cout << "   2D " << file2D << std::endl;
        std::cout << "   3D " << file3D << std::endl;
        std::exit(1);
    }

    for (int entry = 0; entry < entries3D; ++entry) {
        std::cout << "Entry " << entry << std::endl;

        gTree2D->GetEntry(entry);
        gTree3D->GetEntry(entry);

        FillAugmented(*gHits, *gHitsYZ, *gHitsXZ, *gHitsXY);

        // Relax for a very long time.  This could be a lot more efficient,
        // but it's not so slow, so WTH.
        double alpha = 0.1;
        double totalChange = 0.0;
        for (int i=0; i<10000; ++i) {
            double step = 0.5/alpha;
            // Take one step.
            totalChange = EvolveDeposits(step,alpha);
            // Find the total measured charge, and the current sum of the
            // charge distributed to the cubes.  When the measured charge and
            // the distributed charge are close, stop the iterations.
            double measuredCharge = 0.0;
            double fiberCharge = 0.0;
            for (std::vector<struct AugmentedDeposit>::const_iterator
                     d = gAugmentedDeposits.begin();
                 d != gAugmentedDeposits.end(); ++d) {
                measuredCharge += (*d).GetMeasurement();
                fiberCharge += (*d).GetFiberMeasurement()/(*d).GetFiberCubes();
            }
            double diff = std::abs(measuredCharge - fiberCharge);
            double delta = diff/fiberCharge;
#ifdef DEBUG_EVOLUTION
            std::cout << i
                      << " " << totalChange
                      << " " << alpha
                      << " " << step
                      << " " << diff
                      << " " << delta
                      << std::endl;
#endif
            if (diff < 1.0 && delta < 1.0E-4) break;
            // Increase the Lagrange multiplier.
            alpha = std::min(1.001*alpha,100.0);
        };

        // Save the output
        gOutput3D.clear();
        std::copy(gHits->begin(), gHits->end(),
                  std::back_inserter(gOutput3D));
        gOutput2D.clear();
        std::copy(gUnusedHits->begin(), gUnusedHits->end(),
                  std::back_inserter(gOutput2D));
        for (std::vector<struct AugmentedCube>::const_iterator
                 c = gAugmentedCubes.begin();
             c != gAugmentedCubes.end(); ++c) {
            gOutput3D[c->Index].SetCharge(c->GetDeposit());
        }
        gOutputTree->Fill();

        // Count the number of non-overlapped fibers
        int totalFibers = 0;
        int withChargeFibers = 0;
        int multiCubeFibers = 0;
        for (std::map<int,struct AugmentedFiber>::iterator
                 f = gAugmentedFibers.begin();
             f != gAugmentedFibers.end(); ++f) {
            ++totalFibers;
            if (f->second.Deposits.size() < 1) continue;
            ++withChargeFibers;
            if (f->second.Deposits.size() < 2) continue;
            ++multiCubeFibers;
        }

        // Sum up the measured charge, and the sum of the divided deposits.
        // The values should be very close.  This only includes fibers that
        // contribute to a cube.
        double measuredCharge = 0.0;
        double depositedCharge = 0.0;
        double fiberCharge = 0.0;
        for (std::vector<struct AugmentedDeposit>::const_iterator
                 d = gAugmentedDeposits.begin();
             d != gAugmentedDeposits.end(); ++d) {
            measuredCharge += (*d).GetMeasurement();
            depositedCharge += (*d).GetDeposit();
            fiberCharge += (*d).GetFiberMeasurement()/(*d).GetFiberCubes();
        }

        // Sum up the charge in all the cubes.
        double cubeCharge = 0.0;
        for (std::vector<struct AugmentedCube>::const_iterator
                 c = gAugmentedCubes.begin();
             c != gAugmentedCubes.end(); ++c) {
            cubeCharge += c->GetDeposit();
#ifdef DEBUG_CHARGE_SPLITTING
            double qsum = 0.0;
            double qqsum = 0.0;
            double sum = 0.0;
            int sharedFibers = 0;
            for (std::vector<int>::const_iterator d = c->Deposits.begin();
                 d != c->Deposits.end(); ++d) {
                double q =  gAugmentedDeposits[*d].GetDeposit();
                std::cout << q;
                qsum += q;
                qqsum += q*q;
                sum += 1.0;
                if (gAugmentedDeposits[*d].GetFiberCubes() > 1) {
                    ++sharedFibers;
                    std::cout << "+";
                }
                std::cout << " ";
            }
            double qrms = 0.0;
            if (sum > 1.0) {
                qsum /= sum;
                qqsum /= sum;
                qrms = sqrt(qqsum-qsum*qsum);
            }
            std::cout << qrms << " " << sharedFibers << std::endl;
#endif
        }

        std::cout << " Sum of charge in used fibers " << fiberCharge
                  << std::endl;
        std::cout << " Sum of charge distributed to cubes " << measuredCharge
                  << " (should about match the above)."
                  << std::endl;
        std::cout << " Sum of the deposits assigned to cubes "
                  << depositedCharge
                  << std::endl;
        std::cout << " Sum of charge in all the cubes " << cubeCharge
                  << std::endl;
        std::cout << " Fibers --  Total: " << totalFibers
                  << ",  Above Threshold: " << withChargeFibers
                  << ",  Overlaps: " << multiCubeFibers << std::endl;
    }
}
