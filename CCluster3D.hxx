#ifndef CCluster3D_hxx_seen
#define CCluster3D_hxx_seen

#include "TROOT.h"
#include "TObject.h"
#include "TMath.h"
#include <iostream>
#include <vector>
#include "TVector3.h"
#include "CHit3D.hxx"

/// The base class for a clustered hits

class CCluster3D : public TObject {
public:
    CCluster3D() {}
    virtual ~CCluster3D() {}
    
    void SetId(int id){fID=id;}
    
    void SetConstituents(std::vector<CHit3D> hits){
        fConstituents = hits;
    }
   
    int GetId() const{return fID;}
    
    std::vector<CHit3D> GetConstituents() const{
        return fConstituents;
    }
    
    bool operator==(const CCluster3D& rhs) const { return this->GetId() == rhs.GetId();}
    
   
private:
    
    int fID;
    
    std::vector<CHit3D> fConstituents;

    ClassDef(CCluster3D,1);

};
#endif
