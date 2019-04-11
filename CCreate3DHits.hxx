#ifndef CCreate3DHits_hxx_seen
#define CCreate3DHits_hxx_seen

#include "TROOT.h"
#include "TObject.h"
#include "TMath.h"
#include "CHit2D.hxx"
#include "CHit3D.hxx"
#include <iostream>
#include <vector>
#include "TVector3.h"

typedef std::vector<CHit2D>::iterator Chit2DIter;


struct CompareColumn
{
    double asColumn(const CHit2D& h) const{
        return h.GetColumn();
    }
    double asColumn(double c) const{
        return c;
    }
    template<typename T1, typename T2>
    bool operator()(T1 const& lhs, T2 const& rhs) const {
        return asColumn(lhs)<asColumn(rhs);
    }
};
struct CompareRow
{
    double asRow(const CHit2D& h) const{
        return h.GetRow();
    }
    double asRow(double c) const{
        return c;
    }
    template<typename T1, typename T2>
    bool operator()(T1 const& lhs, T2 const& rhs) const {
        return asRow(lhs)<asRow(rhs);
    }
};

bool SortHits(const CHit2D& lhs,  const CHit2D& rhs){
    
        if(lhs.GetColumn()!=rhs.GetColumn()){
            return lhs.GetColumn()<rhs.GetColumn();
            }else if(lhs.GetColumn()==rhs.GetColumn()){
                return lhs.GetRow()<rhs.GetRow();
                }
    
    return false;
    
};

bool SortHits3D(const CHit3D& lhs,  const CHit3D& rhs){
    
    if(lhs.GetPosition().X()!=rhs.GetPosition().X()){
        return lhs.GetPosition().X()<rhs.GetPosition().X();
    }else if(lhs.GetPosition().X()==rhs.GetPosition().X()){
        if(lhs.GetPosition().Z()!=rhs.GetPosition().Z()){
        return lhs.GetPosition().Z()<rhs.GetPosition().Z();
        }else if(lhs.GetPosition().Z()==rhs.GetPosition().Z()){
            return lhs.GetPosition().Y()<rhs.GetPosition().Y();
        }
    }
    return false;
};

std::set<double> UniqueElements(std::vector<CHit2D> tgt,int c=0){
    std::set<double> set;
    if(c==0){
        for(std::size_t i=0;i<tgt.size();++i){
          set.insert(tgt[i].GetColumn());
        }
    }else if(c==1){
        for(std::size_t i=0;i<tgt.size();++i){
            set.insert(tgt[i].GetRow());
        }
    }
    return set;
};

void C3HitCase(std::vector<CHit2D>* hitsXY,std::vector<CHit2D>* hitsXZ,std::vector<CHit2D>* hitsYZ,std::vector<CHit2D>& unusedXY,std::vector<CHit2D>& unusedXZ,std::vector<CHit2D>& unusedYZ,std::vector<CHit3D>& hits3D,double ChargeCut){
    std::sort(hitsXY->begin(),hitsXY->end(),SortHits);
    std::sort(hitsXZ->begin(),hitsXZ->end(),SortHits);
    std::sort(hitsYZ->begin(),hitsYZ->end(),SortHits);
    
    std::set<double> UcolYZ = UniqueElements(*hitsYZ,0);
    
    std::vector<CHit2D> usedHitXY;
    std::vector<CHit2D> usedHitXZ;
    std::vector<CHit2D> usedHitYZ;
    

    
    for(std::set<double>::iterator s=UcolYZ.begin();s!=UcolYZ.end();++s){
        Chit2DIter lowYZ = std::lower_bound(hitsYZ->begin(),hitsYZ->end(),*s,CompareColumn());
        Chit2DIter upYZ = std::upper_bound(hitsYZ->begin(),hitsYZ->end(),*s,CompareColumn());
        

        
    
        if(std::binary_search(hitsXZ->begin(),hitsXZ->end(),*s,CompareColumn())){
            Chit2DIter lowXZ = std::lower_bound(hitsXZ->begin(),hitsXZ->end(),*s,CompareColumn());
            Chit2DIter upXZ = std::upper_bound(hitsXZ->begin(),hitsXZ->end(),*s,CompareColumn());

            for(Chit2DIter hyz=lowYZ;hyz!=upYZ;++hyz){

                if(std::binary_search(hitsXY->begin(),hitsXY->end(),(*hyz).GetRow(),CompareColumn())){

                    Chit2DIter lowXY = std::lower_bound(hitsXY->begin(),hitsXY->end(),(*hyz).GetRow(),CompareColumn());
                    Chit2DIter upXY = std::upper_bound(hitsXY->begin(),hitsXY->end(),(*hyz).GetRow(),CompareColumn());

                    
                    for(Chit2DIter hxy=lowXY;hxy!=upXY;++hxy){

                        if(std::binary_search(lowXZ,upXZ,(*hxy).GetRow(),CompareRow())){
                            Chit2DIter lowXZ_2 = std::lower_bound(lowXZ,upXZ,(*hxy).GetRow(),CompareRow());
                            Chit2DIter upXZ_2 = std::upper_bound(lowXZ,upXZ,(*hxy).GetRow(),CompareRow());
                            for(Chit2DIter hxz=lowXZ_2;hxz!=upXZ_2;++hxz){

                            //CREATE 3D HIT
                            CHit3D hit;
                            double time = std::min((*hxz).GetTime(),(*hxy).GetTime());
                            time=std::min(time,(*hyz).GetTime());
                            hit.SetTime(time);
                            hit.SetFiberCharge((*hxy).GetCharge(),2);
                            hit.SetFiberCharge((*hxz).GetCharge(),1);
                            hit.SetFiberCharge((*hyz).GetCharge(),0);
                                
                                if((*hxy).GetCharge()<ChargeCut || (*hxz).GetCharge()<ChargeCut || (*hyz).GetCharge()<ChargeCut)continue;
                                
                            hit.AddConstituent((*hxy).GetId(),2);
                            hit.AddConstituent((*hxz).GetId(),1);
                            hit.AddConstituent((*hyz).GetId(),0);
                            TVector3 position((*hxy).GetRow(),(*hxy).GetColumn(),(*hxz).GetColumn());
                                hit.SetCharge(-1);
                            hit.SetPosition(position);
                            hits3D.push_back(hit);
                                usedHitXY.push_back(*hxy);
                                usedHitXZ.push_back(*hxz);
                                usedHitYZ.push_back(*hyz);

                            }
                        }
                    }
                }
            }
        }
    }
    
    for(Chit2DIter h=hitsXY->begin();h!=hitsXY->end();++h){
        if(std::find(usedHitXY.begin(),usedHitXY.end(),*h)!=usedHitXY.end())continue;
        unusedXY.push_back(*h);
    }
    for(Chit2DIter h=hitsXZ->begin();h!=hitsXZ->end();++h){
        if(std::find(usedHitXZ.begin(),usedHitXZ.end(),*h)!=usedHitXZ.end())continue;
        unusedXZ.push_back(*h);
    }
    for(Chit2DIter h=hitsYZ->begin();h!=hitsYZ->end();++h){
        if(std::find(usedHitYZ.begin(),usedHitYZ.end(),*h)!=usedHitYZ.end())continue;
        unusedYZ.push_back(*h);
    }
    
};



void CCreate3DHits(TTree* tree2D){
    
    double ChargeCut = 1.5;
    
    std::vector<CHit2D>* hitsXY=0;
    std::vector<CHit2D>* hitsXZ=0;
    std::vector<CHit2D>* hitsYZ=0;

    
    tree2D->SetBranchAddress("2DHitsXY",&hitsXY);
    tree2D->SetBranchAddress("2DHitsXZ",&hitsXZ);
    tree2D->SetBranchAddress("2DHitsYZ",&hitsYZ);
    
    std::vector<CHit3D> hits3D;
    std::vector<CHit2D> unused2D;
    
    TFile* hfile = new TFile("FileWith3DHits.root","RECREATE");
    TTree *newTree = new TTree("treeWith3DHitt","tree with 3D hits");
    
    newTree->Branch("3DHits","vector<CHit3D>",&hits3D,64000,0);
    newTree->Branch("Unused2DHits","vector<CHit2D>",&unused2D,64000,0);
    
    int nentry = tree2D->GetEntries();
    for(int i=0;i<nentry;++i){
        
        std::vector<CHit2D> unusedXY;
        std::vector<CHit2D> unusedXZ;
        std::vector<CHit2D> unusedYZ;
        
        
        std::vector<CHit2D> used_for2hXY;
        std::vector<CHit2D> used_for2hXZ;
        std::vector<CHit2D> used_for2hYZ;
        
        tree2D->GetEntry(i);
        std::cout<<"Number of 2DHits"<<std::endl;
        std::cout<<"NXY="<<hitsXY->size()<<"; XZ="<<hitsXZ->size()<<"; YZ="<<hitsYZ->size()<<std::endl;
        
        if(hitsXY->size()>0,hitsXZ->size()>0,hitsYZ->size()>0){
            C3HitCase(hitsXY,hitsXZ,hitsYZ,unusedXY,unusedXZ,unusedYZ,hits3D,ChargeCut);
        }
        
        std::cout<<"3D Hits Formed From 3 fibers="<<hits3D.size()<<std::endl;
        
        std::cout<<"Number of Unused 2D Hits after 3 Fiber matching"<<std::endl;
        std::cout<<"XY="<<unusedXY.size()<<"; XZ="<<unusedXZ.size()<<"; YZ="<<unusedYZ.size()<<std::endl;
        
        std::sort(unusedXY.begin(),unusedXY.end(),SortHits);
        std::sort(unusedXZ.begin(),unusedXZ.end(),SortHits);
        std::sort(unusedYZ.begin(),unusedYZ.end(),SortHits);
        
        if(unusedYZ.size()>0 && unusedXZ.size()>0){
            for(std::vector<CHit2D>::iterator hyz=unusedYZ.begin();hyz!=unusedYZ.end();++hyz){
                   for(std::vector<CHit2D>::iterator hxz=unusedXZ.begin();hxz!=unusedXZ.end();++hxz){
                       if((*hxz).GetCharge()<ChargeCut || (*hyz).GetCharge()<ChargeCut)continue;
                       if((*hyz).GetColumn()==(*hxz).GetColumn()){
                           //CREATE3DHIT
                           CHit3D hit;
                           double time = std::min((*hxz).GetTime(),(*hyz).GetTime());
                           hit.SetTime(time);
                           hit.SetFiberCharge(-999,2);
                           hit.SetFiberCharge((*hxz).GetCharge(),1);
                           hit.SetFiberCharge((*hyz).GetCharge(),0);                           hit.AddConstituent((*hxz).GetId(),1);
                           hit.AddConstituent((*hyz).GetId(),0);
                           TVector3 position((*hxz).GetRow(),(*hyz).GetRow(),(*hxz).GetColumn());
                           hit.SetPosition(position);
                           hit.SetCharge(-1);
                           hits3D.push_back(hit);
                           used_for2hYZ.push_back(*hyz);
                           used_for2hXZ.push_back(*hxz);
                       }
                   }
            }
        }
        if(unusedXY.size()>0 && unusedYZ.size()>0){
                for(std::vector<CHit2D>::iterator hyz=unusedYZ.begin();hyz!=unusedYZ.end();++hyz){
                        for(std::vector<CHit2D>::iterator hxy=unusedXY.begin();hxy!=unusedXY.end();++hxy){
                             if((*hxy).GetCharge()<ChargeCut || (*hyz).GetCharge()<ChargeCut)continue;
                            if((*hyz).GetRow()==(*hxy).GetColumn()){
                            //CREATE3DHIT
                                CHit3D hit;
                                double time = std::min((*hxy).GetTime(),(*hyz).GetTime());
                                hit.SetTime(time);
                                hit.SetFiberCharge((*hxy).GetCharge(),2);
                                hit.SetFiberCharge(-999,1);
                                hit.SetFiberCharge((*hyz).GetCharge(),0);
                                hit.AddConstituent((*hxy).GetId(),2);
                                hit.AddConstituent((*hyz).GetId(),0);
                                TVector3 position((*hxy).GetRow(),(*hyz).GetRow(),(*hyz).GetColumn());
                                hit.SetPosition(position);
                                hit.SetCharge(-1);
                                hits3D.push_back(hit);
                                used_for2hYZ.push_back(*hyz);
                                used_for2hXY.push_back(*hxy);

                            }
                        }
                }
            }
                      

        
        if(unusedXY.size()>0 && unusedXZ.size()>0){
          for(std::vector<CHit2D>::iterator hxz=unusedXZ.begin();hxz!=unusedXZ.end();++hxz){
              for(std::vector<CHit2D>::iterator hxy=unusedXY.begin();hxy!=unusedXY.end();++hxy){
                  if((*hxz).GetCharge()<ChargeCut || (*hxy).GetCharge()<ChargeCut)continue;
                  if((*hxz).GetRow()==(*hxy).GetRow()){
                      //CREATE3DHIT
                      CHit3D hit;
                      double time = std::min((*hxy).GetTime(),(*hxz).GetTime());
                      hit.SetTime(time);
                      hit.SetFiberCharge((*hxy).GetCharge(),2);
                      hit.SetFiberCharge((*hxz).GetCharge(),1);
                      hit.SetFiberCharge(-999,0);
                      hit.AddConstituent((*hxy).GetId(),2);
                      hit.AddConstituent((*hxz).GetId(),1);
                      TVector3 position((*hxy).GetRow(),(*hxy).GetColumn(),(*hxz).GetColumn());
                      hit.SetPosition(position);
                      hit.SetCharge(-1);
                      hits3D.push_back(hit);
                      used_for2hXZ.push_back(*hxz);
                      used_for2hXY.push_back(*hxy);

                  }
              }
          }
        }
        

        for(Chit2DIter h=unusedXY.begin();h!=unusedXY.end();++h){
            if(std::find(used_for2hXY.begin(),used_for2hXY.end(),*h)!=used_for2hXY.end())continue;
             unused2D.push_back(*h);
        }
        for(Chit2DIter h=unusedXZ.begin();h!=unusedXZ.end();++h){
            if(std::find(used_for2hXZ.begin(),used_for2hXZ.end(),*h)!=used_for2hXZ.end())continue;
             unused2D.push_back(*h);
        }
        for(Chit2DIter h=unusedYZ.begin();h!=unusedYZ.end();++h){
            if(std::find(used_for2hYZ.begin(),used_for2hYZ.end(),*h)!=used_for2hYZ.end())continue;
             unused2D.push_back(*h);
        }
        
        
        //Set Unique IDs for 3D hits
      //  std::sort(hits3D.begin(),hits3D.end(),SortHits3D);
        for(std::size_t h=0;h<hits3D.size();++h){
            hits3D[h].SetId(h);
        }

        
        std::cout<<"3D Hits Formed Total="<<hits3D.size()<<std::endl;
        std::cout<<"Unused Hits Total="<<unused2D.size()<<std::endl;
        newTree->Fill();
       

        hits3D.clear();
        unused2D.clear();
        
    }

    hfile->Write();
    hfile->Close();
    delete hfile;
    
    delete hitsXY;
    delete hitsXZ;
    delete hitsYZ;
};


#endif
