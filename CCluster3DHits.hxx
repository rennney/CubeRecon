#ifndef CCluster3DHits_hxx_seen
#define CCluster3DHits_hxx_seen

#include "TROOT.h"
#include "TObject.h"
#include "TMath.h"
#include "CHit3D.hxx"
#include "CCluster3D.hxx"
#include "CBond3D.hxx"
#include <iostream>
#include <vector>
#include "TVector3.h"

#include "TTmplDensityCluster.hxx"
#include "TTmplMinimalSpanningTree.hxx"

#include "LinearFit.hxx"



template<typename T>
class CPositionMetric{
public:
    double operator()(T lhs,T rhs){
        double x1=lhs.GetPosition().X();
        double y1=lhs.GetPosition().Y();
        double z1=lhs.GetPosition().Z();
        
        double x2=rhs.GetPosition().X();
        double y2=rhs.GetPosition().Y();
        double z2=rhs.GetPosition().Z();
        
        return sqrt(pow(x1-x2,2)+pow(y1-y2,2)+pow(z1-z2,2));
    }
    
};
template<typename T>
class CPositionMetricWithCharge{
public:
    double operator()(T lhs,T rhs){
        double x1=lhs.GetPosition().X();
        double y1=lhs.GetPosition().Y();
        double z1=lhs.GetPosition().Z();
        double charge1=lhs.GetCharge();
        
        double x2=rhs.GetPosition().X();
        double y2=rhs.GetPosition().Y();
        double z2=rhs.GetPosition().Z();
        double charge2=rhs.GetCharge();
        
        double chargeCorr=exp(-(charge2+charge1))*0.01;
        return sqrt(pow(x1-x2,2)+pow(y1-y2,2)+pow(z1-z2,2))+chargeCorr;
    }
    
};
typedef TTmplDensityCluster<CHit3D,CPositionMetric<CHit3D>> AlgDBScan;
typedef TTmplMinimalSpanningTree<CHit3D,CPositionMetricWithCharge<CHit3D>,int> AlgMST;

std::vector<int> indexNParents(const AlgMST::MinimalSpanningTree& tree,int i, int n){
    std::vector<int> indexes;
    int currentIndex=i;
    int parent=0;
    if(tree[currentIndex].Parent==-1){return indexes;}
    while(parent<n && tree[currentIndex].Parent!=-1){
        int parentIndex=tree[currentIndex].Parent;
        
        currentIndex=parentIndex;
   
        indexes.push_back(currentIndex);
        parent++;
    }
    return indexes;
};

std::vector<CHit3D> ConvertIdToHit(const std::vector<int>& IDs,const std::vector<CHit3D>& hits){
    std::vector<CHit3D> actualHits;
    for(std::size_t i=0;i<IDs.size();++i){
        for(std::vector<CHit3D>::const_iterator h=hits.begin();h!=hits.end();++h){
            if(IDs[i]==h->GetId()){actualHits.push_back(*h);break;}
        }
    }
    return actualHits;
};

void CCluster3DHits(TTree* tree3D){
    
    const unsigned int minHits=2;
    const double maxDist=3;
    
#define SPLITTING
#ifdef SPLITTING
    int SPLITSIZE=19; //must be odd
    int NHITSFIT=(SPLITSIZE-1)/2;
    double ANGLE=0.8;
#endif
#ifndef SPLITTING
    int SPLITSIZE=9999;
#endif
    
    std::vector<CHit3D>* hits3D=0;
    std::vector<CHit2D>* hits2D=0;
    
    
    tree3D->SetBranchAddress("3DHits",&hits3D);
    tree3D->SetBranchAddress("Unused2DHits",&hits2D);

    
    std::vector<CCluster3D> clustersFinal;
    std::vector<CHit3D> all3D;
    std::vector<CHit3D> unused3D;
    std::vector<CHit2D> unused2D;
    std::vector<CBond3D> bonds;
   
    
    TFile* hfile = new TFile("FileWith3DClusters.root","RECREATE");
    TTree *newTree = new TTree("treeWith3DClusters","tree with 3D clusters");
    
    newTree->Branch("3DClusters","vector<CCluster3D>",&clustersFinal,64000,1);
    newTree->Branch("Bonds","vector<CBond3D>",&bonds,64000,1);
    newTree->Branch("All3DHits","vector<CHit3D>",&all3D,64000,1);
    newTree->Branch("Unused3DHits","vector<CHit3D>",&unused3D,64000,1);
    newTree->Branch("Unused2DHits","vector<CHit2D>",&unused2D,64000,1);

    int nentry = tree3D->GetEntries();
    for(int i=0;i<nentry;++i){
        std::cout<<"Process event "<<i+1<<"/"<<nentry<<std::endl;
        
        tree3D->GetEntry(i);
        unused2D=(*hits2D);
        all3D=(*hits3D);
        
        std::shared_ptr<AlgDBScan> clus(new AlgDBScan(minHits,maxDist));
        std::vector<CHit3D> hits3D_mod;
        for(std::size_t h=0;h<hits3D->size();++h){
            if((*hits3D)[h].GetCharge()<2)continue;
            hits3D_mod.push_back((*hits3D)[h]);
        }
        
        clus->Cluster(hits3D_mod);
        
        unsigned int nClusters=clus->GetClusterCount();
        
        std::cout<<"Number of Clusters Formed="<<nClusters<<std::endl;
        
        std::vector<CHit3D> used3DHits;
        
        
        //nClusters
        for(std::size_t c=0;c<nClusters;++c){
            std::cout<<"+++++++++++NEW CLUSTER++++++++++++"<<std::endl;
            std::vector<CHit3D> clusteredHits=clus->GetPoints(c);

            std::shared_ptr<AlgMST> mstAlg(new AlgMST);
            std::cout<<"Cluster Size="<<clusteredHits.size()<<std::endl;
            mstAlg->AddVertices(clusteredHits.begin(),clusteredHits.end());
            mstAlg->MakeTree(*clusteredHits.begin());
            
            const AlgMST::MinimalSpanningTree& tree = mstAlg->GetTree();
            
            int placeInTree=-2;
            int topDepth=-2;
            
            for (int t = 0; t < tree.size(); ++t) {
                if(topDepth<tree[t].VertexDepth){
                    topDepth=tree[t].VertexDepth;
                    placeInTree=t;
                }
            }
            
            std::shared_ptr<AlgMST> mstAlg2(new AlgMST);
            mstAlg2->AddVertices(clusteredHits.begin(),clusteredHits.end());
            mstAlg2->MakeTree(*std::find(clusteredHits.begin(),clusteredHits.end(),tree[placeInTree].Object));
            
            
            const AlgMST::MinimalSpanningTree& tree2 = mstAlg2->GetTree();
            
            for(int t = 0; t < tree2.size(); ++t){
                 tree2[t].UserData=0;
               
            }
            
            if((int)tree2.size()<SPLITSIZE){
                CCluster3D cluster;
                for(int t = 0; t < tree2.size(); ++t){
                    if((int)t==0)cluster.SetStartPoint(tree2[t].Object.GetId());
                    if((int)t==(int)tree2.size()-1)cluster.SetEndPoint(tree2[t].Object.GetId());
                    cluster.AddConstituent(tree2[t].Object.GetId());
                }
                clustersFinal.push_back(cluster);
            }else{
#ifdef SPLITTING
                bool allVisited=0;
                std::vector<CCluster3D> clusters;
                std::vector<CHit3D> shortTails;
                std::vector<int> shortTailsInd;
                
                while(!allVisited){
                    int placeInTree2=-2;
                    int topDepth2=-2;
                    for (int t = 0; t < tree2.size(); ++t) {
                        if(topDepth2<tree2[t].VertexDepth && tree2[t].UserData==0){
                            topDepth2=tree2[t].VertexDepth;
                            placeInTree2=t;
                        }
                    }
                   
                    int startIndex=placeInTree2;
                    
                    if(startIndex==-2){allVisited=1;continue;}
                    
                    std::vector<int> parents = indexNParents(tree2,startIndex,SPLITSIZE);
              
                    tree2[startIndex].UserData=1;
                    int distToUsed=0;
                    for(std::size_t p=0;p<parents.size();++p){
                        if(tree2[parents[p]].UserData==0){distToUsed++;}else{break;}
                    }
                    if((int)parents.size()==SPLITSIZE && distToUsed>14){
                        
                        std::vector<std::pair<int,double>> runningAlphaPairs;
                        int currentIndex = parents[6];
                        bool hitroot=0;
                        bool hitUsedHit=0;
                        int step=0;
                        runningAlphaPairs.push_back(std::make_pair(startIndex,100));
                        for(int a=0;a<6;++a){
                            runningAlphaPairs.push_back(std::make_pair(parents[a],100));
                        }
                        
                        while (!hitroot && !hitUsedHit) {
                            
                            std::vector<int> parentsBefore = indexNParents(tree2,startIndex,NHITSFIT);
                            parentsBefore.insert(parentsBefore.begin(),startIndex);
                            std::vector<int> parentsAfter = indexNParents(tree2,currentIndex,NHITSFIT);
                            parentsAfter.insert(parentsAfter.begin(),currentIndex);
                            std::vector<TVector3> pointsBefore;
                            std::vector<TVector3> pointsAfter;
                            for(std::size_t p=0;p<parentsBefore.size();++p){
                                pointsBefore.push_back(tree2[parentsBefore[p]].Object.GetPosition());
                                pointsAfter.push_back(tree2[parentsAfter[p]].Object.GetPosition());
                            }
                            
                            double p0Before[6];
                            double p0After[6];
                            TVector3 unitBefore=(pointsBefore.front()-pointsBefore.back()).Unit();
                            p0Before[0]=unitBefore.X();
                            p0Before[2]=unitBefore.Y();
                            p0Before[4]=unitBefore.Z();
                            p0Before[1]=pointsBefore[NHITSFIT].X();
                            p0Before[3]=pointsBefore[NHITSFIT].Y();
                            p0Before[5]=pointsBefore[NHITSFIT].Z();
                            TVector3 unitAfter=(pointsAfter.front()-pointsAfter.back()).Unit();
                            p0After[0]=unitAfter.X();
                            p0After[2]=unitAfter.Y();
                            p0After[4]=unitAfter.Z();
                            p0After[1]=pointsAfter[0].X();
                            p0After[3]=pointsAfter[0].Y();
                            p0After[5]=pointsAfter[0].Z();
                            
                            std::vector<double> outBefore = LinearFit(pointsBefore,p0Before);
                            std::vector<double> outAfter = LinearFit(pointsAfter,p0After);
                            
                            TVector3 vectBefore(outBefore[0],outBefore[2],outBefore[4]);
                            TVector3 vectAfter(outAfter[0],outAfter[2],outAfter[4]);
                            
                            double alpha =(vectBefore.Unit()).Dot((vectAfter.Unit()));///(unitBefore.Mag()*vectBefore.Mag());
                            
                            if(fabs(alpha)<ANGLE){runningAlphaPairs.push_back(std::make_pair(currentIndex,alpha));}else{runningAlphaPairs.push_back(std::make_pair(currentIndex,100));}
                            tree2[currentIndex].UserData=1;
                            currentIndex=tree2[currentIndex].Parent;
                            startIndex=tree2[startIndex].Parent;
                            if(tree2[parentsAfter.back()].Parent==-1)hitroot=1;
                            if(tree2[parentsAfter.back()].UserData==1)hitUsedHit=1;
                            
                            for(std::size_t p=0;p<parentsBefore.size();++p){
                                tree2[parentsBefore[p]].UserData=1;
                                tree2[parentsAfter[p]].UserData=1;
                            }
                            
                        }
                        
                        std::vector<int> parentsLast = indexNParents(tree2,currentIndex,NHITSFIT);
                        parentsLast.insert(parentsLast.begin(),currentIndex);
                        for(int a=0;a<NHITSFIT;++a){
                            runningAlphaPairs.push_back(std::make_pair(parentsLast[a],100));
                        }

                        
                        std::vector<std::pair<int,double>>::iterator it= runningAlphaPairs.begin();
                        std::vector<std::pair<int,double>>::iterator clbegin= runningAlphaPairs.begin();
                        std::vector<std::pair<int,double>>::iterator clend= runningAlphaPairs.begin();
                        while (it!=runningAlphaPairs.end()) {
                            if(it->second!=100){
                                std::vector<std::pair<int,double>> kinkClust;
                                std::vector<std::pair<int,double>>::iterator it2=it;
                                while(it2->second!=100 && it2!=runningAlphaPairs.end()){
                                    kinkClust.push_back(*it2);
                                    ++it2;
                                    ++it;
                                }
                            
                                int indexAlptha=-1;
                                double minAlpha=9999;
                                for(std::size_t kk=0;kk<kinkClust.size();++kk){
                                    if(kinkClust[kk].second<minAlpha){minAlpha=kinkClust[kk].second;indexAlptha=kk;}
                                }

                                clend=it-(int)kinkClust.size()+indexAlptha;

                                CCluster3D clusterPart;
                                for(std::vector<std::pair<int,double>>::iterator f=clbegin;f!=clend+1;++f){
                                    if(f==clbegin)clusterPart.SetStartPoint(tree2[f->first].Object.GetId());
                                    if(f==clend)clusterPart.SetEndPoint(tree2[f->first].Object.GetId());
                                    clusterPart.AddConstituent(tree2[f->first].Object.GetId());
                                }
                                clusters.push_back(clusterPart);
                                clbegin=clend;
                            }else{
                                ++it;
                            }
                        }

                        clend=runningAlphaPairs.end();
                        CCluster3D clusterPart;
                        for(std::vector<std::pair<int,double>>::iterator f=clbegin;f!=clend;++f){
                            if(f==clbegin)clusterPart.SetStartPoint(tree2[f->first].Object.GetId());
                            if(f==clend-1)clusterPart.SetEndPoint(tree2[f->first].Object.GetId());
                            clusterPart.AddConstituent(tree2[f->first].Object.GetId());
                        }
                        clusters.push_back(clusterPart);
                        
                        tree2[currentIndex].UserData=1;

                    }else{
                        parents.insert(parents.begin(),startIndex);
                        if(distToUsed<5){
                            for(std::size_t p=0;p<distToUsed+1;++p){
                                shortTails.push_back(tree2[parents[p]].Object);
                                shortTailsInd.push_back(parents[p]);
                                tree2[parents[p]].UserData=1;
                            }
                        }else{
                            //CreateCluster
                            CCluster3D clusterPart;
                            int p_final=-1;
                            if(distToUsed+2<(int)parents.size()){p_final=distToUsed+2;}else{p_final=(int)parents.size();}
                            for(int p=0;p<p_final;++p){
                                if(p==0)clusterPart.SetStartPoint(tree2[parents[p]].Object.GetId());
                                if(p==p_final-1)clusterPart.SetEndPoint(tree2[parents[p]].Object.GetId());
                                clusterPart.AddConstituent(tree2[parents[p]].Object.GetId());
                                tree2[parents[p]].UserData=1;
                            }
                            clusters.push_back(clusterPart);
                        }
                    }
                    
                    //Check Visited Hits
                    int nAnused=0;
                    for (int t = 0; t < tree2.size(); ++t) {
                        if(tree2[t].UserData==0)nAnused++;
                        }
                    if(nAnused>0)allVisited=0;
                    
                    }
                std::cout<<"shortTails="<<shortTails.size()<<"; clusterparts="<<clusters.size()<<std::endl;
                int nShortsUsed=0;
                for(std::size_t st=0;st<shortTails.size();++st){
                    CHit3D hit= shortTails[st];
                    int clustnum=-1;
                    double mindist=99999;
                    for(std::size_t cc=0;cc<clusters.size();++cc){
                        std::vector<CHit3D> constituents=ConvertIdToHit(clusters[cc].GetConstituents(),all3D);
                        for(std::size_t hitc=0;hitc<constituents.size();++hitc){
                            if((constituents[hitc].GetPosition()-hit.GetPosition()).Mag()<mindist){mindist=(constituents[hitc].GetPosition()-hit.GetPosition()).Mag();clustnum=cc;}
                        }
                    }
                    if(clustnum>-1){
                        nShortsUsed++;
                        clusters[clustnum].AddConstituent(hit.GetId());
                    }
                }

                
                for(std::size_t cc=0;cc<clusters.size();++cc){
                    clustersFinal.push_back(clusters[cc]);
                }
#endif
                
            }
            
            for(std::size_t h=0;h<clusteredHits.size();++h){
                used3DHits.push_back(clusteredHits[h]);
            }

            
        }
        
        
      //GIVE CLUSTERS IDs and MAKE Bonds
        for(std::size_t cc=0;cc<clustersFinal.size();++cc){
            clustersFinal[cc].SetId((int)cc);
        }
        
        for(std::size_t cc=0;cc<clustersFinal.size();++cc){
            
            int start=clustersFinal[cc].GetStartPoint();
            bool unUsedBond=1;
            for(std::size_t j=0;j<bonds.size();++j){
                if(bonds[j].GetPoint()==start)unUsedBond=0;
            }
            if(unUsedBond){
                CBond3D newBond;
                newBond.SetPoint(start);
                newBond.AddConstituent(clustersFinal[cc].GetId());
            for(std::size_t cc2=cc+1;cc2<clustersFinal.size();++cc2){
                int start2=clustersFinal[cc2].GetStartPoint();
                int end2=clustersFinal[cc2].GetEndPoint();
                if(start==start2 || start==end2)newBond.AddConstituent(clustersFinal[cc2].GetId());
            }
                bonds.push_back(newBond);
            }
            int end=clustersFinal[cc].GetEndPoint();
             unUsedBond=1;
            for(std::size_t j=0;j<bonds.size();++j){
                if(bonds[j].GetPoint()==end)unUsedBond=0;
            }
            if(unUsedBond){
                CBond3D newBond;
                newBond.SetPoint(end);
                newBond.AddConstituent(clustersFinal[cc].GetId());
                for(std::size_t cc2=cc+1;cc2<clustersFinal.size();++cc2){
                    int start2=clustersFinal[cc2].GetStartPoint();
                    int end2=clustersFinal[cc2].GetEndPoint();
                    if(end==start2 || end==end2)newBond.AddConstituent(clustersFinal[cc2].GetId());
                }
                bonds.push_back(newBond);
            }
            
        }
        
        for(std::size_t j=0;j<bonds.size();++j){
            bonds[j].SetId((int)j);
            int pnt=bonds[j].GetPoint();
                for(std::size_t cc=0;cc<clustersFinal.size();++cc){
                    if(pnt==clustersFinal[cc].GetStartPoint() || pnt==clustersFinal[cc].GetEndPoint())clustersFinal[cc].AddBond(bonds[j].GetId());
                }
        }
        
        for(std::vector<CHit3D>::iterator h=hits3D->begin();h!=hits3D->end();++h){
            if(std::find(used3DHits.begin(),used3DHits.end(),*h)!=used3DHits.end())continue;
            unused3D.push_back(*h);
        }
        
        newTree->Fill();
        
        unused3D.clear();
        unused2D.clear();
        all3D.clear();
        clustersFinal.clear();
        bonds.clear();
}
    hfile->Write();
    hfile->Close();
    delete hfile;
};


#endif
