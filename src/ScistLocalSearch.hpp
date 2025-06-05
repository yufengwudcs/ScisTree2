//
//  ScistLocalSearch.hpp
//  
//
//  Created by Yufeng Wu on 12/16/22.
//  For fast local search of MLE
//  For the moment, only support haploid version

#ifndef ScistLocalSearch_hpp
#define ScistLocalSearch_hpp

#include <vector>
#include <map>
#include <set>
#include <string>
#include "Utils4.h"
#include "MarginalTree.h"
#include "ScistPerfPhyUtils.hpp"
class ScistPerfPhyGuideTree;
class ScistGenGenotypeMat;

// **************************************************************
// Fast rSPR local search

class ScistFastSPRLocalSearch
{
public:
    ScistFastSPRLocalSearch(const ScistGenGenotypeMat &genosInputIn, //const ScistPerfPhyGuideTree &treeGuideIn,
                            const string &strTreeInit, int nt);
    void SetNumThreads(int t) { numThreads = t; }
    void SetHeuristicMode(bool f, double frac, int t) { fHeu = f; fracHeuSPRSrc = frac; thresSPRDropStop = t; }
    void SetInitTree(const MarginalTree &mtreeInit);
    double CalcCurrMaxProb() const;
    void GetCurrTree(MarginalTree &treeOut) const;
    double GetOpt1SPRMove3(int &nodeSPRSrc, int &nodeSPRDest, int &nodeLCA, std::vector<std::pair<int,double> > *pMapSPRVal );
    
    // find best SPR moves that are ancestral to pruned subtree
    double GetOpt1SPRMove3Ances( int &nodeSPRSrc, int &nodeSPRAncDest, std::vector<std::pair<int,double> > &listSPRValAnc );
    double GetOpt1SPRMove3AncesMT(int &nodeSPRSrc, int &nodeLCADest,  std::vector<std::pair<int,double> > &listSPRValAnc );
    
    double GetOptSubtreeReroot( int &nodeSTRootNew, int &nodeST );
    int GetNumSites() const;
    int GetNumCells() const;
    double GetCurrMaxProb() const { return this->logProbCurr; }
    void SetCurrMaxProb(double p) { this->logProbCurr = p; }
    double GetSumLogAllele0Prob() const { return this->sumLogprobAllele0; }
    
    // not really interface functions but for multithreading code...
    double CalcMaxQForSPR(int site, int noder, int nodeu, int nodew, int nodev, /*int nodev1,*/ int nodewAnces, int nodevAnces, /*int nodev1Ances,*/ int nodevp, int nodevpAnces, int nodeuAnces) const;
    double CalcMaxQForReroot(int site, int node, int nodeSub, vector<vector<map<int,double> > > & tblM4) const;
    double GetQVal(int site, int node) const { return this->tblQ[site][node]; }
    double GetM0Val(int site, int node) const { return this->tblM0[site][node]; }
    bool IsNodeAncestralTo(int node1, int node2) const;
    int GetNumNodesInTree() const;
    int GetAncesIndex(int nodeDesc, int nodeAnces) const;
    double CalcMaxQPruneSite(int site, int noder, int nodeu) const;
    bool FastEvalSTForSPR(const std::set<int> &setTaxaST, double probMaxCurr ) const;
    double GetOpt1SPRMove3OnePruneST( int nodeuPrune, int &nodeSPRDest, int &nodeLCA );
    double GetOpt1SPRMove3OnePruneSTAnc( int nodeuPrune, int &nodeLCARegraft );
    
private:
    void Init();
    void InitM0Tbl();
    void InitM0TblMT();
    void InitM1Tbl();
    void InitM2Tbl();
    std::string ConsTreeFromSetClusters( const std::set<ScistPerfPhyCluster> &setClusters ) const;
    int GetNumAncesForNode(int node) const { return mapListAnces[node].size(); }
    bool QuickCheckSPRSrc(int noder, int nodeu, int noderChildOtherSide, double maxQRootCur ) const;
    bool QuickCheckSPRHeu(int noder, int nodeChild, int noderChildOtherSide, double maxQRootCur) const;
    double CalcMaxQPrune(int noder, int nodeu, int noderChildOtherSide) const;
    double CalcMaxQPruneHeu(int noder, int nodeChild, int noderChildOtherSide) const;
    double CalcMaxQCurr() const;
    
    // bounding approach
    void PreCalcQBoundForST();
    // bounding the max Q when pruning the subtree rooted at nodeu into a subtree that is rooted at the other child of noder for a site
    double CalcMaxQBoundPrunedST(int site, int nodeu, int noder, int noderPreIndex, int noderSideST, int rootIndex) const;
    // bounding the max Q when pruning the subtree rooted at nodeu into a subtree that is rooted at the other side of noder for a site
    double CalcMaxQBoundPrunedSTGen(int site, int nodeu, int noder, int noderPreIndex, int nodew, int rootIndexu, int nodeAncIndexw) const;
    // estimate qbound for a SPR from an ST into another ST
    double EstQBoundForSPR(int nodeu, int noder, int nodew) const;
    // calc max Q value for a SPR
    double CalcMaxQForSPRImp(int nodeu, int noder, int nodew) const;
    // calc max Q for a SPR where subtree is pruned to be as sibling of an ancestral node (a special case of SPR)
    double CalcMaxQForSPRAncesImp(int nodeu, int noder) const;
    // Recursively search for better SPR within a target tree to regraft
    double FindBestRegraft(int nodeu, int noder, int nodew, double qMaxCur, int &nodeSTRegraftOpt);
    // non-recursive version of FindBestRegraft
    double FindBestRegraftNoRec(int nodeu, int noder, int noderDesc, double qMaxCur, int &nodeSTRegraftOpt);
    double GetOpt1SPRMove3MT(int &nodeSPRSrc, int &nodeSPRDest, int &nodeLCA, std::vector<std::pair<int,double> > *pMapSPRVal );
    
    const ScistGenGenotypeMat &genosInput;
    //const ScistPerfPhyGuideTree &treeGuide;
    const string &strTreeInit;
    MarginalTree mtree;
    int numThreads;
    bool fHeu;
    double fracHeuSPRSrc;
    int thresSPRDropStop;
    double logProbCurr;
    double sumLogprobAllele0;
    std::vector<std::vector<double> > tblQ;
    std::vector<std::vector<double> > tblM0;
    std::vector<std::vector<double> > tblM0a;
    std::vector<std::vector< std::vector<double> > > tblM1;
    std::vector<std::vector< std::vector<double> > > tblM2;
    
    // for bounding
    // [pruned ST][parent node] = max Q if pruning ST into the other subtree under the parent node
    std::vector<std::vector<double> > tblSPRPruneSTQBounds;

    // helper: decide whether a node is ancestral to a node
    std::vector< set<int> > mapListAnces;
    std::vector< map<int,int> > mapListAncesLookup;     // which index of the ancestor for the particular ancestor to look up
};



// **************************************************************
// Fast rSPR local search to find opt tree
class ScistPerfPhyProbOnTree;

class ScistFastSPRLocalSearchLoop
{
public:
    ScistFastSPRLocalSearchLoop(const ScistGenGenotypeMat &genosInputIn, const string &strTreeInit /*const ScistPerfPhyGuideTree &treeGuideIn*/);
    double FindOpt(std::string &strTreeOpt);
    //double FindOptMulti(string &strTreeOpt);
    double FindOpt2(std::string &strTreeOpt);
    double FindOpt2Multi(string &strTreeOpt);
    void SetNumThreads(int t) {numThreads = t;}
    void SetHeuristicMode(bool f, double frac, int t) { fHeu = f; heuFracSPR = frac; thresSPRDrop = t;}
    void SetVerbose(bool f) { fVerbose = f; }
    
private:
    void FindSPRSrcLoserFrom( const std::map<int, double> &mapSrcScore, const std::set<int> &setAvoided, double frac, std::set<int> &setSrcLosers ) const;
    double ScoreTree(MarginalTree &mtree, ScistPerfPhyProbOnTree *pTreeCalc);
    double ScoreTreeMulti(MarginalTree &mtree, ScistPerfPhyProbOnTree *pTreeCalc);
    
    const ScistGenGenotypeMat &genosInput;
    const string &strTreeInit;
    MarginalTree mtreeOpt;
    double minProbInc;
    int maxIters;
    int numThreads;
    bool fHeu;
    double heuFracSPR;
    int thresSPRDrop;
    bool fVerbose;
};

// **************************************************************
// Fast NNI local search

class ScistFastNNILocalSearch
{
public:
    //ScistFastNNILocalSearch(const ScistGenGenotypeMat &genosInputIn, const ScistPerfPhyGuideTree &treeGuideIn, int nt);
    ScistFastNNILocalSearch(const ScistGenGenotypeMat &genosInputIn, const string &treeInitIn, int nt);
    double CalcCurrMaxProb(std::string &strNNIBestTree) const;
    void SetNumThreads(int n) { numThreads = n; }

private:
    double CalcCurrMaxProbMT(std::string &strNNIBestTree) const;
    void Init();
    void InitMT();
    int GetNumSites() const;
    int GetNumCells() const;
    int GetNumNodesInTree() const;
    
    const ScistGenGenotypeMat &genosInput;
    //const ScistPerfPhyGuideTree &treeGuide;
    MarginalTree mtree;
    int numThreads;
    std::vector<std::vector<double> > tblQ;
    //std::vector<double> tblQMax;    // max Q value at sites
    std::vector<std::vector<double> > tblM0;   // tblM0[s][node] = max Q value of nodes in the subtree rooted at node for site s
    std::vector<std::vector<double> > tblM1;   // tblM1[s][node] = max Q value of nodes NOT in the subtree rooted at node for site s
};

// **************************************************************
// Fast re-rooting local search

class ScistFastRerootLocalSearch
{
public:
    ScistFastRerootLocalSearch(const ScistGenGenotypeMat &genosInputIn,
                            const string &strTreeInit, int nt);
    void SetNumThreads(int t) { numThreads = t; }
    void SetInitTree(const MarginalTree &mtreeInit);
    void SetExcludeSTs(const std::set<int> &excSTs) { setSTExcluded = excSTs; }
    //void SetCandSizeLimit(double f) { fracCands = f;}    // set to any size limit f: from 0.0 to number of taxa
    double CalcCurrMaxProb() const;
    void GetCurrTree(MarginalTree &treeOut) const;
    double GetOptSubtreeReroot( int &nodeSTRootNew, int &nodeST );
    double GetOptSubtreeRerootMT( int &nodeSTRootNew, int &nodeST );
    int GetNumSites() const;
    int GetNumCells() const;
    double GetCurrMaxProb() const { return this->logProbCurr; }
    void SetCurrMaxProb(double p) { this->logProbCurr = p; }
    double GetSumLogAllele0Prob() const { return this->sumLogprobAllele0; }
    
    // not really interface functions but for multithreading code...
    double CalcMaxQForReroot(int site, int node, int nodeSub, int rootIndex, int rootSubIndex, int childOtherside, int nodeSubNodeIndex) const;
    double GetQVal(int site, int node) const { return this->tblQ[site][node]; }
    double GetM0Val(int site, int node) const { return this->tblM0[site][node]; }
    bool IsNodeAncestralTo(int node1, int node2) const;
    int GetAncesIndex(int nodeDesc, int nodeAnces) const;
    
private:
    void Init();
    void InitQTbl();
    void InitM0Tbl();
    void InitM1Tbl();
    void InitM1TblMT();
    void InitM2Tbl();
    void InitM2TblMT();
    void InitM4Tbl();
    void InitM4TblMT();
    void InitCandidates();
    int GetNumNodesInTree() const;
    std::string ConsTreeFromSetClusters( const std::set<ScistPerfPhyCluster> &setClusters ) const;
    void CollectRerootCands( );
    int GetNumAncesForNode(int node) const { return mapListAnces[node].size(); }
    
    const ScistGenGenotypeMat &genosInput;
    //const ScistPerfPhyGuideTree &treeGuide;
    const string &strTreeInit;
    MarginalTree mtree;
    int numThreads;
    bool fHeu;
    int thresSPRDropStop;
    double logProbCurr;
    double sumLogprobAllele0;
    std::vector<std::vector<double> > tblQ;
    std::vector<std::vector<double> > tblM0;
    std::vector<std::vector< std::vector<double> > > tblM1;
    std::vector<std::vector< std::vector<double> > > tblM2;
    std::vector<std::vector< std::vector<double> > > tblM4;
    // helper: decide whether a node is ancestral to a node
    std::vector< set<int> > mapListAnces;
    std::vector< map<int,int> > mapListAncesLookup;     // which index of the ancestor for the particular ancestor to look up
    std::vector<std::pair<int,int> > listCandsReroot;
    std::set<int> setSTExcluded;
    //double fracCands;
};




#endif /* ScistLocalSearch_hpp */
