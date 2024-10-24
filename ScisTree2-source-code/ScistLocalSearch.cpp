//
//  ScistLocalSearch.cpp
//  
//
//  Created by Yufeng Wu on 12/16/22.
//

#include "ScistLocalSearch.hpp"
#include "ScistPerfPhyUtils.hpp"
#include "ScistPerfPhyImp.hpp"
#include "ScistGenotype.hpp"
#include "Utils3.h"
#include "Utils4.h"
#include "TreeBuilder.h"
#include "MarginalTree.h"
#include "RBT.h"
#include "PhylogenyTree.h"
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <thread>
#include <mutex>
#include <queue>
#include <stack>
#include "UtilsNumerical.h"

// **************************************************************
// Fast rSPR local search

ScistFastSPRLocalSearch :: ScistFastSPRLocalSearch(const ScistGenGenotypeMat &genosInputIn, /*const ScistPerfPhyGuideTree &treeGuideIn, */ const string &strTreeInitIn, int nt) : genosInput(genosInputIn), //treeGuide(treeGuideIn),
    strTreeInit(strTreeInitIn), numThreads(nt), fHeu(true), fracHeuSPRSrc(1.0), thresSPRDropStop(100)
{
    Init();
//cout << "SPR: initialized.\n";
}

double ScistFastSPRLocalSearch :: CalcCurrMaxProb() const
{
    double probMaxAll = 0.0;
    for(int site=0; site<GetNumSites(); ++site)
    {
        double val = this->tblM0[site][this->mtree.GetRoot()];
        if( val < 0.0 )
        {
            // if no best cut, then we just choose to not to cut
            val = 0.0;
        }
        probMaxAll += val;
    }
    return probMaxAll;
}

void ScistFastSPRLocalSearch :: GetCurrTree(MarginalTree &treeOut) const
{
    treeOut = mtree;
}

double ScistFastSPRLocalSearch :: GetOpt1SPRMove3(int &nodeSPRSrc, int &nodeSPRDest, int &nodeLCA, std::vector<std::pair<int,double> > *pMapSPRVal )
{
    // init to be non-informative
    if( pMapSPRVal != NULL )
    {
        pMapSPRVal->clear();
        (*pMapSPRVal).resize( GetNumNodesInTree() );
        for(int u=0; u<GetNumNodesInTree(); ++u)
        {
            (*pMapSPRVal)[u].first = -1;
            (*pMapSPRVal)[u].second = MAX_NEG_DOUBLE_VAL;
        }
    }
    if( this->numThreads > 1 && GetNumNodesInTree() > this->numThreads  )
    {
        //
        return GetOpt1SPRMove3MT(nodeSPRSrc, nodeSPRDest, nodeLCA, pMapSPRVal);
    }
//cout << "GetOpt1SPRMove: current tree: " << this->mtree.GetNewickNoBrLen() << endl;
//this->mtree.Dump();
    
    // consider all rSPR move
    double probMax = MAX_NEG_DOUBLE_VAL;
    nodeSPRSrc = -1;
    nodeSPRDest = -1;
    nodeLCA = -1;
    
    double maxQSoFar = MAX_NEG_DOUBLE_VAL;
    //int numSkippedSTs = 0;
    //int numEarlyAbortSTs = 0;
    
    // statt from each node
    for( int nodeu = 0; nodeu < GetNumNodesInTree(); ++nodeu )
    {
        int nodeSPRDestStep = -1, nodeLCAStep = -1;
        double maxQStep = GetOpt1SPRMove3OnePruneST(nodeu, nodeSPRDestStep, nodeLCAStep);
//cout << "nodeu:" << nodeu << ", nodeSPRDestStep:" << nodeSPRDestStep << ",nodeLCAStep:" << nodeLCAStep << ", maxQStep:" << maxQStep << ", sumLogprobAllele0:" << this->sumLogprobAllele0 << endl;
        
        // update prune ST result
        if( pMapSPRVal != NULL )
        {
            (*pMapSPRVal)[nodeu].first = nodeSPRDestStep;
            (*pMapSPRVal)[nodeu].second = maxQStep + this->sumLogprobAllele0;
        }
        
        if( maxQStep > maxQSoFar )
        {
            maxQSoFar = maxQStep;
            nodeSPRSrc = nodeu;
            nodeSPRDest = nodeSPRDestStep;
            nodeLCA = nodeLCAStep;
        }
//cout << "--- max Q value over ALL sites: " << probMaxStep << endl;
    }

//cout << "++++ Number of pruned subtrees skipped: " << numSkippedSTs << ", number of early aborted subtrees: " << numEarlyAbortSTs  << " among total " << totPrunedSTs <<" prunable subtrees" << ", probMax: " << probMax << endl;
    
    
    // if fail to find SPR, just return the smallest
    if( nodeSPRSrc < 0 || nodeSPRDest < 0 || nodeLCA < 0 )
    {
        //cout << "Warning: fail to find 1-SPR move. \n";
        return MAX_NEG_DOUBLE_VAL;
    }
    
    probMax = maxQSoFar + this->sumLogprobAllele0;
    
    //YW_ASSERT_INFO(nodeSPRSrc >= 0 && nodeSPRDest >=0 && nodeLCA >=0, "Fail to find optimal SPR");
//cout << "***** Max prob of one SPR move: " << probMax << ", nodeSPRSrc:" << nodeSPRSrc << ", nodeSPRDest:" << nodeSPRDest << endl;
    //
    return probMax;
}

// use mutex to
static void UtilFastSPRSearch4( ScistFastSPRLocalSearch *pthis, int numTotNodes, vector<bool> *pvecNodeProcessed, vector<double> *pListProb, vector<int> *pListDest, vector<int> *pListLCA, double sumLogprobAllele0 )
{
    int nodeu = 0;
    static std::mutex mut;
    while(nodeu < numTotNodes)
    {
        // find a node that is not being processed
        mut.lock();
        while( nodeu < numTotNodes && (*pvecNodeProcessed)[nodeu] == true )
        {
            ++nodeu;
        }
        if( nodeu < numTotNodes )
        {
            (*pvecNodeProcessed)[nodeu] = true;
        }
        mut.unlock();
        if( nodeu >= numTotNodes )
        {
            break;
        }
        //
        int nodeSPRDestStep = -1, nodeLCAStep = -1;
        double maxQStep = pthis->GetOpt1SPRMove3OnePruneST(nodeu, nodeSPRDestStep, nodeLCAStep);
        
        // update prune ST result
        (*pListProb)[nodeu] = maxQStep + sumLogprobAllele0;
        (*pListDest)[nodeu] = nodeSPRDestStep;
        (*pListLCA)[nodeu] = nodeLCAStep;
    }
}

double ScistFastSPRLocalSearch :: GetOpt1SPRMove3MT(int &nodeSPRSrc, int &nodeSPRDest, int &nodeLCA,  std::vector<std::pair<int,double> > *pMapSPRVal )
{
    // do multi-threading
//cout << "Start multi-threading...\n";
    int numThreadsUse = numThreads;
    if( GetNumNodesInTree() < numThreadsUse )
    {
        numThreadsUse = GetNumNodesInTree();
    }
    vector<double> listProbs(GetNumNodesInTree());
    vector<int> listDests(GetNumNodesInTree());
    vector<int> listLCAs(GetNumNodesInTree());
    //int numItems = GetNumNodesInTree();
    // vector for a node being processed
    vector<bool> vecNodeProcessed(GetNumNodesInTree());
    for(int i=0; i<GetNumNodesInTree(); ++i)
    {
        vecNodeProcessed[i] = false;
    }
    vector<thread *> listPtrThreads;
    for(int t=0; t<numThreadsUse; ++t)
    {
        thread *pthr = new thread(UtilFastSPRSearch4, this, GetNumNodesInTree(), &vecNodeProcessed, &listProbs, &listDests, &listLCAs, this->sumLogprobAllele0 );
        listPtrThreads.push_back(pthr);
    }
    
    // wait to join
    for(unsigned int j=0; j<listPtrThreads.size(); ++j)
    {
      listPtrThreads[j]->join();
    }
    for(unsigned int j=0; j<listPtrThreads.size(); ++j)
    {
      delete listPtrThreads[j];
    }
    listPtrThreads.clear();
    
    // consider all rSPR move
    double probMax = MAX_NEG_DOUBLE_VAL;
    nodeSPRSrc = -1;
    nodeSPRDest = -1;
    nodeLCA = -1;
    
    // statt from each node
    for( int nodeu = 0; nodeu < (int)listProbs.size(); ++nodeu )
    {
        int nodeSPRDestStep = listDests[nodeu], nodeLCAStep = listLCAs[nodeu];
        double maxQStep = listProbs[nodeu];
        
        // update prune ST result
        if( pMapSPRVal != NULL )
        {
            (*pMapSPRVal)[nodeu].first = nodeSPRDestStep;
            (*pMapSPRVal)[nodeu].second = maxQStep;
        }
        
        if( maxQStep > probMax )
        {
            probMax = maxQStep;
            nodeSPRSrc = nodeu;
            nodeSPRDest = nodeSPRDestStep;
            nodeLCA = nodeLCAStep;
        }
//cout << "--- max Q value over ALL sites: " << probMaxStep << endl;
    }

//cout << "++++ Number of pruned subtrees skipped: " << numSkippedSTs << ", number of early aborted subtrees: " << numEarlyAbortSTs  << " among total " << totPrunedSTs <<" prunable subtrees" << ", probMax: " << probMax << endl;
    
    
    // if fail to find SPR, just return the smallest
    if( nodeSPRSrc < 0 || nodeSPRDest < 0 || nodeLCA < 0 )
    {
        //cout << "Warning: fail to find 1-SPR move. \n";
        return MAX_NEG_DOUBLE_VAL;
    }
    
    //YW_ASSERT_INFO(nodeSPRSrc >= 0 && nodeSPRDest >=0 && nodeLCA >=0, "Fail to find optimal SPR");
//cout << "***** Max prob of one SPR move: " << probMax << ", nodeSPRSrc:" << nodeSPRSrc << ", nodeSPRDest:" << nodeSPRDest << endl;
    //
    return probMax;
}


// find best SPR moves that are ancestral to pruned subtree
double ScistFastSPRLocalSearch :: GetOpt1SPRMove3Ances(int &nodeSPRSrc, int &nodeSPRAncDest, std::vector<std::pair<int,double> > &listSPRValAnc )
{
//cout << "GetOpt1SPRMove3Ances: sum of allele 0 prob: " << this->sumLogprobAllele0 << endl;
    // init to be non-informative
    nodeSPRSrc = -1;
    nodeSPRAncDest = -1;
    listSPRValAnc.clear();
    listSPRValAnc.resize( GetNumNodesInTree() );
    for(int u=0; u<GetNumNodesInTree(); ++u)
    {
        listSPRValAnc[u].first = -1;
        listSPRValAnc[u].second = MAX_NEG_DOUBLE_VAL;
    }
    if( this->numThreads > 1 && GetNumNodesInTree() > this->numThreads  )
    {
        return GetOpt1SPRMove3AncesMT(nodeSPRSrc, nodeSPRAncDest, listSPRValAnc);
    }
//cout << "GetOpt1SPRMove: current tree: " << this->mtree.GetNewickNoBrLen() << endl;
//this->mtree.Dump();
    
    // statt from each node
    double probMax = MAX_NEG_DOUBLE_VAL;
    for( int nodeu = 0; nodeu < GetNumNodesInTree(); ++nodeu )
    {
        int  nodeLCAStep = -1;
        double maxQStep = GetOpt1SPRMove3OnePruneSTAnc(nodeu, nodeLCAStep);
        
        // update prune ST result
        listSPRValAnc[nodeu].first = nodeLCAStep;
        listSPRValAnc[nodeu].second = maxQStep + this->sumLogprobAllele0;
//cout << "--- max prob value over ALL sites: " << listSPRValAnc[nodeu].second << endl;
        
        if( listSPRValAnc[nodeu].second > probMax )
        {
            probMax = listSPRValAnc[nodeu].second;
            nodeSPRSrc = nodeu;
            nodeSPRAncDest = nodeLCAStep;
        }
    }
    
    return probMax;
}

// use mutex to
static void UtilFastSPRSearchAnces( ScistFastSPRLocalSearch *pthis, int numTotNodes, vector<bool> *pvecNodeProcessed, vector<double> *pListProb, vector<int> *pListLCA, double sumLogprobAllele0 )
{
    int nodeu = 0;
    static std::mutex mut;
    while(nodeu < numTotNodes)
    {
        // find a node that is not being processed
        mut.lock();
        while( nodeu < numTotNodes && (*pvecNodeProcessed)[nodeu] == true )
        {
            ++nodeu;
        }
        if( nodeu < numTotNodes )
        {
            (*pvecNodeProcessed)[nodeu] = true;
        }
        mut.unlock();
        if( nodeu >= numTotNodes )
        {
            break;
        }
        //
        int nodeLCAStep = -1;
        double maxQStep = pthis->GetOpt1SPRMove3OnePruneSTAnc(nodeu, nodeLCAStep);
        
        // update prune ST result
        (*pListProb)[nodeu] = maxQStep + sumLogprobAllele0;
        (*pListLCA)[nodeu] = nodeLCAStep;
    }
}

double ScistFastSPRLocalSearch :: GetOpt1SPRMove3AncesMT(int &nodeSPRSrc, int &nodeLCADest,  std::vector<std::pair<int,double> > &listSPRValAnc )
{
    // do multi-threading
//cout << "Start multi-threading...\n";
    int numThreadsUse = numThreads;
    if( GetNumNodesInTree() < numThreadsUse )
    {
        numThreadsUse = GetNumNodesInTree();
    }
    vector<double> listProbs(GetNumNodesInTree());
    vector<int> listLCAs(GetNumNodesInTree());
    // vector for a node being processed
    vector<bool> vecNodeProcessed(GetNumNodesInTree());
    for(int i=0; i<GetNumNodesInTree(); ++i)
    {
        vecNodeProcessed[i] = false;
    }
    vector<thread *> listPtrThreads;
    for(int t=0; t<numThreadsUse; ++t)
    {
        thread *pthr = new thread(UtilFastSPRSearchAnces, this, GetNumNodesInTree(), &vecNodeProcessed, &listProbs, &listLCAs, this->sumLogprobAllele0 );
        listPtrThreads.push_back(pthr);
    }
    
    // wait to join
    for(unsigned int j=0; j<listPtrThreads.size(); ++j)
    {
      listPtrThreads[j]->join();
    }
    for(unsigned int j=0; j<listPtrThreads.size(); ++j)
    {
      delete listPtrThreads[j];
    }
    listPtrThreads.clear();
    
    // consider all rSPR move
    double probMax = MAX_NEG_DOUBLE_VAL;
    nodeSPRSrc = -1;
    nodeLCADest = -1;
    
    // statt from each node
    for( int nodeu = 0; nodeu < (int)listProbs.size(); ++nodeu )
    {
        int nodeLCAStep = listLCAs[nodeu];
        double maxQStep = listProbs[nodeu];
        
        // update prune ST result
        listSPRValAnc[nodeu].first = nodeLCAStep;
        listSPRValAnc[nodeu].second = maxQStep;
        
        if( maxQStep > probMax )
        {
            probMax = maxQStep;
            nodeSPRSrc = nodeu;
            nodeLCADest = nodeLCAStep;
        }
//cout << "--- max Q value over ALL sites: " << probMaxStep << endl;
    }
    
    // if fail to find SPR, just return the smallest
    if( nodeSPRSrc < 0 || nodeLCADest < 0 )
    {
        //cout << "Warning: fail to find 1-SPR move. \n";
        return MAX_NEG_DOUBLE_VAL;
    }
    
    //YW_ASSERT_INFO(nodeSPRSrc >= 0 && nodeSPRDest >=0 && nodeLCA >=0, "Fail to find optimal SPR");
//cout << "***** Max prob of one SPR move: " << probMax << ", nodeSPRSrc:" << nodeSPRSrc << ", nodeSPRDest:" << nodeSPRDest << endl;
    //
    return probMax;
}

double ScistFastSPRLocalSearch :: GetOpt1SPRMove3OnePruneST( int nodeuPrune, int &nodeSPRDest, int &nodeLCA )
{
    //
    // try all ancestors
    double maxQSoFar = MAX_NEG_DOUBLE_VAL;
    nodeSPRDest = -1;
    nodeLCA = -1;
    int nodeAncesDesc = nodeuPrune;
    int nodeAnces = mtree.GetParent(nodeuPrune);
    if( nodeAnces < 0 )
    {
        return maxQSoFar;
    }
    
    while( nodeAnces >= 0 )
    {
        int nodew = mtree.GetSibling(nodeAncesDesc);
        int nodeRegraft = -1;
        //double qmaxOptStep = FindBestRegraft( nodeuPrune, nodeAnces, nodew, maxQSoFar, nodeRegraft );
        double qmaxOptStep = FindBestRegraftNoRec( nodeuPrune, nodeAnces, nodew, maxQSoFar, nodeRegraft );
        
        if( qmaxOptStep > maxQSoFar )
        {
            nodeSPRDest = nodeRegraft;
            nodeLCA = nodeAnces;
            maxQSoFar = qmaxOptStep;
//cout << "Max prob updated to: " << maxQSoFar + this->sumLogprobAllele0 << endl;
        }
        
        // move up
        nodeAncesDesc = nodeAnces;
        nodeAnces = mtree.GetParent(nodeAnces);
    }
    
    return maxQSoFar;
}

double ScistFastSPRLocalSearch :: GetOpt1SPRMove3OnePruneSTAnc( int nodeuPrune, int &nodeLCARegraft )
{
    // try all ancestors
    double maxQSoFar = MAX_NEG_DOUBLE_VAL;
    nodeLCARegraft = -1;
    int nodeAnces = mtree.GetParent(nodeuPrune);
    if( nodeAnces < 0 )
    {
        return maxQSoFar;
    }
    
    while( nodeAnces >= 0 )
    {
        // if it is root, do nothing
        if( nodeAnces == mtree.GetRoot() )
        {
            break;
        }
        
        double qmaxOptStep = CalcMaxQForSPRAncesImp( nodeuPrune, nodeAnces );
        
        if( qmaxOptStep > maxQSoFar )
        {
            nodeLCARegraft = nodeAnces;
            maxQSoFar = qmaxOptStep;
//cout << "Max prob updated to: " << maxQSoFar + this->sumLogprobAllele0 << endl;
        }
        
        // move up
        nodeAnces = mtree.GetParent(nodeAnces);
    }
    
    return maxQSoFar;
}

// initialization code
void ScistFastSPRLocalSearch :: Init()
{
//cout << "ScistFastSPRLocalSearch :: Init()\n";
    // construct marg tree to work with
    string strTree = this->strTreeInit;
//cout << "Initial tree: " << strTree << endl;
    ReadinMarginalTreesNewickWLenString(strTree, GetNumCells(), mtree);
    //cout << "Marginal tree: ";
    //mtree.Dump();
    
    SetInitTree(mtree);
}

void ScistFastSPRLocalSearch :: SetInitTree(const MarginalTree &mtreeInit)
{
    this->mtree = mtreeInit;
    // initialize ancestral relations
    mapListAnces.clear();
    mapListAncesLookup.clear();
    mapListAnces.resize( mtree.GetTotNodesNum() );
    mapListAncesLookup.resize( mtree.GetTotNodesNum() );
    for(int node=0; node<mtree.GetTotNodesNum(); ++node)
    {
     int nodeAnces = node;
     set<int> ss;
     int ancesIndex = 0;
     while(true)
     {
         ss.insert(nodeAnces);
         mapListAncesLookup[node][ nodeAnces ] = ancesIndex++;
         if( mtree.GetRoot() == nodeAnces)
         {
             break;
         }
         else
         {
             nodeAnces = mtree.GetParent(nodeAnces);
         }
     }
     mapListAnces[node] = ss;
    }
//cout << "Now initialize tables\n";
    //InitQTbl();
    //if( listEdgeCandidacy.size() == 0 )
    //{
    //    InitCandidates();
    //}
//cout << "pass1\n";
    InitM0Tbl();
//cout << "pass2\n";
    InitM1Tbl();
//cout << "pass3\n";
    InitM2Tbl();
//cout << "pass4\n";

    // init prob
    this->sumLogprobAllele0 = this->genosInput.CalcSumLogProbAllele0All();
    this->logProbCurr = CalcCurrMaxProb();
    this->logProbCurr += this->sumLogprobAllele0;
    
    // Testing: 01/25/2024
    // init bounding
    this->tblSPRPruneSTQBounds.clear();
    this->tblSPRPruneSTQBounds.resize( mtree.GetTotNodesNum() );

    for(int node=0; node<mtree.GetTotNodesNum(); ++node)
    {
        this->tblSPRPruneSTQBounds[node].resize( GetNumAncesForNode(node) );
        for(int kk=0; kk<(int)this->tblSPRPruneSTQBounds[node].size(); ++kk)
        {
            this->tblSPRPruneSTQBounds[node][kk] = 0.0;
        }
    }
    // now init q-bounds
    PreCalcQBoundForST();
}

void ScistFastSPRLocalSearch :: InitM0Tbl()
{
    //cout << "InitM0Tbl: " << endl;
    // collect the maximum Q values of each subtree
//#if 0
    this->tblQ.clear();
    this->tblQ.resize( GetNumSites() );
    for(int s=0; s<GetNumSites(); ++s)
    {
        this->tblQ[s].resize( mtree.GetTotNodesNum() );
    }
//#endif
    this->tblM0.clear();
    this->tblM0.resize( GetNumSites() );
    this->tblM0a.clear();
    this->tblM0a.resize(GetNumSites());
    for(int s=0; s<GetNumSites(); ++s)
    {
        this->tblM0[s].resize( mtree.GetTotNodesNum() );
        this->tblM0a[s].resize( mtree.GetTotNodesNum() );
    }
    // also init other tables
    this->tblM1.clear();
    this->tblM1.resize( GetNumSites() );
    for(int s=0; s<GetNumSites(); ++s)
    {
        this->tblM1[s].resize(  GetNumNodesInTree() );
        for(int v=0; v<GetNumNodesInTree(); ++v)
        {
            this->tblM1[s][v].resize( GetNumAncesForNode(v) );
        }
    }
    this->tblM2.clear();
    this->tblM2.resize( GetNumSites() );
    for(int s=0; s<GetNumSites(); ++s)
    {
        this->tblM2[s].resize( GetNumNodesInTree() );
        for(int v=0; v<GetNumNodesInTree(); ++v)
        {
            this->tblM2[s][v].resize( GetNumAncesForNode(v) );
        }
    }
    
    if( this->numThreads > 1 )
    {
        InitM0TblMT();
        return;
    }
//#if 0
    // Init Q first
    for(int site=0; site<GetNumSites(); ++site)
    {
//cout << "InitQTbl: site:" << site << endl;
        // do a bottom up
        vector<double> &listNodeSplitProb = this->tblQ[site];
        // init to be bad
        for(int node=0; node<mtree.GetTotNodesNum(); ++node)
        {
            //listNodeSplitProb.push_back( -1.0*HAP_MAX_INT );
            listNodeSplitProb[node] =  -1.0*HAP_MAX_INT ;
        }
        
    //cout << "CalcProbMaxForSiteHap: mtree: " << mtree.GetNewickSorted(false) << endl;
    //mtree.Dump();
        
        for(int node=0; node<mtree.GetTotNodesNum(); ++node)
        {
            //cout << "node " << node << endl;
            double logpStep = 0.0;
            if( mtree.IsLeaf(node) )
            {
                // a single leaf in the split
                int lvlbl = mtree.GetLabel(node)-1;
                //cout << "Leaf: " << lvlbl << endl;
                double p0 = this->genosInput.GetGenotypeProbAllele0At(lvlbl, site);
                if( p0 < YW_VERY_SMALL_FRACTION)
                {
                    p0 = YW_VERY_SMALL_FRACTION;
                }
                else if( p0 > 1.0-YW_VERY_SMALL_FRACTION)
                {
                    p0 = 1.0-YW_VERY_SMALL_FRACTION;
                }
                logpStep = log((1-p0)/p0);
    //cout << "Set leaf " << node << " log prob to: " << logpStep << ", p0=" << p0 << endl;
            }
            else
            {
                // get the two children and add them up
                int childLeft = mtree.GetLeftDescendant(node);
                int childRight = mtree.GetRightDescendant(node);
    //cout << "node: " << node << ", childLeft: " << childLeft << ", childRight: " << childRight << endl;
                //cout << "childLeft: " << childLeft << ", right: " << childRight << endl;
                
                YW_ASSERT_INFO( listNodeSplitProb[childLeft] > -1.0*HAP_MAX_INT, "Bad left-a" );
                YW_ASSERT_INFO( listNodeSplitProb[childRight] > -1.0*HAP_MAX_INT, "Bad right1-a" );
                logpStep = listNodeSplitProb[childLeft] + listNodeSplitProb[childRight];
            }
//cout << "log prob: " << logpStep << " for node: " << node << endl;
            listNodeSplitProb[node] = logpStep;
        }
        
        // save the probs
        //this->tblQ[site] = listNodeSplitProb;
    }
//#endif
    for(int site=0; site<GetNumSites(); ++site)
    {
//cout << "site: " << site << endl;
        for(int node=0; node<mtree.GetTotNodesNum(); ++node)
        {
            if( mtree.IsLeaf(node) )
            {
                // a single leaf in the split
                tblM0[site][node] = tblQ[site][node];
            }
            else
            {
                // get the two children and add them up
                int childLeft = mtree.GetLeftDescendant(node);
                int childRight = mtree.GetRightDescendant(node);
    //cout << "node: " << node << ", childLeft: " << childLeft << ", childRight: " << childRight << endl;
                //cout << "childLeft: " << childLeft << ", right: " << childRight << endl;
                
                YW_ASSERT_INFO( tblM0[site][childLeft] > -1.0*HAP_MAX_INT, "Bad left-b" );
                YW_ASSERT_INFO( tblM0[site][childRight] > -1.0*HAP_MAX_INT, "Bad right1-b" );
                double probM0 = std::max(tblM0[site][childLeft], tblM0[site][childRight]);
                tblM0[site][node] = std::max(probM0, tblQ[site][node]);
            }
//cout << "Set tblM0: site:" << site << ", node:" << node << ", to: " << tblM0[site][node] << endl;
        }
        
        // init M0a table: M0a[v] = max Q outside T_v (but v is allowed)
        // now init M1 table.

    //cout << "Init tblM1: " << site << endl;
        // start from root down to leaf
        // for root, set to the smallest log-prob
        this->tblM0a[site].resize(GetNumNodesInTree());
        tblM0a[site][mtree.GetRoot()] = this->tblQ[site][mtree.GetRoot()];
        for(int node=GetNumNodesInTree()-2; node >=0; --node )
        {
            int nodeSib = mtree.GetSibling(node);
            int nodePar = mtree.GetParent(node);
            YW_ASSERT_INFO(nodeSib >= 0, "Fail to find sibling");
            YW_ASSERT_INFO(nodePar >= 0, "Fail to find parent");
            tblM0a[site][node] = std::max( this->tblQ[site][node], tblM0[site][nodeSib] );
            tblM0a[site][node] = std::max( this->tblM0a[site][node], tblM0a[site][nodePar]  );
        }
    }
}

// muliti-threading
static void UtilsInitM0Tbl( ScistFastSPRLocalSearch *pthis, MarginalTree *ptree, int posStart, int posEnd, std::vector<std::vector<double> >  *pTblQ, std::vector<std::vector<double> >  *pTblM0, std::vector<std::vector<double> >  *pTblM0a, const ScistGenGenotypeMat *pgenosInput, std::vector<std::vector< vector<double> > > *pTblM1, std::vector<std::vector< vector<double> > > *pTblM2 )
{
//#if 0
    // Init Q first
    for(int node=0; node<ptree->GetTotNodesNum(); ++node)
    {
//cout << "InitQTbl: site:" << site << endl;
        // do a bottom up
        //vector<double> &listNodeSplitProb = (*pTblQ)[site];
        // init to be bad
        if( ptree->IsLeaf(node) )
        {
            // a single leaf in the split
            int lvlbl = ptree->GetLabel(node)-1;
            for(int site=posStart; site<=posEnd; ++site)
            {
                //listNodeSplitProb.push_back( -1.0*HAP_MAX_INT );
                //(*pTblQ)[site][node] =  -1.0*HAP_MAX_INT ;
                
                //cout << "node " << node << endl;
                //double logpStep = 0.0;
                
                //cout << "Leaf: " << lvlbl << endl;
                double p0 = pgenosInput->GetGenotypeProbAllele0At(lvlbl, site);
                if( p0 < YW_VERY_SMALL_FRACTION)
                {
                    p0 = YW_VERY_SMALL_FRACTION;
                }
                else if( p0 > 1.0-YW_VERY_SMALL_FRACTION)
                {
                    p0 = 1.0-YW_VERY_SMALL_FRACTION;
                }
                double logpStep = log((1-p0)/p0);
                (*pTblQ)[site][node] = logpStep;
                //cout << "Set leaf " << node << " log prob to: " << logpStep << ", p0=" << p0 << endl;
            }
        }
        else
        {
            // get the two children and add them up
            int childLeft = ptree->GetLeftDescendant(node);
            int childRight = ptree->GetRightDescendant(node);
    //cout << "node: " << node << ", childLeft: " << childLeft << ", childRight: " << childRight << endl;
                //cout << "childLeft: " << childLeft << ", right: " << childRight << endl;
                
                //YW_ASSERT_INFO( listNodeSplitProb[childLeft] > -1.0*HAP_MAX_INT, "Bad left-a" );
                //YW_ASSERT_INFO( listNodeSplitProb[childRight] > -1.0*HAP_MAX_INT, "Bad right1-a" );
            for(int site=posStart; site<=posEnd; ++site)
            {
                double logpStep = (*pTblQ)[site][childLeft] + (*pTblQ)[site][childRight];
                //cout << "log prob: " << logpStep << " for node: " << node << endl;
                (*pTblQ)[site][node] = logpStep;
            }

        }
        
        // save the probs
        //this->tblQ[site] = listNodeSplitProb;
    }
//#endif
    
    // first iterate over sites
    for(int node=0; node<ptree->GetTotNodesNum(); ++node)
    {
        if( ptree->IsLeaf(node) )
        {
            // a single leaf in the split
            for(int site=posStart; site<=posEnd; ++site)
            {
                (*pTblM0)[site][node] = (*pTblQ)[site][node];
            }
        }
        else
        {
            // get the two children and add them up
            int childLeft = ptree->GetLeftDescendant(node);
            int childRight = ptree->GetRightDescendant(node);
            //cout << "node: " << node << ", childLeft: " << childLeft << ", childRight: " << childRight << endl;
            //cout << "childLeft: " << childLeft << ", right: " << childRight << endl;
            for(int site=posStart; site<=posEnd; ++site)
            {
                //YW_ASSERT_INFO( tblM0[site][childLeft] > -1.0*HAP_MAX_INT, "Bad left-b" );
                //YW_ASSERT_INFO( tblM0[site][childRight] > -1.0*HAP_MAX_INT, "Bad right1-b" );
                double probM0 = std::max((*pTblM0)[site][childLeft], (*pTblM0)[site][childRight]);
                (*pTblM0)[site][node] = std::max(probM0, (*pTblQ)[site][node]);
            }
            //cout << "Set tblM0: site:" << site << ", node:" << node << ", to: " << tblM0[site][node] << endl;
        }
    }
    // init M0a table: M0a[v] = max Q outside T_v (but v is allowed)
    // now init M1 table.

    //cout << "Init tblM1: " << site << endl;
        // start from root down to leaf
        // for root, set to the smallest log-prob
    for(int site=posStart; site<=posEnd; ++site)
    {
        (*pTblM0a)[site][ptree->GetRoot()] = (*pTblQ)[site][ptree->GetRoot()];
    }
    
    for(int node=pthis->GetNumNodesInTree()-2; node >=0; --node )
    {
        int nodeSib = ptree->GetSibling(node);
        int nodePar = ptree->GetParent(node);
        for(int site=posStart; site<=posEnd; ++site)
        {
            //YW_ASSERT_INFO(nodeSib >= 0, "Fail to find sibling");
            //YW_ASSERT_INFO(nodePar >= 0, "Fail to find parent");
            (*pTblM0a)[site][node] = std::max( (*pTblQ)[site][node], (*pTblM0)[site][nodeSib] );
            (*pTblM0a)[site][node] = std::max( (*pTblM0a)[site][node], (*pTblM0a)[site][nodePar]  );
        }
    }
    
    // init M1
    for(int site=posStart; site<=posEnd; ++site)
    {
        // do iteration over each node from bottom up
        for(int node=0; node<ptree->GetTotNodesNum(); ++node)
        {
            int nodeAncesIndex = 0;
            // trace up until reaching the root
            int nodeAnces = node;
            double probMaxPath = MAX_NEG_DOUBLE_VAL;
            while(nodeAnces >= 0)
            {
                //
                (*pTblM1)[site][node][nodeAncesIndex++] = probMaxPath;
//cout << "Set tblM1: site:" << site << ", node:" << node << ", nodeAnces: " << nodeAnces << ", to: " << tblM1[site][node][nodeAnces] << endl;
                double prob2 = pthis->GetQVal(site, nodeAnces);
                if( prob2 > probMaxPath )
                {
                    probMaxPath = prob2;
                }
                if( nodeAnces == ptree->GetRoot() )
                {
                    break;
                }
                else
                {
                    nodeAnces = ptree->GetParent(nodeAnces);
                }
            }
        }
    }
    
    // Init M2
    for(int site=posStart; site<=posEnd; ++site)
    {
//cout << "site: " << site << endl;
        // do iteration over each node from bottom up
        for(int node=0; node<ptree->GetTotNodesNum(); ++node)
        {
            (*pTblM2)[site][node][0] = MAX_NEG_DOUBLE_VAL;
            if( node == ptree->GetRoot() )
            {
                continue;
            }
            
            // trace up until reaching the root
            int nodeAnces = ptree->GetParent( node );
            int nodeAncesIndex = 1;
            int nodeAncesBelow = node;
            double probMaxPath = MAX_NEG_DOUBLE_VAL;
            while(nodeAnces >= 0)
            {
                //
                (*pTblM2)[site][node][nodeAncesIndex++] = probMaxPath;
//cout << "set M2: node:" << node << ", nodeAnces: " << nodeAnces << " to : " << probMaxPath << endl;
                if( nodeAnces == ptree->GetRoot() )
                {
                    break;
                }
                else
                {
                    nodeAncesBelow = nodeAnces;
                    nodeAnces = ptree->GetParent(nodeAnces);
                }
                // get the side chain at nodeAnces
                int nodeSide = ptree->GetLeftDescendant(nodeAncesBelow);
                if( pthis->IsNodeAncestralTo( nodeSide, node) )
                {
                    nodeSide = ptree->GetRightDescendant(nodeAncesBelow);
                }
                //YW_ASSERT_INFO(nodeSide != nodeAncesBelow, "Fail to find side chain node");
//cout << "nodeside:" << nodeSide << ", node: " << node << ", M0: at side node: " << this->tblM0[site][nodeSide] << endl;
                
                double prob2 = pthis->GetM0Val(site, nodeSide);
                if( prob2 > probMaxPath )
                {
                    probMaxPath = prob2;
                }
            }
        }
    }
}
void ScistFastSPRLocalSearch :: InitM0TblMT()
{
//cout << "InitM0TblMT:\n";
    int numItems = GetNumSites();
    int blkSize = numItems/this->numThreads;
    int blkStart = 0;
    int numThreadsUse = this->numThreads;
    if( this->numThreads > GetNumSites() )
    {
        numThreadsUse = GetNumSites();
    }
    vector<thread *> listPtrThreads;
    for(int t=0; t<numThreadsUse; ++t)
    {
        int blkEnd = blkStart + blkSize - 1;
        if( blkEnd >= numItems  || t == numThreadsUse-1)
        {
            blkEnd = numItems-1;
        }
        YW_ASSERT_INFO(blkEnd >= blkStart, "WRONG345");
        // start thread
        //thread *pthr = new thread(UtilsInitM0Tbl, this, &mtree, blkStart, blkEnd, &tblQ, &tblM0, &tblM0a, &(this->genosInput), &tblM1, &tblM2, &tblM10, &tblM8, &tblM9, &tblM11 );
        thread *pthr = new thread(UtilsInitM0Tbl, this, &mtree, blkStart, blkEnd, &tblQ, &tblM0, &tblM0a, &(this->genosInput), &tblM1, &tblM2 );
        listPtrThreads.push_back(pthr);
        blkStart += blkSize;
    }
    
    // wait to join
    for(unsigned int j=0; j<listPtrThreads.size(); ++j)
    {
      listPtrThreads[j]->join();
    }
    for(unsigned int j=0; j<listPtrThreads.size(); ++j)
    {
      delete listPtrThreads[j];
    }
    listPtrThreads.clear();
//cout << "InitM1TblMT: multithreading done\n";
}

void ScistFastSPRLocalSearch :: InitM1Tbl()
{
    //if( this->numThreads > 1 && GetNumSites() > this->numThreads )
    if( this->numThreads > 1 )
    {
        return;
    }
    this->tblM1.clear();
    this->tblM1.resize( GetNumSites() );
//cout << "InitM1Tbl: " << endl;
    // first iterate over sites
    for(int site=0; site<GetNumSites(); ++site)
    {
//cout << "site: " << site << endl;
        // do iteration over each node from bottom up
        this->tblM1[site].resize(mtree.GetTotNodesNum());
        for(int node=0; node<mtree.GetTotNodesNum(); ++node)
        {
            this->tblM1[site][node].resize( GetNumAncesForNode(node) );
            //this->tblM1[site][node].resize( mtree.GetTotNodesNum() );
            // trace up until reaching the root
            int nodeAnces = node;
            int nodeAncesIndex = 0;
            double probMaxPath = MAX_NEG_DOUBLE_VAL;
            while(nodeAnces >= 0)
            {
                //
                this->tblM1[site][node][nodeAncesIndex++] = probMaxPath;
                //this->tblM1[site][node][nodeAnces] = probMaxPath;
//cout << "Set tblM1: site:" << site << ", node:" << node << ", nodeAncesIndex: " << nodeAncesIndex << ", nodeAnces: " << nodeAnces << ", to: " << probMaxPath << endl;
                double prob2 = this->tblQ[site][nodeAnces];
                if( prob2 > probMaxPath )
                {
                    probMaxPath = prob2;
                }
                
                if( nodeAnces == mtree.GetRoot() )
                {
                    break;
                }
                else
                {
                    nodeAnces = mtree.GetParent(nodeAnces);
                }
            }
        }
    }
}

void ScistFastSPRLocalSearch :: InitM2Tbl()
{
    if( this->numThreads > 1 )
    {
        //InitM2TblMT();
        return;
    }
//cout << "^^^^^^^Init M2:\n";
    this->tblM2.clear();
    this->tblM2.resize( GetNumSites() );

    // first iterate over sites
    for(int site=0; site<GetNumSites(); ++site)
    {
        this->tblM2[site].resize(mtree.GetTotNodesNum());
//cout << "site: " << site << endl;
        // do iteration over each node from bottom up
        for(int node=0; node<mtree.GetTotNodesNum(); ++node)
        {
            this->tblM2[site][node].resize( GetNumAncesForNode(node) );
            //this->tblM2[site][node].resize( mtree.GetTotNodesNum() );
            this->tblM2[site][node][0] = MAX_NEG_DOUBLE_VAL;
            if( node == mtree.GetRoot() )
            {
                continue;
            }
            
            // trace up until reaching the root
            int nodeAnces = mtree.GetParent( node );
            int nodeAncesIndex = 1;
            int nodeAncesBelow = node;
            double probMaxPath = MAX_NEG_DOUBLE_VAL;
            while(nodeAnces >= 0)
            {
                //
                this->tblM2[site][node][nodeAncesIndex++] = probMaxPath;
                //this->tblM2[site][node][nodeAnces] = probMaxPath;
//cout << "set M2: node:" << node << ", nodeAncesIndex: " << nodeAncesIndex << ", nodeAnces: " << nodeAnces << " to : " << probMaxPath << endl;
                if( nodeAnces == mtree.GetRoot() )
                {
                    break;
                }
                else
                {
                    nodeAncesBelow = nodeAnces;
                    nodeAnces = mtree.GetParent(nodeAnces);
                }
                // get the side chain at nodeAnces
                int nodeSide = mtree.GetLeftDescendant(nodeAncesBelow);
                if( IsNodeAncestralTo( nodeSide, node) )
                {
                    nodeSide = mtree.GetRightDescendant(nodeAncesBelow);
                }
                //YW_ASSERT_INFO(nodeSide != nodeAncesBelow, "Fail to find side chain node");
//cout << "nodeside:" << nodeSide << ", node: " << node << ", M0: at side node: " << this->tblM0[site][nodeSide] << endl;
                
                double prob2 = this->tblM0[site][nodeSide];
                if( prob2 > probMaxPath )
                {
                    probMaxPath = prob2;
                }
            }
        }
    }
}


bool ScistFastSPRLocalSearch :: IsNodeAncestralTo(int node1, int node2) const
{
    // is node1 ancestral to node2?
    //map<int, set<int> > :: const_iterator it1 = mapListAnces.find(node2);
    //return it1->second.find(node1) != it1->second.end();
    return mapListAnces[node2].find(node1) != mapListAnces[node2].end();
}

int ScistFastSPRLocalSearch :: GetNumSites() const
{
    return genosInput.GetNumSites();
}
int ScistFastSPRLocalSearch :: GetNumCells() const
{
    return genosInput.GetNumHaps();
}
int ScistFastSPRLocalSearch :: GetNumNodesInTree() const
{
    return 2*GetNumCells()-1;
}
std::string ScistFastSPRLocalSearch :: ConsTreeFromSetClusters( const std::set<ScistPerfPhyCluster> &setClusters ) const
{
    //
    // now construct tree
    ScistInfPerfPhyUtils treeBuild;
    map<int, ScistPerfPhyCluster> mapPickedClus;
    int s = 0;
    for(set<ScistPerfPhyCluster>:: iterator it = setClusters.begin(); it != setClusters.end(); ++it)
    {
        mapPickedClus[s] = *it;
        ++s;
    }
    string strTree = treeBuild.ConsTreeWCombDistClus( this->genosInput, mapPickedClus );
    return strTree;
}


double ScistFastSPRLocalSearch :: CalcMaxQCurr() const
{
    //
    double res = 0.0;
    for(int site=0; site<GetNumSites(); ++site )
    {
        res += this->tblM0[site][mtree.GetRoot()];
    }
    return res;
}

bool ScistFastSPRLocalSearch :: QuickCheckSPRSrc(int noder, int nodeu, int noderChildOtherSide, double maxQRootCur) const
{
    double qvalPrune = CalcMaxQPrune( noder, nodeu, noderChildOtherSide);
    return qvalPrune > maxQRootCur;    // || qvalRegraft > maxQRootCur;
}

double ScistFastSPRLocalSearch :: CalcMaxQPrune(int noder, int nodeu, int noderChildOtherSide) const
{
    return MAX_NEG_DOUBLE_VAL;
}

bool ScistFastSPRLocalSearch :: QuickCheckSPRHeu(int noder, int nodeChild, int noderChildOtherSide, double maxQRootCur) const
{
    double qvalPrune = CalcMaxQPruneHeu( noder, nodeChild, noderChildOtherSide);
    return qvalPrune > maxQRootCur;    // || qvalRegraft > maxQRootCur;
}

double ScistFastSPRLocalSearch :: CalcMaxQPruneHeu(int noder, int nodeChild, int noderChildOtherSide) const
{
    return MAX_NEG_DOUBLE_VAL;
}


double ScistFastSPRLocalSearch :: CalcMaxQPruneSite(int site, int noder, int nodeu) const
{
    return MAX_NEG_DOUBLE_VAL;
}


double ScistFastSPRLocalSearch :: CalcMaxQForSPR(int site, int noder, int nodeu, int nodew, int nodev, /*int nodev1,*/ int nodewAnces, int nodevAnces, /*int nodev1Ances, */ int nodevp, int nodevpAnces, int nodeuAnces) const
//double ScistFastSPRLocalSearch :: CalcMaxQForSPR(int site, int noder, int nodeu, int nodew, int nodewAnces, int nodeuAnces) const
{
    double prob2 = this->tblM0[site][nodeu];
    double prob3 = this->tblM0[site][nodew];
    double prob4 = this->tblM2[site][nodeu][nodeuAnces];
    double prob5 = this->tblM2[site][nodew][nodewAnces];
    double prob8 = MAX_NEG_DOUBLE_VAL;
    if( nodev != noder)
    {
        if( nodevp != noder )
        {
            prob8 = this->tblM1[site][nodevp][nodevpAnces] - this->tblQ[site][nodeu];
        }
    }
    double prob9 = MAX_NEG_DOUBLE_VAL;
    if( nodew != noder)
    {
        prob9 = this->tblM1[site][nodew][nodewAnces] + this->tblQ[site][nodeu];
    }
    double prob12 = this->tblM0a[site][noder];
    double res = std::max(prob2, std::max(prob3, std::max( prob4, std::max(prob5,std::max(prob8, std::max( prob9, prob12)) ))));
    if( res < 0.0 )
    {
        res = 0.0;
    }
    return res;
}

int ScistFastSPRLocalSearch :: GetAncesIndex(int nodeDesc, int nodeAnces) const
{
    std::map<int,int> :: const_iterator it = mapListAncesLookup[nodeDesc].find(nodeAnces);
    YW_ASSERT_INFO( it != mapListAncesLookup[nodeDesc].end(), "Fail to find ancestor" );
    return it->second;
}

// bounding approach
void ScistFastSPRLocalSearch :: PreCalcQBoundForST()
{
    //
    for(int node=0; node<mtree.GetTotNodesNum(); ++node)
    {
        // trace up until reaching the root
        int nodeAnces = node;
        int nodeAncesPre = -1;
        int nodeAncesInd = 0;
        while(nodeAnces >= 0)
        {
            //
            if( nodeAnces != node )
            {
                int nodeSideST = mtree.GetSibling(nodeAncesPre);
                int rootIndex = GetAncesIndex(nodeAnces, mtree.GetRoot());
                for(int site=0; site<GetNumSites(); ++site)
                {
                    double mq = CalcMaxQBoundPrunedST(site, node, nodeAnces, nodeAncesInd-1, nodeSideST, rootIndex);
                    this->tblSPRPruneSTQBounds[node][nodeAncesInd] += mq;
                }
            }
            ++nodeAncesInd;
            
            if( nodeAnces == mtree.GetRoot() )
            {
                break;
            }
            else
            {
                nodeAncesPre = nodeAnces;
                nodeAnces = mtree.GetParent(nodeAnces);
            }
        }
    }
#if 0
    // calculate the largest Q for now
    double maxQ0 = 0.0;
    for(int site=0; site<GetNumSites(); ++site)
    {
        maxQ0 += this->tblM0[site][mtree.GetRoot()];
    }
    cout << "[maxQ]: " << maxQ0 << ".  Computed Q bounds: \n";
    for(int node=0; node<mtree.GetTotNodesNum(); ++node)
    {
        cout << "Node " << node << ": ";
        DumpDoubleVec( this->tblSPRPruneSTQBounds[node] );
    }
    //exit(1);
#endif
}

// bounding the max Q when pruning the subtree rooted at nodeu into a subtree that is rooted at the other child of noder for a site, noderPreIndex: index of child of noder that is ancestral to nodeu
double ScistFastSPRLocalSearch :: CalcMaxQBoundPrunedST(int site, int nodeu, int noder, int noderPreIndex, int noderSideST, int rootIndex) const
{
    double res = 0.0;
    YW_ASSERT_INFO( noderPreIndex >=0 && noderSideST >= 0, "Cannot be negative" );
    //
    double qu = tblQ[site][nodeu];
    //double qconst = std::max( this->tblM0[site][nodeu], std::max( tblQ[site][noder], std::max( this->tblM2[site][noder][rootIndex], this->tblM2[site][nodeu][noderPreIndex+1]  ) ) );
    double qconst = std::max( this->tblM0[site][nodeu], std::max( tblM0a[site][noder], this->tblM2[site][nodeu][noderPreIndex+1] ) );
    res = std::max( qconst, std::max(this->tblM1[site][nodeu][noderPreIndex] - qu, std::max(this->tblM0[site][noderSideST] + qu, this->tblM0[site][noderSideST]) ) );
    return max(res, 0.0);
}

// bounding the max Q when pruning the subtree rooted at nodeu into a subtree that is rooted at the other side of noder for a site, noderPreIndex: index of child of noder that is ancestral to nodeu
double ScistFastSPRLocalSearch :: CalcMaxQBoundPrunedSTGen(int site, int nodeu, int noder, int noderPreIndex, int nodew, int rootIndexu, int nodeAncIndexw) const
{
    double res = 0.0;
    YW_ASSERT_INFO( noderPreIndex >=0 && nodeAncIndexw >= 0, "Cannot be negative" );
    //
    double qu = tblQ[site][nodeu];
    //double qconst = std::max( this->tblM0[site][nodeu], std::max( tblQ[site][noder], std::max( this->tblM2[site][noder][rootIndex], this->tblM2[site][nodeu][noderPreIndex+1]  ) ) );
    double qconst = std::max( this->tblM0[site][nodeu], std::max( tblM0a[site][noder], this->tblM2[site][nodeu][noderPreIndex+1] ) );
    res = std::max( qconst, std::max(this->tblM1[site][nodeu][noderPreIndex+1] - qu, std::max(this->tblM0[site][nodew] + qu, this->tblM0[site][nodew]) ) );
    res = std::max( res, this->tblM2[site][nodew][nodeAncIndexw]  );
    res = std::max( res, this->tblM1[site][nodew][nodeAncIndexw] + qu );
    return max(res, 0.0);
}

double ScistFastSPRLocalSearch :: EstQBoundForSPR(int nodeu, int noder, int nodew) const
{
    double res = 0.0;
    int noderAncIndex = GetAncesIndex(nodeu, noder);
    int rootIndex = GetAncesIndex(noder, mtree.GetRoot());
    int nodeAncIndexw = GetAncesIndex(nodew, noder);
    for(int site=0; site<GetNumSites(); ++site)
    {
        res += CalcMaxQBoundPrunedSTGen( site, nodeu, noder, noderAncIndex-1, nodew, rootIndex, nodeAncIndexw );
    }
    return res;
}

// Recursively search for better SPR within a target tree to regraft
double ScistFastSPRLocalSearch :: FindBestRegraft(int nodeu, int noder, int nodew, double qMaxCur, int &nodeSTRegraftOpt)
{
    //
    nodeSTRegraftOpt = -1;      // set to no-result for now
    
    // calc a quick bound; if not better than current one, stop
    double qbound = EstQBoundForSPR(nodeu, noder, nodew);
    //const double MIN_INC = log(1.00001);
    const double MIN_INC = 0.0;
    if( qbound <= qMaxCur + MIN_INC )
    {
        return MAX_NEG_DOUBLE_VAL;
    }
    
    // first calc the regraft to nodew's maxQ
    double maxQw = CalcMaxQForSPRImp(nodeu, noder, nodew);
    nodeSTRegraftOpt = nodew;

    // if it is a leaf, just return its SPR value
    if( mtree.IsLeaf(nodew) == false )
    {
        int cn1 = mtree.GetLeftDescendant(nodew);
        int cn2 = mtree.GetRightDescendant(nodew);
        int pruneST1 = -1, pruneST2 = -1;
        double maxQcn1 = FindBestRegraft(nodeu, noder, cn1, std::max(qMaxCur, maxQw), pruneST1);
        if( maxQcn1 > maxQw )
        {
            maxQw = maxQcn1;
            nodeSTRegraftOpt = pruneST1;
        }
        double maxQcn2 = FindBestRegraft(nodeu, noder, cn2, std::max(qMaxCur, maxQw), pruneST2);
        if( maxQcn2 > maxQw )
        {
            maxQw = maxQcn2;
            nodeSTRegraftOpt = pruneST2;
        }
    }
    //YW_ASSERT_INFO( nodeSTRegraftOpt >=0, "Fail to find regraft" );
    return maxQw;
}

// non-recursive version of FindBestRegraft
double ScistFastSPRLocalSearch :: FindBestRegraftNoRec(int nodeu, int noder, int noderDesc, double qMaxCur, int &nodeSTRegraftOpt)
{
    //
    nodeSTRegraftOpt = -1;      // set to no-result for now
    
    // create a queue of subtrees to evaluate
    stack<int> queueSTsRegraft;
    queueSTsRegraft.push(noderDesc);
    
    double maxQSoFar = qMaxCur;
    
    // process until no STs to process
    while( queueSTsRegraft.empty() == false )
    {
        int nodew = queueSTsRegraft.top();
        queueSTsRegraft.pop();
        
        // calc a quick bound; if not better than current one, stop
        double qbound = EstQBoundForSPR(nodeu, noder, nodew);
        //const double MIN_INC = log(1.00001);
        const double MIN_INC = 0.0;
        if( qbound <= maxQSoFar + MIN_INC )
        {
            // donot go down this subtree
            continue;
        }
        
        // first calc the regraft to nodew's maxQ
        double maxQw = CalcMaxQForSPRImp(nodeu, noder, nodew);
        int nodeSTRegraftOptStep = nodew;
        
        if( maxQw > maxQSoFar )
        //if( maxQw > maxQSoFar + MIN_INC )
        {
            maxQSoFar = maxQw;
            nodeSTRegraftOpt = nodeSTRegraftOptStep;
        }
        // if it is a leaf, just return its SPR value
        if( mtree.IsLeaf(nodew) == false )
        {
            int cn1 = mtree.GetLeftDescendant(nodew);
            int cn2 = mtree.GetRightDescendant(nodew);
            queueSTsRegraft.push(cn1);
            queueSTsRegraft.push(cn2);
        }
    }
    //YW_ASSERT_INFO( nodeSTRegraftOpt >=0, "Fail to find regraft" );
    return maxQSoFar;
}

double ScistFastSPRLocalSearch :: CalcMaxQForSPRImp(int nodeu, int noder, int nodew) const
{
    int nodev = mtree.GetParent(nodeu);
    //int nodev1 = mtree.GetSibling(nodeu);
    int nodewAnces = GetAncesIndex(nodew, noder);
    int nodevAnces = GetAncesIndex(nodev, noder);
    //int nodev1Ances = GetAncesIndex(nodev1, noder);
    int nodeuAnces = GetAncesIndex(nodeu, noder);
    int nodevp = -1;
    int nodevpAnces = -1;
    if( nodev != noder)
    {
        nodevp = mtree.GetParent(nodev);
        nodevpAnces = GetAncesIndex(nodevp, noder);
    }
    double probMaxStep = 0.0;
    for(int site=0; site<GetNumSites(); ++site)
    {
        double probs = CalcMaxQForSPR(site, noder, nodeu, nodew, nodev, /*nodev1,/*/ nodewAnces, nodevAnces, /*nodev1Ances,*/ nodevp, nodevpAnces, nodeuAnces);
        //cout << "Max Q value for SPR: site:" << site << ", noder:" << noder << ", nodeu:" << nodeu << ", nodew:" << nodew << ", v,v1,wa,va,v1a " << nodev << "," << nodev1 << "," << nodewAnces << "," << nodevAnces << "," << nodev1Ances << " is " << probs << endl;
        probMaxStep += probs;
    }
    return probMaxStep;
}

// calc max Q for a SPR where subtree is pruned to be as sibling of an ancestral node (a special case of SPR)
double ScistFastSPRLocalSearch :: CalcMaxQForSPRAncesImp(int nodeu, int noder) const
{
//cout << "CalcMaxQForSPRAncesImp: nodeu: " << nodeu << ", noder: " << noder << endl;
    // assume: noder != nodeu and noder != parent(nodeu)
    int nodev = mtree.GetParent(nodeu);
    
    if( nodeu == noder || nodev == noder)
    {
        return MAX_NEG_DOUBLE_VAL;
    }
    int nodeuAnces = GetAncesIndex(nodeu, noder);
    int nodevp = -1;
    int nodevpAnces = -1;
    nodevp = mtree.GetParent(nodev);
    nodevpAnces = GetAncesIndex(nodevp, noder);
    int noderSib = -1;
    if( noder != mtree.GetRoot() )
    {
        noderSib = mtree.GetSibling(noder);
    }
    // child of noder that is on the other side
    int nodeSide = -1;
    if( mtree.IsLeaf(noder) == false && noder != nodeu)
    {
        nodeSide = mtree.GetLeftDescendant(noder);
        if( IsNodeAncestralTo( nodeSide, nodeu ) )
        {
            nodeSide = mtree.GetRightDescendant(noder);
        }
    }
    double probMaxStep = 0.0;
    for(int site=0; site<GetNumSites(); ++site)
    {
        //
        double maxQSite = MAX_NEG_DOUBLE_VAL;
        if( noderSib >= 0 )
        {
            maxQSite = this->tblM0a[site][noder];
        }
        if( nodeSide >= 0 )
        {
            maxQSite = std::max(maxQSite, this->tblM0[site][nodeSide]);
        }
        
        double prob4 = this->tblM2[site][nodeu][nodeuAnces];

        double prob8 = MAX_NEG_DOUBLE_VAL;
        if( nodevp != noder )
        {
            prob8 = this->tblM1[site][nodevp][nodevpAnces] - this->tblQ[site][nodeu];
        }
        maxQSite = std::max( maxQSite, std::max(prob4, prob8) );
        
        //cout << "Max Q value for SPR: site:" << site << ", noder:" << noder << ", nodeu:" << nodeu << ", nodew:" << nodew << ", v,v1,wa,va,v1a " << nodev << "," << nodev1 << "," << nodewAnces << "," << nodevAnces << "," << nodev1Ances << " is " << probs << endl;
        probMaxStep += maxQSite;
//cout << "site: " << site << ", maxQsite: " << maxQSite << endl;
    }
    return probMaxStep;
}

// **************************************************************
// Fast rSPR local search to find opt tree

// default search opt parameters
const int DEF_MAX_NUM_ITERS = 1000;
//const int DEF_MAX_MULTI_SPRS_MOVES = 0;
const double DEF_THRES_PROB_INC = log(1.00001);

ScistFastSPRLocalSearchLoop :: ScistFastSPRLocalSearchLoop(const ScistGenGenotypeMat &genosInputIn, const string &strTreeInitIn /*const ScistPerfPhyGuideTree &treeGuideIn*/) : genosInput(genosInputIn), strTreeInit(strTreeInitIn), /*treeGuide(treeGuideIn),*/ minProbInc(DEF_THRES_PROB_INC), maxIters(DEF_MAX_NUM_ITERS), numThreads(1), /*maxMultiSPRsMoves(DEF_MAX_MULTI_SPRS_MOVES), */ fHeu(true), heuFracSPR(1.0), thresSPRDrop(100), fVerbose(false)
{
}
double ScistFastSPRLocalSearchLoop :: FindOpt(string &strTreeOpt)
{
    return FindOpt2(strTreeOpt);
}

double ScistFastSPRLocalSearchLoop :: FindOpt2(string &strTreeOpt)
{
    // TEST, don't do multiple moves for now
    return FindOpt2Multi(strTreeOpt);
#if 0
    // now adjust by the overall prob(allele=0)
    // strTreeOpt: already initialzed with the initial tree
    //double sumLogprobAllele0 = this->genosInput.CalcSumLogProbAllele0All();
//cout << "sum of log prob allele 0 = " << sumLogprobAllele0 << endl;
//cout << "FindOpt...\n";
    //
    ScistFastSPRLocalSearch treeSearchCur(this->genosInput, /*this->treeGuide, */ this->strTreeInit, this->numThreads );
    treeSearchCur.SetHeuristicMode( this->fHeu, this->heuFracSPR, this->thresSPRDrop );
    //double probMaxCurr = treeSearchCur.CalcCurrMaxProb();
    double probMaxCurr = treeSearchCur.GetCurrMaxProb();
    //cout << "Initial SPR log-likelihood: " << probMaxCurr << endl;
//cout << "Initial max Q value : " << probMaxCurr << endl;
    
    // keep track for each SPR source, the maximum changed prob it can get
    int numHaps = this->genosInput.GetNumHaps();
    vector<double> mapMaxProbSPRSrc( 2*numHaps ), mapMaxProbSPRDest(2*numHaps);
    for(int i=0; i<2*numHaps; ++i)
    {
        mapMaxProbSPRSrc[i] = MAX_NEG_DOUBLE_VAL;
        mapMaxProbSPRDest[i] = MAX_NEG_DOUBLE_VAL;
    }
    //set<int> setLowPerformerSrc;
    
    //probMaxCurr += sumLogprobAllele0;
    
    MarginalTree mTreeCurr;
    treeSearchCur.GetCurrTree(mTreeCurr);
    
    //double score0 = ScoreTree(mTreeCurr);
    //cout << "(INITIIAL tree: Re-computing log-likelihood: " << score0 << endl;
    
    // now start SPR local search
    //bool fMultiModeEnabled = true;
    int numIters = 0;
    
    while(++numIters <= this->maxIters)
    {
//cout << "**** Current max prob: " << probMaxCurr << "    treeSearchCur: ";
//cout << mTreeCurr.GetNewickNoBrLen() << endl;
        // search for opt
        int nodeSPRSrc, nodeSPRDest, nodeLCA;
        double probMaxStep = treeSearchCur.GetOpt1SPRMove3(nodeSPRSrc, nodeSPRDest, nodeLCA, NULL );
//cout << "Done with 1-SPR neighborhood search\n";
        
        // adjust
        //probMaxStep += sumLogprobAllele0;
        
        if( probMaxStep < probMaxCurr + this->minProbInc)
        {
            //cout << "SPR: cannot improve.\n";
            break;
        }
        
        probMaxCurr = probMaxStep;
        // move to the next
        //cout << "Prior to SPR: tree = " << mTreeCurr.GetNewickNoBrLen() << endl;
        //mTreeCurr.Dump();
        //double scorep = ScoreTree(mTreeCurr);
        //cout << "(Before SPR: Re-computing log-likelihood: " << scorep << endl;
        //mTreeCurr.Dump();
        
        mTreeCurr.PerformSPR(nodeSPRSrc, nodeSPRDest);
//cout << "AFTER performing 1 SPR: tree is: " << mTreeCurr.GetNewickNoBrLen() << endl;
//mTreeCurr.Dump();
        
        // temp code
        //double score = ScoreTree(mTreeCurr);
        //cout << "Re-computing log-likelihood: " << score << endl;
        
        treeSearchCur.SetInitTree(mTreeCurr);
        treeSearchCur.SetCurrMaxProb( probMaxCurr );
        
        if( fVerbose )
        {
            //cout << "++ 1-rSPR local search: likelihood improved to " << probMaxCurr << endl;
        }
    }
    //cout << "Total number of rSPR local search steps: " << numIters << endl;
    //cout << "**** Maximum log-likelihood: " << probMaxCurr << endl;
    //cout << "Constructed single cell phylogeny: " << mTreeCurr.GetNewickNoBrLen() << endl;
    
    strTreeOpt = mTreeCurr.GetNewickNoBrLen();
    
    // now valid to make sure it is indeed the same prob
    //ScistPerfPhyGuideTree treeGuidThis;
    //treeGuidThis.InitDecAll(strTreeOpt);
    //ScistFastSPRLocalSearch treeSearchCur2(this->genosInput, treeGuidThis );
    //double probMaxCurr2 = treeSearchCur.CalcCurrMaxProb();
    //probMaxCurr2 += sumLogprobAllele0;
    //cout << "Recomputed log-likelihood: " << probMaxCurr2 << endl;
    
    return probMaxCurr;
#endif
}

double ScistFastSPRLocalSearchLoop :: FindOpt2Multi(string &strTreeOpt)
{
    //
    ScistFastSPRLocalSearch treeSearchCur(this->genosInput,this->strTreeInit, this->numThreads );
    //treeSearchCur.SetHeuristicMode( this->fHeu, this->heuFracSPR, this->thresSPRDrop );
    //treeSearchCur.SetHeuristicMode( false, this->heuFracSPR, this->thresSPRDrop );
    //double probMaxCurr = treeSearchCur.CalcCurrMaxProb();
    double probMaxCurr = treeSearchCur.GetCurrMaxProb();
    //cout << "Initial SPR log-likelihood: " << probMaxCurr << endl;
//cout << "Initial max Q value : " << probMaxCurr << endl;
    
    // keep track for each SPR source, the maximum changed prob it can get
    int numHaps = this->genosInput.GetNumHaps();
    
    MarginalTree mTreeCurr;
    treeSearchCur.GetCurrTree(mTreeCurr);
    
    //double score0 = ScoreTree(mTreeCurr);
    //cout << "(INITIIAL tree: Re-computing log-likelihood: " << score0 << endl;
    
    // now start SPR local search
    int numIters = 0;
    
    //set<int> setDropped, setDroppedDest;
    //set<set<int> > setDroppedSrcSTsHist, setDroppedDestSTsHist;
    
    //bool fRetry = false;
    while(++numIters <= this->maxIters)
    {
        if( fVerbose )
        {
            //cout << "++FindOptMulti: iteration " << numIters << ", max prob: " << probMaxCurr << endl;
        }
        
        //vector<set<int> > listDescLavesTaxa;
        //mTreeCurr.ConsDecedentLeavesInfoLabels( listDescLavesTaxa );
        
        // find out all the clades of current ree
        vector<set<int> > listDescLaves;
        mTreeCurr.ConsDecedentLeavesInfo( listDescLaves );
        
//cout << "**** Current max prob: " << probMaxCurr << "    treeSearchCur: ";
//cout << mTreeCurr.GetNewickNoBrLen() << endl;
//mTreeCurr.Dump();
        // search for opt
        int nodeSPRSrc, nodeSPRDest, nodeLCA;
        std::vector<std::pair<int,double> > mapSPRProbs(2*numHaps);
        double probMaxStep = treeSearchCur.GetOpt1SPRMove3(nodeSPRSrc, nodeSPRDest, nodeLCA, &mapSPRProbs);
//cout << "Done with 1-SPR neighborhood search\n";
        //cout << "Set of normal SPR: probMaxStep: " << probMaxStep  << "\n";
        //for(unsigned int jjj=0; jjj<mapSPRProbs.size(); ++jjj)
        //{
        //    cout << "node: " << mapSPRProbs[jjj].first << ": likeli: " << mapSPRProbs[jjj].second << endl;
        //}
        
        // also get the ances special cases
        std::vector<std::pair<int,double> > mapSPRProbsAnc;
        int nodeSPRSrcAnc, nodeLCAAnc;
        double probMaxStepAnc = treeSearchCur.GetOpt1SPRMove3Ances( nodeSPRSrcAnc, nodeLCAAnc, mapSPRProbsAnc );
        //cout << "Set of ancestral SPR: probMaxStepAnc: " << probMaxStepAnc << "\n";
        //for(unsigned int jjj=0; jjj<mapSPRProbsAnc.size(); ++jjj)
        //{
        //    cout << "Ancestor: " << mapSPRProbsAnc[jjj].first << ": likeli: " << mapSPRProbsAnc[jjj].second << endl;
        //}
//exit(1);
        
//cout << "probMaxStep: " << probMaxStep << ", probMaxStepAnc: " << probMaxStepAnc << endl;
        
        // pick the overall better one
        double probMaxStepCombo = std::max(probMaxStep, probMaxStepAnc);
        
        
        // adjust
        //probMaxStep += sumLogprobAllele0;
        
        if( probMaxStepCombo < probMaxCurr + this->minProbInc)
        {
            //cout << "SPR: cannot improve.\n";
            break;
        }
        
        double probMaxCurrPre = probMaxCurr;
        probMaxCurr = probMaxStepCombo;
        if( probMaxStepAnc > probMaxStep )
        {
            nodeSPRSrc = nodeSPRSrcAnc;
            nodeSPRDest = nodeLCAAnc;
//cout << "Using the ancesteral version of SPR: nodeSPRSrcAnc=" << nodeSPRSrcAnc << ", nodeLCAAnc:" << nodeLCAAnc << endl;
        }
        // move to the next
        //cout << "Prior to SPR: tree = " << mTreeCurr.GetNewickNoBrLen() << endl;
        //mTreeCurr.Dump();
        //double scorep = ScoreTree(mTreeCurr);
        //cout << "(Before SPR: Re-computing log-likelihood: " << scorep << endl;
        //mTreeCurr.Dump();
        
        // keep track what STs have been pruned; no duplicate prunning of the same ST is allowed
        //set<int> setPrunedSTs;
        
        mTreeCurr.PerformSPR(nodeSPRSrc, nodeSPRDest);
        //setPrunedSTs.insert(nodeSPRSrc);
        //setDropped.insert(nodeSPRSrc);
//cout << "AFTER performing 1 SPR: tree is: " << mTreeCurr.GetNewickNoBrLen() << endl;
//mTreeCurr.Dump();
        //treeSearchCur.SetInitTree(mTreeCurr);
        //treeSearchCur.SetCurrMaxProb( probMaxCurr );
        
        //set< set<int> > setCladesInSPRs;
        //setCladesInSPRs.insert( listDescLaves[nodeSPRSrc] );
        //setCladesInSPRs.insert( listDescLaves[nodeSPRDest] );
        map<set<int>, int> mapSTToLeaves;
        mTreeCurr.ConsDecedentLeavesInfo2(mapSTToLeaves);
        
        // prepare to score tree
        ScistGenGenotypeMat &genosInputUse = const_cast<ScistGenGenotypeMat &>(this->genosInput);
        ScistPerfPhyProbOnTree phTreeCalc( genosInputUse, mTreeCurr );
        phTreeCalc.PrepareProbMaxQ();
        
        // now perform all possible SPR moves that is good
        int numExtraSPR = 0;
        
        //int numSPREval = 0;
        map<double, vector<pair<int,int> > > mapSPRProbs2;
        const double MIN_INC2 = log(1.01);
        for(int nodeaa = 0; nodeaa<(int)mapSPRProbs.size(); ++nodeaa)
        {
            if( mapSPRProbs[nodeaa].first < 0 || mapSPRProbs[nodeaa].second < probMaxCurrPre + MIN_INC2 )
            {
                continue;
            }
            //cout << "extra SPR: max prob " << mapSPRProbs[nodeaa].second << ", nodeaa:" << nodeaa << ", dest:" << mapSPRProbs[nodeaa].first  << endl;
            mapSPRProbs2[mapSPRProbs[nodeaa].second].push_back( std::make_pair( nodeaa, mapSPRProbs[nodeaa].first ) );
            //++numSPREval;
        }
        // add the second version
        //int numAdded2 = 0;
        for(int nodeaa = 0; nodeaa<(int)mapSPRProbsAnc.size(); ++nodeaa)
        {
            if( mapSPRProbsAnc[nodeaa].first < 0 || mapSPRProbsAnc[nodeaa].second < probMaxCurrPre + MIN_INC2 )
            {
                continue;
            }
            //cout << "extra SPR: max prob " << mapSPRProbs[nodeaa].second << ", nodeaa:" << nodeaa << ", dest:" << mapSPRProbs[nodeaa].first  << endl;
            mapSPRProbs2[mapSPRProbsAnc[nodeaa].second].push_back( std::make_pair( nodeaa, mapSPRProbsAnc[nodeaa].first ) );
            //++numSPREval;
            //++numAdded2;
        }
//cout << "Type-2 SPR added: " << numAdded2 << endl;

        if( fVerbose )
        {
            //cout << "Number of extra SPRs to evaluate: " << numSPREval << endl;
        }
//cout << "Current tree: " << mTreeCurr.GetNewickSorted(false) << endl;
//#if 0
        // now work the way back
        //int numExtraCand = 0, numExtraBad = 0;
        for( map<double, vector<pair<int,int> > > :: reverse_iterator it22 = mapSPRProbs2.rbegin(); it22 != mapSPRProbs2.rend(); ++it22 )
        {
            //if( it22->first <= probMaxCurrPre )
            //{
            //    break;
            //}
            
            for(vector<pair<int,int> > :: iterator it23 = it22->second.begin(); it23 != it22->second.end(); ++it23)
            {
                int nodeSrcStep = it23->first, nodeDestStep = it23->second;
                
                //if(setPrunedSTs.find(nodeSrcStep) != setPrunedSTs.end() )
                //{
                //    continue;
                //}
                
                // get their subtrees
                set<int> ssSrc = listDescLaves[nodeSrcStep];
                set<int> ssDest = listDescLaves[nodeDestStep];
                // find the STs for these in the current tree
                map<set<int>, int> :: iterator itt4 = mapSTToLeaves.find(ssSrc);
                if( itt4 == mapSTToLeaves.end() )
                {
                    continue;
                }
                map<set<int>, int> :: iterator itt5 = mapSTToLeaves.find(ssDest);
                if( itt5 == mapSTToLeaves.end() )
                {
                    continue;
                }
                int nodeSrc2 = itt4->second, nodeDest2 = itt5->second;
                if( mTreeCurr.GetSibling(nodeSrc2) == nodeDest2 )
                {
                    continue;
                }
//cout << "Evaluating source: " << nodeSrc2 << ", dest: " << nodeDest2 << ", source subtree: ";
//DumpIntSet(ssSrc);
//cout << "   dest subtree: ";
//DumpIntSet(ssDest);
                // now try to perform the rSPR (on a copy of current tree)
                MarginalTree mTreeCurrCopy = mTreeCurr;
                mTreeCurrCopy.PerformSPR(nodeSrc2, nodeDest2);
                double score2 = ScoreTree(mTreeCurrCopy, &phTreeCalc);
                //++numExtraCand;
                // TEST: use larger threshold
                const double THRES_MIN_INC = log(1.001);
                if( score2 > probMaxCurr + THRES_MIN_INC )
                {
                    // take it
                    if( fVerbose )
                    {
                        //cout << "----- Taking one extra SPR move: src subtree: likelihood improved to: " << score2 << endl;
                    }
                    ++numExtraSPR;
                    probMaxCurr = score2;
                    mTreeCurr.PerformSPR(nodeSrc2, nodeDest2);
                    //treeSearchCur.SetInitTree(mTreeCurr);
                    //treeSearchCur.SetCurrMaxProb( probMaxCurr );
                    mTreeCurr.ConsDecedentLeavesInfo2(mapSTToLeaves);
                    // testing: 01/30/24
                    //phTreeCalc.PrepareProbMaxQ();
                    
                    // remember we have performed SPR on this ST
                    //setPrunedSTs.insert(nodeSrcStep);
                    //setDropped.insert(nodeSrcStep);
#if 0
                    double scoreRecalc = ScoreTree(mTreeCurr);
                    YW_ASSERT_INFO( std::fabs(scoreRecalc-probMaxCurr) < 0.01, "Wrong in additional SPRs" );
                    //cout << "Re-computing log-likelihood: " << score << endl;
#endif
                    
//cout << "After extra SPR move: tree is: " << mTreeCurr.GetNewickSorted(false) << endl;
                }
                else
                {
                    //++numExtraBad;
                }
            }
        }
//cout << "numExtraCand: " << numExtraCand << ", numExtraBad: " << numExtraBad << endl;
        
        // update search
        treeSearchCur.SetInitTree(mTreeCurr);
        treeSearchCur.SetCurrMaxProb( probMaxCurr );
        
        if( fVerbose && numExtraSPR > 0 )
        {
            //cout << "Number of extra SPR moves taken: " << numExtraSPR << endl;
        }
//#endif
        // temp code
        //double score = ScoreTree(mTreeCurr);
        //cout << "Re-computing log-likelihood: " << score << endl;
        
        if( fVerbose )
        {
            cout << "++ 1-rSPR local search: likelihood improved to " << probMaxCurr << endl;
        }

//#endif
    }
    //cout << "Total number of rSPR local search steps: " << numIters << endl;
    //cout << "**** Maximum log-likelihood: " << probMaxCurr << endl;
    //cout << "Constructed single cell phylogeny: " << mTreeCurr.GetNewickNoBrLen() << endl;
    
    strTreeOpt = mTreeCurr.GetNewickNoBrLen();
    
    // now valid to make sure it is indeed the same prob
    //ScistPerfPhyGuideTree treeGuidThis;
    //treeGuidThis.InitDecAll(strTreeOpt);
    //ScistFastSPRLocalSearch treeSearchCur2(this->genosInput, treeGuidThis );
    //double probMaxCurr2 = treeSearchCur.CalcCurrMaxProb();
    //probMaxCurr2 += sumLogprobAllele0;
    //cout << "Recomputed log-likelihood: " << probMaxCurr2 << endl;
    
//cout << "FindOptMulti log-likelihood: " << probMaxCurr << ", tree: " << strTreeOpt << endl;
    
    return probMaxCurr;
}

void ScistFastSPRLocalSearchLoop :: FindSPRSrcLoserFrom( const std::map<int, double> &mapSrcScore, const std::set<int> &setAvoided, double frac, std::set<int> &setSrcLosers ) const
{
    // find the bottom frac (say 0.5) of spr candidates that are not faring well
    setSrcLosers.clear();
    vector<pair<double, int> > listPairs;
    for( std::map<int, double> :: const_iterator it = mapSrcScore.begin(); it != mapSrcScore.end(); ++it )
    {
        if( setAvoided.find(it->first) == setAvoided.end() )
        {
            listPairs.push_back(std::make_pair( it->second, it->first ));
        }
    }
    std::sort( listPairs.begin(), listPairs.end() );
    // put the bottom range
    int posLoser = (int)listPairs.size() * (1.0-frac);
    for(int i=0; i<posLoser; ++i)
    {
        setSrcLosers.insert( listPairs[i].second );
    }
}

double ScistFastSPRLocalSearchLoop :: ScoreTree(MarginalTree &mtree, ScistPerfPhyProbOnTree *pTreeCalc)
{
    if( this->numThreads > 1)
    {
        return ScoreTreeMulti(mtree, pTreeCalc);
    }
//cout << "Score tree: " << treeToScore.GetNewick() << endl;
    //set<pair<ScistPerfPhyCluster,ScistPerfPhyCluster> > setClusDone;
    //map<pair<ScistPerfPhyCluster,ScistPerfPhyCluster>, pair<ScistPerfPhyCluster,ScistPerfPhyCluster> > mapChangedClus;
    double res = 0.0;;
    //ScistGenGenotypeMat &genosInputUse = const_cast<ScistGenGenotypeMat &>(this->genosInput);
    //ScistPerfPhyProbOnTree probTree( genosInputUse, mtree );
    //probTree.PrepareProbMaxQ();
    
    // YW: do we really need the map as above? TBD
    for(int site=0; site<genosInput.GetNumSites(); ++site)
    {
//cout << "ScoreTree: site " << site << endl;
//cout << "Heterozygote clus: ";
//listClusMutsInputHetero[site].Dump();
//cout << "Homozygous clus: ";
//listClusMutsInputHomo[site].Dump();
        //pair<ScistPerfPhyCluster,ScistPerfPhyCluster> clusChanged;
        //double loglikeliSite = probTree.CalcProbMaxForSite( site, clusChanged.first, clusChanged.second );
        //double loglikeliSite = probTree.CalcProbMaxForSite2( site );
        double loglikeliSite = pTreeCalc->CalcProbMaxForSiteForTree( site, mtree );
        //mapChangedClus[ pp0 ] = clusChanged;
        //listChangedCluster.push_back(clusChanged);
        res += loglikeliSite;
        //setClusDone.insert( pp0 );
//cout << "site: " << site << ", loglikesite: " << loglikeliSite << endl;
//cout << "site prob: " << loglikeliSite << ": clusChanged: ";
//clusChanged.first.Dump();
//cout << "  and ";
//clusChanged.second.Dump();
    }
    return res;
}

static void UtilsScoreTree(ScistGenGenotypeMat *pMatIn, MarginalTree *ptree, ScistPerfPhyProbOnTree *pProbTree, int blkStart, int blkEnd, double *pProb)
{
//cout << "Score tree: " << treeToScore.GetNewick() << endl;
    //set<pair<ScistPerfPhyCluster,ScistPerfPhyCluster> > setClusDone;
    //map<pair<ScistPerfPhyCluster,ScistPerfPhyCluster>, pair<ScistPerfPhyCluster,ScistPerfPhyCluster> > mapChangedClus;
    double res = 0.0;;
    //ScistGenGenotypeMat &genosInputUse = *pMatIn;
    
    // YW: do we really need the map as above? TBD
    for(int site=blkStart; site<=blkEnd; ++site)
    {
//cout << "ScoreTree: site " << site << endl;
//cout << "Heterozygote clus: ";
//listClusMutsInputHetero[site].Dump();
//cout << "Homozygous clus: ";
//listClusMutsInputHomo[site].Dump();
        //pair<ScistPerfPhyCluster,ScistPerfPhyCluster> clusChanged;
        //double loglikeliSite = pProbTree->CalcProbMaxForSite( site, clusChanged.first, clusChanged.second );
        double loglikeliSite = pProbTree->CalcProbMaxForSiteForTree( site, *ptree );
        //mapChangedClus[ pp0 ] = clusChanged;
        //listChangedCluster.push_back(clusChanged);
        res += loglikeliSite;
        //setClusDone.insert( pp0 );
//cout << "site: " << site << ", loglikesite: " << loglikeliSite << endl;
//cout << "site prob: " << loglikeliSite << ": clusChanged: ";
//clusChanged.first.Dump();
//cout << "  and ";
//clusChanged.second.Dump();
    }
    *pProb = res;
}

double ScistFastSPRLocalSearchLoop :: ScoreTreeMulti(MarginalTree &mtree, ScistPerfPhyProbOnTree *pTreeCalc)
{
//cout << "ScoreTreeMulti\n";
    ScistGenGenotypeMat &genosInputUse = const_cast<ScistGenGenotypeMat &>(this->genosInput);
    //ScistPerfPhyProbOnTree probTree( genosInputUse, mtree );
    //probTree.PrepareProbMaxQ();
    
    int numItems = genosInputUse.GetNumSites();
    int blkSize = numItems/this->numThreads;
    int blkStart = 0;
    int numThreadsUse = this->numThreads;
    if( this->numThreads > genosInputUse.GetNumSites() )
    {
        numThreadsUse = genosInputUse.GetNumSites();
    }
    vector<thread *> listPtrThreads;
    vector<double> listProbOut( numThreadsUse );
    for(int t=0; t<numThreadsUse; ++t)
    {
        int blkEnd = blkStart + blkSize - 1;
        if( blkEnd >= numItems  || t == numThreadsUse-1)
        {
            blkEnd = numItems-1;
        }
        YW_ASSERT_INFO(blkEnd >= blkStart, "WRONG345");
        // start thread
        //thread *pthr = new thread(UtilsInitM0Tbl, this, &mtree, blkStart, blkEnd, &tblQ, &tblM0, &tblM0a, &(this->genosInput), &tblM1, &tblM2, &tblM10, &tblM8, &tblM9, &tblM11 );
        //thread *pthr = new thread(UtilsScoreTree, &genosInputUse, &probTree, blkStart, blkEnd, &listProbOut[t] );
        thread *pthr = new thread(UtilsScoreTree, &genosInputUse, &mtree, pTreeCalc, blkStart, blkEnd, &listProbOut[t] );
        listPtrThreads.push_back(pthr);
        blkStart += blkSize;
    }
    
    // wait to join
    for(unsigned int j=0; j<listPtrThreads.size(); ++j)
    {
      listPtrThreads[j]->join();
    }
    for(unsigned int j=0; j<listPtrThreads.size(); ++j)
    {
      delete listPtrThreads[j];
    }
    listPtrThreads.clear();
    
    double res = 0.0;
    for(unsigned int i=0; i<listProbOut.size(); ++i)
    {
        res += listProbOut[i];
    }
    return res;
}



// **************************************************************
// Fast NNI local search

ScistFastNNILocalSearch :: ScistFastNNILocalSearch(const ScistGenGenotypeMat &genosInputIn, const string &treeInitIn, int nt) : genosInput(genosInputIn), numThreads(nt)
{
//cout << "**ScistFastNNILocalSearch:Initial tree: " << strTree << endl;
    ReadinMarginalTreesNewickWLenString(treeInitIn, GetNumCells(), mtree);
    mtree.RearrangeParIncOrder();
    mtree.BuildDescendantInfo();
//cout << "Marginal tree: ";
//mtree.Dump();
    
    Init();
}


void ScistFastNNILocalSearch :: Init()
{
    if( this->numThreads > 1)
    {
        InitMT();
        return;
    }
//cout << "ScistFastNNILocalSearch: Init\n";
    // store values of Q table
    this->tblQ.clear();
    this->tblQ.resize( GetNumSites() );
    this->tblM0.resize(GetNumSites() );
    //this->tblQMax.resize( GetNumSites() );
    for(int s=0; s<GetNumSites(); ++s)
    {
        this->tblQ[s].resize( mtree.GetTotNodesNum() );
        this->tblM0[s].resize( mtree.GetTotNodesNum() );
    }
    //
    for(int site=0; site<GetNumSites(); ++site)
    {
//cout << "InitQTbl: site:" << site << endl;
        // do a bottom up
        vector<double> listNodeSplitProb;
        // init to be bad
        for(int node=0; node<mtree.GetTotNodesNum(); ++node)
        {
            listNodeSplitProb.push_back( -1.0*HAP_MAX_INT );
        }
        
    //cout << "CalcProbMaxForSiteHap: mtree: " << mtree.GetNewickSorted(false) << endl;
    //mtree.Dump();
        //this->tblQMax[site] = MAX_NEG_DOUBLE_VAL;
        
        for(int node=0; node<mtree.GetTotNodesNum(); ++node)
        {
//cout << "node " << node << endl;
            double logpStep = 0.0;
            if( mtree.IsLeaf(node) )
            {
//cout << "Leaf label: " << mtree.GetLabel(node) << endl;
                // a single leaf in the split
                int lvlbl = mtree.GetLabel(node)-1;
                //cout << "Leaf: " << lvlbl << endl;
                double p0 = this->genosInput.GetGenotypeProbAllele0At(lvlbl, site);
                if( p0 < YW_VERY_SMALL_FRACTION)
                {
                    p0 = YW_VERY_SMALL_FRACTION;
                }
                else if( p0 > 1.0-YW_VERY_SMALL_FRACTION)
                {
                    p0 = 1.0-YW_VERY_SMALL_FRACTION;
                }
                logpStep = log((1-p0)/p0);
                
                this->tblM0[site][node] = logpStep;
//cout << "Set leaf " << node << " log prob to: " << logpStep << ", p0=" << p0 << endl;
            }
            else
            {
                // get the two children and add them up
                int childLeft = mtree.GetLeftDescendant(node);
                int childRight = mtree.GetRightDescendant(node);
    //cout << "node: " << node << ", childLeft: " << childLeft << ", childRight: " << childRight << endl;
                //cout << "childLeft: " << childLeft << ", right: " << childRight << endl;
                
                YW_ASSERT_INFO( listNodeSplitProb[childLeft] > -1.0*HAP_MAX_INT, "Bad left-c" );
                YW_ASSERT_INFO( listNodeSplitProb[childRight] > -1.0*HAP_MAX_INT, "Bad right1-c" );
                logpStep = listNodeSplitProb[childLeft] + listNodeSplitProb[childRight];
                
                this->tblM0[site][node] = logpStep;
                this->tblM0[site][node] = std::max(this->tblM0[site][node], this->tblM0[site][childLeft] );
                this->tblM0[site][node] = std::max(this->tblM0[site][node], this->tblM0[site][childRight] );
//cout << "node: " << node << ", prob: " << logpStep << ", childLeft: " << childLeft << ", childRight: " << childRight << endl;
            }
//cout << "log prob: " << logpStep << " for node: " << node << endl;
            listNodeSplitProb[node] = logpStep;
            
            //if( logpStep > this->tblQMax[site] )
            //{
            //    this->tblQMax[site] = logpStep;
            //}
        }
        
        // save the probs
        this->tblQ[site] = listNodeSplitProb;
    }
    
    // now init M1 table.
    this->tblM1.resize(GetNumSites());
    for(int site=0; site<GetNumSites(); ++site)
    {
//cout << "Init tblM1: " << site << endl;
        // start from root down to leaf
        // for root, set to the smallest log-prob
        this->tblM1[site].resize(GetNumNodesInTree());
        tblM1[site][mtree.GetRoot()] = MAX_NEG_DOUBLE_VAL;
        for(int node=GetNumNodesInTree()-2; node >=0; --node )
        {
            int nodeSib = mtree.GetSibling(node);
            int nodePar = mtree.GetParent(node);
            YW_ASSERT_INFO(nodeSib >= 0, "Fail to find sibling");
            YW_ASSERT_INFO(nodePar >= 0, "Fail to find parent");
            tblM1[site][node] = std::max( this->tblQ[site][nodePar], tblM0[site][nodeSib] );
            tblM1[site][node] = std::max( this->tblM1[site][node], tblM1[site][nodePar]  );
        }
    }
#if 0
    cout << "****** TblQ: \n";
    for(int site=0; site<GetNumSites(); ++site)
    {
        DumpDoubleVec(this->tblQ[site]);
    }
    cout << "****** TblM0: \n";
    for(int site=0; site<GetNumSites(); ++site)
    {
        DumpDoubleVec(this->tblM0[site]);
    }
    cout << "****** TblM1: \n";
    for(int site=0; site<GetNumSites(); ++site)
    {
        DumpDoubleVec(this->tblM1[site]);
    }
#endif
}

static void UtilsInitMTProc(MarginalTree *ptree, const ScistGenGenotypeMat *pgenosInput, int siteStart, int siteEnd, int numSites, vector<std::vector<double> > *pTblQ, vector<std::vector<double> > *pTblM0, vector<std::vector<double> > *pTblM1 )
{
    //
    for(int site=siteStart; site<=siteEnd; ++site)
    {
//cout << "InitQTbl: site:" << site << endl;
        // do a bottom up
        vector<double> listNodeSplitProb;
        // init to be bad
        for(int node=0; node<ptree->GetTotNodesNum(); ++node)
        {
            listNodeSplitProb.push_back( -1.0*HAP_MAX_INT );
        }
        
    //cout << "CalcProbMaxForSiteHap: mtree: " << mtree.GetNewickSorted(false) << endl;
    //mtree.Dump();
        //this->tblQMax[site] = MAX_NEG_DOUBLE_VAL;
        
        for(int node=0; node<ptree->GetTotNodesNum(); ++node)
        {
//cout << "node " << node << endl;
            double logpStep = 0.0;
            if( ptree->IsLeaf(node) )
            {
//cout << "Leaf label: " << mtree.GetLabel(node) << endl;
                // a single leaf in the split
                int lvlbl = ptree->GetLabel(node)-1;
                //cout << "Leaf: " << lvlbl << endl;
                double p0 = pgenosInput->GetGenotypeProbAllele0At(lvlbl, site);
                if( p0 < YW_VERY_SMALL_FRACTION)
                {
                    p0 = YW_VERY_SMALL_FRACTION;
                }
                else if( p0 > 1.0-YW_VERY_SMALL_FRACTION)
                {
                    p0 = 1.0-YW_VERY_SMALL_FRACTION;
                }
                logpStep = log((1-p0)/p0);
                
                (*pTblM0)[site][node] = logpStep;
//cout << "Set leaf " << node << " log prob to: " << logpStep << ", p0=" << p0 << endl;
            }
            else
            {
                // get the two children and add them up
                int childLeft = ptree->GetLeftDescendant(node);
                int childRight = ptree->GetRightDescendant(node);
    //cout << "node: " << node << ", childLeft: " << childLeft << ", childRight: " << childRight << endl;
                //cout << "childLeft: " << childLeft << ", right: " << childRight << endl;
                
                YW_ASSERT_INFO( listNodeSplitProb[childLeft] > -1.0*HAP_MAX_INT, "Bad left-d" );
                YW_ASSERT_INFO( listNodeSplitProb[childRight] > -1.0*HAP_MAX_INT, "Bad right1-d" );
                logpStep = listNodeSplitProb[childLeft] + listNodeSplitProb[childRight];
                
                (*pTblM0)[site][node] = logpStep;
                (*pTblM0)[site][node] = std::max( (*pTblM0)[site][node], (*pTblM0)[site][childLeft] );
                (*pTblM0)[site][node] = std::max((*pTblM0)[site][node], (*pTblM0)[site][childRight] );
//cout << "node: " << node << ", prob: " << logpStep << ", childLeft: " << childLeft << ", childRight: " << childRight << endl;
            }
//cout << "log prob: " << logpStep << " for node: " << node << endl;
            listNodeSplitProb[node] = logpStep;
            
            //if( logpStep > this->tblQMax[site] )
            //{
            //    this->tblQMax[site] = logpStep;
            //}
        }
        
        // save the probs
        (*pTblQ)[site] = listNodeSplitProb;
    }
    
    // now init M1 table.
    for(int site=siteStart; site<=siteEnd; ++site)
    {
//cout << "Init tblM1: " << site << endl;
        // start from root down to leaf
        // for root, set to the smallest log-prob
        (*pTblM1)[site][ptree->GetRoot()] = MAX_NEG_DOUBLE_VAL;
        for(int node=ptree->GetTotNodesNum()-2; node >=0; --node )
        {
            int nodeSib = ptree->GetSibling(node);
            int nodePar = ptree->GetParent(node);
            YW_ASSERT_INFO(nodeSib >= 0, "Fail to find sibling");
            YW_ASSERT_INFO(nodePar >= 0, "Fail to find parent");
            (*pTblM1)[site][node] = std::max( (*pTblQ)[site][nodePar], (*pTblM0)[site][nodeSib] );
            (*pTblM1)[site][node] = std::max( (*pTblM1)[site][node], (*pTblM1)[site][nodePar]  );
        }
    }
}

void ScistFastNNILocalSearch :: InitMT()
{
//cout << "ScistFastNNILocalSearch: Init\n";
    // store values of Q table
    this->tblQ.clear();
    this->tblQ.resize( GetNumSites() );
    this->tblM0.resize(GetNumSites() );
    this->tblM1.resize(GetNumSites());
    //this->tblQMax.resize( GetNumSites() );
    for(int s=0; s<GetNumSites(); ++s)
    {
        this->tblQ[s].resize( mtree.GetTotNodesNum() );
        this->tblM0[s].resize( mtree.GetTotNodesNum() );
        this->tblM1[s].resize( mtree.GetTotNodesNum() );
    }
    int numThreadsUse = this->numThreads;
    if( numThreadsUse > GetNumSites() )
    {
        numThreadsUse = GetNumSites();
    }
    int blkSize = GetNumSites()/numThreadsUse;
    int blkStart = 0;
    vector<thread *> listPtrThreads;
    for(int t=0; t<numThreadsUse; ++t)
    {
        int blkEnd = blkStart + blkSize - 1;
        if( blkEnd >= GetNumSites()  || t == numThreadsUse-1)
        {
            blkEnd = GetNumSites()-1;
        }
        YW_ASSERT_INFO(blkEnd >= blkStart, "WRONG345");
        // start thread
        thread *pthr = new thread(UtilsInitMTProc, &mtree, &genosInput, blkStart, blkEnd, GetNumSites(), &tblQ, &tblM0, &tblM1);
        listPtrThreads.push_back(pthr);
        blkStart += blkSize;
    }
    
    // wait to join
    for(unsigned int j=0; j<listPtrThreads.size(); ++j)
    {
      listPtrThreads[j]->join();
    }
    for(unsigned int j=0; j<listPtrThreads.size(); ++j)
    {
      delete listPtrThreads[j];
    }
    listPtrThreads.clear();

#if 0
    cout << "****** TblQ: \n";
    for(int site=0; site<GetNumSites(); ++site)
    {
        DumpDoubleVec(this->tblQ[site]);
    }
    cout << "****** TblM0: \n";
    for(int site=0; site<GetNumSites(); ++site)
    {
        DumpDoubleVec(this->tblM0[site]);
    }
    cout << "****** TblM1: \n";
    for(int site=0; site<GetNumSites(); ++site)
    {
        DumpDoubleVec(this->tblM1[site]);
    }
#endif
}


double ScistFastNNILocalSearch :: CalcCurrMaxProb(std::string &strNNIBestTree) const
{
    //if(this->numThreads > 1 && GetNumNodesInTree() > this->numThreads )
    if(this->numThreads > 1 )
    {
        return CalcCurrMaxProbMT(strNNIBestTree);
    }
    
    //
    double logProbMaxNNI = MAX_NEG_DOUBLE_VAL;
    int nodeNNI = -1, nodeNNIDest = -1;;
    for(int node=0; node<GetNumNodesInTree(); ++node)
    {
        if( node == mtree.GetRoot() )
        {
            continue;
        }
        int nodePar = mtree.GetParent(node);
        if( nodePar == mtree.GetRoot() )
        {
            continue;
        }
        int nodeParSib =mtree.GetSibling(nodePar);
        int nodeSib = mtree.GetSibling(node);
        YW_ASSERT_INFO(nodeParSib >= 0, "Fail to find sibling of parent node");
        YW_ASSERT_INFO(nodeSib >=0, "Fail to find sibling of node");
        //
//cout << "--NNI: node:" << node << ", sibliing: " << nodeParSib << endl;
        double logProbNewNode = 0.0;
        for(int site=0; site<GetNumSites(); ++site)
        {
            double probMaxNode = this->tblQ[site][node] + this->tblQ[site][nodeParSib];
            double probSub1 = this->tblM0[site][node];
            double probSub = std::max(probSub1, this->tblM0[site][nodeSib]);
            double maxProbStep = std::max(probMaxNode, probSub);
            //double maxProbStep2 = std::max(maxProbStep, this->tblM0[site][nodeParSib]);
            double probPar = this->tblM1[site][nodePar];
            maxProbStep = std::max(maxProbStep, probPar);
            
            // finally if the result is not better than all-0, use all-0. That is, if max is smaller than 0, use 0
            if( 0.0 > maxProbStep )
            {
//cout << "**** TAKING all-0 option at site: " << site << endl;
                maxProbStep = 0.0;
            }
            
            logProbNewNode += maxProbStep;
//cout << "max prob at site: " << site << ": " << maxProbStep << "  probMaxNode:" << probMaxNode << ", probSub:" << probSub << ", probPar:" << probPar << ", max BEFORE NNI: " << this->tblM0[site][mtree.GetRoot()] << endl;
        }
        if( logProbMaxNNI + DEF_THRES_PROB_INC  < logProbNewNode )
        {
//cout << "Higher NNI: nodeNNI: " << node << ", NNI dest:" << nodeParSib << ", logProbNewNode: " << logProbNewNode << endl;
            logProbMaxNNI = logProbNewNode;
            nodeNNI = node;
            nodeNNIDest = nodeParSib;
        }
    }
    
    // construct current tree
    YW_ASSERT_INFO(nodeNNI >= 0 && nodeNNIDest >= 0, "Fail to find NNI" );
    MarginalTree mTreeCurr = this->mtree;
    mTreeCurr.PerformSPR(nodeNNI, nodeNNIDest);
    strNNIBestTree = mTreeCurr.GetNewickNoBrLen();
    
    // add on the overall prob of zero for all genotypes
    double prob0All = this->genosInput.CalcSumLogProbAllele0All();
//cout << "ProbZeroAll: " << prob0All << endl;
    double logProbMaxNNIRes = logProbMaxNNI + prob0All;
//cout << "logProbMaxNNIRes: " << logProbMaxNNIRes << ", nodeNNI:" << nodeNNI << ", nodeNNIDest:" << nodeNNIDest << endl;
    
// test code
//double probM0Orig = 0.0;
//for(int site=0; site<GetNumSites(); ++site)
//{
//probM0Orig += this->tblM0[site][mtree.GetRoot()];
//}
//probM0Orig += prob0All;
//cout << "CHECKING: log-likelihood computed by M0 table on original tree: " << probM0Orig << endl;
    
    return logProbMaxNNIRes;
}


static void UtilsFindBestNNIMove( const MarginalTree *ptree, int tindex, int siteStart, int siteEnd, int numNodes, const vector<std::vector<double> > *pTblQ, const vector<std::vector<double> > *pTblM0, const vector<std::vector<double> > *pTblM1, vector<double> *pListProb)
{
    //double logProbMaxNNI = MAX_NEG_DOUBLE_VAL;
    //int nodeNNISrc = -1;
    //int nodeNNIDest = -1;
    for(int node=0; node<numNodes; ++node)
    {
        if( node == ptree->GetRoot() )
        {
            continue;
        }
        int nodePar = ptree->GetParent(node);
        if( nodePar == ptree->GetRoot() )
        {
            continue;
        }
        int nodeParSib =ptree->GetSibling(nodePar);
        int nodeSib = ptree->GetSibling(node);
        YW_ASSERT_INFO(nodeParSib >= 0, "Fail to find sibling of parent node");
        YW_ASSERT_INFO(nodeSib >=0, "Fail to find sibling of node");
        //
//cout << "--NNI: node:" << node << ", sibliing: " << nodeParSib << endl;
        double logProbNewNode = 0.0;
        for(int site=siteStart; site<=siteEnd; ++site)
        {
            double probMaxNode = (*pTblQ)[site][node] + (*pTblQ)[site][nodeParSib];
            double probSub1 = (*pTblM0)[site][node];
            double probSub = std::max(probSub1, (*pTblM0)[site][nodeSib]);
            double maxProbStep = std::max(probMaxNode, probSub);
            double probPar = (*pTblM1)[site][nodePar];
            maxProbStep = std::max(maxProbStep, probPar);
            
            // finally if the result is not better than all-0, use all-0. That is, if max is smaller than 0, use 0
            if( 0.0 > maxProbStep )
            {
//cout << "**** TAKING all-0 option at site: " << site << endl;
                maxProbStep = 0.0;
            }
            
            logProbNewNode += maxProbStep;
//cout << "max prob at site: " << site << ": " << maxProbStep << "  probMaxNode:" << probMaxNode << ", probSub:" << probSub << ", probPar:" << probPar << ", max BEFORE NNI: " << this->tblM0[site][mtree.GetRoot()] << endl;
        }
        (*pListProb)[node] = logProbNewNode;
    }
}

double ScistFastNNILocalSearch :: CalcCurrMaxProbMT(std::string &strNNIBestTree) const
{
    // split nodes into buckets, each for a thread
    int numThreadsUse = this->numThreads;
    if( GetNumSites() < numThreadsUse )
    {
        numThreadsUse = GetNumSites();
    }
    vector<vector<double> > listNNIProbs(numThreadsUse);
    for(int i=0; i<numThreadsUse; ++i)
    {
        listNNIProbs[i].resize( GetNumNodesInTree() );
        for(int j=0; j<GetNumNodesInTree(); ++j)
        {
            listNNIProbs[i][j] = 0.0;
        }
    }
    int blkSize = GetNumSites()/numThreadsUse;
    int blkStart = 0;
    vector<thread *> listPtrThreads;
    for(int t=0; t<numThreadsUse; ++t)
    {
        int blkEnd = blkStart + blkSize - 1;
        if( blkEnd >= GetNumSites()  || t == numThreadsUse-1)
        {
            blkEnd = GetNumSites()-1;
        }
        YW_ASSERT_INFO(blkEnd >= blkStart, "WRONG345");
        // start thread
        thread *pthr = new thread(UtilsFindBestNNIMove, &mtree, t, blkStart, blkEnd, GetNumNodesInTree(), &tblQ, &tblM0, &tblM1, &listNNIProbs[t] );
        listPtrThreads.push_back(pthr);
        blkStart += blkSize;
    }
    
    // wait to join
    for(unsigned int j=0; j<listPtrThreads.size(); ++j)
    {
      listPtrThreads[j]->join();
    }
    for(unsigned int j=0; j<listPtrThreads.size(); ++j)
    {
      delete listPtrThreads[j];
    }
    listPtrThreads.clear();
    
    // find the best NNI move
    double logProbMaxNNI = MAX_NEG_DOUBLE_VAL;
    int nodeNNI = -1;
    int nodeNNIDest = -1;
    // sum over the main
    vector<double> listNNIProbsSum( GetNumNodesInTree() );
    for(int node=0; node<GetNumNodesInTree(); ++node)
    {
        if( node == this->mtree.GetRoot() )
        {
            continue;
        }
        int nodePar = this->mtree.GetParent(node);
        if( nodePar == this->mtree.GetRoot() )
        {
            continue;
        }
        int nodeParSib =this->mtree.GetSibling(nodePar);
        
        listNNIProbsSum[node] = 0;
        for(int t=0; t<numThreadsUse; ++t )
        {
            listNNIProbsSum[node] += listNNIProbs[t][node];
        }
        if( listNNIProbsSum[node] > logProbMaxNNI + DEF_THRES_PROB_INC )
        //if( listNNIProbsSum[node] > logProbMaxNNI  )
        {
            logProbMaxNNI = listNNIProbsSum[node];
            nodeNNI = node;
            nodeNNIDest = nodeParSib;
        }
    }
    
    // construct current tree
    YW_ASSERT_INFO(nodeNNI >= 0 && nodeNNIDest >= 0, "Fail to find NNI" );
    MarginalTree mTreeCurr = this->mtree;
    mTreeCurr.PerformSPR(nodeNNI, nodeNNIDest);
    strNNIBestTree = mTreeCurr.GetNewickNoBrLen();
    
    // add on the overall prob of zero for all genotypes
    double prob0All = this->genosInput.CalcSumLogProbAllele0All();
//cout << "ProbZeroAll: " << prob0All << endl;
    double logProbMaxNNIRes = logProbMaxNNI + prob0All;
//cout << "logProbMaxNNIRes: " << logProbMaxNNIRes << ", nodeNNI:" << nodeNNI << ", nodeNNIDest:" << nodeNNIDest << endl;
//cout << "CalcCurrMaxProbMT: \n";
    return logProbMaxNNIRes;
}

int ScistFastNNILocalSearch :: GetNumSites() const
{
    return genosInput.GetNumSites();
}
int ScistFastNNILocalSearch :: GetNumCells() const
{
    return genosInput.GetNumHaps();
}
int ScistFastNNILocalSearch :: GetNumNodesInTree() const
{
    return 2*GetNumCells()-1;
}

// **************************************************************
// Fast re-rooting local search

ScistFastRerootLocalSearch :: ScistFastRerootLocalSearch(const ScistGenGenotypeMat &genosInputIn, const string &strTreeInitIn, int nt) : genosInput(genosInputIn),
    strTreeInit(strTreeInitIn), numThreads(nt)
{
    Init();
}

double ScistFastRerootLocalSearch :: CalcCurrMaxProb() const
{
    double probMaxAll = 0.0;
    for(int site=0; site<GetNumSites(); ++site)
    {
        double val = this->tblM0[site][this->mtree.GetRoot()];
        if( val < 0.0 )
        {
            // if no best cut, then we just choose to not to cut
            val = 0.0;
        }
        probMaxAll += val;
    }
    return probMaxAll;
}

void ScistFastRerootLocalSearch :: GetCurrTree(MarginalTree &treeOut) const
{
    treeOut = mtree;
}


double ScistFastRerootLocalSearch :: GetOptSubtreeReroot( int &nodeSTRootNew, int &nodeST )
{
    // collect candidates
    CollectRerootCands();
    
    if( this->numThreads > 1 && GetNumNodesInTree() > this->numThreads  )
    {
        //
        return GetOptSubtreeRerootMT(nodeSTRootNew, nodeST);
    }
//cout << "GetOptSubtreeReroot: sum of prob 0: " << this->sumLogprobAllele0 << endl;
//cout << "tree: ";
//mtree.Dump();

    // nodeST: the subtree; nodeSTRootNew: a descedant node of nodeST which is the new root
    // best way of rotating the subtree at nodeST so that the resulting Q is max
    vector<pair<int,int> > listRerootCands = this->listCandsReroot;

    vector<double> listProbs( listRerootCands.size() );
    for(unsigned int i=0; i< listRerootCands.size(); ++i)
    {
        //
        int nodeSub = listRerootCands[i].first, node = listRerootCands[i].second;
        int rootIndex = GetAncesIndex(node, mtree.GetRoot());
        int rootSubIndex = GetAncesIndex(nodeSub, mtree.GetRoot());
        int childRootOtherSide = -1;
        if( node != mtree.GetRoot() )
        {
            childRootOtherSide = mtree.GetLeftDescendant(mtree.GetRoot());
            if( this->mapListAnces[node].find(childRootOtherSide) != this->mapListAnces[node].end() )
            {
                childRootOtherSide = mtree.GetRightDescendant(mtree.GetRoot());
            }
        }
        int nodeSubNodeIndex = GetAncesIndex(nodeSub, node);
        
        //
        double probMaxStep = this->sumLogprobAllele0;
        for(int site=0; site<GetNumSites(); ++site)
        {
            double probs = CalcMaxQForReroot(site, node, nodeSub, rootIndex, rootSubIndex, childRootOtherSide, nodeSubNodeIndex);
//cout << "Max Q value for reroot: site:" << site << ", node:" << node << ", nodeSub:" << nodeSub << " is " << probs << endl;
            probMaxStep += probs;
        }
//cout << "noder: " << noder << ", nodeu: " << nodeu << ", nodew: " << nodew << ", probMaxStep: " << probMaxStep << endl;
        listProbs[i] = probMaxStep;
//cout << "Reroot: node:" << node << ", nodeSub:" << nodeSub << " --- max Q value over ALL sites: " << probMaxStep << endl;
    }
    
    // consider all rSPR move
    double probMax = MAX_NEG_DOUBLE_VAL;
    double THRES_MIN_INC = 0.0000001;
    nodeSTRootNew = -1;
    nodeST = -1;
    // find the best move
    for(unsigned int i=0; i<listRerootCands.size(); ++i)
    {
        double probMaxStep = listProbs[i];
        int nodeu = listRerootCands[i].first;
        int noder = listRerootCands[i].second;
        if( probMaxStep > probMax + THRES_MIN_INC )
        {
            probMax = probMaxStep;
            nodeSTRootNew = nodeu;
            nodeST = noder;
        }
    }
//cout << "GetOptSubtreeReroot: listProbs: \n";
//DumpDoubleVec(listProbs);
    // if fail to find SPR, just return the smallest
    if( nodeSTRootNew < 0 || nodeST < 0 )
    {
        cout << "Warning: fail to find re-rootiing move. \n";
        return MAX_NEG_DOUBLE_VAL;
    }
    
//cout << "^^^^^^ Max prob of one rerooting move: " << probMax << ", subtree rooted at:" << nodeST << ", to a new node (within the subtree):" << nodeSTRootNew << endl;
    //
    return probMax;
}
    
static void UtilFastRerootSearch2( ScistFastRerootLocalSearch *pthis, int tindex, const std::vector<std::pair<int, int> > *plistRerootCand, unsigned int posStart, unsigned int posEnd, vector<double> *pListProb, vector<int> *pListNodes, vector<int> *pListNodeSTs, double sumLogprobAllele0, MarginalTree *ptree )
{
    //double probMax = MAX_NEG_DOUBLE_VAL;
    //int nodeBest = -1;
    //int nodeSTBest = -1;
    
    for(unsigned int i=posStart; i<= posEnd; ++i)
    {
        //
        double probMaxStep = sumLogprobAllele0;
        int nodeSubStep = (*plistRerootCand)[i].first, nodeStep = (*plistRerootCand)[i].second;
        int rootIndex = pthis->GetAncesIndex(nodeStep, ptree->GetRoot());
        int rootSubIndex = pthis->GetAncesIndex(nodeSubStep, ptree->GetRoot());
        int childRootOtherSide = -1;
        if( nodeStep != ptree->GetRoot() )
        {
            childRootOtherSide = ptree->GetLeftDescendant(ptree->GetRoot());
            //if( pthis->mapListAnces[node].find(childRootOtherSide) != this->mapListAnces[node].end() )
            if( pthis->IsNodeAncestralTo( childRootOtherSide, nodeStep ) )
            {
                childRootOtherSide = ptree->GetRightDescendant(ptree->GetRoot());
            }
        }
        int nodeSubNodeIndex = pthis->GetAncesIndex(nodeSubStep, nodeStep);
        for(int site=0; site<pthis->GetNumSites(); ++site)
        {
            double probs = pthis->CalcMaxQForReroot(site, nodeStep, nodeSubStep, rootIndex, rootSubIndex, childRootOtherSide, nodeSubNodeIndex );
//cout << "Max Q value for SPR: site:" << site << ", noder:" << noder << ", nodeu:" << nodeu << ", nodew:" << nodew << " is " << probs << endl;
            probMaxStep += probs;
        }
        (*pListProb)[i] = probMaxStep;
    }
}

double ScistFastRerootLocalSearch :: GetOptSubtreeRerootMT( int &nodeSTRootNew, int &nodeST )
{
    // do multi-threading
//cout << "Start multi-threading...\n";
    
    // nodeST: the subtree; nodeSTRootNew: a descedant node of nodeST which is the new root
    // best way of rotating the subtree at nodeST so that the resulting Q is max
    vector<pair<int,int> > listRerootCands = this->listCandsReroot;

    int numThreadsUse = numThreads;
    if( (int)listRerootCands.size() < numThreadsUse )
    {
        numThreadsUse = (int)listRerootCands.size();
    }
    
//cout << "Number of re-rooting candidates: " << listRerootCands.size() << endl;
    
    vector<double> listProbs(listRerootCands.size());
    //vector<double> listProbs(numThreadsUse);
    vector<int> listNodes(numThreadsUse);
    vector<int> listNodeSTs(numThreadsUse);
    int numItems = listRerootCands.size();
    int blkSize = numItems/numThreadsUse;
    int blkStart = 0;
    vector<thread *> listPtrThreads;
    for(int t=0; t<numThreadsUse; ++t)
    {
        int blkEnd = blkStart + blkSize - 1;
        if( blkEnd >= numItems  || t == numThreadsUse-1)
        {
            blkEnd = numItems-1;
        }
        YW_ASSERT_INFO(blkEnd >= blkStart, "WRONG345");
        // start thread
        thread *pthr = new thread(UtilFastRerootSearch2, this, t, &listRerootCands, blkStart, blkEnd, &listProbs, &listNodes, &listNodeSTs, this->sumLogprobAllele0, &mtree );
        listPtrThreads.push_back(pthr);
        blkStart += blkSize;
    }
    
    // wait to join
    for(unsigned int j=0; j<listPtrThreads.size(); ++j)
    {
      listPtrThreads[j]->join();
    }
    for(unsigned int j=0; j<listPtrThreads.size(); ++j)
    {
      delete listPtrThreads[j];
    }
    listPtrThreads.clear();
    
    // consider all rSPR move
    double probMax = MAX_NEG_DOUBLE_VAL;
    double THRES_MIN_INC = 0.0000001;
    nodeSTRootNew = -1;
    nodeST = -1;
    // find the best move
    for(unsigned int i=0; i<listProbs.size(); ++i)
    {
        double probMaxStep = listProbs[i];
        int node = listRerootCands[i].second;
        int nodes = listRerootCands[i].first;
        //int node = listNodes[i];
        //int nodes = listNodeSTs[i];
        if( probMaxStep > probMax + THRES_MIN_INC )
        {
            probMax = probMaxStep;
            nodeSTRootNew = nodes;
            nodeST = node;
        }
    }
//cout << "GetOptSubtreeRerootMT: listProbs: \n";
//DumpDoubleVec(listProbs);
    
    // if fail to find SPR, just return the smallest
    if( nodeSTRootNew < 0 || nodeST < 0 )
    {
        cout << "Warning: fail to find re-rootiing move. \n";
        return MAX_NEG_DOUBLE_VAL;
    }
    
//cout << "^^^^^^ Max prob of one rerooting move: " << probMax << ", subtree rooted at:" << nodeST << ", to a new node (within the subtree):" << nodeSTRootNew << endl;
    //
    return probMax;
}


// initialization code
void ScistFastRerootLocalSearch :: Init()
{
//cout << "ScistFastSPRLocalSearch :: Init()\n";
    string strTree = this->strTreeInit;
//cout << "Initial tree: " << strTree << endl;
    ReadinMarginalTreesNewickWLenString(strTree, GetNumCells(), mtree);
    //cout << "Marginal tree: ";
    //mtree.Dump();
    
    SetInitTree(mtree);
}

void ScistFastRerootLocalSearch :: SetInitTree(const MarginalTree &mtreeInit)
{
    this->mtree = mtreeInit;
    // initialize ancestral relations
    mapListAnces.clear();
    mapListAncesLookup.clear();
    mapListAnces.resize( mtree.GetTotNodesNum() );
    mapListAncesLookup.resize( mtree.GetTotNodesNum() );
    for(int node=0; node<mtree.GetTotNodesNum(); ++node)
    {
     int nodeAnces = node;
     set<int> ss;
     int ancesIndex = 0;
     while(true)
     {
         ss.insert(nodeAnces);
         mapListAncesLookup[node][ nodeAnces ] = ancesIndex++;
         if( mtree.GetRoot() == nodeAnces)
         {
             break;
         }
         else
         {
             nodeAnces = mtree.GetParent(nodeAnces);
         }
     }
     mapListAnces[node] = ss;
    }

//cout << "Now initialize tables\n";
    InitQTbl();
    InitM0Tbl();
    InitM1Tbl();
    InitM2Tbl();
    // Init tblM4 first
    InitM4Tbl();
    
    // init prob
    this->sumLogprobAllele0 = this->genosInput.CalcSumLogProbAllele0All();
    this->logProbCurr = CalcCurrMaxProb();
    this->logProbCurr += this->sumLogprobAllele0;
//cout << "Done ScistFastRerootLocalSearch: SetInitTree\n";
}
void ScistFastRerootLocalSearch :: InitQTbl()
{
    this->tblQ.clear();
    this->tblQ.resize( GetNumSites() );
    for(int s=0; s<GetNumSites(); ++s)
    {
        this->tblQ[s].resize( mtree.GetTotNodesNum() );
    }
    //
    for(int site=0; site<GetNumSites(); ++site)
    {
//cout << "InitQTbl: site:" << site << endl;
        // do a bottom up
        vector<double> listNodeSplitProb;
        // init to be bad
        for(int node=0; node<mtree.GetTotNodesNum(); ++node)
        {
            listNodeSplitProb.push_back( -1.0*HAP_MAX_INT );
        }
        
    //cout << "CalcProbMaxForSiteHap: mtree: " << mtree.GetNewickSorted(false) << endl;
    //mtree.Dump();
        
        for(int node=0; node<mtree.GetTotNodesNum(); ++node)
        {
            //cout << "node " << node << endl;
            double logpStep = 0.0;
            if( mtree.IsLeaf(node) )
            {
                // a single leaf in the split
                int lvlbl = mtree.GetLabel(node)-1;
                //cout << "Leaf: " << lvlbl << endl;
                double p0 = this->genosInput.GetGenotypeProbAllele0At(lvlbl, site);
                if( p0 < YW_VERY_SMALL_FRACTION)
                {
                    p0 = YW_VERY_SMALL_FRACTION;
                }
                else if( p0 > 1.0-YW_VERY_SMALL_FRACTION)
                {
                    p0 = 1.0-YW_VERY_SMALL_FRACTION;
                }
                logpStep = log((1-p0)/p0);
    //cout << "Set leaf " << node << " log prob to: " << logpStep << ", p0=" << p0 << endl;
            }
            else
            {
                // get the two children and add them up
                int childLeft = mtree.GetLeftDescendant(node);
                int childRight = mtree.GetRightDescendant(node);
    //cout << "node: " << node << ", childLeft: " << childLeft << ", childRight: " << childRight << endl;
                //cout << "childLeft: " << childLeft << ", right: " << childRight << endl;
                
                YW_ASSERT_INFO( listNodeSplitProb[childLeft] > -1.0*HAP_MAX_INT, "Bad left-a" );
                YW_ASSERT_INFO( listNodeSplitProb[childRight] > -1.0*HAP_MAX_INT, "Bad right1-a" );
                logpStep = listNodeSplitProb[childLeft] + listNodeSplitProb[childRight];
            }
//cout << "log prob: " << logpStep << " for node: " << node << endl;
            listNodeSplitProb[node] = logpStep;
        }
        
        // save the probs
        this->tblQ[site] = listNodeSplitProb;
    }
}

void ScistFastRerootLocalSearch :: InitM0Tbl()
{
//cout << "InitM0Tbl: " << endl;
    // collect the maximum Q values of each subtree
    this->tblM0.clear();
    this->tblM0.resize( GetNumSites() );
    for(int s=0; s<GetNumSites(); ++s)
    {
        this->tblM0[s].resize( mtree.GetTotNodesNum() );
    }
    for(int site=0; site<GetNumSites(); ++site)
    {
//cout << "site: " << site << endl;
        for(int node=0; node<mtree.GetTotNodesNum(); ++node)
        {
            if( mtree.IsLeaf(node) )
            {
                // a single leaf in the split
                tblM0[site][node] = tblQ[site][node];
            }
            else
            {
                // get the two children and add them up
                int childLeft = mtree.GetLeftDescendant(node);
                int childRight = mtree.GetRightDescendant(node);
    //cout << "node: " << node << ", childLeft: " << childLeft << ", childRight: " << childRight << endl;
                //cout << "childLeft: " << childLeft << ", right: " << childRight << endl;
                
                YW_ASSERT_INFO( tblM0[site][childLeft] > -1.0*HAP_MAX_INT, "Bad left-b" );
                YW_ASSERT_INFO( tblM0[site][childRight] > -1.0*HAP_MAX_INT, "Bad right1-b" );
                double probM0 = std::max(tblM0[site][childLeft], tblM0[site][childRight]);
                tblM0[site][node] = std::max(probM0, tblQ[site][node]);
            }
//cout << "Set tblM0: site:" << site << ", node:" << node << ", to: " << tblM0[site][node] << endl;
        }
    }
}
void ScistFastRerootLocalSearch :: InitM1Tbl()
{
    //if( this->numThreads > 1 && GetNumSites() > this->numThreads )
    if( this->numThreads > 1 )
    {
        InitM1TblMT();
        return;
    }
    this->tblM1.clear();
    this->tblM1.resize( GetNumSites() );
//cout << "InitM1Tbl: " << endl;
    // first iterate over sites
    for(int site=0; site<GetNumSites(); ++site)
    {
//cout << "site: " << site << endl;
        // do iteration over each node from bottom up
        this->tblM1[site].resize(mtree.GetTotNodesNum());
        for(int node=0; node<mtree.GetTotNodesNum(); ++node)
        {
            this->tblM1[site][node].resize( GetNumAncesForNode(node) );
            // trace up until reaching the root
            int nodeAnces = node;
            int nodeAncesInd = 0;
            double probMaxPath = MAX_NEG_DOUBLE_VAL;
            while(nodeAnces >= 0)
            {
                //
                this->tblM1[site][node][nodeAncesInd++] = probMaxPath;
//cout << "Set tblM1: site:" << site << ", node:" << node << ", nodeAnces: " << nodeAnces << ", to: " << tblM1[site][node][nodeAnces] << endl;
                double prob2 = this->tblQ[site][nodeAnces];
                if( prob2 > probMaxPath )
                {
                    probMaxPath = prob2;
                }
                
                if( nodeAnces == mtree.GetRoot() )
                {
                    break;
                }
                else
                {
                    nodeAnces = mtree.GetParent(nodeAnces);
                }
            }
        }
    }
}

// muliti-threading
static void UtilsInitM1Tbl2( ScistFastRerootLocalSearch *pthis, MarginalTree *ptree, int posStart, int posEnd, std::vector<std::vector< std::vector<double> > > *pTblM1 )
{
    // first iterate over sites
    for(int site=posStart; site<=posEnd; ++site)
    {
        // do iteration over each node from bottom up
        for(int node=0; node<ptree->GetTotNodesNum(); ++node)
        {
            // trace up until reaching the root
            int nodeAnces = node;
            int nodeAncesIndex = 0;
            double probMaxPath = MAX_NEG_DOUBLE_VAL;
            while(nodeAnces >= 0)
            {
                //
                (*pTblM1)[site][node][nodeAncesIndex++] = probMaxPath;
//cout << "Set tblM1: site:" << site << ", node:" << node << ", nodeAnces: " << nodeAnces << ", to: " << tblM1[site][node][nodeAnces] << endl;
                double prob2 = pthis->GetQVal(site, nodeAnces);
                if( prob2 > probMaxPath )
                {
                    probMaxPath = prob2;
                }
                if( nodeAnces == ptree->GetRoot() )
                {
                    break;
                }
                else
                {
                    nodeAnces = ptree->GetParent(nodeAnces);
                }
            }
        }
    }
}
void ScistFastRerootLocalSearch :: InitM1TblMT()
{
//cout << "InitM1TblMT:\n";
    this->tblM1.clear();
    this->tblM1.resize( GetNumSites() );
    for(int s=0; s<GetNumSites(); ++s)
    {
        this->tblM1[s].resize(  GetNumNodesInTree() );
        for(int node=0; node<GetNumNodesInTree(); ++node)
        {
            this->tblM1[s][node].resize( GetNumAncesForNode(node) );
        }
    }

    int numItems = GetNumSites();
    int blkSize = numItems/this->numThreads;
    int blkStart = 0;
    int numThreadsUse = this->numThreads;
    if( this->numThreads > GetNumSites() )
    {
        numThreadsUse = GetNumSites();
    }
    vector<thread *> listPtrThreads;
    for(int t=0; t<numThreadsUse; ++t)
    {
        int blkEnd = blkStart + blkSize - 1;
        if( blkEnd >= numItems  || t == numThreadsUse-1)
        {
            blkEnd = numItems-1;
        }
        YW_ASSERT_INFO(blkEnd >= blkStart, "WRONG345");
        // start thread
        thread *pthr = new thread(UtilsInitM1Tbl2, this, &mtree, blkStart, blkEnd, &tblM1);
        listPtrThreads.push_back(pthr);
        blkStart += blkSize;
    }
    
    // wait to join
    for(unsigned int j=0; j<listPtrThreads.size(); ++j)
    {
      listPtrThreads[j]->join();
    }
    for(unsigned int j=0; j<listPtrThreads.size(); ++j)
    {
      delete listPtrThreads[j];
    }
    listPtrThreads.clear();
//cout << "InitM1TblMT: multithreading done\n";
}

void ScistFastRerootLocalSearch :: InitM2Tbl()
{
    if( this->numThreads > 1 )
    {
        InitM2TblMT();
        return;
    }
//cout << "^^^^^^^Init M2:\n";
    this->tblM2.clear();
    this->tblM2.resize( GetNumSites() );

    // first iterate over sites
    for(int site=0; site<GetNumSites(); ++site)
    {
        this->tblM2[site].resize(mtree.GetTotNodesNum());
//cout << "site: " << site << endl;
        // do iteration over each node from bottom up
        for(int node=0; node<mtree.GetTotNodesNum(); ++node)
        {
            this->tblM2[site][node].resize( GetNumAncesForNode(node) );
            this->tblM2[site][node][0] = MAX_NEG_DOUBLE_VAL;
            
            if( node == mtree.GetRoot() )
            {
                continue;
            }
            
            // trace up until reaching the root
            int nodeAnces = mtree.GetParent( node );
            int nodeAncesInd = 1;
            int nodeAncesBelow = node;
            double probMaxPath = MAX_NEG_DOUBLE_VAL;
            while(nodeAnces >= 0)
            {
                //
                this->tblM2[site][node][nodeAncesInd++] = probMaxPath;
//cout << "set M2: node:" << node << ", nodeAnces: " << nodeAnces << " to : " << probMaxPath << endl;
                if( nodeAnces == mtree.GetRoot() )
                {
                    break;
                }
                else
                {
                    nodeAncesBelow = nodeAnces;
                    nodeAnces = mtree.GetParent(nodeAnces);
                }
                // get the side chain at nodeAnces
                int nodeSide = mtree.GetLeftDescendant(nodeAncesBelow);
                if( IsNodeAncestralTo( nodeSide, node) )
                {
                    nodeSide = mtree.GetRightDescendant(nodeAncesBelow);
                }
                //YW_ASSERT_INFO(nodeSide != nodeAncesBelow, "Fail to find side chain node");
//cout << "nodeside:" << nodeSide << ", node: " << node << ", M0: at side node: " << this->tblM0[site][nodeSide] << endl;
                
                double prob2 = this->tblM0[site][nodeSide];
                if( prob2 > probMaxPath )
                {
                    probMaxPath = prob2;
                }
            }
        }
    }
}

// multi-threading
static void UtilsInitM2Tbl2(ScistFastRerootLocalSearch *pthis, MarginalTree *ptree, int posStart, int posEnd, std::vector<std::vector< std::vector<double> > > *pTblM2)
{
    // first iterate over sites
    for(int site=posStart; site<=posEnd; ++site)
    {
//cout << "site: " << site << endl;
        // do iteration over each node from bottom up
        for(int node=0; node<ptree->GetTotNodesNum(); ++node)
        {
            (*pTblM2)[site][node][0] = MAX_NEG_DOUBLE_VAL;
            if( node == ptree->GetRoot() )
            {
                continue;
            }
            
            // trace up until reaching the root
            int nodeAnces = ptree->GetParent( node );
            int nodeAncesInd = 1;
            int nodeAncesBelow = node;
            double probMaxPath = MAX_NEG_DOUBLE_VAL;
            while(nodeAnces >= 0)
            {
                //
                (*pTblM2)[site][node][nodeAncesInd++] = probMaxPath;
//cout << "set M2: node:" << node << ", nodeAnces: " << nodeAnces << " to : " << probMaxPath << endl;
                if( nodeAnces == ptree->GetRoot() )
                {
                    break;
                }
                else
                {
                    nodeAncesBelow = nodeAnces;
                    nodeAnces = ptree->GetParent(nodeAnces);
                }
                // get the side chain at nodeAnces
                int nodeSide = ptree->GetLeftDescendant(nodeAncesBelow);
                if( pthis->IsNodeAncestralTo( nodeSide, node) )
                {
                    nodeSide = ptree->GetRightDescendant(nodeAncesBelow);
                }
                //YW_ASSERT_INFO(nodeSide != nodeAncesBelow, "Fail to find side chain node");
//cout << "nodeside:" << nodeSide << ", node: " << node << ", M0: at side node: " << this->tblM0[site][nodeSide] << endl;
                
                double prob2 = pthis->GetM0Val(site, nodeSide);
                if( prob2 > probMaxPath )
                {
                    probMaxPath = prob2;
                }
            }
        }
    }
}

void ScistFastRerootLocalSearch :: InitM2TblMT()
{
//cout << "^^^^^^^Init M2:\n";
    this->tblM2.clear();
    this->tblM2.resize( GetNumSites() );
    for(int s=0; s<GetNumSites(); ++s)
    {
        this->tblM2[s].resize( GetNumNodesInTree() );
        for(int node=0; node<GetNumNodesInTree(); ++node)
        {
            this->tblM2[s][node].resize( GetNumAncesForNode(node) );
        }
    }

    int numItems = GetNumSites();
    int numThreadsUse = this->numThreads;
    if( numThreadsUse > numItems )
    {
        numThreadsUse = numItems;
    }
    int blkSize = numItems/numThreadsUse;
    int blkStart = 0;
    vector<thread *> listPtrThreads;
    for(int t=0; t<numThreadsUse; ++t)
    {
        int blkEnd = blkStart + blkSize - 1;
        if( blkEnd >= numItems  || t == numThreadsUse-1)
        {
            blkEnd = numItems-1;
        }
        YW_ASSERT_INFO(blkEnd >= blkStart, "WRONG345");
        // start thread
        thread *pthr = new thread(UtilsInitM2Tbl2, this, &mtree, blkStart, blkEnd, &tblM2);
        listPtrThreads.push_back(pthr);
        blkStart += blkSize;
    }
    
    // wait to join
    for(unsigned int j=0; j<listPtrThreads.size(); ++j)
    {
      listPtrThreads[j]->join();
    }
    for(unsigned int j=0; j<listPtrThreads.size(); ++j)
    {
      delete listPtrThreads[j];
    }
    listPtrThreads.clear();
//cout << "InitM2TblMT: multithreading done\n";
}

void ScistFastRerootLocalSearch :: InitM4Tbl()
{
    tblM4.clear();
    this->tblM4.resize( GetNumSites() );
    for(int s=0; s<GetNumSites(); ++s)
    {
        this->tblM4[s].resize( GetNumNodesInTree() );
        for(int node=0; node<GetNumNodesInTree(); ++node)
        {
            this->tblM4[s][node].resize(GetNumAncesForNode(node));
        }
    }
//return;
    /*tblM5.clear();
    this->tblM5.resize( GetNumSites() );
    for(int s=0; s<GetNumSites(); ++s)
    {
        this->tblM5[s].resize( GetNumNodesInTree() );
    }*/
    
    if( this->numThreads > 1 )
    {
        InitM4TblMT();
        return;
    }
    
    // collect M4[s][nodeu][noder] = max partial sum of off-path [noder->nodeu] Q values
    // now collect M4
    for(int site=0; site<GetNumSites(); ++site)
    {
//cout << "site: " << site << endl;
        // do iteration over each node from bottom up
        for(int node=0; node<mtree.GetTotNodesNum(); ++node)
        {
            //tblM4[site][node][0] = 0;
            // trace up until reaching the root
            int nodeCurr = node;
            int nodeAnces = mtree.GetParent(node);
            
            if( nodeAnces < 0 )
            {
                continue;
            }
            int nodeAncesInd = 1;
            tblM4[site][node][nodeAncesInd] = 0;
            while(nodeAnces >= 0)
            {
                //
                int nodeAncesAnces = mtree.GetParent(nodeAnces);
                int nodeAncesAncesInd = nodeAncesInd;
                ++nodeAncesAncesInd;
                if( nodeAncesAnces >= 0 )
                {
                    int nodeSib = mtree.GetSibling(nodeCurr);
                    //tblM4[site][node][nodeAncesAnces] = std::max( tblM4[site][node][nodeAnces] + this->tblQ[site][nodeSib], this->tblQ[site][nodeSib]);
                    tblM4[site][node][nodeAncesAncesInd] = std::max( tblM4[site][node][nodeAncesInd] + this->tblQ[site][nodeSib], this->tblQ[site][nodeSib]);
//cout << "Set tblM4: site:" << site << ", node:" << node << ", nodeAncesAnces: " << nodeAncesAnces << ", to: " << tblM4[site][node][nodeAncesAnces] << "      tblQ[sib] = " << this->tblQ[site][nodeSib] << endl;
                }
                if( nodeAnces == mtree.GetRoot() )
                {
                    break;
                }
                else
                {
                    nodeCurr = nodeAnces;
                    nodeAnces = mtree.GetParent(nodeAnces);
                    ++nodeAncesInd;
                }
            }
        }
    }
    //cout << "Done with tblM4 init.\n";
}

// multi-threading
static void UtilsInitM4Tbl2(ScistFastRerootLocalSearch *pthis, MarginalTree *ptree, int posStart, int posEnd, std::vector<std::vector< std::vector<double> > > *pTblM4)
{
    // first iterate over sites
    for(int site=posStart; site<=posEnd; ++site)
    {
//cout << "site: " << site << endl;
        // do iteration over each node from bottom up
        for(int node=0; node<ptree->GetTotNodesNum(); ++node)
        {
            //(*pTblM4)[site][node][node] = 0;
            // trace up until reaching the root
            int nodeCurr = node;
            int nodeAnces = ptree->GetParent(node);
            if( nodeAnces < 0 )
            {
                continue;
            }
            int nodeAncesIndex = 1;
            (*pTblM4)[site][node][nodeAncesIndex] = 0;
            while(nodeAnces >= 0)
            {
                //
                int nodeAncesAnces = ptree->GetParent(nodeAnces);
                int nodeAncesAncesIndex = nodeAncesIndex+1;
                if( nodeAncesAnces >= 0 )
                {
                    int nodeSib = ptree->GetSibling(nodeCurr);
                    double prob2 = pthis->GetQVal(site, nodeSib);
                    //(*pTblM4)[site][node][nodeAncesAnces] = std::max( (*pTblM4)[site][node][nodeAnces] + prob2, prob2 );
                    (*pTblM4)[site][node][nodeAncesAncesIndex++] = std::max( (*pTblM4)[site][node][nodeAncesIndex] + prob2, prob2 );
                    //cout << "Set tblM4: site:" << site << ", node:" << node << ", nodeAncesAnces: " << nodeAncesAnces << ", to: " << tblM4[site][node][nodeAncesAnces] << "      tblQ[sib] = " << this->tblQ[site][nodeSib] << endl;
                }
                if( nodeAnces == ptree->GetRoot() )
                {
                    break;
                }
                else
                {
                    nodeCurr = nodeAnces;
                    nodeAnces = ptree->GetParent(nodeAnces);
                    ++ nodeAncesIndex;
                }
            }
        }
    }
}


void ScistFastRerootLocalSearch :: InitM4TblMT()
{
    int numItems = GetNumSites();
    int numThreadsUse = this->numThreads;
    if( numThreadsUse > numItems )
    {
        numThreadsUse = numItems;
    }
    int blkSize = numItems/numThreadsUse;
    int blkStart = 0;
//cout << "InitM4TblMT : numThresuse: " << numThreadsUse << ", blk size: " << blkSize << endl;
    vector<thread *> listPtrThreads;
    for(int t=0; t<numThreadsUse; ++t)
    {
        int blkEnd = blkStart + blkSize - 1;
        if( blkEnd >= numItems  || t == numThreadsUse-1)
        {
            blkEnd = numItems-1;
        }
        YW_ASSERT_INFO(blkEnd >= blkStart, "WRONG345");
        // start thread
        thread *pthr = new thread(UtilsInitM4Tbl2, this, &mtree, blkStart, blkEnd, &tblM4);
        listPtrThreads.push_back(pthr);
        blkStart += blkSize;
    }
    
    // wait to join
    for(unsigned int j=0; j<listPtrThreads.size(); ++j)
    {
      listPtrThreads[j]->join();
    }
    for(unsigned int j=0; j<listPtrThreads.size(); ++j)
    {
      delete listPtrThreads[j];
    }
    listPtrThreads.clear();
//cout << "Done with InitM4TblMT\n";
}

//*************************************************************
bool ScistFastRerootLocalSearch :: IsNodeAncestralTo(int node1, int node2) const
{
    // is node1 ancestral to node2?
    //map<int, set<int> > :: const_iterator it1 = mapListAnces.find(node2);
    //return it1->second.find(node1) != it1->second.end();
    return mapListAnces[node2].find(node1) != mapListAnces[node2].end();
}

int ScistFastRerootLocalSearch :: GetNumSites() const
{
    return genosInput.GetNumSites();
}
int ScistFastRerootLocalSearch :: GetNumCells() const
{
    return genosInput.GetNumHaps();
}
int ScistFastRerootLocalSearch :: GetNumNodesInTree() const
{
    return 2*GetNumCells()-1;
}
std::string ScistFastRerootLocalSearch :: ConsTreeFromSetClusters( const std::set<ScistPerfPhyCluster> &setClusters ) const
{
    //
    // now construct tree
    ScistInfPerfPhyUtils treeBuild;
    map<int, ScistPerfPhyCluster> mapPickedClus;
    int s = 0;
    for(set<ScistPerfPhyCluster>:: iterator it = setClusters.begin(); it != setClusters.end(); ++it)
    {
        mapPickedClus[s] = *it;
        ++s;
    }
    string strTree = treeBuild.ConsTreeWCombDistClus( this->genosInput, mapPickedClus );
    return strTree;
}

double ScistFastRerootLocalSearch :: CalcMaxQForReroot(int site, int node, int nodeSub, int rootIndex, int rootSubIndex, int childRootOtherSide, int nodeSubNodeIndex) const
{
//cout << "CalcMaxQForReroot: site:" << site << ", node:" << node << ", nodeSub:" << nodeSub << ", rootIndex:" << rootIndex << ", rootSubIndex:" << rootSubIndex << ", childRootOtherSide=" << childRootOtherSide << ", nodeSubNodeIndex: " << nodeSubNodeIndex << endl;

    // get the child of node that is not ancestral nodeSub (i.e., on the other side)
    int nodeChildOther = mtree.GetLeftDescendant(node);
    if( IsNodeAncestralTo(nodeChildOther, nodeSub) )
    {
        nodeChildOther = mtree.GetRightDescendant(node);
    }
    
    // prob1: max within subtree nodeSub
    double prob1 = this->tblM0[site][nodeSub];
    // prob2: max within subtree nodeOtherChild
    double prob2 = this->tblM0[site][nodeChildOther];
    // prob3: the original Q at node
    double prob3 = this->tblQ[site][node];
    // prob5: max along path from node to root
    double prob5 = MAX_NEG_DOUBLE_VAL;
    if( node != mtree.GetRoot() )
    {
        prob5 = this->tblM1[site][node][rootIndex];
    }
    // prob6: max along nodes from
    double prob6 = MAX_NEG_DOUBLE_VAL;
    //it2 = this->tblM2[site][nodeSub].find(mtree.GetRoot());
    //YW_ASSERT_INFO(it2 != this->tblM2[site][nodeSub].end(), "Fail112-5");
    //prob6 = it2->second;
    prob6 = this->tblM2[site][nodeSub][rootSubIndex];
    //}
    // prob7: max over the other side of subtree
    double prob7 = MAX_NEG_DOUBLE_VAL;
    if( node != mtree.GetRoot() )
    {
        prob7 = this->tblM0[site][childRootOtherSide];
    }
    double prob8 = MAX_NEG_DOUBLE_VAL;
    int nodev1 = mtree.GetParent(nodeSub);
    if( nodev1 != node )
    {
        prob8 = this->tblM4[site][nodeSub][nodeSubNodeIndex] + this->tblQ[site][nodeChildOther];
    }
    
    // finally allow the root Q value to be included
    double prob9 = this->tblQ[site][mtree.GetRoot()];
//cout << "CalcMaxQReroot: p1:" << prob1 << ", p2:" << prob2 << ", p3:" << prob3 << ", p5:" << prob5 << ", p6:" << prob6 << ", p7:" << prob7 << ", p8:" << prob8 << ", p9:" << prob9 << endl;
    double res = std::max(prob1, std::max(prob2, std::max(prob3, std::max(prob5, std::max(prob6, std::max(prob7, std::max(prob8, prob9 )))))));
    if( res < 0.0 )
    {
        res = 0.0;
    }
//cout << "maxQ: " << res << endl;
    return res;
}

void ScistFastRerootLocalSearch :: CollectRerootCands( )
{
    //
    listCandsReroot.clear();
    //int numExc = 0;
    for(unsigned int i=0; i<mapListAnces.size(); ++i)
    {
        int nodeSub = i;
        for(set<int> :: iterator it = mapListAnces[i].begin(); it != mapListAnces[i].end(); ++it)
        {
            int node = *it;
            
            if( this->setSTExcluded.find(node) != this->setSTExcluded.end() )
            {
                //++numExc;
                continue;
            }
            
            // skip if nodeSub is direct desendant of node
            if( node == nodeSub || mtree.GetParent(nodeSub) == node )
            {
                continue;
            }
            pair<int,int> pp(nodeSub, node);
            listCandsReroot.push_back(pp);
        }
    }
//cout << "Number of excluded rooting: " << numExc << endl;
}

int ScistFastRerootLocalSearch :: GetAncesIndex(int nodeDesc, int nodeAnces) const
{
    std::map<int,int> :: const_iterator it = mapListAncesLookup[nodeDesc].find(nodeAnces);
    YW_ASSERT_INFO( it != mapListAncesLookup[nodeDesc].end(), "Fail to find ancestor" );
    return it->second;
}
