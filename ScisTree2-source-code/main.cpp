#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <cstring>
#include <vector>
#include <map>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

using namespace std;

#include "Utils2.h"
#include "Utils3.h"
#include "Utils4.h"
#include "ScistPerfPhyUtils.hpp"
#include "ScistPerfPhyImp.hpp"
#include "ScistGenotype.hpp"
#include "ScistDoublet.hpp"
#include "ScistErrRateInf.hpp"

//*****************************************************************************
// Main driving functions
// ***************************************************************************


// ***************************************************************************
// Main for computing lower bound
// ***************************************************************************
static void Usage()
{
	cout << "Usage: ./scistree <options> <input file> " << endl;
    cout << "Options:\n";
    //cout << "\t -d <dn>           dn: number of allowed doublet genotypes; dc: cost of having a doublet\n";
#if 0   // not currently supported
    cout << "\t -d <dn>           dn: number of allowed doublet genotypes\n";
#endif
    cout << "\t -v                Turn on verbose mode  \n";
    cout << "\t -T <nthread>      Set number of threads to nthread\n";
    //cout << "\t -p                Find optimal false positive rate and false negative rate\n";
    //cout << "\t -l                Find cell tree with branch length (by default, constructed cell trees don't have branch length\n";
    cout << "\t -n                Only build simple neighbor joining tree (may be useful for very large data)\n";
    cout << "\t -e                Output mutation tree (may not be binary tree) from called genotypes branch labels.\n";
    cout << "\t -e0               Output mutation tree but don't output labels (for visualizing large trees).\n";
    //cout << "\t -s <level>        Use SPR tree search (this will be slower); level: # of SPRs to allow (default is 1)\n";
    //cout << "\t -s                Use SPR local tree search (this is the default)\n";
    cout << "\t -q                Use NNI local tree search (NNI is faster but less accurate)\n";
    //cout << "\t -S                Turn on exhaustive SPR local tree search (this will be even slower)\n";
//    cout << "\t -r  <f> <d>       Config how SPR runs: f between 0 and 1 (default: 0.5); the smaller, the faster but less accurate); dropStop: an integer (default: 100; the larger, the more accurate but slower)\n";
    //cout << "\t -r  <f> <d>       Config SPR: f from 0 and 1 (default: 0.5); d: integer (default: 100)\n";
    //cout << "\t -u                turn on subtree-rerooting search (this can be a little slower)\n";
    //cout << "\t -I                run tree search iteratively (this can be a little slower)\n";
    //cout << "\t -x                run tree search faster (this may somewhat reduce the accuracy)\n";
    cout << "\t -o <output-file>  Set output file (used for mutation tree output (in GML) format; should have suffix .gml (default: mutation-tree.gml)\n";
    //cout << "\t -t <threshold>    Discard somewhat ambigous genotyeps when constructing intial trees: \n\t\t\t genotypes discarded if the prob. of alternative genotypes is less than <threshold> \n\t\t\t(default is 0, i.e. use all genotypes)\n";
#if 0  // not supported in this version
    cout << "\t -f <initial-trees-file>     Perform local search from the list of trees (in Newick format, one tree per line) instead of constructing initial trees by ScisTree\n";
#endif
	exit(1);
}

// settings
static int fileInArgIndex = 1;
static int numDoublets = 0;
static double costDoublet = 0.0;
static bool fVerbose = false;
static bool fOptParam = false;
static bool fOptBrLen = false;
static bool fNJOnly = false;
static bool fSPR = true;
//static int numHeuSPRs = 0;
static bool fExactSPRSearch = false;
static bool fReroot = false;
static double fracSPRSrc = 0.5;
static int thresSPRDropStop = 100;
//static int numSPR = 1;
//static bool fTest = false;
static double thresProbSignificance = 0.0;
static vector<string> listCellNames;
static vector<string> listSiteNames;
static int numSites = 0;
static int numSCs = 0;
static int numThreads = 1;
static string strMutTreeOutFile = "mutation-tree.gml";
static bool fOutPPEdgeLabel = false;
static bool fOutputLabel=true;
static bool fFastTree = false;
static bool fIterative = false;
//static vector<string> listInitTrees;
// GLobal variables



// Local functions
static bool CheckArguments(int argc, char **argv) 
{
    if( argc <= 1  )
    {
        return false;
    }

    // Check argument one by one
    //int argpos = 1;
    for(int i = 1; i< argc; ++i)
    {
        if( argv[i][0] == '-' && argv[i][1] == 'l' )
        {
            cout << "Branch length optimization is not yet supported by Scistree2. \n";
            exit(1);
            
            YW_ASSERT_INFO( i <argc-1, "Check input" );
            fOptBrLen = true;
            cout << "Turn on branch optimization. " << endl;
        }
#if 0
        else if( argv[i][0] == '-' && argv[i][1] == 'd' )
        {
            cout << "Doublet is not yet supported by Scistree2. \n";
            exit(1);
            
            YW_ASSERT_INFO( i <argc-1, "Check input" );
            ++i;
            sscanf(argv[i], "%d", &numDoublets);
            //YW_ASSERT_INFO( i <argc-1, "Check input" );
            //++i;
            //float costDoubletThis;
            //sscanf(argv[i], "%f", &costDoubletThis);
            //costDoublet = costDoubletThis;
            //cout << "Setting doublet number to " << numDoublets << ", and doublet cost to " << costDoublet << endl;
            cout << "Setting doublet number to " << numDoublets << endl;
        }
#endif
        else if( argv[i][0] == '-' && argv[i][1] == 'v' )
        {
            YW_ASSERT_INFO( i <argc-1, "Check input" );
            fVerbose = true;
            cout << "Turn on verbose mode" << endl;
        }
        else if( argv[i][0] == '-' && argv[i][1] == 'n' )
        {
            YW_ASSERT_INFO( i <argc-1, "Check input" );
            fNJOnly = true;
            cout << "Only build initial tree (no local search). Note: turn on -v to see the initial tree." << endl;
        }
#if 0
        else if( argv[i][0] == '-' && argv[i][1] == 'p' )
        {
            YW_ASSERT_INFO( i <argc-1, "Check input" );
            fOptParam = true;
            cout << "Search for optimal genotype error rates" << endl;
        }
#endif
        else if( argv[i][0] == '-' && argv[i][1] == 'e' )
        {
            YW_ASSERT_INFO( i <argc-1, "Check input" );
            fOutPPEdgeLabel = true;
            cout << "Output perfect phylogeny with edge labels" << endl;
            
            string strOpt=argv[i];
            if(strOpt.length()>=3 && strOpt[2] == '0')
            {
                cout << "  -- no labels in mutation tree\n";
                fOutputLabel=false;
            }
        }
#if 0
        else if( argv[i][0] == '-' && argv[i][1] == 's' )
        {
            YW_ASSERT_INFO( i <argc-1, "Check input" );
            fSPR = true;
            //++i;
            //sscanf(argv[i], "%d", &numSPR);
            //cout << "Use SPR tree search: level set to " << numSPR << endl;
            //cout << "Turn on rSPR local search. Note: genotype calling under SPR search is not implemented yet. \n";
            cout << "Turn on rSPR local search (this is the default). \n";
        }
#endif
        else if( argv[i][0] == '-' && argv[i][1] == 'q' )
        {
            YW_ASSERT_INFO( i <argc-1, "Check input" );
            fSPR = false;
            //++i;
            //sscanf(argv[i], "%d", &numSPR);
            //cout << "Use SPR tree search: level set to " << numSPR << endl;
            //cout << "Turn on rSPR local search. Note: genotype calling under SPR search is not implemented yet. \n";
            cout << "Turn off rSPR local search and use NNI search instead. \n";
        }
#if 0
        else if( argv[i][0] == '-' && argv[i][1] == 'S' )
        {
            YW_ASSERT_INFO( i <argc-1, "Check input" );
            //++i;
            //sscanf(argv[i], "%d", &numHeuSPRs);
            //cout << "Number of heuristic SPRs (multiple SPRs per step): set to " << numHeuSPRs << endl;
            fExactSPRSearch = true;
            // also turn on SPR
            fSPR = true;
            cout << "Turn on exaustive 1-rSPR local search mode. Note: this can be slow for large trees...\n";
        }
        else if( argv[i][0] == '-' && argv[i][1] == 'r' )
        {
            YW_ASSERT_INFO( i <argc-2, "Check input" );
            ++i;
            float frac;
            sscanf(argv[i], "%f", &frac);
            fracSPRSrc = frac;
            // now read in stopping
            ++i;
            sscanf(argv[i], "%d", &thresSPRDropStop);
            //cout << "Use SPR tree search: level set to " << numSPR << endl;
            //cout << "Turn on rSPR local search. Note: genotype calling under SPR search is not implemented yet. \n";
            cout << "Set rSPR local search parameter to: fraction=" << fracSPRSrc << ", stopping= " << thresSPRDropStop << "\n";
        }
#endif
#if 0
        else if( argv[i][0] == '-' && argv[i][1] == 't' )
        {
            YW_ASSERT_INFO( i <argc-1, "Check input" );
            ++i;
            float thresUse = 0.0;
            sscanf(argv[i], "%f", &thresUse);
            thresProbSignificance = thresUse;
            cout << "Threshold for probability significance: set to " << thresProbSignificance << endl;
        }
#endif
        else if( argv[i][0] == '-' && argv[i][1] == 'o' )
        {
            YW_ASSERT_INFO( i <argc-1, "Check input" );
            ++i;
            strMutTreeOutFile = argv[i];
            cout << "Use mutation tree file name to " << strMutTreeOutFile << endl;
        }
        else if( argv[i][0] == '-' && argv[i][1] == 'T' )
        {
            //YW_ASSERT_INFO( i <argc-1, "Test code" );
            //fTest = true;
            //cout << "Test SPR local search " << endl;
            YW_ASSERT_INFO( i <argc-1, "Check input" );
            ++i;
            sscanf(argv[i], "%d", &numThreads);
            cout << "Number of threads: set to " << numThreads << endl;
        }
#if 0
        else if( argv[i][0] == '-' && argv[i][1] == 'I' )
        {
            YW_ASSERT_INFO( i <argc-1, "Check input" );
            fIterative = true;
            cout << "Turn on iterative tree search mode...\n";
        }
        else if( argv[i][0] == '-' && argv[i][1] == 'u' )
        {
            YW_ASSERT_INFO( i <argc-1, "Check input" );
            fReroot = true;
            cout << "Turn on subtree re-rooting search mode...\n";
        }
#endif
#if 0
        else if( argv[i][0] == '-' && argv[i][1] == 'x' )
        {
            YW_ASSERT_INFO( i <argc-1, "Check input" );
            // turn off iterative mode
            fIterative = false;
            // turn off exhaustive search
            fExactSPRSearch = false;
            cout << "Turn on fast search mode...\n";
        }
#endif
#if 0
        else if( argv[i][0] == '-' && argv[i][1] == 'f' )
        {
            fFastTree = true;
            cout << "Construct initial tree using fast heuristics. " << endl;
        }
#endif
        /*else if( argv[i][0] == '-' && argv[i][1] == 'f' )
        {
            YW_ASSERT_INFO( i <argc-1, "Check input" );
            ++i;
            string strInitTreeFile = argv[i];
            cout << "Read inital tree from file:" << strInitTreeFile << endl;
            if( ReadLinesFromFile (strInitTreeFile.c_str(), listInitTrees) == false )
            {
                YW_ASSERT_INFO(false, "Fail to read the initial tree file");
            }
            cout << "Number of initial trees: " << listInitTrees.size() << endl;
        } */
        else if( argv[i][0] != '-' )
        {
            // not an option one. Right now the only one is file
            fileInArgIndex = i;
            //filenameGMLPrefix = argv[i];
        }
        else
        {
            return false;
        }
    }


    return true;
}

// input handling
static ScistGenGenotypeMat* ReadsInput(const char *filename )
{
    //
    ifstream inFile(filename);
    if(!inFile)
    {
        cout << "Can not open "<< filename <<endl;
        YW_ASSERT_INFO( false, "Stop");
    }
    ScistGenGenotypeMat *pMatIn = NULL;
    while( inFile.eof() == false )
    {
        const int BUF_SZ = 102400;
        char buffer[BUF_SZ];
        inFile.getline(buffer, BUF_SZ);
        if( strlen(buffer) > 0 )
        {
//cout << "read one line: " << buffer << endl;
            // now try to read alleles
            std::istringstream is( buffer );
            
            // looking for keyword
            string strKey;
            is >> strKey;
            if( strKey == "HAPLOTYPES" || strKey == "HAPLOID" )
            {
                is >> numSites >> numSCs;
//cout << "numSites: " << numSites << ", numSCs: " << numSCs << endl;
                YW_ASSERT_INFO(numSites >0 && numSCs > 0, "Site and single cells numbers: Cannot be zeros");
                
                // read in names if specified
                while( is.eof() == false )
                {
                    string strName;
                    is >> strName;
                    if( strName.length()>0)
                    {
                        listCellNames.push_back(strName);
//cout << "One lineage name: " << strName << endl;
                    }
                    if( (int)listCellNames.size() > numSCs )
                    {
                        break;
                    }
                }
                if( listCellNames.size() > 0 && (int)listCellNames.size() != numSCs)
                {
                    YW_ASSERT_INFO(false, "Fatal error: you must provide names for each lineage");
                }
                bool fSiteName=false;
                if(listCellNames.size()>0)
                {
                    fSiteName=true;;
                }
                
                pMatIn = new ScistHaplotypeMat;
                for(int i=0; i<(int)listCellNames.size(); ++i)
                {
                    pMatIn->AddGenotypeName( listCellNames[i] );
                }
                
                pMatIn->ReadFromFile(inFile, numSites, numSCs, fSiteName);
 
#if 0
if( fSiteName )
{
cout << "List of site names: ";
for(int i=0; i<numSites; ++i)
{
cout << pMatIn->GetSiteName(i) << " ";
}
cout << endl;
}
#endif
                
                break;
            }
            else if( strKey == "TERNARY" )
            {
                cout << "TERNARY data is not supported yet by Scistree2. Please consider using the HAPLOID data for the moment.\n";
                exit(1);
                
                is >> numSites >> numSCs;
                //cout << "numSites: " << numSites << ", numSCs: " << numSCs << endl;
                YW_ASSERT_INFO(numSites >0 && numSCs > 0, "Site and single cells numbers: Cannot be zeros");
                
                // read in names if specified
                while( is.eof() == false )
                {
                    string strName;
                    is >> strName;
                    if( strName.length()>0)
                    {
                        listCellNames.push_back(strName);
                    //cout << "One lineage name: " << strName << endl;
                    }
                    if( (int)listCellNames.size() > numSCs )
                    {
                        break;
                    }
                }
                if( listCellNames.size() > 0 && (int)listCellNames.size() != numSCs)
                {
                    YW_ASSERT_INFO(false, "Fatal error: you must provide names for each lineage");
                }
                bool fSiteName=false;
                if(listCellNames.size()>0)
                {
                    fSiteName=true;;
                }
                
                pMatIn = new ScistTernaryMat;
                for(int i=0; i<(int)listCellNames.size(); ++i)
                {
                    pMatIn->AddGenotypeName( listCellNames[i] );
                }
                
                pMatIn->ReadFromFile(inFile, numSites, numSCs, fSiteName);
                
                break;
            }
        }
    }
    pMatIn->SetSignificantThres(thresProbSignificance);
    
    // initialize cell names to plain 1, 2, ... if not specified
    if( listCellNames.size() == 0 )
    {
        YW_ASSERT_INFO(numSCs > 0, "Number of SCs: not intiialized");
        for(int c=1; c<=numSCs; ++c)
        {
            string str = std::to_string(c);
            listCellNames.push_back(str);
        }
    }
    pMatIn->GetSiteNamesAll(listSiteNames);
    
    // preprocess
    pMatIn->Preprocess();
    
    return pMatIn;
}

// test code
static void TestCode( const char *filename )
{
    ScistGenGenotypeMat *pMatInput = ReadsInput( filename );
    string filenameUse = filename;
    pMatInput->SetFileName(filenameUse);
    
    // set num of threads in case we want to speed up
    pMatInput->SetNumThreads(numThreads);
     
    if( fOptParam )
    {
        cout << "Now searching for optimal genotype error rates...\n";
        ScistErrRateInf serInf(*pMatInput);
        serInf.SetVerbose(fVerbose);
        serInf.Infer();
    }
    else
    {
        //string treeNJ = pMatInput->ConsNJTree();
        string treeNJ;
        if( fFastTree )
        {
            treeNJ = pMatInput->ConsPerfPhylogenyHeu();
            //set<set<int> > setCladesMust;
            //treeNJ = pMatInput->ConsConstrainedUPGMATreeZeroRoot( setCladesMust );
        }
        else
        {
            // use neighbor joining
            treeNJ = pMatInput->ConsNJTreeZeroRoot();
        }
        if( fVerbose )
        {
            cout << "Initial tree from noisy genotypes: " << treeNJ << endl;
        }
        if(fNJOnly)
        {
            delete pMatInput;
            return;
        }
        
        //ScistInfPerfPhyTest();
        // plain mode if no double is allowed
        if( numDoublets == 0 )
        {
            ScistPerfPhyMLE ppInfHeu(*pMatInput, &treeNJ);
            ppInfHeu.SetBrOpt(fOptBrLen);
            ppInfHeu.SetVerbose(fVerbose);
            ppInfHeu.SetPPOut(fOutPPEdgeLabel);
            ppInfHeu.SetPPOutLabel(fOutputLabel);
            ppInfHeu.SetSPR(fSPR);
            ppInfHeu.SetExactSPR(fExactSPRSearch);
            ppInfHeu.SetReroot(fReroot);
            ppInfHeu.SetSPRSrcFrac(fracSPRSrc, thresSPRDropStop);
            //ppInfHeu.SetHeuSPRNum(numHeuSPRs);
            ppInfHeu.SetNumThreads(numThreads);
            ppInfHeu.SetCellNames(listCellNames);
            ppInfHeu.SetSiteNames(listSiteNames);
            ppInfHeu.SetMutTreeFileName(strMutTreeOutFile);
            ppInfHeu.SetIterateMode(fIterative);
            std::set< std::pair<std::pair<int,int>, int> > listChangedPlaces;
            std::string strTreeNW;
            if( fVerbose )
            {
                cout << "Now searching for maximum likelihood tree... " << endl;
            }
            ppInfHeu.Infer(&listChangedPlaces, &strTreeNW);
            
#if 0
            //cout << "Matrix: ";
             //pMatInput->Dump();
             // do it again
             ScistPerfPhyMLE ppInfHeu1(*pMatInput, &strTreeNW);
             ppInfHeu1.SetBrOpt(fOptBrLen);
             ppInfHeu1.SetVerbose(fVerbose);
             ppInfHeu1.SetPPOut(fOutPPEdgeLabel);
             ppInfHeu1.SetPPOutLabel(fOutputLabel);
             ppInfHeu1.SetSPR(fSPR);
             ppInfHeu1.SetExactSPR(fExactSPRSearch);
             ppInfHeu1.SetReroot(fReroot);
             ppInfHeu1.SetSPRSrcFrac(fracSPRSrc, thresSPRDropStop);
             //ppInfHeu.SetHeuSPRNum(numHeuSPRs);
             ppInfHeu1.SetNumThreads(numThreads);
             ppInfHeu1.SetCellNames(listCellNames);
             ppInfHeu1.SetSiteNames(listSiteNames);
             ppInfHeu1.SetMutTreeFileName(strMutTreeOutFile);
             ppInfHeu1.SetIterateMode(fIterative);
             std::set< std::pair<std::pair<int,int>, int> > listChangedPlaces2;
             std::string strTreeNW2;
             if( fVerbose )
             {
                 cout << "######### Now re-searching for maximum likelihood tree... " << endl;
             }
             ppInfHeu1.Infer(&listChangedPlaces2, &strTreeNW2);
#endif
        }
        else
        {
            // right now only work with haplotype matrix
            ScistHaplotypeMat *pMatInputUse = dynamic_cast<ScistHaplotypeMat *>(pMatInput);
            YW_ASSERT_INFO( pMatInputUse != NULL, "At present, doublet feature only works for binary genotype matrix." );
            
            cout << "SEARCHING FOR DOUBLETS...\n";
            ScistDoubletSearch sds( *pMatInput, numDoublets );
            sds.SetVerbose(fVerbose);
            sds.SetDouletCost( costDoublet );
            sds.SetMutTreeOut(fOutPPEdgeLabel);
            sds.SetCellNames(listCellNames);
            sds.SetSiteNames(listSiteNames);
            sds.SetMutTreeFileName(strMutTreeOutFile);
            sds.SearchInc();
        }
    }
    
    
    delete pMatInput;
    
}


////////////////////////////////////////////////////////////////////////////////////////

//const char *CODE_VER_INFO ="*** SCISTREE ver. 1.2.0.6, May 19, 2019 ***";
//const char *CODE_VER_INFO ="*** SCISTREE ver. 1.3.0.2, January 16, 2023 ***";
//const char *CODE_VER_INFO ="*** SCISTREE ver. 1.3.1.1, January 30, 2023 ***";
//const char *CODE_VER_INFO ="*** SCISTREE ver. 1.4.1.0, Feburary 16, 2023 ***";
//const char *CODE_VER_INFO ="*** SCISTREE ver. 1.5.1.0, October 13, 2023 ***";
//const char *CODE_VER_INFO ="*** SCISTREE ver. 1.5.2.0, October 14, 2023 ***";
//const char *CODE_VER_INFO ="*** SCISTREE ver. 1.5.2.6, October 30, 2023 ***";
//const char *CODE_VER_INFO ="*** SCISTREE ver. 1.6.1.0, November 3, 2023 ***";
//const char *CODE_VER_INFO ="*** SCISTREE ver. 2.0.5.4, December 7, 2023 ***";
//const char *CODE_VER_INFO ="*** SCISTREE ver. 2.1.0.0, Janurary 22, 2024 ***";
//const char *CODE_VER_INFO ="*** SCISTREE ver. 2.1.1.0, Feburary 2, 2024 ***";
//const char *CODE_VER_INFO ="*** SCISTREE ver. 2.1.1.1, Feburary 5, 2024 ***";
const char *CODE_VER_INFO ="*** SCISTREE ver. 2.2.0.0, October 24, 2024 ***";

//******************************************************************
int main(int argc, char **argv)
{
//    int seq = 0x001;
//    int seqMut;
//    MutateHCSeqAt(seq, seqMut, 4, 2);
//cout << "mutated seq = " << seqMut << endl;
    
	cout << CODE_VER_INFO << endl << endl;

	// first verify usage
	if( CheckArguments(  argc, argv  ) == false)
	{
		Usage(); 
	}
    
//cout << "here0\n";
    long tstart1 = GetCurrentTimeTick();

    TestCode(argv[ fileInArgIndex ]);

    cout << "Elapsed time = " << GetElapseTime( tstart1 ) << " seconds." << endl;
    
    // dump out stats
    //ApproxGTPStats::Instance().DumpStats();

    return 0;

}

