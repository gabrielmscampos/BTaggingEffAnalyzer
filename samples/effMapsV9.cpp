#include "HEPAnalysis.h"

//-------------------------------------------------------------------------------------------------
// Description:
//-------------------------------------------------------------------------------------------------

//-------------------------------------------------------------------------------------------------
// Define output variables
//-------------------------------------------------------------------------------------------------
namespace effMapsV9{
    int Jet_jetId;
    int Jet_hadronFlavour;
    float Jet_pt;
    float Jet_eta;
    float Jet_btagDeepB;
    float Jet_btagDeepFlavB;
}


//-------------------------------------------------------------------------------------------------
// Define output derivatives
//-------------------------------------------------------------------------------------------------
void HEPAnalysis::SetupeffMapsV9() {

    //======SETUP CUTFLOW==========================================================================
    _cutFlow.insert(pair<string,double>("00_TwoLepOS", 0) );
    _cutFlow.insert(pair<string,double>("01_MET", 0) );
    _cutFlow.insert(pair<string,double>("02_LeadingLep_Pt", 0) );
    _cutFlow.insert(pair<string,double>("03_LepLep_DM", 0) );
    _cutFlow.insert(pair<string,double>("04_LepLep_Pt", 0) );
    _cutFlow.insert(pair<string,double>("05_LepLep_DR", 0) );
    _cutFlow.insert(pair<string,double>("06_Selected", 0) );

    //======SETUP HISTOGRAMS=======================================================================
    //makeHist( "histogram1DName", 40, 0., 40., "xlabel", "ylabel" );   [example]
    //makeHist( "histogram2DName", 40, 0., 40., 100, 0., 50., "xlabel",  "ylabel", "zlabel", "COLZ" );   [example]

    //======SETUP SYSTEMATIC HISTOGRAMS============================================================
    //sys_regions = { 0, 1, 2 }; [example] // Choose regions as defined in RegionID. Empty vector means that all events will be used.
    //makeSysHist( "histogram1DSysName", 40, 0., 40., "xlabel", "ylabel" );   [example]
    //makeSysHist( "histogram2DSysName", 40, 0., 40., 100, 0., 50., "xlabel",  "ylabel", "zlabel", "COLZ" );   [example]

    //======SETUP OUTPUT BRANCHES==================================================================
    //_outputTree->Branch("variable1NameInTheTree", &effMapsV9::variable1Name );  [example]

    //======SETUP INFORMATION IN OUTPUT HDF5 FILE==================================================
    HDF_insert("EventPosition", &_EventPosition );
    HDF_insert("Jet_jetId", &effMapsV9::Jet_jetId );
    HDF_insert("Jet_hadronFlavour", &effMapsV9::Jet_hadronFlavour );
    HDF_insert("Jet_pt", &effMapsV9::Jet_pt );
    HDF_insert("Jet_eta", &effMapsV9::Jet_eta );
    HDF_insert("Jet_btagDeepB", &effMapsV9::Jet_btagDeepB );
    HDF_insert("Jet_btagDeepFlavB", &effMapsV9::Jet_btagDeepFlavB );

    return;
}


//-------------------------------------------------------------------------------------------------
// Define the selection region
//-------------------------------------------------------------------------------------------------
bool HEPAnalysis::effMapsV9Region() {

    LeptonSelection();

    if( !(RecoLepID > 0) ) return false;                                        // Has two reconstructed leptons with opposite signal
    _cutFlow.at("00_TwoLepOS") += evtWeight;

    // JetSelection();
    METCorrection();

    if( !(MET_pt > MET_CUT) ) return false;                                 // MET > BASE_CUT
    _cutFlow.at("01_MET") += evtWeight;

    Get_Leptonic_Info(true, true);

    if( !(LeadingLep_pt > LEADING_LEP_PT_CUT) ) return false;               // Leading lepton pt > BASE_CUT
    _cutFlow.at("02_LeadingLep_Pt") += evtWeight;

    Get_LepLep_Variables(true, true);

    if( !(LepLep_deltaM < LEPLEP_DM_CUT) ) return false;                        // Difference between Z boson mass and the inv. mass of two leptons < CUT
    _cutFlow.at("03_LepLep_DM") += evtWeight;

    if( !(LepLep_pt > LEPLEP_PT_CUT) ) return false;                       // Two leptons system pt > BASE_CUT
    _cutFlow.at("04_LepLep_Pt") += evtWeight;

    if( !(LepLep_deltaR < LEPLEP_DR_CUT) ) return false;                   // Upper cut in LepLep Delta R
    _cutFlow.at("05_LepLep_DR") += evtWeight;

    bool GoodEvent = lumi_certificate.GoodLumiSection( _datasetName, run, luminosityBlock );
    if( !GoodEvent ) return false;                                              // Select only certified data events

    if( !METFilters() ) return false;                                           // Selected by MET filters

    if( !Trigger() ) return false;                                              // Selected by triggers
    _cutFlow.at("06_Selected") += evtWeight;

    Weight_corrections();

    return true;
}


//-------------------------------------------------------------------------------------------------
// Write your analysis code here
//-------------------------------------------------------------------------------------------------
void HEPAnalysis::effMapsV9Selection() {


    for (unsigned int ijet = 0; ijet < nJet; ++ijet) {

        if( Jet_pt[ijet] <= JET_PT_CUT ) continue;
        if( Jet_jetId[ijet] < JET_ID_WP ) continue;
        if( abs(Jet_eta[ijet]) >= JET_ETA_CUT ) continue;

        effMapsV9::Jet_jetId = Jet_jetId[ijet];
        effMapsV9::Jet_hadronFlavour = Jet_hadronFlavour[ijet];
        effMapsV9::Jet_pt = Jet_pt[ijet];
        effMapsV9::Jet_eta = Jet_eta[ijet];
        effMapsV9::Jet_btagDeepB = Jet_btagDeepB[ijet]; // DeepCSV
        effMapsV9::Jet_btagDeepFlavB = Jet_btagDeepFlavB[ijet]; // DeepJet
        HDF_fill();
    }

    return;
}


//-------------------------------------------------------------------------------------------------
// Produce systematic histograms
//-------------------------------------------------------------------------------------------------
void HEPAnalysis::effMapsV9Systematic() {

    //FillSystematic( "histogram1DSysName", var, evtWeight );  [Example]
    //FillSystematic( "histogram2DSysName", var1, var2, evtWeight );  [Example]
}


//-------------------------------------------------------------------------------------------------
// Make efficiency plots
//-------------------------------------------------------------------------------------------------
void HEPAnalysis::FinisheffMapsV9() {

    //MakeEfficiencyPlot( _histograms1D.at("Matched_pt"), _histograms1D.at("all_pt"), "Match_pt" );   [example]

    return;
}
