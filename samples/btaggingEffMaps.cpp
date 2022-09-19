#include "HEPAnalysis.h"

//-------------------------------------------------------------------------------------------------
// Description:
//-------------------------------------------------------------------------------------------------

//-------------------------------------------------------------------------------------------------
// Define output variables
//-------------------------------------------------------------------------------------------------
namespace btaggingEffMaps{
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
void HEPAnalysis::SetupbtaggingEffMaps() {

    // Setup cutflow
    _cutFlow.insert(pair<string,double>("00_TwoLepOS", 0) );
    _cutFlow.insert(pair<string,double>("01_MET", 0) );
    _cutFlow.insert(pair<string,double>("02_LeadingLep_Pt", 0) );
    _cutFlow.insert(pair<string,double>("03_LepLep_DM", 0) );
    _cutFlow.insert(pair<string,double>("04_LepLep_Pt", 0) );
    _cutFlow.insert(pair<string,double>("05_LepLep_DR", 0) );
    _cutFlow.insert(pair<string,double>("06_Selected", 0) );

    // Identifiers
    _outputTree->Branch("EventPosition", &_EventPosition );

    // Variables
    _outputTree->Branch("Jet_jetId", &btaggingEffMaps::Jet_jetId);
    _outputTree->Branch("Jet_hadronFlavour", &btaggingEffMaps::Jet_hadronFlavour);
    _outputTree->Branch("Jet_pt", &btaggingEffMaps::Jet_pt);
    _outputTree->Branch("Jet_eta", &btaggingEffMaps::Jet_eta);
    _outputTree->Branch("Jet_btagDeepB", &btaggingEffMaps::Jet_btagDeepB);
    _outputTree->Branch("Jet_btagDeepFlavB", &btaggingEffMaps::Jet_btagDeepFlavB);

    return;
}


//-------------------------------------------------------------------------------------------------
// Define the selection region
//-------------------------------------------------------------------------------------------------
bool HEPAnalysis::btaggingEffMapsRegion() {

    LeptonSelection();
    
    if( !(RecoLepID > 0) ) return false;                                        // Has two reconstructed leptons with opposite signal
    _cutFlow.at("00_TwoLepOS") += evtWeight;
    
    JetSelection();
    METCorrection();
    if( !(MET_pt > MET_BASE_CUT) ) return false;                                 // MET > BASE_CUT
    _cutFlow.at("01_MET") += evtWeight; 
    
    Get_LeadingAndTrailing_Lepton_Variables();
    
    if( !(LeadingLep_pt > LEADING_LEP_PT_BASE_CUT) ) return false;               // Leading lepton pt > BASE_CUT
    _cutFlow.at("02_LeadingLep_Pt") += evtWeight; 
    
    if( !DYJetsToLL_processing() ) return false;
    
    Get_LepLep_Variables();
    
    if( !(LepLep_deltaM < LEPLEP_DM_CUT) ) return false;                        // Difference between Z boson mass and the inv. mass of two leptons < CUT
    _cutFlow.at("03_LepLep_DM") += evtWeight;  
    
    if( !(LepLep_pt > LEPLEP_PT_BASE_CUT) ) return false;                       // Two leptons system pt > BASE_CUT
    _cutFlow.at("04_LepLep_Pt") += evtWeight;
    
    if( !(LepLep_deltaR < LEPLEP_DR_BASE_CUT) ) return false;                   // Upper cut in LepLep Delta R 
    _cutFlow.at("05_LepLep_DR") += evtWeight; 
    
    bool GoodEvent = lumi_certificate.GoodLumiSection( _datasetName, run, luminosityBlock );
    if( !GoodEvent ) return false;                                              // Select only certified data events
    
    if( !METFilters() ) return false;                                           // Selected by MET filters
    
    if( !Trigger() ) return false;                                              // Selected by triggers
    _cutFlow.at("06_Selected") += evtWeight; 
    
    //Get_Jet_Angular_Variables( );
    Weight_corrections();

    return true;
}


//-------------------------------------------------------------------------------------------------
// Write your analysis code here
//-------------------------------------------------------------------------------------------------
void HEPAnalysis::btaggingEffMapsSelection() {

    // Only selected jets
    // * Jet_pt > 20
    // * Jet_jetId >= 6 (TightLepVeto)
    // * |Jet_eta| < 2.4
    for (unsigned int ijet = 0; ijet < nJet; ++ijet) {
        if (
            Jet_pt[ijet] > JET_PT_CUT
            && Jet_jetId[ijet] >= JET_ID_WP
            && abs(Jet_eta[ijet]) < JET_ETA_CUT
        ) {
            btaggingEffMaps::Jet_jetId = Jet_jetId[ijet];
            btaggingEffMaps::Jet_hadronFlavour = Jet_hadronFlavour[ijet];
            btaggingEffMaps::Jet_pt = Jet_pt[ijet];
            btaggingEffMaps::Jet_eta = Jet_eta[ijet];
            btaggingEffMaps::Jet_btagDeepB = Jet_btagDeepB[ijet]; // DeepCSV
            btaggingEffMaps::Jet_btagDeepFlavB = Jet_btagDeepFlavB[ijet]; // DeepJet
            _outputTree->Fill();
        }
    }

    return;
}


//-------------------------------------------------------------------------------------------------
// Produce systematic histograms
//-------------------------------------------------------------------------------------------------
void HEPAnalysis::btaggingEffMapsSystematic() { }


//-------------------------------------------------------------------------------------------------
// Make efficiency plots
//-------------------------------------------------------------------------------------------------
void HEPAnalysis::FinishbtaggingEffMaps() {
    return;
}