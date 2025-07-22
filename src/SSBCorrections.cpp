#include "../interface/SSBCorrections.h"
#include "../TextReader/TextReader.hpp"
#include "../CorrectionFiles/METXY/XYMETCorrection_withUL17andUL18andUL16.h"
#include "../CorrectionFiles/Rochester/RoccoR.cc"
#include <fstream>      
#include <cstdlib>      
#include "correction.h"
#include "TRandom3.h"
#include <cmath>
#include <iostream>
#include <filesystem>
#include <vector>
#include "TLorentzVector.h"



SSBCorrections::SSBCorrections(TextReader* reader, const std::string inputfileName) {
    std::cout << "TextReader in SSBCorrections ! " << std::endl;
    std::cout << "Current directory: " << std::filesystem::current_path() << std::endl;
    reader->PrintoutVariables();


    std::string jsonDir;// = std::filesystem::current_path().string() + "/jsonpog-integration/POG/";
    const std::string cvmfsPath = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/";
    if (std::filesystem::exists(cvmfsPath)) {
        std::cout << "[INFO] Using CVMFS path for JSONs: " << cvmfsPath << std::endl;
        jsonDir = cvmfsPath;
    } else {
        std::string localPath = std::filesystem::current_path().string() + "/jsonpog-integration/POG/";
        std::cout << "[WARNING] CVMFS path not found. Falling back to local path: " << localPath << std::endl;
        jsonDir = localPath;
    }

//    std::string jsonDir = std::filesystem::current_path().string() + "/jsonpog-integration/POG/";

    std::string puw_path     = reader->GetText("PUWeightPath");
    std::string jec_path     = reader->GetText("JECPath");
    std::string jer_path     = reader->GetText("JERPath");
    std::string jmar_path     = reader->GetText("JMARPath");
    std::string jveto_path    = reader->GetText("JetVetoPath");
    std::string jveto_key    = reader->GetText("JetVetoName");
    std::string jer_sf_path  = reader->GetText("JERSFPath");
    std::string jec_name     = reader->GetText("JECName");
    std::string jer_name     = reader->GetText("JERName");
    std::string jer_res_name = reader->GetText("JERResName");
    std::string jet_btag_conf = reader->GetText("Jet_btag");  // e.g., "deepCSVL", "deepJetM"
    std::string muon_path    = reader->GetText("MuonSFPath");
    std::string elec_path    = reader->GetText("ElecSFPath");
    std::string RunPeriod    = reader->GetText("RunRange");
    std::string puJson  =  "";
    btag_sf_type_            = reader->GetText("BTagSFType");

    // Muon Infor. ID  ISO // 
    std::string muon_id_corName  = reader->GetText( "MuonIDSFName"  );
    std::string muon_iso_corName = reader->GetText( "MuonIsoSFName" );
    
    // Electron ID ISO // 
    std::string ele_sf_name_      = reader->GetText( "ElecIDSFName"  );  
    //std::string ele_reco_sf_name_    = reader->GetText( "ElecRecoSFName" );

    // Trigger SF // 
    std::string Trig_sf_name_      = reader->GetText( "TrigSFFile"  );  
    std::string Trig_sf_histname_  = reader->GetText( "TrigSFHist"  );  



    //if inputfileName    

    if (RunPeriod.find("2016PreVFP") != std::string::npos) {
        year_ = "2016preVFP";
        puJson = "Collisions16_UltraLegacy_goldenJSON";
        
    } else if (RunPeriod.find("2016PostVFP") != std::string::npos) {
        year_ = "2016postVFP";
        puJson = "Collisions16_UltraLegacy_goldenJSON";
    } else if (RunPeriod.find("2017") != std::string::npos) {
        year_ = "2017";
        puJson = "Collisions17_UltraLegacy_goldenJSON";
    } else if (RunPeriod.find("2018") != std::string::npos) {
        year_ = "2018";
        puJson = "Collisions18_UltraLegacy_goldenJSON";
    } 
    else {
        std::cerr << "[SSBCorrections] Unknown RunPeriod: " << RunPeriod << std::endl;
        year_ = "2018"; // fallback or throw error
    }


    bool is_data = inputfileName.find("Data") != std::string::npos;

    std::string era = "Unknown";

    if (is_data) {
   	 // For data: extract year and era from the input file name, e.g. "Run2016F-HIPM"
    	size_t run_pos = inputfileName.find("Run");
    	if (run_pos != std::string::npos && run_pos + 7 <= inputfileName.size()) {
        	// Extract year: 4 digits following "Run"
        	std::string year_digits = inputfileName.substr(run_pos + 3, 4); // e.g., "2016", "2017"

        	// Extract era: characters immediately after "RunYYYY"
        	era = inputfileName.substr(run_pos + 7); // e.g., "B", "F-HIPM"
        	size_t delim = era.find_first_of("._/");
        	if (delim != std::string::npos) {
            		era = era.substr(0, delim); // Trim suffix after era
        	}
    	}
    } else {
        // For MC: use the RunPeriod directly (e.g., "2016PreVFP", "2016PostVFP", "2017", "2018")
       era = "MC";  // Dummy placeholder; era is not needed for MC in most JEC logic
    }

    std::cout << "era : " << era <<  std::endl;
    std::cout << "jec_name: " << jec_name << std::endl; 
    jec_name = ExpandJECName(jec_name, RunPeriod, era, is_data); 
    std::cout << "after ExpandJECName jec_name: " << jec_name << std::endl; 
    


    std::cout << "jsonDir + puw_path : " << jsonDir + puw_path << std::endl;
    auto puset = correction::CorrectionSet::from_file(jsonDir + puw_path);
    
    pu_weight_ = puset->at(puJson);


    auto jec_set = correction::CorrectionSet::from_file(jsonDir + jec_path);
    jec_ = jec_set->compound().at(jec_name);

    auto jer_set = correction::CorrectionSet::from_file(jsonDir + jer_path);
    jer_ = jer_set->at(jer_res_name);

    auto jer_sf_set = correction::CorrectionSet::from_file(jsonDir + jer_sf_path);
    jer_sf_ = jer_sf_set->at(jer_name);

    auto jmar_sf_set = correction::CorrectionSet::from_file(jsonDir + jmar_path);
    pujetid_sf_ = jmar_sf_set->at("PUJetID_eff");

//    auto jetveto_set = correction::CorrectionSet::from_file(jsonDir + jveto_path);
//    jetvetomap_ = jmar_sf_set->at("Summer19UL18_V1");

      
    if (btag_sf_type_ != "comb" && btag_sf_type_ != "mujets") {
        std::cerr << "[WARNING] Invalid BTagSFType: " << btag_sf_type_ 
                  << ". Using default 'comb'." << std::endl;
        btag_sf_type_ = "comb";
    }

    // Analysis type specific recommendation
    if (btag_sf_type_ == "comb") {
        std::cout << "[INFO] Using 'comb' SF (QCD + ttbar enriched regions)" << std::endl;
    } else {
        std::cout << "[INFO] Using 'mujets' SF (QCD enriched regions, bias avoidance)" << std::endl;
    }


    std::string btag_algo = "";
    std::string btag_wp = "";

    if (jet_btag_conf.find("deepCSV") != std::string::npos) {
        btag_algo = "DeepCSV";
    } else if (jet_btag_conf.find("deepJet") != std::string::npos) {
        btag_algo = "DeepJet";
    } else if (jet_btag_conf.find("pfCSVV2") != std::string::npos) {
        btag_algo = "CSVv2";
    } else {
        std::cerr << "[WARNING] Unknown b-tag algorithm in Jet_btag: " << jet_btag_conf << std::endl;
    }

    char last = jet_btag_conf.back();
    if (last == 'L' || last == 'l') btag_wp = "Loose";
    else if (last == 'M' || last == 'm') btag_wp = "Medium";
    else if (last == 'T' || last == 't') btag_wp = "Tight";
    else {
        std::cerr << "[WARNING] Unknown WP in Jet_btag: " << jet_btag_conf << std::endl;
    }

    std::string btag_sf_json = "BTV/" + year_ + "_UL/btagging.json.gz";
    std::string btag_tagger = (btag_algo == "DeepJet") ? "deepJet_comb" : "deepCSV_comb";
    std::string btag_eff_path = "CorrectionFiles/BTag/UL" + year_ + "/btagEff_" + btag_algo + ".root";
    std::cout << "btag_eff_path : " << btag_eff_path << std::endl;

    InitBtagSFCorrection(jsonDir + btag_sf_json, btag_tagger);
    //LoadMCBtagEfficiencies(btag_eff_root, btag_algo);
    LoadMCBtagEfficiencies(btag_eff_path, btag_algo);

    // Load muon SF
    //auto muon_set = CorrectionSet::from_file(jsonDir+muon_path);
    std::cout << "jsonDir+muon_path " << jsonDir+muon_path << std::endl; 
    auto muon_set = correction::CorrectionSet::from_file(jsonDir + muon_path);

    std::cout << "muon_id_corName : " << muon_id_corName << std::endl;
    std::cout << "muon_iso_corName : " << muon_iso_corName << std::endl;

    muon_id_   = muon_set->at(muon_id_corName);
    muon_iso_  = muon_set->at(muon_iso_corName);

    auto cset = correction::CorrectionSet::from_file(jsonDir + elec_path);
    //std::cout << "sk ele 1 ele_sf_name_ : " << ele_sf_name_ << std::endl;
    ele_sf_ = cset->at(ele_sf_name_);
    //std::cout << "sk ele 2 ele_reco_sf_name_: " << ele_reco_sf_name_ << std::endl;
    std::string TrigSFPath = std::filesystem::current_path().string() + "/CorrectionFiles/Trig/";
    TFile *f_trg     = new TFile((TrigSFPath+Trig_sf_name_).c_str());
    H_trig = (TH2D*) f_trg->Get(Trig_sf_histname_.c_str()); 
}

double SSBCorrections::GetCorrectedJetPt(double raw_pt, double eta, double area, double rho) const {
    //std::cout << "in GetCorrectedJetPt raw_pt : " << raw_pt
    //          << " eta " << eta << " area : " << area << " rho : " << rho << std::endl;

    double sf = jec_->evaluate({area, eta, raw_pt, rho});  
    //std::cout << "sf in GetCorrectedJetPt : " << sf << std::endl;

    return raw_pt * sf;
}

double SSBCorrections::GetCorrectedJetMass(double raw_mass, double raw_pt, double eta, double area, double rho) const {
    //double sf = jec_->evaluate({eta, raw_pt, area});
    double sf = jec_->evaluate({area, eta, raw_pt, rho});
    return raw_mass * sf;
}

double SSBCorrections::GetJER(double eta, double pt) const {
    return jer_->evaluate({eta, pt});
}

double SSBCorrections::SmearJER(double reco_pt, double gen_pt, double eta, double rho, const std::string& jer_tag) const {
    double sf = jer_sf_->evaluate({eta, jer_tag});
    double resolution = jer_->evaluate({eta, reco_pt, rho});

    if (gen_pt > 0.0) {
        double delta_pt = reco_pt - gen_pt;
        double smeared_pt = gen_pt + sf * delta_pt;
        return std::max(0.0, smeared_pt);
    }

    double smear_factor = 1.0 + std::sqrt(sf * sf - 1.0) * gRandom->Gaus(0, 1) * resolution;
    return std::max(0.0, reco_pt * smear_factor);
}

TLorentzVector SSBCorrections::RecomputeMET(double raw_met_pt, double raw_met_phi,
                                            const std::vector<TLorentzVector>& rawJets,
                                            const std::vector<TLorentzVector>& corrJets) const {
    // Convert raw MET to px, py components
    double met_px = raw_met_pt * std::cos(raw_met_phi);
    double met_py = raw_met_pt * std::sin(raw_met_phi);

    // Calculate correction in px, py only
    for (size_t i = 0; i < rawJets.size(); ++i) {
        // Add back the difference in transverse momentum
        // (raw jet was subtracted from original MET, now we subtract corrected jet)
        met_px += (rawJets[i].Px() - corrJets[i].Px());
        met_py += (rawJets[i].Py() - corrJets[i].Py());
    }

    // Create corrected MET with mass = 0
    double corrected_met_pt = std::sqrt(met_px * met_px + met_py * met_py);
    double corrected_met_phi = std::atan2(met_py, met_px);

    TLorentzVector correctedMET;
    correctedMET.SetPtEtaPhiM(corrected_met_pt, 0, corrected_met_phi, 0);

    return correctedMET;
}

float SSBCorrections::GetPUJetIDSFAndEff(float pt, float eta, bool passPU, bool genMatched, const std::string& wp, const std::string& syst, bool getEff) const {
    if (!pujetid_sf_) {
        std::cerr << "[SSBCorrections::GetPUJetIDSFAndEff] PUJetID correction not loaded." << std::endl;
        return 1.0;
    }
    
    if (pt >= 50.0) return 1.0;
    
    try {
        std::string eval_type = getEff ? "MCEff" : syst;  // Efficiency or Scale Factor
        std::variant<double, std::vector<double>> val = pujetid_sf_->evaluate({eta, pt, eval_type, wp});
        return std::get<double>(val);
    } catch (const std::exception& e) {
        std::cerr << "[SSBCorrections::GetPUJetIDSFAndEff] Evaluation failed: " << e.what() << std::endl;
        return 1.0;
    }
}


double SSBCorrections::GetMuonRecoSF(double pt, double eta) const {
    //auto val = muon_reco_->evaluate({pt, eta});
    std::variant<double, std::vector<double>> val = muon_reco_->evaluate({pt, eta});
    //return std::get<double>(val);
    return std::get<double>(val);
}

double SSBCorrections::GetMuonIDSF(double pt, double eta, const std::string& syst_tag) const {
    //auto val = muon_id_->evaluate({eta, pt, syst_tag});
    std::variant<double, std::vector<double>> val = muon_id_->evaluate({eta, pt, syst_tag});
    return std::get<double>(val);
}

double SSBCorrections::GetMuonIsoSF(double pt, double eta, const std::string& syst_tag) const {
    //auto val = muon_iso_->evaluate({eta, pt, syst_tag});
    std::variant<double, std::vector<double>> val = muon_iso_->evaluate({eta, pt, syst_tag});
    return std::get<double>(val);
}


double SSBCorrections::DoubleMuon_IDIsoEff(TLorentzVector lep1, TLorentzVector lep2,
                                           TString muidsys, TString muisosys, TString tracksys) const {
    float pt1 = std::min(lep1.Pt(), 119.999);
    float pt2 = std::min(lep2.Pt(), 119.999);
    float abseta1 = std::abs(lep1.Eta());
    float abseta2 = std::abs(lep2.Eta());

    //std::cout << "muidsys : " << muidsys << std::endl;
    //std::cout << "muisosys : " << muisosys << std::endl;

    std::string IDSyst = "nominal", IsoSyst = "nominal";

    if (muidsys.Contains("up", TString::kIgnoreCase)) IDSyst = "up";
    else if (muidsys.Contains("down", TString::kIgnoreCase)) IDSyst = "down";

    if (muisosys.Contains("up", TString::kIgnoreCase)) IsoSyst = "up";
    else if (muisosys.Contains("down", TString::kIgnoreCase)) IsoSyst = "down";
    //std::cout << "IDSyst : " << IDSyst << std::endl;
    double mu1id  = GetMuonIDSF(pt1, abseta1, IDSyst);
    //std::cout << "After mu1id  " << mu1id << std::endl;
    double mu2id  = GetMuonIDSF(pt2, abseta2, IDSyst);
    double mu1iso = GetMuonIsoSF(pt1, abseta1, IsoSyst);
    double mu2iso = GetMuonIsoSF(pt2, abseta2, IsoSyst);

    // track SF is 1.0, so you don't need to apply them... 
    double mu1trk =1.0;// TrackSF(lep1->Eta());
    double mu2trk =1.0;// TrackSF(lep2->Eta());
/*
    if (tracksys.Contains("up", TString::kIgnoreCase)) {
        mu1trk += TrackSFErr(lep1->Eta(), tracksys);
        mu2trk += TrackSFErr(lep2->Eta(), tracksys);
    } else if (tracksys.Contains("down", TString::kIgnoreCase)) {
        mu1trk -= TrackSFErr(lep1->Eta(), tracksys);
        mu2trk -= TrackSFErr(lep2->Eta(), tracksys);
    }
*/
    return mu1id * mu2id * mu1iso * mu2iso * mu1trk * mu2trk;
}

float SSBCorrections::GetElectronIDSF(float pt, float eta, const std::string& wp) const {
    if (!ele_sf_) {
        std::cerr << "[GetElectronIDSF] SF not loaded! Call LoadElectronSF() first.\n";
        return 1.0;
    }

    //float id_sf = ele_sf_->evaluate({ year_, wp, "nominal", eta, pt });
    //float reco_sf = ele_reco_sf_->evaluate({ year_, "nominal", eta, pt });
    return ele_sf_->evaluate({ year_, wp, "nominal", eta, pt });;
}

float SSBCorrections::GetElectronRecoSF(float pt, float eta) const {
    if (!ele_reco_sf_) {
        std::cerr << "[GetElectronRecoSF] SF not loaded! Call LoadElectronSF() first.\n";
        return 1.0;
    }

    //float id_sf = ele_sf_->evaluate({ year_, wp, "nominal", eta, pt });
    //float reco_sf = ele_reco_sf_->evaluate({ year_, "nominal", eta, pt });
    return ele_reco_sf_->evaluate({ year_, "nominal", eta, pt });
}


float SSBCorrections::GetPUWeight(float nTrueInt, const std::string& systTag) const {
    std::string variation = systTag;
    
    if (!pu_weight_) {
        std::cerr << "[SSBCorrections::GetPUWeight] PU correction not loaded.\n";
        return 1.0;
    }
    
    // If systTag is "Central" or empty, treat as "nominal"
    if (systTag == "Central" || systTag == "") {
        variation = "nominal";
    }
    
    //std::cout << "variation in PU " << variation << std::endl;

    // If the variation is not valid, fallback to "nominal"
    if (variation != "nominal" && variation != "up" && variation != "down") {
        std::cerr << "PileUp sys Error ... Defaulting to Weight_PileUp ... Original value: " << variation << std::endl;
        variation = "nominal";
    }
    
    //std::cout << "nTrueInt : " << nTrueInt << std::endl;
    
    try {
        std::variant<double, std::vector<double>> val = pu_weight_->evaluate({ nTrueInt, variation });
        double weight = std::get<double>(val);
        //std::cout << "PileUp weight: " << weight << std::endl;  
        return weight;
    } catch (const std::exception& e) {
        std::cerr << "[SSBCorrections::GetPUWeight] Evaluation failed: " << e.what() << std::endl;
        return 1.0;
    }
}

float SSBCorrections::GetElectronSF(const std::string& sf_type, float eta, float pt, const std::string& syst) const {
    std::string recoPtThreshold = pt >= 20 ? "RecoAbove20" : "RecoBelow20";
    std::string type = sf_type;

    if (sf_type == "Reco") type = recoPtThreshold;

    try {
        return ele_sf_->evaluate({year_, syst, type, eta, pt});
    } catch (const std::exception& e) {
        std::cerr << "[GetElectronSF] Failed to evaluate: " << e.what() << std::endl;
        return 1.0;
    }
}

float SSBCorrections::MatchGenPt(const TLorentzVector& reco_jet,
                                  const std::vector<TLorentzVector>& gen_jets,
                                  float maxDR) const {
    float minDR = maxDR;
    float matched_genpt = -1.0;

    for (const auto& gen_jet : gen_jets) {
        float dR = reco_jet.DeltaR(gen_jet);
        if (dR < minDR) {
            minDR = dR;
            matched_genpt = gen_jet.Pt();
        }
    }
    return matched_genpt;
}

JetCorrectionOutput SSBCorrections::ApplyJetCorrectionsWithMET(
    const std::vector<TLorentzVector>& rawJets,
    const std::vector<float>& rawFactors,
    const std::vector<float>& areas,
    float rho,
    bool isData,
    bool applyJES,
    bool applyJER,
    double raw_met_pt,
    double raw_met_phi,
    const std::vector<TLorentzVector>& genJets,
    const std::vector<int>& genJetIndices
) const {
    std::vector<TLorentzVector> correctedJets;
    std::vector<TLorentzVector> rawJetsRebuilt;

    correctedJets.reserve(rawJets.size());
    rawJetsRebuilt.reserve(rawJets.size());

    for (size_t i = 0; i < rawJets.size(); ++i) {
        double eta  = rawJets[i].Eta();
        double phi  = rawJets[i].Phi();

        double raw_pt   = rawJets[i].Pt() * (1.0 - rawFactors[i]);
        double raw_mass = rawJets[i].M()  * (1.0 - rawFactors[i]);

        TLorentzVector raw_jet;
        raw_jet.SetPtEtaPhiM(raw_pt, eta, phi, raw_mass);
        rawJetsRebuilt.push_back(raw_jet);

        double corrected_pt = raw_pt;
        double corrected_mass = raw_mass;

        if (applyJES) {
            corrected_pt   = GetCorrectedJetPt(raw_pt, eta, areas[i], rho);
            corrected_mass = GetCorrectedJetMass(raw_mass, raw_pt, eta, areas[i], rho);
        }

        if (!isData && applyJER && genJetIndices.size() > i && genJetIndices[i] >= 0 && genJetIndices[i] < genJets.size()) {
            float matched_genpt = genJets[genJetIndices[i]].Pt();
            corrected_pt = SmearJER(corrected_pt, matched_genpt, eta, rho, "nominal");
        }

        if (raw_pt > 0) {
            corrected_mass *= (corrected_pt / raw_pt);
        }

        TLorentzVector corr_jet;
        corr_jet.SetPtEtaPhiM(corrected_pt, eta, phi, corrected_mass);
        correctedJets.push_back(corr_jet);
    }

    TLorentzVector correctedMET = RecomputeMET(raw_met_pt, raw_met_phi, rawJetsRebuilt, correctedJets);
    JetCorrectionOutput result;

    result.corrected_jets = correctedJets;
    result.corrected_met  = correctedMET;

    return result;
}

/*void SSBCorrections::InitBtagSFCorrection(const std::string& json_path, const std::string& tagger_name) {
    std::cout << "[SSBCorrections] Loading b-tagging SF from JSON: " << json_path << std::endl;
    
    auto cset = correction::CorrectionSet::from_file(json_path);
    
    // Create correction names
    std::string comb_name = tagger_name; // e.g., "deepCSV_comb"
    std::string incl_name = tagger_name;
    incl_name.replace(incl_name.find("comb"), 4, "incl"); // "deepCSV_incl"
    
    // Store corrections in map
    btag_corrections_["comb"] = cset->at(comb_name);
    btag_corrections_["incl"] = cset->at(incl_name);
    
    std::cout << "[SSBCorrections] Loaded corrections: " << comb_name 
              << " and " << incl_name << std::endl;
}*/
void SSBCorrections::InitBtagSFCorrection(const std::string& json_path, 
                                          const std::string& tagger_name) {
    std::cout << "[SSBCorrections] Loading b-tagging SF from JSON: " << json_path << std::endl;
    
    auto cset = correction::CorrectionSet::from_file(json_path);
    
    // Heavy flavor correction name (config-based selection)
    std::string heavy_flavor_name = tagger_name;  // e.g., "deepJet_comb"
    
    // Replace "comb" with config-specified type
    size_t pos = heavy_flavor_name.find("comb");
    if (pos != std::string::npos) {
        heavy_flavor_name.replace(pos, 4, btag_sf_type_);
    } else {
        std::cerr << "[ERROR] Expected 'comb' in tagger_name: " << tagger_name << std::endl;
        return;
    }
    
    // Light flavor correction name (always "incl")
    std::string light_flavor_name = tagger_name;
    pos = light_flavor_name.find("comb");
    if (pos != std::string::npos) {
        light_flavor_name.replace(pos, 4, "incl");
    }
    
    // Load corrections with generic keys
    btag_corrections_["heavy"] = cset->at(heavy_flavor_name);
    btag_corrections_["light"] = cset->at(light_flavor_name);
    
    std::cout << "[SSBCorrections] Loaded corrections: " 
              << heavy_flavor_name << " and " << light_flavor_name << std::endl;
}


/*std::string SSBCorrections::getBtagCorrectionName(int flavor) const {
    return (flavor == 0) ? "incl" : "comb";
}*/
std::string SSBCorrections::getBtagCorrectionName(int flavor) const {
    // 0 = light flavor (u,d,s,g), others = heavy flavor (b,c)
    return (flavor == 0) ? "light" : "heavy";
}

/*float SSBCorrections::GetBtagSF(float pt, float eta, int flav, const std::string& wp, const std::string& syst) const {
    float pt_clamped = std::clamp(pt, 20.0f, 1000.0f);
    float eta_abs = std::fabs(eta);

    std::string corr_name = getBtagCorrectionName(flav);
    auto it = btag_corrections_.find(corr_name);
    
    if (it == btag_corrections_.end()) {
        std::cerr << "[ERROR] B-tag correction not found: " << corr_name << std::endl;
        return 1.0;
    }

    try {
        //std::cout << "[DEBUG] Using " << corr_name << " correction for flavor " << flav << std::endl;
        return it->second->evaluate({syst, wp, flav, eta_abs, pt_clamped});
    } catch (const std::exception& e) {
        std::cerr << "[WARNING] GetBtagSF failed: " << e.what() << std::endl;
        return 1.0;
    }
}*/

float SSBCorrections::GetBtagSF(float pt, float eta, int flav, 
                                const std::string& wp, const std::string& syst) const {
    float pt_clamped = std::clamp(pt, 20.0f, 1000.0f);
    float eta_abs = std::fabs(eta);

    std::string corr_name = getBtagCorrectionName(flav);  // "heavy" or "light"
    auto it = btag_corrections_.find(corr_name);
    
    if (it == btag_corrections_.end()) {
        std::cerr << "[ERROR] B-tag correction not found: " << corr_name << std::endl;
        return 1.0;
    }

    try {
        return it->second->evaluate({syst, wp, flav, eta_abs, pt_clamped});
    } catch (const std::exception& e) {
        std::cerr << "[WARNING] GetBtagSF failed: " << e.what() << std::endl;
        return 1.0;
    }
}

float SSBCorrections::ComputeBTagEventWeight(const std::vector<float>& pts,
                                             const std::vector<float>& etas,
                                             const std::vector<int>& flavs,
                                             const std::vector<bool>& isTagged,
                                             const std::string& algo,
                                             const std::string& wp,
                                             const std::string& syst) const {
    
    // Check input vector sizes
    if (pts.size() != etas.size() || 
        pts.size() != flavs.size() || 
        pts.size() != isTagged.size()) {
        std::cout << "Input vectors have different sizes" << std::endl;
        return 1.0;
    }
    //std::cout << "start ! " << std::endl; 
    // Initialize probabilities
    float p_mc = 1.0;
    float p_data = 1.0;
    
    // Loop over all jets
    for (size_t i = 0; i < pts.size(); ++i) {
        float pt   = pts[i];
        float eta  = etas[i];
        int flav   = flavs[i];
        bool tagged = isTagged[i];
         
        // Get MC efficiency
        float eff = GetMCBtagEfficiency(pt, eta, flav, algo, wp);
        
        // Get scale factor
        float sf = GetBtagSF(pt, eta, flav, wp, syst);
        
        // Avoid division by zero
        if (eff < 1e-5) eff = 1e-5;
        if (eff > 1.0 - 1e-5) eff = 1.0 - 1e-5;
        
        // Calculate probabilities
        if (tagged) {
            // Tagged jet contribution
            p_mc *= eff;
            p_data *= (eff * sf);
        } else {
            // Not tagged jet contribution
            p_mc *= (1.0f - eff);
            p_data *= (1.0f - eff * sf);
        }
        
        /*std::cout << "  Jet " << i << ": pt=" << pt << ", eta=" << eta 
                  << ", flav=" << flav << ", tagged=" << tagged 
                  << ", SF=" << sf << ", Eff=" << eff 
                  << ", p_mc=" << p_mc << ", p_data=" << p_data << std::endl;*/
    }
    
    // Calculate final event weight
    float btag_evt_weight = (p_mc > 1e-10) ? p_data / p_mc : 1.0;
    
    return btag_evt_weight;
}

void SSBCorrections::LoadMCBtagEfficiencies(const std::string& filepath, const std::string& algo) {
    TFile* f = TFile::Open(filepath.c_str(), "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "[ERROR] Failed to open efficiency file: " << filepath << std::endl;
        return;
    }

    for (const std::string& flav : {"b", "c", "l"}) {
        for (const std::string& wp : {"Loose", "Medium", "Tight"}) {
            std::string name = "eff_" + algo + "_" + flav + "_" + wp;
	    //std::cout << "name in LoadMCBtagEfficiencies : " << name << std::endl;
            TH2D* hist = (TH2D*)f->Get(name.c_str());
            if (hist) {
                hist->SetDirectory(0); // Detach histogram from file
                eff_histograms_[algo + "_" + flav + "_" + wp] = hist;
                //std::cout << "[INFO] Loaded hist: " << name << std::endl;
            } else {
                std::cerr << "[WARNING] Histogram not found: " << name << std::endl;
            }
        }
    }
    f->Close();
}

float SSBCorrections::GetMCBtagEfficiency(float pt, float eta, int flav, const std::string& algo, const std::string& wp) const {
    std::string flav_str = "l";
    if (flav == 5) flav_str = "b";
    else if (flav == 4) flav_str = "c";

    // Convert single character WP to full name for efficiency lookup
    std::string wp_full = wp;
    if (wp == "L") wp_full = "Loose";
    else if (wp == "M") wp_full = "Medium";
    else if (wp == "T") wp_full = "Tight";

    std::string key = algo + "_" + flav_str + "_" + wp_full;
    //std::cout << "[Info] key in GetMCBtagEfficiency " << key << std::endl;
    auto it = eff_histograms_.find(key);
    if (it == eff_histograms_.end()) {
        std::cerr << "[WARNING] Efficiency hist not found: " << key << std::endl;
        return 1.0;
    }

    TH2D* hist = it->second;
    int bin_x = hist->GetXaxis()->FindBin(pt);
    int bin_y = hist->GetYaxis()->FindBin(eta);
    float eff = hist->GetBinContent(bin_x, bin_y);

    return std::clamp(eff, 0.0f, 1.0f);
}

double SSBCorrections::RochesterCorrectionData(TString year, int Q, double pt, double eta, double phi, int s,int m) const{

    RoccoR rc;
    std::string correctionFile;
    if (year.Contains("2016Pre"))  correctionFile =  "./CorrectionFiles/Rochester/RoccoR2016aUL.txt";
    if (year.Contains("2016Post")) correctionFile =  "./CorrectionFiles/Rochester/RoccoR2016bUL.txt";  
    if (year.Contains("2017"))     correctionFile =  "./CorrectionFiles/Rochester/RoccoR2017UL.txt";
    if (year.Contains("2018"))     correctionFile =  "./CorrectionFiles/Rochester/RoccoR2018UL.txt";

    std::ifstream file(correctionFile);
    if(!file.good()){
	    std::cerr << "ERROR:  Rochester file is not found:" << correctionFile << std::endl; 
    	    std::exit(1);
    }    

    rc.init(correctionFile);

    double correction;
    correction = rc.kScaleDT(Q,pt,eta,phi,s,m); return correction;
}

double SSBCorrections::RochesterCorrectionMC(TString year, int Q, double pt, double eta,double phi,int genID,double genPt,int nl, int s,int m) const{

    RoccoR rc;
    std::string correctionFile;
    if (year.Contains("2016Pre"))  correctionFile = "./CorrectionFiles/Rochester/RoccoR2016aUL.txt";
    if (year.Contains("2016Post")) correctionFile = "./CorrectionFiles/Rochester/RoccoR2016bUL.txt";  
    if (year.Contains("2017"))     correctionFile = "./CorrectionFiles/Rochester/RoccoR2017UL.txt";
    if (year.Contains("2018"))     correctionFile = "./CorrectionFiles/Rochester/RoccoR2018UL.txt";

    std::ifstream file(correctionFile);
    if(!file.good()){
            std::cerr << "ERROR:  Rochester file is not found:" << correctionFile << std::endl;
            std::exit(1);
    }

    rc.init(correctionFile);

    double correction; double u;

    bool genMatch = genID != -1;

    if(genMatch) correction = rc.kSpreadMC(Q,pt,eta,phi,genPt,s,m);
    else if(!genMatch){u = gRandom->Rndm(); correction = rc.kSmearMC(Q,pt,eta,phi,nl,u,s,m);} //Random number is needed when gen-mathcing is failed
    return correction;
}


TLorentzVector SSBCorrections::METXYCorrection(const TLorentzVector& type1_met,
                                               int runnb, TString year, bool isMC, int npv, bool isUL, bool ispuppi
                                               ) const {


    std::pair<double, double> correctedMET = METXYCorr_Met_MetPhi(type1_met.Pt(),type1_met.Phi(),runnb,year,isMC,npv,isUL,ispuppi);

    double met_pt  = correctedMET.first;
    double met_phi = correctedMET.second;

    TLorentzVector corrected_met;
    corrected_met.SetPtEtaPhiM(met_pt,0,met_phi,0);
    return corrected_met;
}


TLorentzVector SSBCorrections::METXYCorrection_corrlib(const TLorentzVector& type1_met,
                                               const std::string& era,
                                               bool isData,
                                               int npv) const {
    if (!metphi_corr_) {
        std::cerr << "[SSBCorrections::METXYCorrection] MET correction object not loaded.\n";
        return type1_met;
    }

    std::string data_tag = isData ? "data" : "mc";

    double corr_x = 0.0;
    double corr_y = 0.0;

    try {
        std::variant<double, std::vector<double>> val_x = metphi_corr_->evaluate({era, data_tag, npv, "x"});
        std::variant<double, std::vector<double>> val_y = metphi_corr_->evaluate({era, data_tag, npv, "y"});

        if (std::holds_alternative<double>(val_x)) {
            corr_x = std::get<double>(val_x);
        } else {
            std::cerr << "[METXYCorrection] Warning: unexpected type for x correction\n";
        }

        if (std::holds_alternative<double>(val_y)) {
            corr_y = std::get<double>(val_y);
        } else {
            std::cerr << "[METXYCorrection] Warning: unexpected type for y correction\n";
        }
    } catch (const std::exception& e) {
        std::cerr << "[SSBCorrections::METXYCorrection] Evaluation failed: " << e.what() << std::endl;
        return type1_met;
    }

    double met_x = type1_met.Px() - corr_x;
    double met_y = type1_met.Py() - corr_y;

    TLorentzVector corrected_met;
    corrected_met.SetPxPyPzE(met_x, met_y, 0, std::sqrt(met_x * met_x + met_y * met_y));
    return corrected_met;
}
/// Trigger SF 
double SSBCorrections::TrigDiMuon_Eff(TLorentzVector lep1, TLorentzVector lep2, TString Sys_) {
    double pt1 = lep1.Pt();
    double pt2 = lep2.Pt();

    double leading_pt = std::max(pt1, pt2);
    double subleading_pt = std::min(pt1, pt2);

    return GetTrgEff(leading_pt, subleading_pt, Sys_);
}

double SSBCorrections::TrigDiElec_Eff(TLorentzVector lep1, TLorentzVector lep2, TString Sys_) {
    double pt1 = lep1.Pt();
    double pt2 = lep2.Pt();

    double leading_pt = std::max(pt1, pt2);
    double subleading_pt = std::min(pt1, pt2);

    return GetTrgEff(leading_pt, subleading_pt, Sys_);
}

double SSBCorrections::TrigMuElec_Eff(TLorentzVector lep1, TLorentzVector lep2, TString Sys_) {
    double pt_mu = lep1.Pt();   // Assuming lep1 is muon
    double pt_ele = lep2.Pt();  // Assuming lep2 is electron

    double leading_pt = std::max(pt_mu, pt_ele);
    double subleading_pt = std::min(pt_mu, pt_ele);

    return GetTrgEff(leading_pt, subleading_pt, Sys_);
}

double SSBCorrections::GetTrgEff(double pt1, double pt2, TString Sys_) {
    // Clamp values below the histogram max range (assumed 500)
    double max_pt = 499.999;

    pt1 = std::min(pt1, max_pt);
    pt2 = std::min(pt2, max_pt);

    int xbin = H_trig->GetXaxis()->FindBin(pt1);
    int ybin = H_trig->GetYaxis()->FindBin(pt2);

    double trgsf = H_trig->GetBinContent(xbin, ybin);
    double trgsferr = 0.0;

    if (Sys_ == "nominal" || Sys_ == "Central") {trgsferr = 0.0;}
    else if (Sys_ == "Up" || Sys_ == "up")
        trgsferr = H_trig->GetBinError(xbin, ybin);
    else if (Sys_ == "Down" || Sys_ == "down")
        trgsferr = -H_trig->GetBinError(xbin, ybin);

    return trgsf + trgsferr;
}

std::string SSBCorrections::ExpandJECName(const std::string& base_jec_name, const std::string runPeriod, const std::string& era, bool is_data) {
    //  "Summer19UL16_v7" -> prefix = "Summer19UL16", version = "7"
    size_t pos = base_jec_name.find("_V");
    if (pos == std::string::npos) {
        std::cerr << "[ERROR] Invalid jec_name format: " << base_jec_name << std::endl;
        return base_jec_name; // fallback
    }

    std::string prefix = base_jec_name.substr(0, pos);             // "Summer19UL16"
    std::string version = base_jec_name.substr(pos + 2);           // "7"
    std::string suffix = "_L1L2L3Res_AK4PFchs";

    std::string expanded = prefix;  // prefix

    // Year-specific logic
    if (runPeriod.find("16Pre") != std::string::npos) {
        if (is_data) {
            if (era == "B" || era == "C" || era == "D") {
                expanded += "_RunBCD";
            } else if (era == "E") {
                expanded += "_RunEF";
            } else if (era.find("F") == 0) {
                expanded += "_RunEF";
            } else if (era == "G" || era == "H") {
                expanded += "_RunGH";
            } else {
                expanded += "_RunFGH";  // fallback
            }
        } else {
            expanded += (era == "B" || era == "C" || era == "D" || era == "E") ? "APV" : "";
        }
    } else if (runPeriod.find("16Post") != std::string::npos){
        if (is_data) {
	    if (era.find("F") == 0) {
                expanded += "_RunFGH";
            } else if (era == "G" || era == "H") {
                expanded += "_RunFGH";
            } else {
                expanded += "_RunFGH";  // fallback
            }
        } else {
            expanded += (era == "F" || era == "G" || era == "H") ? "APV" : "";
        }

    } else if (runPeriod.find("2017") != std::string::npos) {
        if (is_data) {
            if (era == "B" ) {
                expanded += "_RunB";
            } else if (era == "C") {
                expanded += "_RunC";
            } else if (era == "D" ) {
                expanded += "_RunD";
            } else if (era == "E") {
                expanded += "_RunE";
            } else if (era == "F") {
                expanded += "_RunF"; 
            } else {
                expanded += "_RunB"; 
            }
        
        }
    } else if (runPeriod.find("2018") != std::string::npos) {
        if (is_data) {
            if (era == "A" ) {
                expanded += "_RunA";
            } else if (era == "B") {
                expanded += "_RunB";
            } else if (era == "C" ) {
                expanded += "_RunC";
            } else if (era == "D") {
                expanded += "_RunD";
            } else {
                expanded += "_RunB"; 
            }
        }
    }

    // Add version and data/MC tag
    expanded += "_V" + version + "_";
    expanded += is_data ? "DATA" : "MC";
    expanded += suffix;

    return expanded;
}
