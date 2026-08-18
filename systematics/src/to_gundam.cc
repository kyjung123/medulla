#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TH1D.h>
#include <TGraph.h>
#include <TClonesArray.h>
#include <TChain.h>

#include <vector>
#include <string>
#include <iostream>
#include <unordered_set>
#include <algorithm>
#include <map>
#include <cmath>
#include <sstream>

#include "configuration.h"
#include "trees.h"

static TDirectory* ensure_dir_path(TFile* f, const std::string& path){
  if(path.empty()) return (TDirectory*)f;

  TDirectory* dir = (TDirectory*)f;
  std::string token;
  std::stringstream ss(path);

  while(std::getline(ss, token, '/')){
    if(token.empty()) continue;

    TDirectory* next = dir->GetDirectory(token.c_str());
    if(!next) next = dir->mkdir(token.c_str());
    dir = next;
  }
  return dir;
}

static TDirectory* get_dir_path(TFile* f, const std::string& path){
  if(path.empty()) return (TDirectory*)f;
  return f->GetDirectory(path.c_str());
}

static bool starts_with(const std::string& s, const std::string& prefix){
  return s.rfind(prefix, 0) == 0;
}

// ------------------------------------------------------------
// tree name -> cut_type
// 여기만 바꾸면 allowed nominal tree + cut_type 둘 다 같이 바뀜
// ------------------------------------------------------------
static const std::map<std::string, int> kTreeCutType = {
//  {"signal", 0},
//  {"signal_multiple", 1},
//  {"selected", 0},
//  {"selected_multiple", 1},
//  {"selected_renzo_low_threshold", 2},
//  {"selected_michel", 3}
//  {"signal_one_proton", 0},
//  {"signal", 0},
  {"signal_N_proton", 0},
  {"selected_one_proton_michel", 0},
  {"selected_renzo_low_threshold", 1},
  {"selected_multiple", 2},
};

static bool is_allowed_selected_tree(const std::string& tname){
  if(!starts_with(tname, "selected")) return false;
  return kTreeCutType.count(tname) > 0;
}

static bool is_allowed_signal_tree(const std::string& tname){
  if(!starts_with(tname, "signal")) return false;
  return kTreeCutType.count(tname) > 0;
}

static bool is_allowed_nominal_tree(const std::string& tname, const std::string& prefix){
  if(prefix == "selected") return is_allowed_selected_tree(tname);
  if(prefix == "signal")   return is_allowed_signal_tree(tname);
  return false;
}

static bool is_syst_tree_name(const std::string& tname){
  return (tname.find("_multisigmaTree")  != std::string::npos) ||
         (tname.find("_variationTree")   != std::string::npos) ||
         (tname.find("_NuMIfluxsimTree") != std::string::npos) ||
         (tname.find("_fluxsimTree")     != std::string::npos) ||
         (tname.find("_multisimTree")    != std::string::npos);
//         (tname.find("_fluxsimTree")     != std::string::npos) ;
}

// ------------------------------------------------------------
// systematic tree finder
// nominal_name + suffix 후보들을 순서대로 확인
// 있으면 return, 없으면 nullptr
// ------------------------------------------------------------
static TTree* find_syst_tree(TDirectory* dir,
                             const std::string& nominal_name,
                             const std::vector<std::string>& candidates)
{
  if(!dir) return nullptr;

  for(const auto& suf : candidates){
    std::string name = nominal_name + suf;
    TObject* obj = dir->Get(name.c_str());
    if(obj && obj->InheritsFrom(TTree::Class())){
      return (TTree*)obj;
    }
  }
  return nullptr;
}

static const std::vector<std::string> kMultisigmaSuffixes = {
  "_multisigmaTree"
};

static const std::vector<std::string> kVariationSuffixes = {
  "_variationTree"
};

static const std::vector<std::string> kFluxSuffixes = {
  "_NuMIfluxsimTree",
  "_fluxsimTree"
};

static const std::vector<std::string> kMultisimSuffixes = {
  "_multisimTree"
//  "_asdfasdf"
};

// ------------------------------------------------------------
// cut_type mapping
// ------------------------------------------------------------
static int cut_type_from_tree_name_safe(const std::string& tname){
  auto it = kTreeCutType.find(tname);
  if(it != kTreeCutType.end()) return it->second;

  std::cout << "WARN: cut_type_from_tree_name_safe: unknown tree name '"
            << tname << "' -> cut_type=99\n";
  return 99;
}

// ------------------------------------------------------------
// Merge + syst를 한 이벤트 루프에서 out_tree->Fill()
// ------------------------------------------------------------
struct MergeWriter {
  int run=0, subrun=0, event=0;
  int cut_type=0, is_nu=0, is_data=0, category=0;

  double _category = -5;

  bool use_pick_existing_hadron_branch = true;

  std::vector<std::string> br_names;
  std::vector<double>      br_vals;

  double reco_had_cos_src=-9999, true_had_cos_src=-9999, reco_had_ke_src=-9999, true_had_ke_src=-9999;
  float ppfx_cv_weight=1.0;

  double reco_leading_hadron_cos_wrt_beam=-9999;
  double true_leading_hadron_cos_wrt_beam=-9999;
  double reco_leading_hadron_ke=-9999;
  double true_leading_hadron_ke=-9999;

  bool booked_nominal=false;

  struct SystParam {
    std::string full;
    std::string out;
    TClonesArray* arr=nullptr;
    bool export_flat_sigma=false;
    double flat_m1 = 1.0;
    double flat_p1 = 1.0;
    bool export_single_weight=false;
    std::string single_weight_name;
    double single_weight = 1.0;

    std::vector<Double_t>* wD = nullptr;
    std::vector<Float_t>*  wF = nullptr;

    std::vector<Double_t>* xD = nullptr;
    std::vector<Float_t>*  xF = nullptr;

    std::string wType;
    std::string xType;

    bool present_in_this_tree=false;
  };

  struct SystTypeBlock {
    std::string type;
    std::vector<SystParam> params;
    bool booked=false;
  };

  std::map<std::string, SystTypeBlock> syst_blocks;

  std::unordered_set<std::string> flat_sigma_names = {
    "ZExpA1CCQE",
    "ZExpA2CCQE",
    "ZExpA3CCQE",
    "ZExpA4CCQE",
    "VecFFCCQEshape",
    "RPA_CCQE",
    "CoulombCCQE",
    "NormCCMEC",
    "NormNCMEC",
    "DecayAngMEC",
    "MaNCEL",
    "EtaNCEL",
    "MaCCRES",
    "MvCCRES",
    "MaNCRES",
    "MvNCRES"
  };

  const std::map<std::string, std::string> single_weight_names = {
    {
      "GENIEReWeight_SBNNuSyst_LQCDZExpFit_correction_ZExpAVariationResponse",
      "LQCDZExpFitCorrection_weight"
    },
    {
      "GENIEReWeight_SBNNuSyst_MINERvAZExpFit_correction_ZExpAVariationResponse",
      "MINERvAZExpFitCorrection_weight"
    }
  };

  std::vector<std::string> morph_names = {
    "GENIEReWeight_SBN_v1_multisigma_VecFFCCQEshape",
    "GENIEReWeight_SBN_v1_multisigma_DecayAngMEC",
    "GENIEReWeight_SBN_v1_multisigma_Theta_Delta2Npi",
    "GENIEReWeight_SBN_v1_multisigma_ThetaDelta2NRad"
  };

  static const char* pick_existing(TTree* in_tree, const char* pion, const char* proton){
    bool hasPion   = (in_tree->GetBranch(pion)   != nullptr);
    bool hasProton = (in_tree->GetBranch(proton) != nullptr);
    if (hasPion && !hasProton) return pion;
    if (!hasPion && hasProton) return proton;
    if (hasPion && hasProton)  return pion;
    return nullptr;
  }

  void book_nominal_once(TTree* out_tree, TTree* first_in_tree){
    if(booked_nominal) return;

    std::unordered_set<std::string> skip = {
      "Run","Subrun","Evt",
      "reco_cut_type","reco_is_nu","reco_is_data","true_category"
    };

    if(use_pick_existing_hadron_branch){
      skip.insert("reco_leading_pion_cos_wrt_beam");
      skip.insert("reco_leading_proton_cos_wrt_beam");
      skip.insert("true_leading_pion_cos_wrt_beam");
      skip.insert("true_leading_proton_cos_wrt_beam");
      skip.insert("reco_leading_pion_ke");
      skip.insert("reco_leading_proton_ke");
      skip.insert("true_leading_pion_ke");
      skip.insert("true_leading_proton_ke");
    }

    br_names.clear();
    TObjArray* blist = first_in_tree->GetListOfBranches();
    br_names.reserve(blist->GetEntries());
    for(int i=0;i<blist->GetEntries();++i){
      auto* b = (TBranch*)blist->At(i);
      std::string name = b->GetName();
      if(skip.count(name)) continue;
      br_names.push_back(name);
    }
    br_vals.assign(br_names.size(), 0.0);

    for(size_t i=0;i<br_names.size();++i){
      out_tree->Branch(br_names[i].c_str(), &br_vals[i]);
    }

    out_tree->Branch("Run", &run);
    out_tree->Branch("Subrun", &subrun);
    out_tree->Branch("Evt", &event);

    out_tree->Branch("cut_type", &cut_type, "cut_type/I");
    out_tree->Branch("is_nu", &is_nu, "is_nu/I");
    out_tree->Branch("is_data", &is_data, "is_data/I");
    out_tree->Branch("category", &category, "category/I");
    out_tree->Branch("ppfx_cv_weight", &ppfx_cv_weight, "ppfx_cv_weight/F");

    if(use_pick_existing_hadron_branch){
      out_tree->Branch("reco_leading_hadron_cos_wrt_beam", &reco_leading_hadron_cos_wrt_beam);
      out_tree->Branch("true_leading_hadron_cos_wrt_beam", &true_leading_hadron_cos_wrt_beam);
      out_tree->Branch("reco_leading_hadron_ke",           &reco_leading_hadron_ke);
      out_tree->Branch("true_leading_hadron_ke",           &true_leading_hadron_ke);
    }

    booked_nominal = true;
  }

  void book_syst_once(cfg::ConfigurationTable& config,
                      TTree* out_tree,
                      const std::string& syst_type)
  {
    auto& block = syst_blocks[syst_type];
    block.type = syst_type;
    if(block.booked) return;

    size_t n = 0;
    for(cfg::ConfigurationTable& t : config.get_subtables("sys")){
      if(t.get_string_field("type") == syst_type) n++;
    }
    block.params.reserve(n);

    for(cfg::ConfigurationTable& t : config.get_subtables("sys")){
      if(t.get_string_field("type") != syst_type) continue;

      block.params.emplace_back();
      SystParam& p = block.params.back();

      p.full = t.get_string_field("name");
      p.out  = t.has_field("name_short") ? t.get_string_field("name_short") : p.full;

      p.arr = new TClonesArray("TGraph", 1);
      out_tree->Branch(p.out.c_str(), &p.arr, 8000, 0);

      if(block.type == "multisigma" && flat_sigma_names.count(p.out) > 0){
        p.export_flat_sigma = true;
        out_tree->Branch((p.out + "_m1").c_str(), &p.flat_m1);
        out_tree->Branch((p.out + "_p1").c_str(), &p.flat_p1);
      }

      auto single_weight_it = single_weight_names.find(p.full);
      if(block.type == "multisim" && single_weight_it != single_weight_names.end()){
        p.export_single_weight = true;
        p.single_weight_name = single_weight_it->second;
        out_tree->Branch(p.single_weight_name.c_str(), &p.single_weight);
      }
    }

    block.booked = true;
  }

  void bind_nominal_input(TTree* in_tree){
    in_tree->SetBranchAddress("Run", &run);
    in_tree->SetBranchAddress("Subrun", &subrun);
    in_tree->SetBranchAddress("Evt", &event);

    _category = -5;
    if(in_tree->GetBranch("true_category")){
      in_tree->SetBranchAddress("true_category", &_category);
    }

    for(size_t i=0;i<br_names.size();++i){
      br_vals[i] = 0.0;
      if(in_tree->GetBranch(br_names[i].c_str())){
        in_tree->SetBranchAddress(br_names[i].c_str(), &br_vals[i]);
      }
    }

    reco_had_cos_src=-9999;
    true_had_cos_src=-9999;
    reco_had_ke_src=-9999;
    true_had_ke_src=-9999;
    ppfx_cv_weight=1.0;

    if(use_pick_existing_hadron_branch){
      const char* reco_cos_br = pick_existing(in_tree,
        "reco_leading_pion_cos_wrt_beam","reco_leading_proton_cos_wrt_beam");
      const char* true_cos_br = pick_existing(in_tree,
        "true_leading_pion_cos_wrt_beam","true_leading_proton_cos_wrt_beam");
      const char* reco_ke_br  = pick_existing(in_tree,
        "reco_leading_pion_ke","reco_leading_proton_ke");
      const char* true_ke_br  = pick_existing(in_tree,
        "true_leading_pion_ke","true_leading_proton_ke");

      if(reco_cos_br) in_tree->SetBranchAddress(reco_cos_br, &reco_had_cos_src);
      if(true_cos_br) in_tree->SetBranchAddress(true_cos_br, &true_had_cos_src);
      if(reco_ke_br)  in_tree->SetBranchAddress(reco_ke_br,  &reco_had_ke_src);
      if(true_ke_br)  in_tree->SetBranchAddress(true_ke_br,  &true_had_ke_src);
    }
  }

  void bind_ppfx_input(TTree* numi_tree){
    ppfx_cv_weight = 1.0;
    if(!numi_tree) return;
    if(numi_tree->GetBranch("ppfx_cv_weight"))
      numi_tree->SetBranchAddress("ppfx_cv_weight", &ppfx_cv_weight);
  }

  void bind_syst_input(TTree* syst_tree, const std::string& syst_type){
    auto it = syst_blocks.find(syst_type);
    if(it == syst_blocks.end()) return;
    auto& block = it->second;

    for(auto& p : block.params){
      p.present_in_this_tree = false;

      p.wD = nullptr; p.wF = nullptr;
      p.xD = nullptr; p.xF = nullptr;
      p.wType.clear();
      p.xType.clear();
    }

    if(!syst_tree) return;

    for(auto& p : block.params){
      TBranch* bw = syst_tree->GetBranch(p.full.c_str());
      if(!bw) continue;

      p.present_in_this_tree = true;

      std::string cls = bw->GetClassName();
      std::string wt;

      if(!cls.empty()){
        wt = cls;
      } else {
        TLeaf* lw = bw->GetLeaf(p.full.c_str());
        wt = (lw ? lw->GetTypeName() : "");
      }

      auto has = [&](const std::string& needle){
        return wt.find(needle) != std::string::npos;
      };

      if(has("vector<double>") || has("std::vector<double>")){
        p.wD = nullptr; p.wF = nullptr;
        syst_tree->SetBranchAddress(p.full.c_str(), &p.wD);
        p.wType = "double";
      }
      else if(has("vector<float>") || has("std::vector<float>")){
        p.wF = nullptr; p.wD = nullptr;
        syst_tree->SetBranchAddress(p.full.c_str(), &p.wF);
        p.wType = "float";
      }
      else{
        std::cout << "WARN: unexpected syst weight type: " << p.full
                  << " class/type=" << wt << " -> disable knob\n";
        p.present_in_this_tree = false;
        continue;
      }

      std::string ns = p.full + "_sigma";
      TBranch* bx = syst_tree->GetBranch(ns.c_str());
      if(!bx){ p.xType=""; continue; }

      std::string clsx = bx->GetClassName();
      std::string xt;
      if(!clsx.empty()) xt = clsx;
      else {
        TLeaf* lx = bx->GetLeaf(ns.c_str());
        xt = (lx ? lx->GetTypeName() : "");
      }

      if(xt.find("vector<double>")!=std::string::npos || xt.find("std::vector<double>")!=std::string::npos){
        p.xD=nullptr; p.xF=nullptr;
        syst_tree->SetBranchAddress(ns.c_str(), &p.xD);
        p.xType="double";
      }
      else if(xt.find("vector<float>")!=std::string::npos || xt.find("std::vector<float>")!=std::string::npos){
        p.xF=nullptr; p.xD=nullptr;
        syst_tree->SetBranchAddress(ns.c_str(), &p.xF);
        p.xType="float";
      }
      else{
        std::cout << "WARN: unexpected syst sigma type: " << ns
                  << " class/type=" << xt << " -> ignore sigma\n";
        p.xType="";
      }
    }
  }

  void fill_one_syst_event(const std::string& syst_type, bool table_is_nu){
    auto it = syst_blocks.find(syst_type);
    if(it == syst_blocks.end()) return;
    auto& block = it->second;

    auto copy_vecD = [&](const std::vector<Double_t>* v, std::vector<Double_t>& out){
      if(!v) return;
      out.reserve(out.size() + v->size());
      for(auto val : *v) out.push_back(val);
    };
    auto copy_vecF = [&](const std::vector<Float_t>* v, std::vector<Double_t>& out){
      if(!v) return;
      out.reserve(out.size() + v->size());
      for(auto val : *v) out.push_back((Double_t)val);
    };

    for(auto& p : block.params){
      std::vector<Double_t> nsigmas;
      std::vector<Double_t> weights;

      p.flat_m1 = 1.0;
      p.flat_p1 = 1.0;
      p.single_weight = 1.0;

      if(p.arr == nullptr){
        std::cout << "WARN: p.arr is null for " << p.full << "\n";
        continue;
      }

      bool is_morph = (syst_type=="multisigma") &&
        (std::find(morph_names.begin(), morph_names.end(), p.full) != morph_names.end());

      if(!p.present_in_this_tree){
        if(syst_type=="multisim"){
          nsigmas = {0.0};
          weights = {1.0};
        } else {
          nsigmas = {-1.0, 0.0, 1.0};
          weights = { 1.0, 1.0, 1.0};
        }
      }
      else if(is_morph){
        Double_t w1 = 1.0;
        if(p.wType=="double" && p.wD && !p.wD->empty()) w1 = p.wD->front();
        if(p.wType=="float"  && p.wF && !p.wF->empty()) w1 = (Double_t)p.wF->front();
        if(!std::isfinite(w1)) w1 = 1.0;

        if(!table_is_nu) w1 = 1.0;

        Double_t w05 = 1.0 + 0.5*(w1 - 1.0);

        nsigmas = {-1.0, -0.5, 0.0, 0.5, 1.0};
        weights = { w1,   w05,  1.0, w05,  w1};
      }
      else{
        if(p.xType=="double") copy_vecD(p.xD, nsigmas);
        else if(p.xType=="float") copy_vecF(p.xF, nsigmas);
        else nsigmas = {-1.0, 0.0, 1.0};

        if(p.wType=="double") copy_vecD(p.wD, weights);
        else if(p.wType=="float") copy_vecF(p.wF, weights);
        else weights.assign(nsigmas.size(), 1.0);

        if(!table_is_nu){
          for(auto& wv : weights) wv = -5.0;
        }

        if(syst_type=="multisigma"){
            bool has_nominal = std::any_of(nsigmas.begin(),nsigmas.end(),[](double x){ return std::abs(x) < 1e-12; }); //Check there is 0 sigma or not

            if(!has_nominal){
                nsigmas.push_back(0.0);
                weights.push_back(table_is_nu ? 1.0 : -5.0);
            }
        }

        // For multisim: weights vector is the universe weights (no sigma axis).
        // Store as TGraph with x=universe_index, y=weight so GUNDAM can use
        // ByUniverseReweight. Override the fallback nsigmas built above.
        if(syst_type=="multisim" && p.present_in_this_tree){
          nsigmas.clear();
          // weights already filled from wD/wF above
          for(std::size_t iu=0; iu<weights.size(); ++iu)
            nsigmas.push_back((Double_t)iu);
        }

        if(weights.size() != nsigmas.size()){
          std::cout << "WARN: size mismatch for knob " << p.full
                    << " nsigmas=" << nsigmas.size()
                    << " weights=" << weights.size()
                    << " -> force neutral\n";
          nsigmas = {-1.0, 0.0, 1.0};
          weights = { 1.0, 1.0, 1.0};
        }
      }

      TGraph g((int)nsigmas.size(), nsigmas.data(), weights.data());
      g.Sort();

      if(p.export_flat_sigma){
        p.flat_m1 = g.Eval(-1.0);
        p.flat_p1 = g.Eval( 1.0);
        if(!std::isfinite(p.flat_m1)) p.flat_m1 = 1.0;
        if(!std::isfinite(p.flat_p1)) p.flat_p1 = 1.0;
      }

      if(p.export_single_weight && p.present_in_this_tree && !weights.empty()){
        p.single_weight = weights.front();
        if(!std::isfinite(p.single_weight)) p.single_weight = 1.0;
      }

      p.arr->Clear("C");
      new((*p.arr)[0]) TGraph(g.GetN(), g.GetX(), g.GetY());
    }
  }
};

static void merge_prefix_under_directory(
    TFile* input,
    TFile* output,
    cfg::ConfigurationTable& config,
    cfg::ConfigurationTable& table,
    const std::string& in_dir_path,
    const std::string& out_dir_path,
    const std::string& out_tree_name,
    const std::string& prefix
){
  auto* indir = get_dir_path(input, in_dir_path);
  if(!indir){
    std::cout << "WARN: missing input directory " << in_dir_path << std::endl;
    return;
  }

  TDirectory* outdir = ensure_dir_path(output, out_dir_path);
  outdir->cd();
  outdir->Delete((out_tree_name + ";*").c_str());

//  TDirectory* parent = (TDirectory*)input;
//  parent = get_parent_directory(parent, in_dir_path.c_str());
//  if(parent){
//    TH1D* pot = (TH1D*) parent->Get("POT");
//    TH1D* livetime = (TH1D*) parent->Get("Livetime");
//    if(pot)      outdir->WriteObject(pot, "POT");
//    if(livetime) outdir->WriteObject(livetime, "Livetime");
//  }
  TH1D* pot = (TH1D*) indir->Get("POT");
  TH1D* livetime = (TH1D*) indir->Get("Livetime");
  if(pot)      outdir->WriteObject(pot, "POT");
  if(livetime) outdir->WriteObject(livetime, "Livetime");

  TTree* out_tree = new TTree(out_tree_name.c_str(), out_tree_name.c_str());
  out_tree->SetAutoFlush(500);       // 500 entries마다 basket flush -> 메모리 해제
  out_tree->SetAutoSave(50000000);   // 50MB마다 디스크에 save

  MergeWriter w;

  if(table.has_field("use_pick_existing_hadron_branch")){
    w.use_pick_existing_hadron_branch =
      table.get_bool_field("use_pick_existing_hadron_branch");
  }
  else if(config.has_field("use_pick_existing_hadron_branch")){
    w.use_pick_existing_hadron_branch =
      config.get_bool_field("use_pick_existing_hadron_branch");
  }

  bool want_syst = table.has_field("gundam_store_syst") ? table.get_bool_field("gundam_store_syst") : false;
  if(want_syst){
    w.book_syst_once(config, out_tree, "multisigma");
    w.book_syst_once(config, out_tree, "variation");
    w.book_syst_once(config, out_tree, "NuMIfluxsim");
  }
  w.book_syst_once(config, out_tree, "multisim");

  bool booked = false;

  const bool table_is_nu   = table.get_bool_field("is_nu");
  const bool table_is_data = table.get_bool_field("is_data");
  w.is_nu = table_is_nu ? 1 : 0;

  TIter nextkey(indir->GetListOfKeys());
  TKey* key;
  while((key=(TKey*)nextkey())){
    TObject* obj = key->ReadObj();
    if(!obj || !obj->InheritsFrom(TTree::Class())) continue;

    TTree* in_tree = (TTree*)obj;
    std::string tname = in_tree->GetName();

    if(!starts_with(tname, prefix)) continue;
    if(is_syst_tree_name(tname)) continue;
    if(!is_allowed_nominal_tree(tname, prefix)) {
      std::cout << "SKIP nominal tree (not allowed): " << tname << "\n";
      continue;
    }

    int forced_cut = cut_type_from_tree_name_safe(tname);

    TTree* ms_tree   = want_syst ? find_syst_tree(indir, tname, kMultisigmaSuffixes) : nullptr;
    TTree* var_tree  = want_syst ? find_syst_tree(indir, tname, kVariationSuffixes)  : nullptr;
    TTree* numi_tree = find_syst_tree(indir, tname, kFluxSuffixes);
    TTree* msim_tree = find_syst_tree(indir, tname, kMultisimSuffixes);

    std::cout << "MERGE " << in_dir_path << "/" << tname
              << " -> " << out_dir_path << "/" << out_tree_name
              << " (cut_type=" << forced_cut << ")"
              << "  ms=" << ms_tree
              << " var=" << var_tree
              << " numi=" << numi_tree
              << " msim=" << msim_tree
              << "  use_pick_existing_hadron_branch=" << w.use_pick_existing_hadron_branch
              << std::endl;

    if(!booked){
      w.book_nominal_once(out_tree, in_tree);
      booked = true;
    }

    w.bind_nominal_input(in_tree);
    if(want_syst){
      w.bind_syst_input(ms_tree,   "multisigma");
      w.bind_syst_input(var_tree,  "variation");
      w.bind_syst_input(numi_tree, "NuMIfluxsim");
    }
    w.bind_ppfx_input(numi_tree);
    w.bind_syst_input(msim_tree, "multisim");

    // msim_tree: vector<double>(1000) branch -> ROOT read cache가 수십GB를 잡아먹음
    // cache 완전히 끄기
    if(msim_tree){
      msim_tree->SetCacheSize(0);
      msim_tree->SetMaxVirtualSize(0);
    }

    Long64_t n_nom = in_tree->GetEntries();
    Long64_t n = n_nom;
    if(ms_tree)   n = std::min(n, ms_tree->GetEntries());
    if(var_tree)  n = std::min(n, var_tree->GetEntries());
    if(numi_tree) n = std::min(n, numi_tree->GetEntries());
    if(msim_tree) n = std::min(n, msim_tree->GetEntries());

    if(want_syst){
      if(ms_tree && ms_tree->GetEntries() != n_nom)
        std::cout << "WARN: entries mismatch: " << in_dir_path << "/" << tname
                  << " nominal=" << n_nom
                  << " multisigma=" << ms_tree->GetEntries()
                  << " -> using min=" << n << "\n";
      if(var_tree && var_tree->GetEntries() != n_nom)
        std::cout << "WARN: entries mismatch: " << in_dir_path << "/" << tname
                  << " nominal=" << n_nom
                  << " variation=" << var_tree->GetEntries()
                  << " -> using min=" << n << "\n";
    }
    if(numi_tree && numi_tree->GetEntries() != n_nom)
      std::cout << "WARN: entries mismatch: " << in_dir_path << "/" << tname
                << " nominal=" << n_nom
                << " numi=" << numi_tree->GetEntries()
                << " -> using min=" << n << "\n";
    if(msim_tree && msim_tree->GetEntries() != n_nom)
      std::cout << "WARN: entries mismatch: " << in_dir_path << "/" << tname
                << " nominal=" << n_nom
                << " multisim=" << msim_tree->GetEntries()
                << " -> using min=" << n << "\n";

    for(Long64_t i=0;i<n;++i){
      in_tree->GetEntry(i);
      if(ms_tree)   ms_tree->GetEntry(i);
      if(var_tree)  var_tree->GetEntry(i);
      if(numi_tree) numi_tree->GetEntry(i);
      if(msim_tree) msim_tree->GetEntry(i);

      w.cut_type = forced_cut;

      if(table_is_data){
        w.is_data  = 1;
        w.category = 10;
      }else{
        w.is_data  = 0;
        w.category = (int)w._category;
      }

      if(w.use_pick_existing_hadron_branch){
        w.reco_leading_hadron_cos_wrt_beam = w.reco_had_cos_src;
        w.true_leading_hadron_cos_wrt_beam = w.true_had_cos_src;
        w.reco_leading_hadron_ke           = w.reco_had_ke_src;
        w.true_leading_hadron_ke           = w.true_had_ke_src;
      }

      if(want_syst){
        w.fill_one_syst_event("multisigma", table_is_nu);
        w.fill_one_syst_event("variation",  table_is_nu);
        w.fill_one_syst_event("NuMIfluxsim", table_is_nu);
      }
      w.fill_one_syst_event("multisim", table_is_nu);

      if(std::isnan(w.ppfx_cv_weight) || w.ppfx_cv_weight <= 0)
          w.ppfx_cv_weight = 1.0;


      out_tree->Fill();

      // 5000 entries마다 강제 flush로 메모리 상한 방지
      if(i > 0 && i % 5000 == 0){
        out_tree->FlushBaskets();
        gDirectory->SaveSelf();
      }
    }
  }

  outdir->cd();
  out_tree->Write("", TObject::kOverwrite);
}

//static void merge_multisim_trees_under_directory(
//    TFile* input,
//    TFile* output,
//    const std::string& in_dir_path,
//    const std::string& out_dir_path,
//    const std::string& out_tree_name,
//    const std::string& prefix
//){
//  auto* indir = get_dir_path(input, in_dir_path);
//  if(!indir){
//    std::cout << "WARN: missing input directory " << in_dir_path << std::endl;
//    return;
//  }
//
//  TDirectory* outdir = ensure_dir_path(output, out_dir_path);
//  outdir->cd();
//  outdir->Delete((out_tree_name + ";*").c_str());
//
//  TChain chain("dummy");
//  bool added = false;
//
//  TIter nextkey(indir->GetListOfKeys());
//  TKey* key;
//  while((key=(TKey*)nextkey())){
//    TObject* obj = key->ReadObj();
//    if(!obj || !obj->InheritsFrom(TTree::Class())) continue;
//
//    std::string tname = obj->GetName();
//
//    if(!starts_with(tname, prefix)) continue;
//    if(is_syst_tree_name(tname)) continue;
//    if(!is_allowed_nominal_tree(tname, prefix)) {
//      std::cout << "SKIP multisim source tree (not allowed): " << tname << "\n";
//      continue;
//    }
//
//    TTree* multisim_tree = find_syst_tree(indir, tname, kMultisimSuffixes);
//    if(!multisim_tree) continue;
//
//    std::string fullpath = input->GetName();
//    fullpath += "/";
//    fullpath += in_dir_path;
//    fullpath += "/";
//    fullpath += multisim_tree->GetName();
//
//    chain.Add(fullpath.c_str());
//    added = true;
//
//    std::cout << "ADD multisim tree: " << in_dir_path << "/" << multisim_tree->GetName()
//              << " -> " << out_dir_path << "/" << out_tree_name << "\n";
//  }
//
//  if(!added){
//    std::cout << "WARN: no multisim trees found for " << in_dir_path
//              << " prefix=" << prefix << std::endl;
//    return;
//  }
//
//  outdir->cd();
//  TTree* merged = chain.CloneTree(-1, "fast");
//  if(!merged){
//    std::cout << "WARN: failed to merge multisim trees into " << out_tree_name << std::endl;
//    return;
//  }
//
//  merged->SetName(out_tree_name.c_str());
//  merged->SetTitle(out_tree_name.c_str());
//  merged->Write("", TObject::kOverwrite);
//
//  std::cout << "WROTE merged multisim tree: "
//            << out_dir_path << "/" << out_tree_name
//            << " entries=" << merged->GetEntries() << std::endl;
//}

int main(int argc, char* argv[])
{
  cfg::ConfigurationTable config;
  config.set_config(argv[1]);
  std::vector<cfg::ConfigurationTable> tables = config.get_subtables("tree");

  std::string input_filename;
  if(argc > 2) input_filename = argv[2];
  else         input_filename = config.get_string_field("input.path");
  TFile* input = TFile::Open(input_filename.c_str(), "READ");

  std::string output_filename;
  if(config.has_field("output.path")) output_filename = config.get_string_field("output.path");
  else                                output_filename = "output.root";
  TFile* output = TFile::Open(output_filename.c_str(), "RECREATE");

  bool did_nominal=false;
  bool did_onbeam=false;

  for(cfg::ConfigurationTable& table : tables){
    std::string origin = table.get_string_field("origin");
    std::string dest   = table.get_string_field("destination");
    while(!dest.empty() && dest.back()=='/') dest.pop_back();

    if(!did_nominal && origin.find("events/nominal/") == 0){
      merge_prefix_under_directory(input, output, config, table,
                                   "events/nominal",
                                   "events/nominal",
                                   "selected_merged",
                                   "selected");

      merge_prefix_under_directory(input, output, config, table,
                                   "events/nominal",
                                   "events/nominal",
                                   "signal_merged",
                                   "signal");

//      merge_multisim_trees_under_directory(input, output,
//              "events/nominal",
//              "events/nominal",
//              "selected_merged_multisimTree",
//              "selected");
//
//      merge_multisim_trees_under_directory(input, output,
//              "events/nominal",
//              "events/nominal",
//              "signal_merged_multisimTree",
//              "signal");
//
      did_nominal = true;
    }

    if(!did_onbeam && origin.find("events/onbeam/") == 0){
      merge_prefix_under_directory(input, output, config, table,
                                   "events/onbeam",
                                   "events/onbeam",
                                   "selected_merged",
                                   "selected");
      did_onbeam = true;
    }
  }

  output->Close();
  input->Close();
  return 0;
}


