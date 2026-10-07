#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLatex.h"
#include "Math/Vector4D.h"
#include "TStyle.h"

using namespace ROOT;
using namespace ROOT::VecOps;
using RNode = ROOT::RDF::RNode;

using Vec_t = const ROOT::RVec<float>&;
using Vec_i = const ROOT::RVec<int>&;

using cRVecF = const ROOT::RVecF &;
using cRVecI = const ROOT::RVecI &;
using cRVecC = const ROOT::RVecC &;
using cRVecU = const ROOT::RVecU &;

float deltaR(float eta_1, float eta_2, float phi_1, float phi_2){
   const float deta = eta_1 - eta_2;
   const float dphi = ROOT::Math::VectorUtil::Phi_mpi_pi(phi_1 - phi_2);
   const float dRsq = std::pow(deta,2) + std::pow(dphi,2);

   return sqrt(dRsq);
}

// Monitoring triggers used in Early Run 3
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauTrigger?extralog=-%20caching%20topic

bool PassMuTauTrig2022(UInt_t ntrig, Vec_i trig_id, Vec_i trig_bits, Vec_t trig_pt,
                       Vec_t trig_eta, Vec_t trig_phi, float tau_pt, float tau_eta,
                       float tau_phi)
{
  if (tau_pt <= 0)
    return false;
  for(int it=0; it < ntrig; it++){
    const ROOT::Math::PtEtaPhiMVector trig(trig_pt[it],trig_eta[it],trig_phi[it],0);
    float dR = deltaR(trig.Eta(),tau_eta,trig.Phi(),tau_phi);
    if (dR < 0.5){ //dR < 0.5
      if ( (trig_bits[it] & (1<<8)) != 0 && (trig_bits[it] & (1<<27)) != 0 && trig_id[it] == 15 ){ 
          return true;
      }
    }
  }
  return false;
}

bool PassEleTauTrig2022(UInt_t ntrig, Vec_t trig_l1pt, Vec_i trig_l1iso, Vec_i trig_id,
                        Vec_i trig_bits, Vec_t trig_pt, Vec_t trig_eta, Vec_t trig_phi,
                        float tau_pt, float tau_eta, float tau_phi)
{
  if (tau_pt <= 0)
    return false;
  for(int it=0; it < ntrig; it++){
    const ROOT::Math::PtEtaPhiMVector trig(trig_pt[it],trig_eta[it],trig_phi[it],0);
    float dR = deltaR(trig.Eta(),tau_eta,trig.Phi(),tau_phi);
    if (dR < 0.5){ //dR < 0.5
      if ( (trig_bits[it] & (1<<8)) != 0 && (trig_bits[it] & (1<<27)) != 0 && trig_id[it] == 15 && trig_pt[it] > 30 ){
        if ( trig_l1iso[it] > 0 && trig_l1pt[it] > 26 ){ // not sure why "|| (trig_l1pt[it] >= 70)" is missing
          return true;
        }
      }
    }
  }
  return false;
}

bool PassDiTauTrig2022(UInt_t ntrig, Vec_t trig_l1pt, Vec_i trig_l1iso, Vec_i trig_id,
                       Vec_i trig_bits, Vec_t trig_pt, Vec_t trig_eta, Vec_t trig_phi,
                       float tau_pt, float tau_eta, float tau_phi)
{
  if (tau_pt <= 0)
    return false;
  for(int it=0; it < ntrig; it++){
    const ROOT::Math::PtEtaPhiMVector trig(trig_pt[it],trig_eta[it],trig_phi[it],0);
    float dR = deltaR(trig.Eta(),tau_eta,trig.Phi(),tau_phi);
    if (dR < 0.5){ //dR < 0.5, 1 => Medium, 17 => Monitoring, 18 => MonitoringForVBFIsoTau, bit1 && bit17 && !bit18
      // drop !bit18 cut
      if ( (trig_bits[it] & (1<<1)) != 0 && (trig_bits[it] & (1<<18)) != 0 && trig_id[it] == 15 && trig_pt[it] > 35 ){ 
        if ( (trig_l1iso[it] > 0 && trig_l1pt[it] >= 32) || (trig_l1pt[it] > 70) ) {
          return true;
        }
      }
    }
  }
  return false;
}

bool PassDiTauJetTrig2022(UInt_t ntrig,Vec_t trig_l1pt, Vec_i trig_l1iso, Vec_i trig_id,
                          Vec_i trig_bits, Vec_t trig_pt, Vec_t trig_eta, Vec_t trig_phi,
                          float tau_pt, float tau_eta, float tau_phi)
{
  if (tau_pt <= 0)
    return false;
  for(int it=0; it < ntrig; it++){
    const ROOT::Math::PtEtaPhiMVector trig(trig_pt[it],trig_eta[it],trig_phi[it],0);
    float dR = deltaR(trig.Eta(),tau_eta,trig.Phi(),tau_phi);
    if (dR < 0.5){ //dR < 0.5, 1 => Medium, 17 => Monitoring, 18 => MonitoringForVBFIsoTau, bit1 && bit17 && !bit18
      // drop !bit18 cut
      if((trig_bits[it] & (1<<1)) != 0 && (trig_bits[it] & (1<<18)) != 0 && trig_id[it] == 15 && trig_pt[it] > 30 ){ 
        if ( (trig_l1iso[it] > 0 && trig_l1pt[it] > 26) )  
          return true;
      }
    }
  }
  return false;
}

// Monitoring triggers used in Late Run 3
// Important to note that in 2022 the lowest pT L1 seed for the ditau trigger was L1_DoubleIsoTau32er2p1,
// while in 2024 it is L1_DoubleIsoTau34er2p1.

bool PassMuTauTrig2024(UInt_t ntrig, Vec_i trig_id, Vec_i trig_bits, Vec_t trig_pt,
                       Vec_t trig_eta, Vec_t trig_phi, float tau_pt, float tau_eta,
                       float tau_phi)
{
  if (tau_pt <= 0)
    return false;
  for(int it=0; it < ntrig; it++){
    const ROOT::Math::PtEtaPhiMVector trig(trig_pt[it],trig_eta[it],trig_phi[it],0);
    float dR = deltaR(trig.Eta(),tau_eta,trig.Phi(),tau_phi);
    if (dR < 0.5){ // HLT_IsoMu20_eta2p1_PNetTauhPFJet27_Medium_eta2p3_CrossL1 trigger bits are 1, 4, 13
      if ( (trig_bits[it] & (1<<1)) != 0 && (trig_bits[it] & (1<<4)) != 0 && (trig_bits[it] & (1<<13)) != 0 && trig_id[it] == 15 ){ 
          return true;
      }
    }
  }
  return false;
}

// Same monitoring path as mutau, just as was done in 22
bool PassEleTauTrig2024(UInt_t ntrig, Vec_t trig_l1pt, Vec_i trig_l1iso, Vec_i trig_id,
                        Vec_i trig_bits, Vec_t trig_pt, Vec_t trig_eta, Vec_t trig_phi,
                        float tau_pt, float tau_eta, float tau_phi)
{
  if (tau_pt <= 0)
    return false;
  for(int it=0; it < ntrig; it++){
    const ROOT::Math::PtEtaPhiMVector trig(trig_pt[it],trig_eta[it],trig_phi[it],0);
    float dR = deltaR(trig.Eta(),tau_eta,trig.Phi(),tau_phi);
    if (dR < 0.5){ // HLT_IsoMu20_eta2p1_PNetTauhPFJet27_Medium_eta2p3_CrossL1 trigger bits are 1, 4, 13
      if ( (trig_bits[it] & (1<<1)) != 0 && (trig_bits[it] & (1<<4)) != 0 && (trig_bits[it] & (1<<13)) != 0 && trig_id[it] == 15 && trig_pt[it] >= 30 ){
        if ( (trig_l1iso[it] > 0 && trig_l1pt[it] >= 26) || (trig_l1pt[it] >= 70) ) { 
          return true;
        }
      }
    }
  }
  return false;
}

bool PassDiTauTrig2024(UInt_t ntrig, Vec_t trig_l1pt, Vec_i trig_l1iso, Vec_i trig_id,
                       Vec_i trig_bits, Vec_t trig_pt, Vec_t trig_eta, Vec_t trig_phi,
                       float tau_pt, float tau_eta, float tau_phi)
{
  if (tau_pt <= 0)
    return false;
  for(int it=0; it < ntrig; it++){
    const ROOT::Math::PtEtaPhiMVector trig(trig_pt[it],trig_eta[it],trig_phi[it],0);
    float dR = deltaR(trig.Eta(),tau_eta,trig.Phi(),tau_phi);
    if (dR < 0.5){ // HLT_IsoMu24_eta2p1_PNetTauhPFJet30_Medium_L2NN_eta2p3_CrossL1 trigger bits are 1, 4, 23
      if((trig_bits[it] & (1<<1)) != 0 && (trig_bits[it] & (1<<4)) != 0 && (trig_bits[it] & (1<<23)) != 0 && trig_id[it] == 15 && trig_pt[it] >= 30){ 
        if ( (trig_l1iso[it] > 0 && trig_l1pt[it] >= 34) || (trig_l1pt[it] >= 70) ) {
          return true;
        }
      }
    }
  }
  return false;
}

bool PassDiTauJetTrig2024(UInt_t ntrig, Vec_t trig_l1pt, Vec_i trig_l1iso, Vec_i trig_id,
                          Vec_i trig_bits, Vec_t trig_pt, Vec_t trig_eta, Vec_t trig_phi,
                          float tau_pt, float tau_eta, float tau_phi)
{
  if (tau_pt <= 0)
    return false;
  for(int it=0; it < ntrig; it++){
    const ROOT::Math::PtEtaPhiMVector trig(trig_pt[it],trig_eta[it],trig_phi[it],0);
    float dR = deltaR(trig.Eta(),tau_eta,trig.Phi(),tau_phi);
    if (dR < 0.5){ // HLT_IsoMu24_eta2p1_PNetTauhPFJet26_L2NN_eta2p3_CrossL1 trigger bits are 4, 20 
      if((trig_bits[it] & (1<<4)) != 0 && (trig_bits[it] & (1<<20)) != 0 && trig_id[it] == 15 && trig_pt[it] >= 26){ 
        if ( (trig_l1iso[it] > 0 && trig_l1pt[it] >= 26) )  
          return true;
      }
    }
  }
  return false;
}

