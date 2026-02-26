#include "ConfigManager.hpp"
#include "CutFlowManager.hpp"
#include "EventReader.hpp"
#include "ExtensionsHelpers.hpp"
#include "HistogramsFiller.hpp"
#include "HistogramsHandler.hpp"
#include "LbLHistogramsFiller.hpp"
#include "LbLObjectsManager.hpp"
#include "Logger.hpp"

using namespace std;

tuple<float, float> findRefEtandPhi(const vector<pair<float, float>> &phiAndEt) {
  float maxTowerEt = 0;
  float refTowerPhi = 0;
  for (auto &[phi, et] : phiAndEt) {
    if (et > maxTowerEt) {
      maxTowerEt = et;
      refTowerPhi = phi;
    }
  }
  return {maxTowerEt, refTowerPhi};
}

tuple<array<float, 5>, array<float, 5>> getVertices(float circleRadius, float trapezoidInnerWidth, float trapezoidOuterWidth,
                                                    float trapezoidHeight, float angle_rad) {
  array<float, 5> verticesX, verticesY;

  double baseX = circleRadius * TMath::Cos(angle_rad);
  double baseY = circleRadius * TMath::Sin(angle_rad);

  verticesX[0] = baseX - (trapezoidInnerWidth / 2.0) * TMath::Sin(angle_rad);
  verticesY[0] = baseY + (trapezoidInnerWidth / 2.0) * TMath::Cos(angle_rad);
  verticesX[1] = baseX + (trapezoidInnerWidth / 2.0) * TMath::Sin(angle_rad);
  verticesY[1] = baseY - (trapezoidInnerWidth / 2.0) * TMath::Cos(angle_rad);

  // Outer vertices (wider end)
  verticesX[2] =
      verticesX[1] + trapezoidHeight * TMath::Cos(angle_rad) + (trapezoidOuterWidth - trapezoidInnerWidth) / 2.0 * TMath::Sin(angle_rad);
  verticesY[2] =
      verticesY[1] + trapezoidHeight * TMath::Sin(angle_rad) - (trapezoidOuterWidth - trapezoidInnerWidth) / 2.0 * TMath::Cos(angle_rad);
  verticesX[3] =
      verticesX[0] + trapezoidHeight * TMath::Cos(angle_rad) - (trapezoidOuterWidth - trapezoidInnerWidth) / 2.0 * TMath::Sin(angle_rad);
  verticesY[3] =
      verticesY[0] + trapezoidHeight * TMath::Sin(angle_rad) + (trapezoidOuterWidth - trapezoidInnerWidth) / 2.0 * TMath::Cos(angle_rad);

  // Closing the trapezoid
  verticesX[4] = verticesX[0];
  verticesY[4] = verticesY[0];

  return {verticesX, verticesY};
}

vector<unique_ptr<TGraph>> getTrapezoids(const vector<pair<float, float>> &phiAndEt, const map<string, float> &visualizationParams, float maxTowerEt, float refTowerPhi) {
  
  float circleRadius = visualizationParams.at("circleRadius");
  float trapezoidInnerWidth = visualizationParams.at("towerInnerWidth");
  float maxTowerHeight = visualizationParams.at("maxTowerHeight");
  
  vector<unique_ptr<TGraph>> trapezoids;
  for (auto &[phi, et] : phiAndEt) {
    double trapezoidHeight = et;

    if (visualizationParams.at("normalizeMomenta")) {
      trapezoidHeight /= maxTowerEt;
      trapezoidHeight *= maxTowerHeight;
    }

    double trapezoidOuterWidth = trapezoidInnerWidth;
    trapezoidOuterWidth += visualizationParams.at("towerOuterWidthScale") * trapezoidHeight;

    double angle_rad = phi;

    if (visualizationParams.at("normalizeRotation")) {
      angle_rad -= refTowerPhi;
    }

    auto [verticesX, verticesY] = getVertices(circleRadius, trapezoidInnerWidth, trapezoidOuterWidth, trapezoidHeight, angle_rad);

    trapezoids.push_back(make_unique<TGraph>(5, verticesX.data(), verticesY.data()));
  }
  return trapezoids;
}

int getColor(const map<string, float> &visualizationParams, bool track = false) {
  string object = track ? "track" : "tower";

  if (visualizationParams.at(object+"Color") >= 0) return visualizationParams.at(object+"Color");
  int colorIndex = 10;
  auto color = gROOT->GetColor(colorIndex);
  color->SetRGB(visualizationParams.at(object+"ColorR"), visualizationParams.at(object+"ColorG"), visualizationParams.at(object+"ColorB"));
  return colorIndex;
}

void addShapesToCanvas(const vector<pair<float, float>> &towersPhiAndEt, TCanvas *canvas, const map<string, float> &visualizationParams) {
  canvas->cd();

  
  auto [maxTowerEt, refTowerPhi] = findRefEtandPhi(towersPhiAndEt);
  auto towersTrapezoids = getTrapezoids(towersPhiAndEt, visualizationParams, maxTowerEt, refTowerPhi);

  int towerColor = getColor(visualizationParams);
  for(auto &trapezoid : towersTrapezoids) {
    trapezoid->SetFillColorAlpha(towerColor, visualizationParams.at("towerAlpha"));
    trapezoid->SetFillStyle(visualizationParams.at("towerFillStyle"));
    trapezoid->DrawClone("f");
  }

  canvas->Update();
}

void addTracksToCanvas(const vector<tuple<float, float, int>> &tracksPhiPtCharge, TCanvas *canvas, const map<string, float> &visualizationParams) {
  canvas->cd();

  for(auto &[phi, pt, charge] : tracksPhiPtCharge){
    float circleRadius = visualizationParams.at("circleRadius");
    // float trackRadius = pt / 3.8;

    // auto track = new TArc(x, y, trackRadius, phiStart, phiEnd);
    float x = charge * circleRadius/sqrt(1+pow(TMath::Tan(phi),2));
    float y = x * TMath::Tan(phi);

    auto track = new TLine(0, 0, x, y);
    track->SetLineColor(getColor(visualizationParams, true));
    track->SetLineWidth(visualizationParams.at("trackWidth"));
    track->DrawClone("same");
  }

  canvas->Update();
}

void drawShapes(const vector<pair<float, float>> &towersPhiAndEt, string outputFileName, const map<string, float> &visualizationParams) {
  auto canvas = new TCanvas(outputFileName.c_str(), outputFileName.c_str(), 800, 800);
  canvas->SetFillColor(visualizationParams.at("backgroundColor"));
  float canvasSize = visualizationParams.at("canvasSize");
  canvas->Range(-canvasSize, -canvasSize, canvasSize, canvasSize);
  double circleRadius = visualizationParams.at("circleRadius");
  TEllipse *circle = new TEllipse(0, 0, circleRadius);
  circle->SetFillColor(visualizationParams.at("backgroundColor"));
  circle->SetFillStyle(1001);
  circle->SetLineColor(visualizationParams.at("circleColor"));
  circle->Draw();

  addShapesToCanvas(towersPhiAndEt, canvas, visualizationParams);

  canvas->Update();
  string outputDir = outputFileName.substr(0, outputFileName.find_last_of("/"));
  gSystem->mkdir(outputDir.c_str());
  canvas->SaveAs(outputFileName.c_str());
}

float getDiphotonAcoplanarity(const shared_ptr<PhysicsObjects> &photons) {
  if (photons->size() < 2) return 1.0;
  auto photon1 = photons->at(0);
  auto photon2 = photons->at(1);
  TLorentzVector photon1vec, photon2vec;
  photon1vec.SetPtEtaPhiM(photon1->Get("et"), photon1->Get("eta"), photon1->Get("phi"), 0);
  photon2vec.SetPtEtaPhiM(photon2->Get("et"), photon2->Get("eta"), photon2->Get("phi"), 0);

  auto diphoton = photon1vec + photon2vec;
  
  double deltaPhi = photon1vec.DeltaPhi(photon2vec);
  double acoplanarity = 1 - (fabs(deltaPhi) / TMath::Pi());

  return acoplanarity;
}

void CheckArgs(int argc, char **argv) {
  if (argc != 2 && argc != 6) {
    fatal() << "Usage: " << argv[0] << " config_path" << endl;
    fatal() << "or" << endl;
    fatal() << argv[0] << " config_path input_path output_path nEvents eventsOffset" << endl;
    exit(1);
  }
}

int main(int argc, char **argv) {
  gROOT->SetBatch(true);

  vector<string> requiredArgs = {"config"};
  vector<string> optionalArgs = {};
  auto args = make_unique<ArgsManager>(argc, argv, requiredArgs, optionalArgs);
  ConfigManager::Initialize(args);

  auto &config = ConfigManager::GetInstance();

  int eventsOffset, nEvents;
  config.GetValue("eventsOffset", eventsOffset);
  config.GetValue("nEvents", nEvents);
  string outputPath = "";

  
   if (argc == 6) {
    config.SetInputPath(argv[2]);
    outputPath = argv[3];
    nEvents = stoi(argv[4]);
    eventsOffset = stoi(argv[5]);
  }

  if(outputPath == ""){
    outputPath = "../plots/visualizations/event_display_"+to_string(eventsOffset)+".pdf";
  }

  map<string, float> visualizationParams;
  config.GetMap("visualizationParams", visualizationParams);

  info() << "Creating objects" << endl;
  auto eventReader = make_shared<EventReader>();
  auto lblObjectsManager = make_unique<LbLObjectsManager>();

  auto overlapCanvas = new TCanvas("overlapCanvas", "overlapCanvas", 2000, 2000);
  overlapCanvas->cd();
  overlapCanvas->SetFillColor(visualizationParams.at("backgroundColor"));
  float canvasSize = visualizationParams.at("canvasSize");
  overlapCanvas->Range(-canvasSize, -canvasSize, canvasSize, canvasSize);
  double circleRadius = visualizationParams.at("circleRadius");
  TEllipse *circle = new TEllipse(0, 0, circleRadius);
  circle->SetFillColor(visualizationParams.at("backgroundColor"));
  circle->SetFillStyle(1001);
  circle->SetLineColor(visualizationParams.at("circleColor"));
  circle->Draw();

  

  info() << "Starting event loop" << endl;
  for (int iEvent = eventsOffset; iEvent < nEvents+eventsOffset; iEvent++) {
    if(iEvent > eventReader->GetNevents()) break;
    auto event = eventReader->GetEvent(iEvent);

    lblObjectsManager->InsertGoodPhotonsCollection(event);
    lblObjectsManager->InsertGoodTracksCollection(event);
    
    
    auto photons = event->GetCollection("goodPhoton");
    float acoplanarity = getDiphotonAcoplanarity(photons);
    // if (acoplanarity > 0.01) continue;

    // vector<pair<float, float>> photonsPhiAndEt;

    // for (auto physicsObject : *photons) {
    //   auto photon = asPhoton(physicsObject);
    //   photonsPhiAndEt.push_back({(float)photon->Get("phi"), (float)photon->Get("et")});
    // }

    auto towers = event->GetCollection("CaloTower");
    vector<pair<float, float>> towersPhiAndEt;

    for (auto physicsObject : *towers) {
      auto tower = asCaloTower(physicsObject);
      if (tower->IsDead()) continue;
      if (tower->IsEtaAboveLimit()) continue;
      if (tower->IsInHEM()) continue;
      if (tower->IsInElectromagneticCrack()) continue;
      //   if (tower->OverlapsWithOtherObjects(event->GetCollection("goodPhoton"))) continue;
      //   if (tower->OverlapsWithOtherObjects(event->GetCollection("goodElectron"))) continue;
      //   if (!tower->IsElectromagneticEnergyAboveNoiseThreshold()) continue;

      float phi = tower->Get("phi");
      float et = tower->Get("et");

      towersPhiAndEt.push_back({phi, et});
    }

    auto [maxTowerEt, refTowerPhi] = findRefEtandPhi(towersPhiAndEt);

    vector<tuple<float, float, int>> tracksPhiPtCharge;
    auto tracks = event->GetCollection("goodTrack");

    for(auto track : *tracks){
      tracksPhiPtCharge.push_back({(float)track->Get("phi")-refTowerPhi, (float)track->Get("pt"), (int)track->Get("charge")});
    }

    // drawShapes(towersPhiAndEt, "../plots/visualizations/event_" + to_string(iEvent) + ".pdf", visualizationParams);
    
    // int runNumber = event->Get("runNumber");
    // int lumiSection = event->Get("lumiSection");
    // int eventNumber = event->Get("eventNumber");
    // info() << runNumber << ":" << lumiSection << ":" << eventNumber << endl;

    addShapesToCanvas(towersPhiAndEt, overlapCanvas, visualizationParams);
    // addShapesToCanvas(photonsPhiAndEt, overlapCanvas, visualizationParams);
    addTracksToCanvas(tracksPhiPtCharge, overlapCanvas, visualizationParams);

  }

  overlapCanvas->Update();
  gSystem->mkdir("../plots/visualizations/");
  overlapCanvas->SaveAs(outputPath.c_str());

  gStyle->SetOptStat(0);
  
  info() << "Finishing up" << endl;
  auto &logger = Logger::GetInstance();
  logger.Print();
  return 0;
}