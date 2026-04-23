#include "LbLEvent.hpp"

#include <cctype>
#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <map>
#include <optional>
#include <sstream>
#include <set>
#include <vector>

#include "ExtensionsHelpers.hpp"
#include "UserExtensionsHelpers.hpp"

using namespace std;

namespace {
  constexpr int kBunchesPerOrbit = 3564;
  const string kMonoSchemesBaseUrl = "https://lpc.web.cern.ch/fillingSchemes/2018";

  string Trim(string value) {
    value.erase(value.begin(), find_if(value.begin(), value.end(), [](unsigned char c) { return !isspace(c); }));
    value.erase(find_if(value.rbegin(), value.rend(), [](unsigned char c) { return !isspace(c); }).base(), value.end());
    return value;
  }

  bool SplitIntPair(const string& line, int& first, int& second) {
    stringstream stream(line);
    string firstToken;
    string secondToken;

    if (!getline(stream, firstToken, ';')) return false;
    if (!getline(stream, secondToken, ';')) return false;

    firstToken = Trim(firstToken);
    secondToken = Trim(secondToken);

    try {
      first = stoi(firstToken);
      second = stoi(secondToken);
    } catch (...) {
      return false;
    }

    return true;
  }

  bool SplitIntStringPair(const string& line, int& first, string& second) {
    stringstream stream(line);
    string firstToken;

    if (!getline(stream, firstToken, ';')) return false;
    if (!getline(stream, second, ';')) return false;

    firstToken = Trim(firstToken);
    second = Trim(second);

    try {
      first = stoi(firstToken);
    } catch (...) {
      return false;
    }

    return !second.empty();
  }

  bool EnsureSchemeFileExists(const string& utilsDir, const string& schemeName, string& schemePath) {
    const filesystem::path schemesDir = filesystem::path(utilsDir) / "schemes";
    schemePath = (schemesDir / (schemeName + ".csv")).string();

    if (filesystem::exists(schemePath)) return true;

    info() << "Scheme file not found locally, will try to download: " << schemePath << endl;

    error_code errorCode;
    filesystem::create_directories(schemesDir, errorCode);
    if (errorCode) {
      error() << "Failed to create schemes directory " << schemesDir.string() << ": " << errorCode.message() << endl;
      return false;
    }

    const string url = kMonoSchemesBaseUrl + "/" + schemeName + ".csv";
    const string command = "curl -fsSL \"" + url + "\" -o \"" + schemePath + "\"";

    info() << "Downloading mono scheme from " << url << endl;
    const int downloadStatus = system(command.c_str());

    if (downloadStatus != 0 || !filesystem::exists(schemePath)) {
      error() << "Failed to download mono scheme file from " << url << endl;
      return false;
    }

    info() << "Downloaded mono scheme to " << schemePath << endl;
    return true;
  }
}

struct LbLEvent::MonoCollisionCache {
  map<int, int> runToFill;
  map<int, string> fillToScheme;
  map<string, optional<set<int>>> collidingBXsByScheme;
};

string LbLEvent::GetMonoUtilsDir() {
  const string sourcePath = __FILE__;
  const string marker = "libs/user_extensions/src/LbLEvent.cpp";
  const auto markerPos = sourcePath.rfind(marker);

  if (markerPos == string::npos) return "utils";

  return sourcePath.substr(0, markerPos) + "utils";
}

LbLEvent::MonoCollisionCache LbLEvent::LoadMonoCollisionCache() {
  MonoCollisionCache cache;
  const string utilsDir = GetMonoUtilsDir();

  info() << "Loading mono collision cache from " << utilsDir << endl;

  {
    const string runToFillPath = utilsDir + "/mono_run_to_fill.csv";
    ifstream file(runToFillPath);
    string line;
    int loadedRows = 0;

    if (!file.is_open()) {
      fatal() << "Failed to open mono run-to-fill mapping: " << runToFillPath << endl;
      throw runtime_error("Failed to open mono_run_to_fill.csv");
    }

    while (getline(file, line)) {
      int run = 0;
      int fill = 0;
      if (SplitIntPair(line, run, fill)) {
        cache.runToFill[run] = fill;
        loadedRows++;
      }
    }

    info() << "Loaded " << loadedRows << " run-to-fill entries" << endl;
  }

  {
    const string fillToSchemePath = utilsDir + "/mono_fill_to_scheme.csv";
    ifstream file(fillToSchemePath);
    string line;
    int loadedRows = 0;

    if (!file.is_open()) {
      fatal() << "Failed to open mono fill-to-scheme mapping: " << fillToSchemePath << endl;
      throw runtime_error("Failed to open mono_fill_to_scheme.csv");
    }

    while (getline(file, line)) {
      int fill = 0;
      string scheme;
      if (SplitIntStringPair(line, fill, scheme)) {
        cache.fillToScheme[fill] = scheme;
        loadedRows++;
      }
    }

    info() << "Loaded " << loadedRows << " fill-to-scheme entries" << endl;
  }

  for (const auto& [fill, scheme] : cache.fillToScheme) {
    if (cache.collidingBXsByScheme.count(scheme) > 0) continue;

    string schemePath;
    if (!EnsureSchemeFileExists(utilsDir, scheme, schemePath)) {
      warn() << "Mono scheme file is unavailable for fill " << fill << ": " << scheme
             << ". Collision-in-previous-BX information will be unavailable for this scheme." << endl;
      cache.collidingBXsByScheme[scheme] = nullopt;
      continue;
    }

    ifstream file(schemePath);
    string line;
    bool isFirstLine = true;
    set<int> beam1;
    set<int> beam2;

    if (!file.is_open()) {
      warn() << "Failed to open mono scheme file for fill " << fill << ": " << schemePath
             << ". Collision-in-previous-BX information will be unavailable for this scheme." << endl;
      cache.collidingBXsByScheme[scheme] = nullopt;
      continue;
    }

    while (getline(file, line)) {
      if (isFirstLine) {
        isFirstLine = false;
        continue;
      }

      stringstream stream(line);
      string token;
      vector<string> parts;

      while (getline(stream, token, ',')) parts.push_back(Trim(token));
      if (parts.size() < 4) continue;

      int rfBucket = 0;

      try {
        rfBucket = stoi(parts[3]);
      } catch (...) {
        continue;
      }

      const int bx = (rfBucket - 1) / 10;

      if (parts[2] == "ring_1") {
        beam1.insert(bx);
      } else if (parts[2] == "ring_2") {
        beam2.insert(bx);
      }
    }

    set<int> colliding;
    for (const auto& bx : beam1) {
      if (beam2.count(bx) > 0) colliding.insert(bx);
    }

    cache.collidingBXsByScheme[scheme] = colliding;
    info() << "Loaded scheme " << scheme << " with " << colliding.size() << " colliding BXs" << endl;
  }

  info() << "Mono collision cache ready with " << cache.collidingBXsByScheme.size() << " schemes" << endl;
  return cache;
}

const LbLEvent::MonoCollisionCache& LbLEvent::GetMonoCollisionCache() {
  static const MonoCollisionCache cache = LoadMonoCollisionCache();
  return cache;
}

optional<bool> LbLEvent::HasCollisionInPreviousBXs(int nBXs) {
  const auto& cache = GetMonoCollisionCache();
  const int runNumber = GetAs<int>("runNumber");
  const int bx = GetAs<int>("bunchNumber");

  if (cache.runToFill.count(runNumber) == 0) {
    warn() << "No mono fill mapping found for run " << runNumber << endl;
    return nullopt;
  }

  const int fill = cache.runToFill.at(runNumber);
  if (cache.fillToScheme.count(fill) == 0) {
    warn() << "No mono scheme mapping found for fill " << fill << " (run " << runNumber << ")" << endl;
    return nullopt;
  }

  const string& scheme = cache.fillToScheme.at(fill);
  if (cache.collidingBXsByScheme.count(scheme) == 0) {
    warn() << "No cached colliding-BX info found for scheme " << scheme << endl;
    return nullopt;
  }

  const auto& collidingBXs = cache.collidingBXsByScheme.at(scheme);
  if (!collidingBXs.has_value()) {
    warn() << "Could not determine colliding-BX info for scheme " << scheme << endl;
    return nullopt;
  }

  const int normalizedBX = ((bx % kBunchesPerOrbit) + kBunchesPerOrbit) % kBunchesPerOrbit;

  for (int distance = 1; distance <= nBXs; distance++) {
    const int previousBX = (normalizedBX - distance + kBunchesPerOrbit) % kBunchesPerOrbit;
    if (collidingBXs->count(previousBX) > 0) return true;
  }

  return false;
}

float LbLEvent::GetDeltaEt() {
  auto photons = GetCollection("goodPhoton");

  if (photons->size() != 2) {
    warn() << "Couldn't calculate deltaEt -- number of photons in the event != 2" << endl;
    return -1;
  }

  auto towers = GetCollection("CaloTower");

  auto photon1 = photons->at(0);
  auto photon2 = photons->at(1);

  float photon1Et = photon1->Get("et");
  float photon2Et = photon2->Get("et");

  float photon1Phi = photon1->Get("phi");
  float photon2Phi = photon2->Get("phi");

  float photon1Eta = photon1->Get("eta");
  float photon2Eta = photon2->Get("eta");

  float maxDeltaR = 1.5;
  float highestTowerEt1 = 0;
  float highestTowerEt2 = 0;

  for (auto physicsObject : *towers) {
    auto tower = asCaloTower(physicsObject);
    if (tower->IsDead()) continue;
    if (tower->IsEtaAboveLimit()) continue;
    if (tower->IsInHEM()) continue;
    if (tower->IsInElectromagneticCrack()) continue;

    float phi = tower->Get("phi");
    float eta = tower->Get("eta");
    float et = tower->Get("et");

    if (sqrt(pow(phi - photon1Phi, 2) + pow(eta - photon1Eta, 2)) < maxDeltaR) {
      if (et > highestTowerEt1) highestTowerEt1 = et;
    }
    if (sqrt(pow(phi - photon2Phi, 2) + pow(eta - photon2Eta, 2)) < maxDeltaR) {
      if (et > highestTowerEt2) highestTowerEt2 = et;
    }
  }

  float deltaEt1 = fabs(photon1Et - highestTowerEt1) / photon1Et;
  float deltaEt2 = fabs(photon2Et - highestTowerEt2) / photon2Et;

  float maxDelta = max(deltaEt1, deltaEt2);
  return maxDelta;
}

float LbLEvent::GetCosThetaStar(bool doElectrons) {
  auto objects = GetCollection(doElectrons ? "goodElectron" : "goodPhoton");

  if (objects->size() != 2) {
    warn() << "Couldn't calculate cosThetaStar -- number of objects in the event != 2" << endl;
    return -1;
  }

  TLorentzVector object1, object2;

  if (doElectrons) {
    object1 = asElectron(objects->at(0))->GetFourMomentum();
    object2 = asElectron(objects->at(1))->GetFourMomentum();
  } else {
    object1 = asPhoton(objects->at(0))->GetFourMomentum();
    object2 = asPhoton(objects->at(1))->GetFourMomentum();
  }
  auto objectsSum = object1 + object2;

  float costhetastarCS = 2. * (objectsSum.E() * object1.Pz() - objectsSum.Pz() * object1.E()) / (objectsSum.M() * objectsSum.Mt());
  return costhetastarCS;
}

float LbLEvent::GetDiphotonAcoplanarity() {
  auto photons = event->GetCollection("goodPhoton");
  if (photons->size() != 2) return -1;

  double deltaPhi = asPhoton(photons->at(0))->GetFourMomentum().DeltaPhi(asPhoton(photons->at(1))->GetFourMomentum());
  double acoplanarity = 1 - (fabs(deltaPhi) / TMath::Pi());

  return acoplanarity;
}

vector<shared_ptr<PhysicsObject>> LbLEvent::GetGenPhotons() {
  vector<shared_ptr<PhysicsObject>> genPhotons;

  auto genParticles = GetCollection("genParticle");
  
  for (auto physicsObject : *genParticles) {
    
    int particlePid = physicsObject->Get("pid");
    float *floatPtr = reinterpret_cast<float *>(&particlePid);
    float floatValue = *floatPtr;
    particlePid = round(floatValue);
    
    if (abs(particlePid) == 22) {
      genPhotons.push_back(physicsObject);
    }
  }

  return genPhotons;
}

vector<shared_ptr<PhysicsObject>> LbLEvent::GetGenMatchedRecoPhotons() {
  vector<shared_ptr<PhysicsObject>> matchedPhotons;

  auto genPhotons = GetGenPhotons();
  auto recoPhotons = GetCollection("photon");

  for (auto genPhoton : genPhotons) {
    float genEta = genPhoton->Get("eta");
    float genPhi = genPhoton->Get("phi");

    shared_ptr<PhysicsObject> bestMatch = nullptr;
    float bestDeltaR = 0.2; // matching cone

    for (auto recoPhoton : *recoPhotons) {
      float recoEta = recoPhoton->Get("eta");
      float recoPhi = recoPhoton->Get("phi");

      float deltaR = sqrt(pow(genEta - recoEta, 2) + pow(genPhi - recoPhi, 2));
      if (deltaR < bestDeltaR) {
        bestDeltaR = deltaR;
        bestMatch = recoPhoton;
      }
    }

    if (bestMatch != nullptr) {
      matchedPhotons.push_back(bestMatch);
    }
  }

  return matchedPhotons;
}
