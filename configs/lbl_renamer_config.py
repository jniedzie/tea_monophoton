nEvents = -1
printEveryNevents = 1000

branchesNames = {
  "run"                   : ("runNumber"                  , "runNumber/I"  ),
  "lumis"                 : ("lumiSection"                , "lumiSection/I"),
  "event"                 : ("eventNumber"                , "eventNumber/I"),

  "xVtx"                  : ("vertex_x"                   , "vector<float>"),
  "yVtx"                  : ("vertex_y"                   , "vector<float>"),
  "zVtx"                  : ("vertex_z"                   , "vector<float>"),

  "mcEta"                 : ("genParticle_eta"            , "vector<float>"),
  "mcPhi"                 : ("genParticle_phi"            , "vector<float>"),
  "mcEt"                  : ("genParticle_et"             , "vector<float>"),
  "mcPID"                 : ("genParticle_pid"            , "vector<int>"),

  "phoE"                  : ("photon_energy"              , "vector<float>"),
  "phoEt"                 : ("photon_et"                  , "vector<float>"),
  "phoEta"                : ("photon_eta"                 , "vector<float>"),
  "phoPhi"                : ("photon_phi"                 , "vector<float>"),
  "phoHoverE"             : ("photon_hOverE"              , "vector<float>"),
  "phoSigmaEtaEta_2012"   : ("photon_sigmaEta2012"        , "vector<float>"),
  "phoSigmaIEtaIEta_2012" : ("photon_sigmaIEtaIEta2012"   , "vector<float>"),
  "phoMaxEnergyXtal"      : ("photon_maxEnergyCrystal"    , "vector<float>"),
  "phoETop"               : ("photon_energyTop"           , "vector<float>"),
  "phoEBottom"            : ("photon_energyBottom"        , "vector<float>"),
  "phoELeft"              : ("photon_energyLeft"          , "vector<float>"),
  "phoERight"             : ("photon_energyRight"         , "vector<float>"),
  "phoHasConversionTracks": ("photon_hasConversionTracks" , "vector<int>"  ),
  "pho_seedTime"          : ("photon_seedTime"            , "vector<float>"),
  "phoSCE"                : ("photon_SCEnergy"            , "vector<float>"),
  "phoSCEt"               : ("photon_SCEt"                , "vector<float>"),
  "phoSCEta"              : ("photon_SCEta"               , "vector<float>"),
  "phoSCPhi"              : ("photon_SCPhi"               , "vector<float>"),
  "phoSCEtaWidth"         : ("photon_SCEtaWidth"          , "vector<float>"),
  "phoSCPhiWidth"         : ("photon_SCPhiWidth"          , "vector<float>"),

  "pho_ecalClusterIsoR2"  : ("photon_ecalClusterIsoR2"    , "vector<float>"),
  "pho_ecalClusterIsoR3"  : ("photon_ecalClusterIsoR3"    , "vector<float>"),
  "pho_ecalClusterIsoR4"  : ("photon_ecalClusterIsoR4"    , "vector<float>"),
  "pho_ecalClusterIsoR5"  : ("photon_ecalClusterIsoR5"    , "vector<float>"),
  "pho_hcalRechitIsoR1"   : ("photon_hcalRechitIsoR1"     , "vector<float>"),
  "pho_hcalRechitIsoR2"   : ("photon_hcalRechitIsoR2"     , "vector<float>"),
  "pho_hcalRechitIsoR3"   : ("photon_hcalRechitIsoR3"     , "vector<float>"),
  "pho_hcalRechitIsoR4"   : ("photon_hcalRechitIsoR4"     , "vector<float>"),
  "pho_hcalRechitIsoR5"   : ("photon_hcalRechitIsoR5"     , "vector<float>"),
  # "pho_trackIsoR1PtCut20" : ("photon_trackIsoR1PtCut20"   , "vector<float>"),
  # "pho_trackIsoR2PtCut20" : ("photon_trackIsoR2PtCut20"   , "vector<float>"),
  # "pho_trackIsoR3PtCut20" : ("photon_trackIsoR3PtCut20"   , "vector<float>"),
  # "pho_trackIsoR4PtCut20" : ("photon_trackIsoR4PtCut20"   , "vector<float>"),
  # "pho_trackIsoR5PtCut20" : ("photon_trackIsoR5PtCut20"   , "vector<float>"),

  "CaloTower_e"           : ("CaloTower_energy"           , "vector<float>"),
  "CaloTower_et"          : ("CaloTower_et"               , "vector<float>"),
  "CaloTower_eta"         : ("CaloTower_eta"              , "vector<float>"),
  "CaloTower_phi"         : ("CaloTower_phi"              , "vector<float>"),
  "CaloTower_hadE"        : ("CaloTower_hadE"             , "vector<float>"),
  "CaloTower_emE"         : ("CaloTower_emE"              , "vector<float>"),

  "trkPt"                 : ("track_pt"                   , "vector<float>"),
  "trkP"                  : ("track_momentum"             , "vector<float>"),
  "trkEta"                : ("track_eta"                  , "vector<float>"),
  "trkPhi"                : ("track_phi"                  , "vector<float>"),
  "trkcharge"             : ("track_charge"               , "vector<int>"  ),
  "trkValidHits"          : ("track_nValidHits"           , "vector<int>"  ),
  "trkMissHits"           : ("track_nMissHits"            , "vector<int>"  ),
  "trkPurity"             : ("track_purity"               , "vector<float>"),
  "trknormchi2"           : ("track_normalizedChi2"       , "vector<float>"),
  "trkd0"                 : ("track_d0"                   , "vector<float>"),
  "trkdxy"                : ("track_dxy"                  , "vector<float>"),
  "trkdz"                 : ("track_dz"                   , "vector<float>"),
  "trkdxyError"           : ("track_dxyError"             , "vector<float>"),
  "trkdzError"            : ("track_dzError"              , "vector<float>"),
  "trkvx"                 : ("track_vx"                   , "vector<float>"),
  "trkvy"                 : ("track_vy"                   , "vector<float>"),
  "trkvz"                 : ("track_vz"                   , "vector<float>"),
  
  "gentrkPt"                 : ("track_pt"                   , "vector<float>"),
  "gentrkP"                  : ("track_momentum"             , "vector<float>"),
  "gentrkEta"                : ("track_eta"                  , "vector<float>"),
  "gentrkPhi"                : ("track_phi"                  , "vector<float>"),
  "gentrkcharge"             : ("track_charge"               , "vector<int>"  ),
  "gentrkValidHits"          : ("track_nValidHits"           , "vector<int>"  ),
  "gentrkMissHits"           : ("track_nMissHits"            , "vector<int>"  ),
  "gentrkPurity"             : ("track_purity"               , "vector<float>"),
  "gentrknormchi2"           : ("track_normalizedChi2"       , "vector<float>"),
  "gentrkd0"                 : ("track_d0"                   , "vector<float>"),
  "gentrkdxy"                : ("track_dxy"                  , "vector<float>"),
  "gentrkdz"                 : ("track_dz"                   , "vector<float>"),
  "gentrkdxyError"           : ("track_dxyError"             , "vector<float>"),
  "gentrkdzError"            : ("track_dzError"              , "vector<float>"),
  "gentrkvx"                 : ("track_vx"                   , "vector<float>"),
  "gentrkvy"                 : ("track_vy"                   , "vector<float>"),
  "gentrkvz"                 : ("track_vz"                   , "vector<float>"),

  "eleEn"                 : ("electron_energy"            , "vector<float>"),
  "elePt"                 : ("electron_pt"                , "vector<float>"),
  "eleEta"                : ("electron_eta"               , "vector<float>"),
  "elePhi"                : ("electron_phi"               , "vector<float>"),
  "eleCharge"             : ("electron_charge"            , "vector<int>"  ),
  "eleMissHits"           : ("electron_nMissHits"         , "vector<int>"  ),
  "eleHoverE"             : ("electron_hOverE"            , "vector<float>"),
  "eleEoverP"             : ("electron_eOverP"            , "vector<float>"),
  "elePFRelIsoWithEA"     : ("electron_PFRelIsoWithEA"    , "vector<float>"),
  "eledEtaAtVtx"          : ("electron_deltaEtaAtVertex"  , "vector<float>"),
  "elePFChIso"            : ("electron_PFChIso"           , "vector<float>"),
  "elePFPhoIso"           : ("electron_PFPhoIso"          , "vector<float>"),
  "elePFNeuIso"           : ("electron_PFNeuIso"          , "vector<float>"),
  "eleConvVeto"           : ("electron_conversionVeto"    , "vector<int>"  ),
  "eleSCEt"               : ("electron_SCEt"              , "vector<float>"),
  "eleSCEta"              : ("electron_SCEta"             , "vector<float>"),
  "eleSCPhi"              : ("electron_SCPhi"             , "vector<float>"),
  "eleSCEn"               : ("electron_SCEnergy"          , "vector<float>"),  

  "muPt"                  : ("muon_pt"                    , "vector<float>"),
  "muEta"                 : ("muon_eta"                   , "vector<float>"),
  "muPhi"                 : ("muon_phi"                   , "vector<float>"),
  "muCharge"              : ("muon_charge"                , "vector<int>"  ),
  "muHoverE"              : ("muon_hOverE"                , "vector<float>"),
  "muMissHits"            : ("muon_nMissHits"             , "vector<int>"  ),
  "muPFRelIsoWithEA"      : ("muon_PFRelIsoWithEA"        , "vector<float>"),
  "mudEtaAtVtx"           : ("muon_deltaEtaAtVertex"      , "vector<float>"),
  "muPFChIso"             : ("muon_PFChIso"               , "vector<float>"),
  "muPFPhoIso"            : ("muon_PFPhoIso"              , "vector<float>"),
  "muPFNeuIso"            : ("muon_PFNeuIso"              , "vector<float>"),
  "muSCEta"               : ("muon_SCEta"                 , "vector<float>"),
  "muSCEt"                : ("muon_SCEt"                  , "vector<float>"),
  "muSCPhi"               : ("muon_SCPhi"                 , "vector<float>"),
  "muSCEn"                : ("muon_SCEn"                  , "vector<float>"),

  "egEt"                  : ("egamma_et"                  , "vector<float>"),
  "egEta"                 : ("egamma_eta"                 , "vector<float>"),
  "egPhi"                 : ("egamma_phi"                 , "vector<float>"),

  "nDisplacedTracks"      : ("nDisplacedTracks"           , "nDisplacedTracks/I"),
  "nPixelClusters"        : ("nPixelClusters"             , "nPixelClusters/I"),
  "nPixelRecHits"         : ("nPixelRecHits"              , "nPixelRecHits/I"),
  "nDedxHits"             : ("nDedxHits"                  , "nDedxHits/I"),

  "n"                     : ("nZDC"                       , "nZDC/I"),
  "e"                     : ("ZDC_energy"                 , "ZDC_energy[nZDC]/F"),
  "saturation"            : ("ZDC_saturation"             , "ZDC_saturation[nZDC]/F"),
  "zside"                 : ("ZDC_zSide"                  , "ZDC_zSide[nZDC]/I"),
  "section"               : ("ZDC_section"                , "ZDC_section[nZDC]/I"),
  "channel"               : ("ZDC_channel"                , "ZDC_channel[nZDC]/I"),
  
  "HLT_HIUPC_SingleEG5_NotMBHF2AND_v1": ("SingleEG5"      , "SingleEG5/I"),
  "HLT_HIUPC_DoubleEG2_NotMBHF2AND_v1": ("DoubleEG2"      , "DoubleEG2/I"),
  "HLT_HIZeroBias_v1"                 : ("ZeroBias"       , "ZeroBias/I"),

  # "nCastorTower": ("nCastorTower", "vector<float>"),
  # "CastorTower_hadE": ("CastorTower_hadE", "vector<float>"),
  # "CastorTower_emE": ("CastorTower_emE", "vector<float>"),
  # "CastorTower_p4": ("CastorTower_p4", "vector<float>"),
  # //"C: (//"C, "vector<float>"),
 
  # pixelTree->SetBranchAddress("nPix"            , &nPhysObjects.at(EPhysObjType::kPixelTrack));
  # pixelTree->SetBranchAddress("pixPt"           , &pixelTrackPt);
  # pixelTree->SetBranchAddress("pixP"            , &pixelTrackP);
  # pixelTree->SetBranchAddress("pixEta"          , &pixelTrackEta);
  # pixelTree->SetBranchAddress("pixPhi"          , &pixelTrackPhi);
  # pixelTree->SetBranchAddress("pixcharge"       , &pixelTrackCharge);
  # pixelTree->SetBranchAddress("pixValidHits"    , &pixelTrackValidHits);
  # pixelTree->SetBranchAddress("pixMissHits"     , &pixelTrackMissingHits);
  # pixelTree->SetBranchAddress("pixPurity"       , &pixelTrackPurity);
  # pixelTree->SetBranchAddress("pixnormchi2"     , &pixelTrackChi2);
  # pixelTree->SetBranchAddress("pixdxy"          , &pixelTrackDxy);
  # pixelTree->SetBranchAddress("pixdz"           , &pixelTrackDz);
  # pixelTree->SetBranchAddress("pixdxyError"     , &pixelTrackDxyErr);
  # pixelTree->SetBranchAddress("pixdzError"      , &pixelTrackDzErr);
  # pixelTree->SetBranchAddress("pixvx"           , &pixelTrackVertexX);
  # pixelTree->SetBranchAddress("pixvy"           , &pixelTrackVertexY);
  # pixelTree->SetBranchAddress("pixvz"           , &pixelTrackVertexZ);  
}
