import FWCore.ParameterSet.Config as cms
from DQMServices.Core.DQMEDHarvester import DQMEDHarvester
from Validation.RecoEgamma.egammaValidator_cff import egammaValidator

def _r(name, title, xt, yt, num, den):
    return f"{name} '{title};{xt};{yt}' {num} {den}"

print(egammaValidator.dirName, "egamma dir name")
def _harvester(subdir, ratios, verbose=0):
    return DQMEDHarvester('DQMGenericClient',
        subDirs = cms.untracked.vstring(egammaValidator.dirName.value() + subdir),
        efficiency = cms.vstring(ratios),
        resolution = cms.vstring(),
        cumulativeDists = cms.untracked.vstring(),
        noFlowDists = cms.untracked.vstring(),
        outputFileName = cms.untracked.string(""),
        verbose = cms.untracked.uint32(verbose))

# ---- electronSeeds
_seed = []
for g in ["Ecal", "Trk", "Both", "Merged"]:
    for v, xt in {"eta": "#eta", "pt": "p_{T} [GeV]"}.items():
        _seed.append(_r(f"eff_{g}Seed_{v}", f"{g} seed efficiency", xt, "efficiency",
                        f"num{g}Seed_vs_{v}", f"den{g}Seed_vs_{v}"))
        if g != "Both":
            _seed.append(_r(f"dup_{g}Seed_{v}", f"{g} seed duplicate rate", xt, "duplicate rate",
                            f"num{g}SeedDup_vs_{v}", f"den{g}Seed_vs_{v}"))

# ---- gsfElectrons (one set per collection label)
_labels = egammaValidator.electronCores + egammaValidator.electrons
print(_labels, "labels")
_ele = []
for L in _labels:
    for v, xt in {"eta": "#eta", "energy": "E [GeV]"}.items():
        _ele += [
            _r(f"eff_ele_{v}_{L}", "Supercluster efficiency", xt, "efficiency", f"num_ele_vs_{v}_{L}", f"den_ele_vs_{v}_{L}"),
            _r(f"dup_ele_{v}_{L}", "Supercluster duplicate rate", xt, "duplicate rate", f"numDuplicate_ele_vs_{v}_{L}", f"den_ele_vs_{v}_{L}"),
            _r(f"fake_ele_{v}_{L}", "Supercluster fake rate", xt, "fake rate", f"numFake_ele_vs_{v}_{L}", f"denFake_ele_vs_{v}_{L}"),
            _r(f"mismatch_ele_{v}_{L}", "Supercluster/GSF mismatch rate", xt, "mismatch", f"mismatchSCvsGsf_vs_{v}_{L}", f"denFake_ele_vs_{v}_{L}"),
        ]
    for v, xt in {"eta": "#eta", "pt": "p_{T} [GeV]"}.items():
        _ele += [
            _r(f"eff_track_{v}_{L}", "GSF track efficiency", xt, "efficiency", f"num_track_vs_{v}_{L}", f"den_track_vs_{v}_{L}"),
            _r(f"fake_track_{v}_{L}", "GSF track fake rate", xt, "fake rate", f"numFake_track_vs_{v}_{L}", f"denFake_track_vs_{v}_{L}"),
        ]

# ---- tracksAndClusters
_tc = []
for v, xt in {"eta": "#eta", "ene": "E [GeV]"}.items():
    _tc += [
        _r(f"eff_SC_{v}", "Supercluster efficiency", xt, "efficiency", f"num_SuperCluster_vs_{v}", f"den_SuperCluster_vs_{v}"),
        _r(f"dup_SC_{v}", "Supercluster duplicate rate", xt, "duplicate rate", f"num_SuperClusterDup_vs_{v}", f"den_SuperCluster_vs_{v}"),
        _r(f"fake_SC_{v}", "Supercluster fake rate", xt, "fake rate", f"num_SuperClusterFake_vs_{v}", f"den_SuperClusterFake_vs_{v}"),
    ]
for v, xt in {"eta": "#eta", "pt": "p_{T} [GeV]"}.items():
    _tc += [
        _r(f"eff_GsfTrack_{v}", "GSF track efficiency", xt, "efficiency", f"num_GsfTrack_vs_{v}", f"den_GsfTrack_vs_{v}"),
        _r(f"fake_GsfTrack_{v}", "GSF track fake rate", xt, "fake rate", f"num_GsfTrackFake_vs_{v}", f"den_GsfTrackFake_vs_{v}"),
    ]

postProcessorEGammaSeeds = _harvester("electronSeeds", _seed)
postProcessorEGammaElectrons = _harvester("gsfElectrons", _ele)
postProcessorEGammaTracksAndClusters = _harvester("tracksAndClusters", _tc)

postProcessorEGamma = cms.Sequence(
    postProcessorEGammaSeeds + postProcessorEGammaElectrons + postProcessorEGammaTracksAndClusters)
