
#include "Framework/EventProcessor.h"

#include "DetDescr/EcalHexReadout.h"

#include "Ecal/Event/EcalHit.h"

// ROOT (MIP tracking)
#include "TMatrixD.h"
#include "TDecompSVD.h"
#include "TVector3.h"

namespace ecal {

class MIPTracking : public framework::Producer {
  // bla
  static float distTwoLines(TVector3 v1, TVector3 v2, TVector3 w1, TVector3 w2);
  static float distPtToLine(TVector3 h1, TVector3 p1, TVector3 p2);
  struct HitData {
    TVector3 pos;
    int layer;
  };
 public:
  MIPTracking(const std::string& name, framework::Process& p)
    : framework::Producer(name,p) {}
  virtual void configure(const framework::config::Parameters& ps) final override;
  virtual void produce(framework::Event& event) final override;
}; // MIPTracking

float MIPTracking::distTwoLines(TVector3 v1, TVector3 v2, TVector3 w1, TVector3 w2) {
  TVector3 e1 = v1 - v2;
  TVector3 e2 = w1 - w2;
  TVector3 crs = e1.Cross(e2);
  if (crs.Mag() == 0) {
    return 100.0;  // arbitrary large number; edge case that shouldn't cause
                   // problems.
  } else {
    return std::abs(crs.Dot(v1 - w1) / crs.Mag());
  }
}

float MIPTracking::distPtToLine(TVector3 h1, TVector3 p1, TVector3 p2) {
  return ((h1 - p1).Cross(h1 - p2)).Mag() / (p1 - p2).Mag();
}

void MIPTracking::configure(const framework::config::Parameters& ps) {
}

void MIPTracking::produce(framework::Event& event) {
  const auto& hexReadout{getCondition<ldmx::EcalHexReadout>(ldmx::EcalHexReadout::CONDITIONS_OBJECT_NAME)};
  // MIP tracking:  vector of hits to be used in the MIP tracking algorithm.  All hits inside the
  // electron ROC (or all hits in the ECal if the event is missing an electron) will be included.
  std::vector<HitData> trackingHitList;
  for (const ldmx::EcalHit& hit : ecalRecHits) {
      // MIP tracking:  Decide whether hit should be added to trackingHitList
      // If inside e- RoC or if etraj is missing, use the hit for tracking:
      if(distance_ele_trajectory >= ele_radii[id.layer()] || distance_ele_trajectory == -1.0) {
        HitData hd;
        hd.pos = TVector3(xy_pair.first, xy_pair.second, hexReadout.getZPosition(id.layer()));
        hd.layer = id.layer();
        trackingHitList.push_back(hd);
      }
  }

  // MIP tracking starts here

  /* Goal:  Calculate 
   *  nStraightTracks (self-explanatory), 
   *  nLinregTracks (tracks found by linreg algorithm),
   */

  // Find epAng and epSep, and prepare EP trajectory vectors:
  TVector3 e_traj_start;
  TVector3 e_traj_end;
  TVector3 p_traj_start;
  TVector3 p_traj_end;
  if(ele_trajectory.size() > 0 && photon_trajectory.size() > 0) {
    // Create TVector3s marking the start and endpoints of each projected trajectory
    e_traj_start.SetXYZ(ele_trajectory[0].first,    ele_trajectory[0].second,   hexReadout.getZPosition(0));
    e_traj_end.SetXYZ(  ele_trajectory[33].first,    ele_trajectory[33].second, hexReadout.getZPosition(33));
    p_traj_start.SetXYZ(photon_trajectory[0].first, photon_trajectory[0].second, hexReadout.getZPosition(0));
    p_traj_end.SetXYZ(  photon_trajectory[33].first, photon_trajectory[33].second, hexReadout.getZPosition(33));

    TVector3 evec   = e_traj_end - e_traj_start;
    TVector3 e_norm = evec.Unit();
    TVector3 pvec   = p_traj_end - p_traj_start;
    TVector3 p_norm = pvec.Unit();
    // Separation variables are currently unused due to pT bias concerns and low efficiency
    // May be removed after a more careful MIP tracking study
    float epDot = e_norm.Dot(p_norm);
    epAng_ = acos(epDot) * 180.0 / M_PI;
    epSep_ = sqrt( pow(e_traj_start.X() - p_traj_start.X(), 2) +
                   pow(e_traj_start.Y() - p_traj_start.Y(), 2) );
  } else {
    // Electron trajectory is missing, so all hits in the Ecal are fair game.
    // Pick e/ptraj so that they won't restrict the tracking algorithm (place them far outside the ECal).
    e_traj_start = TVector3(999,999,hexReadout.getZPosition(0)); //0);
    e_traj_end = TVector3(999,999,hexReadout.getZPosition(33)); // 999);
    p_traj_start = TVector3(1000,1000,hexReadout.getZPosition(0)); //0);
    p_traj_end = TVector3(1000,1000,hexReadout.getZPosition(33)); //1000);
    epAng_ = 3.0 + 1.0;  /*ensures event will not be vetoed by angle/separation cut  */
    epSep_ = 10.0 + 1.0;
  }

  // Near photon step:  Find the first layer of the ECal where a hit near the projected
  // photon trajectory is found
  // Currently unusued pending further study; performance has dropped between v9 and v12.
  firstNearPhLayer_ = hexReadout.getZPosition(33);

  if(photon_trajectory.size() != 0) {  //If no photon trajectory, leave this at the default (ECal back)
    for(std::vector<HitData>::iterator it = trackingHitList.begin(); it != trackingHitList.end(); ++it) {
      float ehDist = sqrt( pow((*it).pos.X() - photon_trajectory[(*it).layer].first,  2)
                         + pow((*it).pos.Y() - photon_trajectory[(*it).layer].second, 2));
      if(ehDist < 8.7 && (*it).layer < firstNearPhLayer_) {
        firstNearPhLayer_ = (*it).layer;
      }
    }
  }

  // Find straight MIP tracks:

  std::sort(trackingHitList.begin(), trackingHitList.end(), [](HitData ha, HitData hb) {return ha.layer > hb.layer;});
  // For merging tracks:  Need to keep track of existing tracks
  // Candidate tracks to merge in will always be in front of the current track (lower z), so only store the last hit
  // 3-layer vector:  each track = vector of 3-tuples (xy+layer).
  std::vector<std::vector<HitData>> track_list;
  
  float cellWidth = 8.7;
  for (int iHit = 0; iHit < trackingHitList.size(); iHit++) {
    int track[34];  // list of hit numbers in track (34 = maximum theoretical length)
    int currenthit;
    int trackLen;

    track[0] = iHit;
    currenthit = iHit;
    trackLen = 1;

    // Search for hits to add to the track:  if hit is in the next two layers behind the current hit,
    // consider adding.
    for (int jHit = 0; jHit < trackingHitList.size(); jHit++) {
      if (trackingHitList[jHit].layer == trackingHitList[currenthit].layer ||
          trackingHitList[jHit].layer > trackingHitList[currenthit].layer + 3) {
        continue;  // if not in the next two layers, continue
      }
      //if it is:  add to track if new hit is directly behind the current hit.
      if (trackingHitList[jHit].pos.X() == trackingHitList[currenthit].pos.X() &&
          trackingHitList[jHit].pos.Y() == trackingHitList[currenthit].pos.Y()) {
        track[trackLen] = jHit;
        trackLen++;
      }
    }
    // Confirm that the track is valid:
    if (trackLen < 2) continue;   // Track must contain at least 2 hits
      float closest_e = distTwoLines(trackingHitList[track[0]].pos, trackingHitList[track[trackLen-1]].pos,
                                     e_traj_start, e_traj_end);
      float closest_p = distTwoLines(trackingHitList[track[0]].pos, trackingHitList[track[trackLen-1]].pos,
                                     p_traj_start, p_traj_end);
      // Make sure that the track is near the photon trajectory and away from the electron trajectory
      // Details of these constraints may be revised
      if (closest_p > cellWidth and closest_e < 2*cellWidth) continue;
      if (trackLen < 4 and closest_e > closest_p) continue;

      //if track found, increment nStraightTracks and remove all hits in track from future consideration
      if (trackLen >= 2) {
        std::vector<HitData> temp_track_list;
        for (int kHit = 0; kHit < trackLen; kHit++) {
          temp_track_list.push_back(trackingHitList[track[kHit]]);
          trackingHitList.erase(trackingHitList.begin() + track[kHit]);
        }
        track_list.push_back(temp_track_list);
        //The *current* hit will have been removed, so iHit is currently pointing to the next hit.
        iHit--;  //Decrement iHit so no hits will get skipped by iHit++
        //nStraightTracks_++; // moved to post-merging
      }
  }

  //Optional addition:  Merge nearby straight tracks.  Not necessary for veto.
  // Criteria:  consider tail of track.  Merge if head of next track is 1/2 layers behind, within 1 cell of xy position.
  //std::cout << "BEGINNING TRACK MERGING" << std::endl;

  for (int track_i = 0; track_i < track_list.size(); track_i++) {
    // for each track, check the remainder of the track list for compatible tracks
    std::vector<HitData> base_track = track_list[track_i];
    HitData tail_hitdata = base_track.back(); // xylayer of last hit in track
    //std::cout << "  Considering track " << track_i << std::endl;
    for (int track_j = track_i+1; track_j < track_list.size(); track_j++) {
      //std::cout << "    Checking for compatibility: " << track_j << std::endl;
      std::vector<HitData> checking_track = track_list[track_j];
      HitData head_hitdata = checking_track.front();
      // if 1-2 layers behind, and xy within one cell...
      if ((head_hitdata.layer == tail_hitdata.layer+1 || head_hitdata.layer == tail_hitdata.layer+2)
          && pow(pow(head_hitdata.pos.X()-tail_hitdata.pos.X(),2)
               + pow(head_hitdata.pos.Y()-tail_hitdata.pos.Y(),2), 0.5) <= cellWidth) {
        // ...then append the second track to the first one and delete it
        // NOTE:  TO ADD:  (trackingHitList[iHit].pos - trackingHitList[jHit].pos).Mag()
        /*
        std::cout << "     **Compatible track found!  Adding track, deleting stuff..." << std::endl;
        std::cout << "     Tail xylayer: " << head_hitdata.pos.X() << "," << head_hitdata.pos.Y() << "," << head_hitdata.layer << std::endl;
        std::cout << "     Head xylayer: " << tail_hitdata.pos.X() << "," << tail_hitdata.pos.Y() << "," << tail_hitdata.layer << std::endl;
        */
        for (int hit_k = 0; hit_k < checking_track.size(); hit_k++) {
          base_track.push_back(track_list[track_j][hit_k]);
        }
        track_list.erase(track_list.begin() + track_j);
        break;
      }
    }
  }
  nStraightTracks_ = track_list.size();

  // Linreg tracking:

  //std::cout << "Finding linreg tracks..." << std::endl;

  for (int iHit = 0; iHit < trackingHitList.size(); iHit++) {
    //std::cout << "hit " << iHit << std::endl;
    int track[34];
    int trackLen;
    int currenthit;
    int hitsInRegion[50]; // Hits being considered at one time
    int nHitsInRegion;    // Number of hits under consideration
    TMatrixD svdMatrix(3,3);
    TMatrixD Vm(3,3);
    TMatrixD hdt(3,3);
    TVector3 slopeVec;
    TVector3 hmean;
    TVector3 hpoint;
    float r_corr_best;
    int hitNums_best[3]; // Hit numbers of current best track candidate
    int hitNums[3];

    trackLen = 0;
    nHitsInRegion = 1;
    currenthit = iHit;
    hitsInRegion[0] = iHit;
    //std::cout << "filling hitsInRegion, tracking list size=" << trackingHitList.size() << std::endl;
    // Find all hits within 2 cells of the primary hit:
    for (int jHit = 0; jHit < trackingHitList.size(); jHit++) {
      //std::cout << "iHit, jHit:  " << iHit << ", " << jHit << std::endl;
      float dstToHit = (trackingHitList[iHit].pos - trackingHitList[jHit].pos).Mag();
      //std::cout << "  Filling hitsInRegion, tracking list len=" << trackingHitList.size() << std::endl;
      if (dstToHit <= 2*cellWidth) {
        hitsInRegion[nHitsInRegion] = jHit;
        nHitsInRegion++;
      }
    }

    //std::cout << "...ready" << std::endl;

    // Look at combinations of hits within the region (do not consider the same combination twice):
    hitNums[0] = iHit;
    for (int jHit = 1; jHit < nHitsInRegion - 1; jHit++) {
      if (trackingHitList.size() < 3) break;
      hitNums[1] = jHit;
      for (int kHit = jHit + 1; kHit < nHitsInRegion; kHit++) {
        hitNums[2] = kHit;
        for (int hInd = 0; hInd < 3; hInd++) {
          //std::cout << "    inner loop" << std::endl;
          // hmean = geometric mean, subtract off from hits to improve SVD performance
          hmean(hInd) = (trackingHitList[hitNums[0]].pos(hInd) +
                         trackingHitList[hitNums[1]].pos(hInd) +
                         trackingHitList[hitNums[2]].pos(hInd))/3.0;
        }
        for (int hInd = 0; hInd < 3; hInd++) {
          for (int lInd = 0; lInd < 3; lInd++) {
            hdt(hInd,lInd) = trackingHitList[hitNums[hInd]].pos(lInd) - hmean(lInd);
          }
        }

        // Perform "linreg" on selected points:
        TDecompSVD svdObj = TDecompSVD(hdt);
        bool decomposed = svdObj.Decompose();
        //std::cout << "    Decomposed" << std::endl;
        if (!decomposed) continue;

        Vm = svdObj.GetV();  // First col of V matrix is the slope of the best-fit line
        for (int hInd = 0; hInd < 3; hInd++) {
          slopeVec(hInd) = Vm[0][hInd];
        }
        hpoint = slopeVec + hmean;  // hmean, hpoint are points on the best-fit line
        //linreg complete:  Now have best-fit line for 3 hits under consideration
        //Check whether the track is valid:  r^2 must be high, and the track must plausibly originate from the photon
        float closest_e = distTwoLines(hmean, hpoint, e_traj_start, e_traj_end);
        float closest_p = distTwoLines(hmean, hpoint, p_traj_start, p_traj_end);
        // Projected track must be close to the photon; details may change after future study.
        if (closest_p > cellWidth or closest_e < 1.5*cellWidth) continue;
        //find r^2
        float vrnc = (trackingHitList[hitNums[0]].pos - hmean).Mag() +
                     (trackingHitList[hitNums[1]].pos - hmean).Mag() +
                     (trackingHitList[hitNums[2]].pos - hmean).Mag();  // ~variance
        float sumerr = distPtToLine(trackingHitList[hitNums[0]].pos, hmean, hpoint) +
                       distPtToLine(trackingHitList[hitNums[1]].pos, hmean, hpoint) +
                       distPtToLine(trackingHitList[hitNums[2]].pos, hmean, hpoint); // sum of |errors|
        float r_corr = 1 - sumerr/vrnc;
        // Check whether r^2 exceeds a low minimum r_corr:  "Fake" tracks are still much more common in background, so making the algorithm
        // oversensitive doesn't lower performance significantly
        if (r_corr > r_corr_best and r_corr > .6) {
          r_corr_best = r_corr;
          trackLen = 0;
          for (int k=0; k<3; k++) { // Only looking for 3-hit tracks currently
            track[k] = hitNums[k];
            trackLen++;
          }
        }
        //std::cout << "    inner loop done!" << std::endl;
      }
    }
    // Ordinarily, additional hits in line w/ track would be added here.  However, this doesn't affect the results of the simple veto.
    // Exclude all hits in a found track from further consideration:
    if (trackLen >= 2) {
      nLinregTracks_++;
      //std::cout << "  Hitlist size = " << trackingHitList.size() << std::endl;
      for (int kHit = 0; kHit < trackLen; kHit++) {
        trackingHitList.erase(trackingHitList.begin() + track[kHit]);
        //std::cout << "  Hitlist size = " << trackingHitList.size() << std::endl;
      }
      iHit--;
    }
  }

  //std::cout << "  MIP tracking completed" << std::endl;
}

}

DECLARE_PRODUCER_NS(ecal,MIPTracking);
