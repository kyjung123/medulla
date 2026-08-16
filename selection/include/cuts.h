/**
 * @file cuts.h
 * @brief Header file for definitions of analysis cuts.
 * @details This file contains definitions of analysis cuts which can be used
 * to select interactions. Each cut is implemented as a function which takes an
 * interaction object as an argument and returns a boolean. These are the
 * building blocks for defining more complex selections.
 * @author mueller@fnal.gov
*/
#ifndef CUTS_H
#define CUTS_H
#include <vector>
#include <numeric>
#include <cmath>
#include <algorithm>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <optional>
#include <limits>


#include "utilities.h"
#include "framework.h"
#include "selectors.h"
#include "include/bivariables.h"
#include "include/biselectors.h"

/**
 * @namespace cuts
 * @brief Namespace for organizing generic cuts which act on interactions.
 * @details This namespace is intended to be used for organizing cuts which act
 * on interactions. Each cut is implemented as a function which takes an
 * interaction object as an argument and returns a boolean. The function should
 * be templated on the type of interaction object if the cut is intended to be
 * used on both true and reconstructed interactions.
 */
namespace cuts
{   
    /**
     * @brief Apply a cut on the validity of the flash match.
     * @details A "valid" flash match is defined as a flash-interaction
     * association with a flash time that is not NaN and a flash match
     * status of 1. The upstream flash matching algorithm (OpT0Finder) has a
     * flash filter that restricts candidate flashes to near the beam window,
     * which means that the majority of cosmogenic interactions are not
     * flash matched. If no flash match is found, the flash time is NaN. This
     * cut is intended to be applied as a preselection cut to reduce comparisons
     * to NaN values, which tend to be noisy on stderr.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction on which to place the flash validity cut.
     * @return true if the interaction is flash matched and the time is valid.
     */
    template<class T>
    bool valid_flashmatch(const T & obj)
    {
        return obj.flash_times.size() > 0 && obj.is_flash_matched == 1 && !std::isnan(obj.flash_times[0]);
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, valid_flashmatch, valid_flashmatch);

    /**
     * @brief Apply no cut; all interactions passed.
     * @details This is a placeholder function for a cut which does not apply
     * any selection criteria. It is intended to be used in cases where a cut
     * function is required, but no selection is desired.
     * @tparam T the type of object (true or reco).
     * @param obj the interaction to select on.
     * @return true (always).
     */
    template<class T>
    bool no_cut(const T & obj) { return true; }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, no_cut, no_cut);

    /**
     * @brief Apply a cut to select neutrinos.
     * @details This function applies a cut to select neutrinos. This cut
     * makes use of the is_neutrino flag in the true interaction object and is
     * intended to be used to identify signal neutrinos.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction is a neutrino.
     * @note This cut is intended to be used for identifying neutrinos in
     * truth, which is useful for making signal definitions.
     */
    template<class T>
    bool neutrino(const T & obj) { return obj.nu_id >= 0; }
    REGISTER_CUT_SCOPE(RegistrationScope::True, neutrino, neutrino);

    /**
     * @brief Apply a cut to select cosmogenic interactions.
     * @details This function applies a cut to select cosmogenic interactions.
     * This cut makes use of the is_neutrino flag in the true interaction
     * object and is intended to be used to identify cosmogenic interactions.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction is a cosmogenic interaction.
     * @note This cut is intended to be used for identifying cosmogenic
     * interactions in truth, which is useful for making background definitions.
     */
    template<class T>
    bool cosmic(const T & obj) { return !neutrino(obj); }
    REGISTER_CUT_SCOPE(RegistrationScope::True, cosmic, cosmic);

    /**
     * @brief Apply a cut to select charged current interactions.
     * @details This function applies a cut to select charged current
     * interactions. This cut makes use of the `current_type` attribute in the
     * true interaction object and is intended to be used to identify charged
     * current interactions.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction is a charged current interaction.
     */
    template<class T>
    bool iscc(const T & obj) { return obj.current_type == 0; }
    REGISTER_CUT_SCOPE(RegistrationScope::True, iscc, iscc);

    /**
     * @brief Apply a cut on the neutrino pdg.
     * @details This function applies a cut to select interactions based on
     * the neutrino pdg
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @param params the parameters for the cut. In this case, this is a vector
     * of neutrino pdg codes accepted
     * @return true if the neutrino pdg is one of the specified pdgs.
     */
    template<class T>
    bool is_neutrino_pdg(const T & obj, std::vector<double> params={})
    {
        if(params.empty())
            return true; // No cut applied if no parameters are given.
        return std::find(params.begin(), params.end(), obj.pdg_code) != params.end();
    }
    REGISTER_CUT_SCOPE(RegistrationScope::True, is_neutrino_pdg, is_neutrino_pdg);
  
    /**
     * @brief Apply a cut on the interaction mode.
     * @details This function applies a cut to select interactions based on
     * the interaction mode. The interaction mode is stored by Genie as an
     * enumerated category.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @param params the parameters for the cut. In this case, this is a vector
     * of interaction modes to select on.
     * @return true if the interaction mode is one of the specified modes.
     */
    template<class T>
    bool is_interaction_mode(const T & obj, std::vector<double> params={})
    {
        if(params.empty())
            return true; // No cut applied if no parameters are given.
        return std::find(params.begin(), params.end(), obj.interaction_mode) != params.end();
    }
    REGISTER_CUT_SCOPE(RegistrationScope::True, is_interaction_mode, is_interaction_mode);

    /**
     * @brief Apply a cut on the interaction type.
     * @details This function applies a cut to select interactions based on
     * the interaction type. The interaction type is stored by Genie as an
     * enumerated category.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @param params the parameters for the cut. In this case, this is a vector
     * of interaction types to select on.
     * @return true if the interaction type is one of the specified types.
     */
    template<class T>
    bool is_interaction_type(const T & obj, std::vector<double> params={})
    {
        if(params.empty())
            return true; // No cut applied if no parameters are given.
        return std::find(params.begin(), params.end(), obj.interaction_type) != params.end();
    }
    REGISTER_CUT_SCOPE(RegistrationScope::True, is_interaction_type, is_interaction_type);

    /**
     * @brief Checks if the primary lepton is mu+ or e+
     * @details This function checks if the primary lepton comes from an
     * anti-neutrino. This is necessary because the current files do not have
     * mctruth neutrino info, only particle truth info. 
     * @param obj the interaction to select on.
     * @return true if the interaction has a primary lepton from an
     * anti-neutrino.
     */
    template<class T>
    bool primary_lepton_from_antineutrino(const T & obj)
    {
        for(const auto & p : obj.particles)
        {
            if((p.pdg_code == -13.0 || p.pdg_code == -11.0) && pvars::primary_classification(p))
                return true;
        }
        return false;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::True, primary_lepton_from_antineutrino, primary_lepton_from_antineutrino);

    /**
     * @brief Apply a fiducial volume cut; the interaction vertex must be
     * reconstructed within the fiducial volume.
     * @details The fiducial volume cut is applied on the reconstructed
     * interaction vertex upstream in SPINE. The fiducial volume is defined
     * (in a SPINE post-processor) as a 25 cm border around the x and y
     * detector faces, a 50 cm border around the downstream (+) z face, and a
     * 30 cm border around the upstream (-) z face. The fiducial volume is
     * intended to reduce the impact of detector edge effects on the analysis.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the vertex is in the fiducial volume.
     */
    //ICARUS volume;: 
    //Positive-x TPC:
//  61.94 + 25 < x < 358.49 - 25
//  86.94 < x < 333.49 cm
//
//y:
//  -181.86 + 25 < y < 134.96 - 25
//  -156.86 < y < 109.96 cm
//
//z:
//  -894.95 + 30 < z < 894.95 - 50
//  -864.95 < z < 844.95 cm
    template<class T>
    bool fiducial_cut(const T & obj)
    {
        return obj.is_fiducial && !(obj.vertex[0] > 210.215 && obj.vertex[1] > 60 && (obj.vertex[2] > 290 && obj.vertex[2] < 390));
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, fiducial_cut, fiducial_cut);

    template<class T>
    bool fiducial_cut_tmp(const T & obj)
    {
        return (
            (abs(obj.vertex[0]) > 10) &&
            (abs(obj.vertex[0]) < 190) &&
            (obj.vertex[2] > 10) &&
            (obj.vertex[2] < 450) &&
            (
                ((obj.vertex[2] > 250) && (obj.vertex[1] > -190) && (obj.vertex[1] < 100)) ||
                ((obj.vertex[2] < 250) && (abs(obj.vertex[1]) < 190))
            )

        );
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, fiducial_cut_tmp, fiducial_cut_tmp);

    /**
     * @brief Veto interactions whose vertex falls in the ICARUS z-gap region.
     * @details A region near |z| < 100 cm in ICARUS exhibits anomalously high
     * rates of poorly-reconstructed interactions ("mystery z-gap"). This cut
     * rejects any interaction whose vertex z-coordinate falls in the range
     * (-100, 100) cm, which is applied on top of the standard fiducial cut
     * when running the pi0 selection.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the vertex z-coordinate is outside (-100, 100) cm.
     */
    template<class T>
    bool avoid_icarus_mystery_zgap(const T & obj)
    {
        return !(obj.vertex[2] > -100 && obj.vertex[2] < 100);
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, avoid_icarus_mystery_zgap, avoid_icarus_mystery_zgap);

    /**
     * @brief Apply a containment cut on the entire interaction.
     * @details The containment cut is applied on the entire interaction. The
     * interaction is considered contained if all particles and all spacepoints
     * are contained within 5cm of the detector edges (configured in a SPINE 
     * post-processor). Additionally, no spacepoints are allowed to be
     * reconstructed in a TPC that did not create it. This is an unphysical
     * condition that can occur when a cosmic muon is moved according to an
     * assumed t0 that is very out-of-time.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction is contained.
     */
    template<class T>
    bool containment_cut(const T & obj) { return obj.is_contained; }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, containment_cut, containment_cut);

    /**
     * @brief Apply a cut to select cathode-crossing interactions.
     * @details This cut is intended to be used in analyses that wish to select
     * (or deselect) interactions that cross the cathode. The cathode-crossing 
     * status is determined by checking the sign of the product of particle
     * extrema x-positions for all particles in the interaction.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction crosses the cathode.
    */
    template<class T>
    bool cathode_crosser(const T & obj)
    {
        for(const auto & p : obj.particles)
        {
            if(!std::isnan(pvars::start_x(p)) && !std::isnan(pvars::end_x(p)))
            {
                if(context::current_detector == caf::Det_t::kSBND)
                {
                    if(pvars::start_x(p) * pvars::end_x(p) < 0)
                        return true;
                }
                else if(context::current_detector == caf::Det_t::kICARUS)
                {
                    if((pvars::start_x(p) < 0 && (pvars::start_x(p) + 210.215) * (pvars::end_x(p) + 210.215) < 0)
                    || (pvars::start_x(p) > 0 && (pvars::start_x(p) - 210.215) * (pvars::end_x(p) - 210.215) < 0))
                        return true;
                }
            }
        }
        return false;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, cathode_crosser, cathode_crosser);

    /**
     * @brief Apply a cut to fiducialize the region around the cathode.
     * @details This cut is intended to be used in analyses that wish to select
     * (or deselect) interactions that occur near the cathode. The
     * fiducialization is applied by checking the x-position of the interaction
     * vertex relative to the cathode position.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction is fiducialized.
     */
    template<class T>
    bool fiducialize_cathode(const T & obj)
    {
        if(context::current_detector == caf::Det_t::kSBND)
        {
            // Apply a cut to fiducialize 5 cm around the cathode of SBND.
            return std::abs(obj.vertex[0]) > 5.0;
        }
        else if(context::current_detector == caf::Det_t::kICARUS)
        {
            // Apply a cut to fiducialize 5 cm around the cathode of ICARUS.
            return std::abs(std::abs(obj.vertex[0]) - 210.215) > 5.0;
        }
        else
        {
            // If the detector is not recognized, do not apply any cut.
            return true;
        }
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, fiducialize_cathode, fiducialize_cathode);

    /**
     * @brief Apply a cut to veto the high-y, high-z region of SBND.
     * @details This cut is intended to be used in analyses that wish to veto
     * the high-y, high-z region of SBND. The high-y, high-z region is the
     * region in positive x where y > 100 cm and z > 250 cm. This region is the
     * subject of some unusual detector effect that hasn't been fully diagnosed
     * at the current date. We apply a cut to veto any activity (has any
     * particle terminating in this region) in this region to mitigate the
     * impact of this effect on analyses.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if no particles in the interaction start or end in the
     * high-y, high-z region.
    */
    template<class T>
    bool veto_sbnd_highy_highz(const T & obj)
    {
        if(context::current_detector == caf::Det_t::kSBND)
        {
            // Apply a cut to veto the high-y, high-z region of SBND.
            bool in_veto_region = false;
            for(const auto & p : obj.particles)
            {
                if((pvars::start_y(p) > 100 && pvars::start_z(p) > 250 && pvars::start_x(p) > 0)
                    || (pvars::end_y(p) > 100 && pvars::end_z(p) > 250 && pvars::end_x(p) > 0))
                {
                    in_veto_region = true;
                    break;
                }
            }
            return !in_veto_region;
        }
        else
        {
            // If the detector is not recognized, do not apply any cut.
            return true;
        }
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, veto_sbnd_highy_highz, veto_sbnd_highy_highz);

    template<class T>
    bool neutrino2026_veto(const T & obj)
    {
        if(context::current_detector == caf::Det_t::kSBND)
        {
            // For SBND, we apply the fiducialize_cathode cut, the
            // veto_sbnd_highy_highz cut, and the cathode-crossing cut to veto
            // the regions of the detector that are not well-modeled and not
            // covered by a validated systematic uncertainty. This is intended
            // to be used for Neutrino 2026.
            return fiducialize_cathode(obj) && veto_sbnd_highy_highz(obj) && !cathode_crosser(obj);
        }
        else if(context::current_detector == caf::Det_t::kICARUS)
        {
            // For ICARUS, we apply the fiducialize_cathode cut to veto the
            // region around the cathode that is not well-modeled and not
            // covered by a validated systematic uncertainty. This is intended
            // to be used for Neutrino 2026.
            return fiducialize_cathode(obj);
        }
        else
        {
            // If the detector is not recognized, do not apply any cut.
            return true;
        }
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, neutrino2026_veto, neutrino2026_veto);

    /**
     * @brief Apply a cut to reject events that have a non-electron particle that
     * is not contained.
     * @details This cut is intended to be used in analyses that select electrons
     * in the final state and wish to allow for electrons that exit the detector.
     * All other particles in the interaction must be contained.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if all non-electron particles are contained.
     */
    template<class T>
    bool nonelectron_containment_cut(const T & obj)
    {
        for(const auto & p : obj.particles)
        {
            if(pvars::pid(p) != pvars::kElectron && !pcuts::containment_cut(p))
                return false;
        }
        return true;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, nonelectron_containment_cut, nonelectron_containment_cut);

    /**
     * @brief Apply a cut to reject events that have a non-muon particle that
     * is not contained.
     * @details This cut is intended to be used in analyses that select muons
     * in the final state and wish to allow for muons that exit the detector.
     * All other particles in the interaction must be contained, which is
     * consistent with our ability to reconstruct exiting muons using the MCS
     * technique.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if all non-muon particles are contained.
     */
    template<class T>
    bool nonmuon_containment_cut(const T & obj)
    {
        for(const auto & p : obj.particles)
        {
            if(pvars::pid(p) != pvars::kMuon && !pcuts::containment_cut(p))
                return false;
        }
        return true;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, nonmuon_containment_cut, nonmuon_containment_cut);
    
    template<class T>
    bool muon_containment_cut(const T & obj,std::vector<double> params={25.0,})
    {
        for(const auto & p : obj.particles)
        {
            if(pvars::pid(p) == pvars::kMuon && pvars::primary_classification(p)  && pcuts::containment_cut(p) && pvars::ke(p) >= params[0])
                return true;
        }
        return false;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, muon_containment_cut, muon_containment_cut);

    template<class T>
    bool pion_containment_cut(const T & obj,std::vector<double> params={25.0,})
    {
        for(const auto & p : obj.particles)
        {
            if(pvars::pid(p) == pvars::kPion && pvars::primary_classification(p) && pcuts::containment_cut(p) && pvars::ke(p) >= params[0])
                return true;
        }
        return false;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, pion_containment_cut, pion_containment_cut);

    /**
     * @brief Apply a cut on the "time containment" of the interaction.
     * @details The time containment cut applies additional restriction that
     * stipulate that all spacepoints must be reconstructed in a feasible TPC.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction is time-contained.
     */
    template<class T>
    bool time_containment_cut(const T & obj) { return obj.is_time_contained; }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, time_containment_cut, time_containment_cut);

    /**
     * @brief Apply a flash time cut on the interaction.
     * @details The flash time cut is applied on the interaction. The flash time
     * is required to be within the beam window, which is expected to be
     * [0 us, 1.6 us] for BNB and [0 us, 9.6 us] for NuMI. This cut is intended
     * to reduce the impact of cosmogenic interactions on analyses.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @param params the parameters for the cut.
     * @return true if the interaction has been matched to an in-time flash.
     * @note The switch to the NuMI beam window is applied by the definition of
     * a preprocessor macro (BEAM_IS_NUMI).
     * @note The cut window has been widened to reconcile the beam window as
     * observed in data and simulation.
     */
    template<class T>
    bool flash_cut(const T & obj, std::vector<double> params={})
    {
        if(!valid_flashmatch(obj))
            return false;
        else if(params.size() == 2 && obj.flash_times[0] >= params[0] && obj.flash_times[0] <= params[1])
            return true;
        else if(params.size() !=2)
            return true;
        else
            return false;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, flash_cut, flash_cut);

    /**
     * @brief Base particle multiplicity for a specific multiplicity.
     * @details This function calculates the multiplicity of a specific
     * particle species in an interaction. The particle species is specified by
     * its SPINE PID index. The function counts the number of primary particles
     * of the specified species with a kinetic energy above a given threshold.
     * @tparam obj the interaction to select on.
     * @param mult the desired multiplicity for the specified particle species.
     * @param particle_species the index of the particle species to count.
     * @param params the parameters for the cut. In this case, this sets the
     * kinetic energy threshold for the particle to count towards the
     * multiplicity. The first element of the vector is used for this purpose.
     * @return the multiplicity of the specified particle species terminated at
     * some maximum value (the desired multiplicity + 1).
     */
    template<class T>
    size_t particle_multiplicity(const T & obj, size_t mult, size_t particle_species, std::vector<double> params={})
    {
        // Default to a kinetic energy threshold of 0 MeV if no parameters are
        // given.
        if(params.empty())
            params.push_back(0.0);

        size_t count(0);
        for(const auto & p : obj.particles)
        {
            if(pvars::pid(p) == particle_species && pvars::primary_classification(p) && pvars::ke(p) >= params[0])
                ++count;
            if(count > mult)
                break; // No need to count further.
        }
        return count;
    }

    /**
     * @brief Binding for a single particle photon multiplicity cut.
     * @details This function binds the single particle multiplicity cut for
     * photons, which corresponds to the index 0 in the
     * @ref utilities::count_primaries function.
     * @param obj the interaction to select on.
     * @param params the parameters for the cut. In this case, this sets the
     * kinetic energy threshold for a photon to count towards the
     * multiplicity. Defaults to 25 MeV.
     * @return true if the interaction has a single primary photon.
     */
    template<class T>
    bool single_photon(const T & obj, std::vector<double> params={25.0,})
    {
        return particle_multiplicity(obj, 1, 0, params) == 1;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, single_photon, single_photon);

    /**
     * @brief Require exactly two above-threshold primary photons in the interaction.
     * @details Counts primary photons whose kinetic energy meets or exceeds the
     * supplied threshold and requires the total to be exactly two. Used as the
     * primary pi0 topology cut to select CC&pi;0 interactions where both photon
     * daughters are reconstructed above threshold. Defaults to 25 MeV.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @param params params[0] sets the photon KE threshold in MeV. Defaults to 25 MeV.
     * @return true if the interaction contains exactly two primary photons above threshold.
     */
    template<class T>
    bool two_photons(const T & obj, std::vector<double> params={25.0,})
    {
        return particle_multiplicity(obj, 2, 0, params) == 2;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, two_photons, two_photons);

    /**
     * @brief Require the diphoton invariant mass to lie within a specified window.
     * @details Selects the best photon pair using the pi0_photon_pair biselector,
     * then computes the diphoton invariant mass as sqrt(2 ke0 ke1 (1 - cos theta)),
     * where the opening angle is determined from the vertex-to-shower-start unit
     * vectors for reconstructed particles and from momentum unit vectors for true
     * particles. The interaction fails the cut if no valid photon pair exists or
     * if the invariant mass falls outside [params[0], params[1]).
     * The KE estimator respects the current pvars::calofn<T> setting.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @param params params[0] lower mass bound (MeV), params[1] upper mass bound (MeV).
     *               Defaults to [60, 300) MeV.
     * @return true if a valid photon pair exists and its invariant mass is in [params[0], params[1]).
     */
    template<class T>
    bool valid_pi0_mass_cut(const T & obj, std::vector<double> params={60.0, 300.0})
    {
        auto [i0, i1] = biselectors::pi0_photon_pair(obj);
        if(i0 == kNoMatch || i1 == kNoMatch) return false;
        const auto & p0 = obj.particles[i0];
        const auto & p1 = obj.particles[i1];
        double vx = obj.vertex[0], vy = obj.vertex[1], vz = obj.vertex[2];
        double ke0 = pvars::calo_ke(p0), ke1 = pvars::calo_ke(p1);
        double ct   = bvars::pi0_opening_costheta_impl(p0, p1, vx, vy, vz);
        double mass = std::sqrt(2.0 * ke0 * ke1 * (1.0 - ct));
        return mass >= params[0] && mass < params[1];
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, valid_pi0_mass_cut, valid_pi0_mass_cut);

    /**
     * @brief Binding for a single particle electron multiplicity cut.
     * @details This function binds the single particle multiplicity cut for
     * electrons, which corresponds to the index 1 in the
     * @ref utilities::count_primaries function.
     * @param obj the interaction to select on.
     * @param params the parameters for the cut. In this case, this sets the
     * kinetic energy threshold for an electron to count towards the
     * multiplicity. Defaults to 25 MeV.
     * @return true if the interaction has a single primary electron.
     */
    template<class T>
    bool single_electron(const T & obj, std::vector<double> params={25.0,})
    {
        return particle_multiplicity(obj, 1, 1, params) == 1;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, single_electron, single_electron);

    /**
     * @brief Binding for a single particle muon multiplicity cut.
     * @details This function binds the single particle multiplicity cut for
     * muons, which corresponds to the index 2 in the
     * @ref utilities::count_primaries function.
     * @param obj the interaction to select on.
     * @param params the parameters for the cut. In this case, this sets the
     * kinetic energy threshold for a muon to count towards the multiplicity.
     * Defaults to 143.425 MeV, which corresponds to a muon of length 50 cm
     * (assuming the muon stops).
     * @return true if the interaction has a single primary muon.
     */
    template<class T>
    bool single_muon(const T & obj, std::vector<double> params={143.425,})
    {
        return particle_multiplicity(obj, 1, 2, params) == 1;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, single_muon, single_muon);

    /**
     * @brief Binding for a single particle pion multiplicity cut.
     * @details This function binds the single particle multiplicity cut for
     * charged pions, which corresponds to the index 3 in the
     * @ref utilities::count_primaries function.
     * @param obj the interaction to select on.
     * @param params the parameters for the cut. In this case, this sets the
     * kinetic energy threshold for a pion to count towards the multiplicity.
     * Defaults to 25 MeV.
     * @return true if the interaction has a single primary charged pion.
     */
    template<class T>
    bool single_pion(const T & obj, std::vector<double> params={25.0,})
    {
        return particle_multiplicity(obj, 1, 3, params) == 1;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, single_pion, single_pion);

    /**
     * @brief Binding for a single particle proton multiplicity cut.
     * @details This function binds the single particle multiplicity cut for
     * protons, which corresponds to the index 4 in the
     * @ref utilities::count_primaries function.
     * @param obj the interaction to select on.
     * @param params the parameters for the cut. In this case, this sets the
     * kinetic energy threshold for a proton to count towards the multiplicity.
     * Defaults to 50 MeV.
     * @return true if the interaction has a single primary proton.
     */
    template<class T>
    bool single_proton(const T & obj, std::vector<double> params={50.0,})
    {
        return particle_multiplicity(obj, 1, 4, params) == 1;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, single_proton, single_proton);

    /**
     * @brief Binding for zero particle photon multiplicity cut (negation of
     * nonzero_particle_multiplicity).
     * @details This function binds the nonzero particle multiplicity cut for
     * photons, which corresponds to the index 0 in the
     * @ref utilities::count_primaries function. The negation of this
     * function is used to select interactions with no primary photons.
     * @param obj the interaction to select on.
     * @param params the parameters for the cut. In this case, this sets the
     * kinetic energy threshold for a photon to count towards the
     * multiplicity. Defaults to 25 MeV.
     * @return true if the interaction has a nonzero primary photon.
     */
    template<class T>
    bool no_photons(const T & obj, std::vector<double> params={25.0,})
    {
        return particle_multiplicity(obj, 0, 0, params) == 0;
    }

    REGISTER_CUT_SCOPE(RegistrationScope::Both, no_photons, no_photons);

    /**
     * @brief Binding for zero particle electron multiplicity cut (negation of
     * nonzero_particle_multiplicity).
     * @details This function binds the nonzero particle multiplicity cut for
     * electrons, which corresponds to the index 1 in the
     * @ref utilities::count_primaries function. The negation of this
     * function is used to select interactions with no primary electrons.
     * @param obj the interaction to select on.
     * @param params the parameters for the cut. In this case, this sets the
     * kinetic energy threshold for an electron to count towards the
     * multiplicity. Defaults to 25 MeV.
     * @return true if the interaction has a nonzero primary electron.
     */
    template<class T>
    bool no_electrons(const T & obj, std::vector<double> params={25.0,})
    {
        return particle_multiplicity(obj, 0, 1, params) == 0;
    }

    REGISTER_CUT_SCOPE(RegistrationScope::Both, no_electrons, no_electrons);
    /**
     * @brief Binding for zero particle muon multiplicity cut (negation of
     * nonzero_particle_multiplicity).
     * @details This function binds the nonzero particle multiplicity cut for
     * muons, which corresponds to the index 2 in the
     * @ref utilities::count_primaries function. The negation of this
     * function is used to select interactions with no primary muons.
     * @param obj the interaction to select on.
     * @param params the parameters for the cut. In this case, this sets the
     * kinetic energy threshold for a muon to count towards the
     * multiplicity. Defaults to 143.425 MeV, which corresponds to a muon of
     * length 50 cm (assuming the muon stops).
     * @return true if the interaction has a nonzero primary muon.
     */

    template<class T>
    bool no_muons(const T & obj, std::vector<double> params={143.425,})
    {
        return particle_multiplicity(obj, 0, 2, params) == 0;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, no_muons, no_muons);

    /**
     * @brief Binding for zero particle pion multiplicity cut (negation of
     * nonzero_particle_multiplicity).
     * @details This function binds the nonzero particle multiplicity cut for
     * charged pions, which corresponds to the index 3 in the
     * @ref utilities::count_primaries function. The negation of this
     * function is used to select interactions with no primary charged pions.
     * @param obj the interaction to select on.
     * @param params the parameters for the cut. In this case, this sets the
     * kinetic energy threshold for a pion to count towards the
     * multiplicity. Defaults to 25 MeV.
     * @return true if the interaction has a nonzero primary charged pion.
     */
    template<class T>
    bool no_charged_pions(const T & obj, std::vector<double> params={25.0,})
    {
        return particle_multiplicity(obj, 0, 3, params) == 0;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, no_charged_pions, no_charged_pions);

    /**
     * @brief Binding for zero particle proton multiplicity cut (negation of
     * nonzero_particle_multiplicity).
     * @details This function binds the nonzero particle multiplicity cut for
     * protons, which corresponds to the index 4 in the
     * @ref utilities::count_primaries function. The negation of this
     * function is used to select interactions with no primary protons.
     * @param obj the interaction to select on.
     * @param params the parameters for the cut. In this case, this sets the
     * kinetic energy threshold for a proton to count towards the
     * multiplicity. Defaults to 50 MeV.
     * @return true if the interaction has a nonzero primary proton.
     */
    template<class T>
    bool no_protons(const T & obj, std::vector<double> params={50.0,})
    {
        return particle_multiplicity(obj, 0, 4, params) == 0;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, no_protons, no_protons);

    /**
     * @brief Cut to select interactions with more than one proton.
     * @details This function applies a cut to select interactions with
     * more than one proton (N > 1). This is complementary to the single_proton
     * cut.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction has more than one proton.
     */
    template<class T>
    bool multiproton(const T & obj, std::vector<double> params={50.0,})
    {
        return particle_multiplicity(obj, 1, 4, params) > 1;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, multiproton, multiproton);

    /**
     * @brief Cut to select interactions with a single Michel electron.
     * @details This function applies a cut to select interactions with a
     * single Michel electron.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @param params the parameters for the cut. In this case, this sets the
     * number of depositions for a Michel to count towards the multiplicity. 
     * Defaults to 10.
     * @return true if the interaction has a single Michel electron.
     */
    template<class T>
    bool single_michel(const T & obj, std::vector<double> params={10.0,})
    {
        size_t count(0);
        for(const auto & p : obj.particles)
        {
            if(pvars::semantic_type(p) == 2 && p.size > params[0])
                ++count;
            if(count > 1)
                break; // No need to count further.
        }
        return count == 1;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, single_michel, single_michel);





//
//    template<class T>
//        bool pion_michel_tag_cut(const T& obj, std::vector<double> params={50.0,20.0})
//        {
//            bool findpion=false;
//            bool michel_tagged = false;
//            size_t pionindex= -99;//gaurantee it is reco pion candidate
//            utilities::three_vector pion_end = {0,0,0};
//            for(size_t i(0); i < obj.particles.size(); ++i)
//            {
//                const auto & p = obj.particles[i];
//                if(pvars::pid(p) == pvars::kPion && pvars::primary_classification(p) && pvars::ke(p) >= params[0]) //same condition of single_pion
//                {
//                    pionindex=p.id;
//                    auto & pion(obj.particles[i]);
//                    pion_end = {pvars::end_x(pion), pvars::end_y(pion), pvars::end_z(pion)};
//                    findpion=true;
//                }
//            }
//            if (findpion)
//            {
//                for(size_t i(0); i < obj.particles.size(); ++i)
//                {
//                    const auto & p = obj.particles[i];
//                    if (p.id == pionindex) continue;
//                    utilities::three_vector particle_vtx = {pvars::start_x(p), pvars::start_y(p), pvars::start_z(p)};
//                    double Atslc = utilities::magnitude(utilities::subtract(pion_end, particle_vtx));
//                    if (p.shape==2 && Atslc<params[1])  
//                    {
//                        michel_tagged = true;
//                    }
//
//                }
//            }
//            else
//            {
//                return michel_tagged;
//            }
//            return michel_tagged;
//        }
//    REGISTER_CUT_SCOPE(RegistrationScope::Both, pion_michel_tag_cut, pion_michel_tag_cut);
    template<class T>
        bool pion_michel_tag_cut(const T& obj, std::vector<double> params={20.0})
        {
            size_t pion_index = selectors::leading_pion(obj);
            if (pion_index == kNoMatch) return false;
            const auto& pion = obj.particles[pion_index];
            // KE threshold: default 50 MeV
            utilities::three_vector pion_end = {
                pvars::end_x(pion),
                pvars::end_y(pion),
                pvars::end_z(pion)
            };
            int64_t pion_id = pion.id;
            bool michel_tagged = false;
            for(size_t i = 0; i < obj.particles.size(); ++i)
            {
                const auto& p = obj.particles[i];
                if (p.id == pion_id) continue;
                utilities::three_vector particle_vtx = {
                    pvars::start_x(p),
                    pvars::start_y(p),
                    pvars::start_z(p)
                };
                double Atslc =utilities::magnitude(utilities::subtract(pion_end, particle_vtx));
                // Michel candidate within 20 cm
                if (p.shape == 2 && Atslc < params[0])
                {
                    michel_tagged = true;
                    break;
                }
            }
            return michel_tagged;
        }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, pion_michel_tag_cut, pion_michel_tag_cut);



    

    template <class T>
        std::vector<int64_t> get_pion_children_pdgs(const T& obj, std::vector<double> params={25.0,})
        {
            std::vector<int64_t> v;
            if (params.empty()) return v;

            size_t pi_idx = kNoMatch;

            for (size_t i = 0; i < obj.particles.size(); ++i) {
                const auto& p = obj.particles[i];
                if (pvars::pid(p) == pvars::kPion && pvars::primary_classification(p) && pvars::ke(p) >= params[0]) {
                    pi_idx = i;                            // <-- keep the INDEX
                }
            }
            if (pi_idx == kNoMatch) return v;

            const auto& pion = obj.particles[pi_idx];

            std::unordered_map<int64_t, size_t> id2i;
            id2i.reserve(obj.particles.size());
            for (size_t i = 0; i < obj.particles.size(); ++i)
                id2i.emplace((int64_t)obj.particles[i].id, i);

            for (auto cid_raw : pion.children_id) {
                auto it = id2i.find((int64_t)cid_raw);
                if (it == id2i.end()) continue;
                v.push_back((int64_t)obj.particles[it->second].pdg_code);
            }
            return v;
        }
    template <class T>
        static inline std::unordered_map<int64_t, std::vector<size_t>>
        build_parent_children_map(const T& obj)
        {
            std::unordered_map<int64_t, std::vector<size_t>> pc;
            pc.reserve(obj.particles.size());
            for (size_t i = 0; i < obj.particles.size(); ++i)
            {
                if (obj.particles[i].parent_id == obj.particles[i].id) continue; // <-- skip self-parented entries
                pc[obj.particles[i].parent_id].push_back(i);
            }
            return pc;
        }

    template <class T>
        struct DescendantFlags {
            bool hasMuon        = false;
            bool hasMichelStrict= false; // (shape==2 && |pdg|==11)
            bool hasChargedPi   = false; // |pdg|==211
            bool hasPi0         = false; // 111
            bool hasNuclear     = false; // p/n or PDG >= 1e9
        };

    template <class T>
        static inline DescendantFlags<T>
        scan_direct_children(const T& obj, size_t pion_idx)
        {
            DescendantFlags<T> f;  // same fields: hasMuon, hasMichelStrict, hasChargedPi, hasPi0, hasNuclear
            if (pion_idx == kNoMatch || pion_idx >= obj.particles.size()) return f;

            const auto pc = build_parent_children_map(obj);
            const int64_t root = obj.particles[pion_idx].id;

            auto it = pc.find(root);
            if (it == pc.end()) return f; // no direct children

            for (size_t idx : it->second) {
                const auto& ch = obj.particles[idx];
                if (ch.id == root) continue; // <-- ignore self as a "child"
                const int pdg = ch.pdg_code;
                const int ap  = std::abs(pdg);

                if (ap == 13)                     f.hasMuon = true;
                if (ch.shape == 2 && ap == 11)    f.hasMichelStrict = true; // Michel-like electron as a *direct* daughter only
                if (ap == 211)                    f.hasChargedPi = true;
                if (pdg == 111)                   f.hasPi0 = true;
                if (ap == 2212 || ap == 2112 || ap >= 1000000000) f.hasNuclear = true; // p/n or ion
            }
            return f;
        }


    template <class T>
        bool pion_michel_decay(const T& obj, std::vector<double> params={25.0,})
        {
            if (params.empty()) return false;

            size_t pi_idx = kNoMatch;
            for (size_t i = 0; i < obj.particles.size(); ++i) {
                const auto& p = obj.particles[i];
                if (pvars::pid(p) == pvars::kPion && pvars::primary_classification(p) && pvars::ke(p) >= params[0]) {
                    pi_idx = i;                               // <-- keep INDEX
                }
            }
            if (pi_idx == kNoMatch) return false;
            const auto ch = get_pion_children_pdgs(obj);
            auto cntN = [&](){
                int c = 0;
                for (auto x: ch) if (x==2212 || x==2112 || std::abs(x)>=1000000000) ++c;
                return c;
            };
            auto onlyGammas = [&](){
                if (ch.empty()) return false;
                for (auto x: ch) if (x != 22) return false;
                return true;
            };

            const int  Nnuc_dir = cntN();
            const bool anyN_dir = (Nnuc_dir > 0);

            // *** Direct-only scan ***
            const auto df = scan_direct_children(obj, pi_idx);

            // -------- DECAY (direct-only) ----------
            // Require a direct muon to call DECAY.
            // (Optional) Allow a direct Michel-shaped e± ONLY if there are no direct pions/pi0/nuclear products.
            if (df.hasMuon)
                return true;

            if (df.hasMichelStrict && !df.hasChargedPi && !df.hasPi0 && !anyN_dir)
                return true;

            // -------- CAPTURE (direct-only) --------
            // π- capture often shows photons only, or nucleons with no pions.
            if (onlyGammas())
                return false;

            if (anyN_dir && !df.hasChargedPi && !df.hasPi0)
                return false;

            // -------- INELASTIC / ELASTIC ----------
            // π0 among direct daughters → inelastic
            if (df.hasPi0)
                return false;

            // Charged π with nucleons → inelastic
            if (df.hasChargedPi && anyN_dir)
                return false;

            // Exactly one direct charged π and nothing else → elastic-like
            {
                int nCpi = 0, nOther = 0;
                for (auto x: ch) {
                    if (std::abs(x) == 211) ++nCpi;
                    else if (x != 0) ++nOther; // ignore 0 if it can appear
                }
                if (nCpi == 1 && nOther == 0)
                    return false; // or INELASTIC if you don't keep ELASTIC separate
            }

            // Any charged π at all (and no stronger signature above) → inelastic bucket
            if (df.hasChargedPi)
                return false;

            // -------- Fallback ----------
            return false;
        }
    REGISTER_CUT_SCOPE(RegistrationScope::True, pion_michel_decay, pion_michel_decay);


    template<class T>
    bool multiple_pion_michel_tag_cut(const T& obj, std::vector<double> params={25.0,25.0})
    {
        // params[0] = pion KE threshold
        // params[1] = distance threshold for Michel proximity
    
        std::vector<size_t> pion_ids;
        std::vector<utilities::three_vector> pion_ends;
    
        // 1) collect ALL reco pion candidates (same condition as your single_pion)
        for (size_t i = 0; i < obj.particles.size(); ++i)
        {
            const auto& p = obj.particles[i];
    
            if (pvars::pid(p) == pvars::kPion && pvars::primary_classification(p) && pvars::ke(p) >= params[0])
            {
                pion_ids.push_back(p.id);
                pion_ends.push_back({pvars::end_x(p), pvars::end_y(p), pvars::end_z(p)});
            }
        }
    
        const int numpion = static_cast<int>(pion_ids.size());
        if (numpion <= 1) return false; // require multiple pions
    
        // 2) count Michel tags (how many particles satisfy your Michel condition
        //    AND are close to ANY pion end)
        int nmichel = 0;
    
        for (size_t i = 0; i < obj.particles.size(); ++i)
        {
            const auto& p = obj.particles[i];
    
            // skip the pion candidates themselves
            if (std::find(pion_ids.begin(), pion_ids.end(), p.id) != pion_ids.end())
                continue;
    
            // your Michel-like condition
            if (p.shape != 2) continue;
    
            utilities::three_vector particle_vtx = {pvars::start_x(p), pvars::start_y(p), pvars::start_z(p)};
    
            bool close_to_any_pion = false;
            for (const auto& pe : pion_ends)
            {
                const double d = utilities::magnitude(utilities::subtract(pe, particle_vtx));
                if (d < params[1]) { close_to_any_pion = true; break; }
            }
    
            if (close_to_any_pion) nmichel++;
        }
    
        // require multiple Michel tags too
        return (nmichel > 1);
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, multiple_pion_michel_tag_cut, multiple_pion_michel_tag_cut);
    
    /**
     * @brief Cut to select interactions with a Michel electron attached to the
     * end of a muon track.
     * @details This function applies a cut to select interactions with a
     * single Michel electron that has a start point within some distance 
     * threshold from the endpoints of a muon track. Checks both endpoint and
     * startpoint, in case of track flipping.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @param params the parameters for the cut. In this case, this sets the
     * distance between the Michel and muon, the Michel's min # of depositions, 
     * and the muon KE threshold.
     * @return true if the interaction has a Michel electron attached to the
     * "end" of the selected muon track.
    */
    template<class T>
    bool michel_attached_muon(const T & obj, std::vector<double> params={})
    {
        // Check that the parameters are given; if not, apply default values.
        if(params.size() < 3)
        {
            params.resize(3);
            params[0] = 10.0;    // Minimum number of depositions for the Michel electron
            params[1] = 143.425; // Muon KE threshold (corresponds to 50 cm track length)
            params[2] = 10.0;    // Distance threshold between Michel and muon start/end points
        }

        for(const auto & p : obj.particles)
        {
            if(pvars::semantic_type(p) != 2 || p.size < params[0])
                continue; // Not target Michel

            for(const auto & p2 : obj.particles)
            {
                if(pvars::pid(p2) != 2 || !pvars::primary_classification(p2) || pvars::ke(p2) < params[1])
                    continue; // Not target muon

                float dx = pvars::start_x(p) - pvars::end_x(p2);
                float dy = pvars::start_y(p) - pvars::end_y(p2);
                float dz = pvars::start_z(p) - pvars::end_z(p2);
                float dist = std::sqrt(dx * dx + dy * dy + dz * dz);

                float dx_flip = pvars::start_x(p) - pvars::start_x(p2);
                float dy_flip = pvars::start_y(p) - pvars::start_y(p2);
                float dz_flip = pvars::start_z(p) - pvars::start_z(p2);
                float dist_flip = std::sqrt(dx_flip * dx_flip + dy_flip * dy_flip + dz_flip * dz_flip);

                if(dist < params[2] || dist_flip < params[2])
                    return true; // Michel is attached to muon
            }
        }

        return false;
    }

    REGISTER_CUT_SCOPE(RegistrationScope::Reco, michel_attached_muon, michel_attached_muon);    

    /**
    * @brief Apply a cut on particle multiplicity with configurable species and count.
    * @details This function applies a cut based on the multiplicity of a specific
    * particle species. The particle species and desired multiplicity are specified
    * via the parameters vector. This allows for flexible multiplicity cuts without
    * needing to define separate functions for each combination.
    * @tparam T the type of interaction (true or reco).
    * @param obj the interaction to select on.
    * @param params the parameters for the cut:
    *   - params[0]: minimum desired multiplicity (converted to size_t)
    *   - params[1]: particle species index (0=photon, 1=electron, 2=muon, 3=pion, 4=proton)
    *   - params[2]: kinetic energy threshold in MeV (optional, defaults to 0.0)
    * @return true if the interaction has exactly the specified multiplicity of the
    * specified particle species above the energy threshold.
    */
    template<class T>
    bool has_particle_multiplicity(const T & obj, std::vector<double> params={1.0, 2.0, 0.0})
    {
        if(params.size() < 2)
        {
            throw std::runtime_error("has_particle_multiplicity requires at least 2 parameters: [multiplicity, particle_species, (optional) ke_threshold]");
        }
        
        size_t desired_mult = static_cast<size_t>(params[0]);
        size_t particle_species = static_cast<size_t>(params[1]);
        double ke_threshold = (params.size() >= 3) ? params[2] : 0.0;
        
        std::vector<double> threshold_params = {ke_threshold};
        size_t actual_mult = particle_multiplicity(obj, 1, particle_species, threshold_params);
        
        return actual_mult >= desired_mult;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, has_particle_multiplicity, has_particle_multiplicity);

    /**
     * @brief Apply a cut on the number of identified primary tracks.
     * @details Counts primary particles with a semantic type of 1 (track).
     * This is a PID-species-agnostic analogue of @ref has_particle_multiplicity,
     * intended for pre-selections that study the impact of PID-based cuts
     * downstream. No requirement is made on a particle's proximity to the
     * interaction vertex.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @param params the parameters for the cut:
     *   - params[0]: minimum desired track count (default 2.0)
     * @return true if the interaction has at least params[0] identified
     * primary tracks.
     */
    template<class T>
    bool track_multiplicity(const T & obj, std::vector<double> params={2.0})
    {
        size_t desired_mult = static_cast<size_t>(params[0]);

        size_t count(0);
        for(size_t i(0); i < obj.particles.size(); ++i)
        {
            const auto & p = obj.particles[i];
            if(pvars::semantic_type(p) != 1 || !pvars::primary_classification(p))
                continue;
            ++count;
        }
        return count >= desired_mult;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, track_multiplicity, track_multiplicity);

    /**
     * @brief Cut to select interactions with an identified longest track.
     * @details Checks that @ref selectors::longest_track finds a valid
     * particle (a particle with a semantic type of 1, i.e. a track).
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction has an identified longest track.
     */
    template<class T>
    bool is_longest_track(const T & obj)
    {
        return selectors::longest_track(obj) != kNoMatch;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, is_longest_track, is_longest_track);

    /**
     * @brief Cut to select interactions with an identified second longest track.
     * @details Checks that @ref selectors::second_longest_track finds a
     * valid particle (a particle with a semantic type of 1, i.e. a track).
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction has an identified second longest track.
     */
    template<class T>
    bool is_second_longest_track(const T & obj)
    {
        return selectors::second_longest_track(obj) != kNoMatch;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, is_second_longest_track, is_second_longest_track);

    /**
     * @brief Cut to select interactions with an identified third longest track.
     * @details Checks that @ref selectors::third_longest_track finds a
     * valid particle (a particle with a semantic type of 1, i.e. a track).
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction has an identified third longest track.
     */
    template<class T>
    bool is_third_longest_track(const T & obj)
    {
        return selectors::third_longest_track(obj) != kNoMatch;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, is_third_longest_track, is_third_longest_track);

    /**
     * @brief Cut to select longest track length below a threshold.
     * @details This function applies a cut to select interactions with a
     * longest track length below a specified threshold.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction's longest track length is below the threshold.
    */
    template<class T>
    bool track_length(const T & obj, std::vector<double> params={25.0})
    {

        if(params.size() < 1)
        {
            throw std::runtime_error("track_length requires at least 1 parameter: [threshold]");
        }

        float maxLen = -1;
        size_t longest_track_idx = kNoMatch;
        for (size_t i = 0; i < obj.particles.size(); ++i) {
            const auto & p = obj.particles[i];
            if (pvars::semantic_type(p) != 1) continue;
            if (pvars::length(p) > maxLen) {
                maxLen = pvars::length(p);
                longest_track_idx = i;
            }
        }
    
        if(longest_track_idx == kNoMatch) return true;
        
        return pvars::length(obj.particles[longest_track_idx]) < params[0];
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, track_length, track_length);

    /**
     * @brief Cut to select interactions with leading electron dE/dx below a threshold.
     * @details This function applies a cut to select interactions with a
     * leading electron dE/dx below a specified threshold.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction has a leading electron with dE/dx below the threshold.
    */
    template<class T>
    bool particle_dedx(const T & obj, std::vector<double> params={0.0,})
    {   
        size_t i = selectors::leading_electron(obj);
        if (i == kNoMatch) return false;
        const auto & p = obj.particles[i];
        bool pass = (double)p.start_dedx < params[0] && (double)p.start_dedx >= 0;
        return pass;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Reco, particle_dedx, particle_dedx);

    /**
     * @brief Cut to select interactions with a leading electron vertex distance below a threshold.
     * @details This function applies a cut to select interactions with a
     * leading electron vertex distance below a specified threshold.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to select on.
     * @return true if the interaction has a leading electron with vertex distance below the threshold.
    */
    template<class T>
    bool vertex_distance_cut(const T & obj, std::vector<double> params={0.0,})
    {   
        size_t i = selectors::leading_electron(obj);
        if (i == kNoMatch) return false;
        const auto & p = obj.particles[i];
        return p.vertex_distance >= 0 ? p.vertex_distance < params[0] : false;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Reco, vertex_distance_cut, vertex_distance_cut);


}
#endif
