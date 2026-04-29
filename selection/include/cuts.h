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
     * @return true if the interaction has a single Michel electron.
    */
    template<class T>
    bool single_michel(const T & obj)
    {
        size_t count(0);
        for(const auto & p : obj.particles)
        {
            if(pvars::semantic_type(p) == 2)
                ++count;
            if(count > 1)
                break; // No need to count further, we only care about multiplicity of 1.
        }
        return count == 1;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::Both, single_michel, single_michel);






    template<class T>
        bool pion_michel_tag_cut(const T& obj, std::vector<double> params={25.0,20.0})
        {
            bool findpion=false;
            bool michel_tagged = false;
            size_t pionindex= -99;//gaurantee it is reco pion candidate
            utilities::three_vector pion_end = {0,0,0};
            for(size_t i(0); i < obj.particles.size(); ++i)
            {
                const auto & p = obj.particles[i];
                if(pvars::pid(p) == pvars::kPion && pvars::primary_classification(p) && pvars::ke(p) >= params[0]) //same condition of single_pion
                {
                    pionindex=p.id;
                    auto & pion(obj.particles[i]);
                    pion_end = {pvars::end_x(pion), pvars::end_y(pion), pvars::end_z(pion)};
                    findpion=true;
                }
            }
            if (findpion)
            {
                for(size_t i(0); i < obj.particles.size(); ++i)
                {
                    const auto & p = obj.particles[i];
                    if (p.id == pionindex) continue;
                    utilities::three_vector particle_vtx = {pvars::start_x(p), pvars::start_y(p), pvars::start_z(p)};
                    double Atslc = utilities::magnitude(utilities::subtract(pion_end, particle_vtx));
                    if (p.shape==2 && Atslc<params[1])  
                    {
                        michel_tagged = true;
                    }

                }
            }
            else
            {
                return michel_tagged;
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

}
#endif
