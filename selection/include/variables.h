/**
 * @file variables.h
 * @brief Header file for definitions of analysis variables.
 *
 * @details This file contains definitions of analysis variables which can be
 * used to extract information from interactions. Each variable is implemented
 * as a function which takes an interaction object as an argument and returns a
 * double. These are the building blocks for producing high-level plots of the
 * selected interactions.
 * @author mueller@fnal.gov
*/
#ifndef VARIABLES_H
#define VARIABLES_H
#define ELECTRON_MASS 0.5109989461
#define MUON_MASS 105.6583745
#define PION_MASS 139.57039
#define PROTON_MASS 938.2720813
#define NEUTRON_MASS 939.56542194

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "sbnanaobj/StandardRecord/SRInteractionDLP.h"
#include "sbnanaobj/StandardRecord/SRInteractionTruthDLP.h"
#include "sbnanaobj/StandardRecord/Proxy/EpilogFwd.h"

#include "include/particle_variables.h"
#include "include/particle_cuts.h"
#include "include/cuts.h"
#include "include/utilities.h"
#include "include/particle_utilities.h"
#include "include/selectors.h"
#include "framework.h"


#include <string>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <optional>
#include <cmath>
#include <limits>
#include <iostream>
#include <cstdlib>   // for std::abs
#include <random>
#include <set>



/**
 * @namespace vars
 * @brief Namespace for organizing generic variables which act on interactions.
 * @details This namespace is intended to be used for organizing variables which
 * act on interactions. Each variable is implemented as a function which takes
 * an interaction object as an argument and returns a double. The function
 * should be templated on the type of interaction object if the variable is
 * intended to be used on both true and reconstructed interactions.
 * @note The namespace is intended to be used in conjunction with the
 * pvars namespace, which is used for organizing variables which act on single
 * particles.
 */
namespace vars
{
    /**
     * @brief Variable for the neutrino ID of the interaction.
     * @details This variable is intended to provide a unique identifier for
     * each parent neutrino within the event record. This number is assigned
     * starting at 0 for the first neutrino in the event and is incremented
     * for each subsequent neutrino. Non-neutrino interactions are assigned
     * a value of -1.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the neutrino ID.
     */
    template<class T>
    double neutrino_id(const T & obj) { return obj.nu_id; }
    REGISTER_VAR_SCOPE(RegistrationScope::True, neutrino_id, neutrino_id);

    template<class T>
    double neutrino_pdg(const T & obj) { return obj.pdg_code; }
    REGISTER_VAR_SCOPE(RegistrationScope::True, neutrino_pdg, neutrino_pdg);

    /**
     * @brief Variable for the interaction ID.
     * @details This variable is intended to provide a unique identifier for
     * each interaction within the event record. This number is assigned
     * starting at 0 for the first interaction in the event and is incremented
     * for each subsequent interaction. This assignment is done upstream in the
     * SPINE reconstruction.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the interaction ID.
     */
    template<class T>
    double interaction_id(const T & obj) { return obj.id; }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, interaction_id, interaction_id);

    /**
     * @brief Variable for the best-match IoU of the interaction.
     * @details The best-match IoU is the intersection over union of the
     * points belonging to a pair of reconstructed and true interactions. The
     * best-match IoU is calculated upstream in the SPINE reconstruction.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the best-match IoU of the interaction.
     */
    template<class T>
    double iou(const T & obj)
    {
        if(obj.match_ids.size() > 0)
            return obj.match_overlaps[0];
        else 
            return PLACEHOLDERVALUE;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, iou, iou);

    /**
     * @brief Variable for the containment status of the interaction.
     * @details The containment status is determined upstream in the SPINE
     * reconstruction and is based on the set of all points in the interaction,
     * which must be contained within the volume of the TPC that created them.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the containment status of the interaction.
     */
    template<class T>
    double containment(const T & obj) { return cuts::containment_cut(obj); }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, containment, containment);

    /**
     * @brief Variable for the fiducial volume status of the interaction.
     * @details The fiducial volume status is determined upstream in the SPINE
     * reconstruction and is a requirement that the interaction vertex is within
     * the fiducial volume of the TPC.
     */
    template<class T>
    double fiducial(const T & obj) { return cuts::fiducial_cut(obj); }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, fiducial, fiducial);

    /**
     * @brief Variable for total visible energy of interaction.
     * @details This function calculates the total visible energy of the
     * interaction by summing the energy of all particles that are identified
     * as counting towards the final state of the interaction.
     * @tparam T the type of interaction (true or reco).
     * @param obj interaction to apply the variable on.
     * @return the total visible energy of the interaction.
     */
    template<class T>
    double visible_energy(const T & obj)
    {
        double energy(0);
        for(const auto & p : obj.particles)
        {
            if(pcuts::final_state_signal(p))
            {
                energy += pvars::energy(p);
                if(pvars::pid(p) == pvars::kProton) energy -= pvars::mass(p) - PROTON_BINDING_ENERGY;
            }
        }
        return energy/1000.0;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, visible_energy, visible_energy);

    /**
     * @brief Variable for energy reconstruction assuming CCQE kinematics using
     * the lepton.
     * @details This function calculates the neutrino energy assuming CCQE
     * kinematics using the leading lepton in the interaction. The leading
     * lepton is defined as the highest kinetic energy electron or muon in the
     * interaction. If no electron or muon is found, the function returns
     * PLACEHOLDERVALUE. This does not check that the interaction is actually
     * QE-like.
     * @tparam T the type of interaction (true or reco).
     * @param obj interaction to apply the variable on.
     * @return the reconstructed neutrino energy assuming CCQE kinematics.
     */
    template<class T>
    double energy_qel(const T & obj)
    {
        size_t ei = selectors::leading_electron(obj);
        size_t mi = selectors::leading_muon(obj);
        size_t li = kNoMatch;

        if(ei != kNoMatch && mi != kNoMatch)
            li = (pvars::ke(obj.particles[ei]) > pvars::ke(obj.particles[mi])) ? ei : mi;
        else if(ei != kNoMatch)
            li = ei;
        else if(mi != kNoMatch)
            li = mi;
        else
            return PLACEHOLDERVALUE;

        double Mn = 939.565;
        double Mp = 938.272;
        double Ml = (li == ei) ? ELECTRON_MASS : MUON_MASS;
        double EB = PROTON_BINDING_ENERGY;

        double El = pvars::energy(obj.particles[li]);
        double pz = pvars::pz(obj.particles[li]);

        double numerator   = 2*(Mn - EB)*El - ((Mn - EB)*(Mn - EB) + Ml*Ml - Mp*Mp);
        double denominator = 2*((Mn - EB) - El + pz);

        return (numerator / denominator) / 1000.0;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, energy_qel, energy_qel);

    /**
     * @brief Variable for energy reconstruction assuming CCQE kinematics using
     * the proton.
     * @details This function calculates the neutrino energy assuming CCQE
     * kinematics using the leading proton in the interaction. The leading
     * proton is defined as the highest kinetic energy proton in the
     * interaction. If no proton is found, the function returns 
     * PLACEHOLDERVALUE. This does not check that the interaction is actually
     * QE-like.
     * @tparam T the type of interaction (true or reco).
     * @param obj interaction to apply the variable on.
     * @return the reconstructed neutrino energy assuming CCQE kinematics.
     */
    template<class T>
    double energy_qep(const T & obj)
    {
        size_t ei = selectors::leading_electron(obj);
        size_t mi = selectors::leading_muon(obj);
        size_t pi = selectors::leading_proton(obj);
        size_t li = kNoMatch;

        if(ei != kNoMatch && mi != kNoMatch)
            li = (pvars::ke(obj.particles[ei]) > pvars::ke(obj.particles[mi])) ? ei : mi;
        else if(ei != kNoMatch)
            li = ei;
        else if(mi != kNoMatch)
            li = mi;
        else
            return PLACEHOLDERVALUE;

        if(pi == kNoMatch)
            return PLACEHOLDERVALUE;

        double Mn = 939.565;
        double Mp = 938.272;
        double Ml = (li == ei) ? ELECTRON_MASS : MUON_MASS;
        double EB = PROTON_BINDING_ENERGY;

        double Ep = pvars::energy(obj.particles[pi]);
        double pz = pvars::pz(obj.particles[pi]);

        double numerator   = 2*(Mn - EB)*Ep - ((Mn - EB)*(Mn - EB) + Mp*Mp - Ml*Ml);
        double denominator = 2*((Mn - EB) - Ep + pz);

        return (numerator / denominator) / 1000.0;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, energy_qep, energy_qep);

    /**
     * @brief Variable for total visible energy of interaction, including
     * sub-threshold particles.
     * @details This function calculates the total visible energy of the
     * interaction by summing the energy of all particles that are identified
     * as counting towards the final state of the interaction. Sub-threshold
     * particles are included calorimetrically.
     * @tparam T the type of interaction (true or reco).
     * @param obj interaction to apply the variable on.
     * @return the total visible energy of the interaction.
     */
    template<class T>
    double visible_energy_calosub(const T & obj)
    {
        double energy(0);
        for(const auto & p : obj.particles)
        {
            if(pcuts::final_state_signal(p))
            {
                energy += pvars::energy(p);
                if(pvars::pid(p) == pvars::kProton) energy -= PROTON_MASS - PROTON_BINDING_ENERGY;
            }
            else if(pcuts::is_primary(p))
                energy += p.calo_ke;
        }
        return energy/1000.0;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, visible_energy_calosub, visible_energy_calosub);

    /**
     * @brief Variable for the flash time of the interaction.
     * @details The flash time is the time of the flash observed in the PMTs
     * and associated with the charge deposition in the interaction using the
     * OpT0Finder likelihood method.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the flash time of the interaction.
     */
    template<class T>
    double flash_time(const T & obj)
    {
        if(obj.flash_times.size() > 0)
            return obj.flash_times[0];
        return PLACEHOLDERVALUE;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Reco, flash_time, flash_time);

    template<class T>
    double flash_time_test(const T & obj)
    {
        if(obj.flash_times.size() > 0)
            return obj.flash_times[0];
        return PLACEHOLDERVALUE;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::True, flash_time_test, flash_time_test);
    /**
     * @brief Variable for the flash score of the interaction.
     * @details The flash score is the likelihood score of the flash observed
     * in the PMTs and associated with the charge deposition in the interaction
     * by the OpT0Finder likelihood method.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the flash score of the interaction.
     */
    template<class T>
    double flash_score(const T & obj)
    {
        if(obj.flash_scores.size() > 0)
            return obj.flash_scores[0];
        return PLACEHOLDERVALUE;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Reco, flash_score, flash_score);

    /**
     * @brief Variable for the flash total photoelectron count of the
     * interaction.
     * @details The flash total photoelectron count is the total number of
     * photoelectrons observed in the PMTs and associated with the charge
     * deposition in the interaction using the OpT0Finder likelihood method.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the flash total photoelectron count of the interaction.
     */
    template<class T>
    double flash_total_pe(const T & obj) { return obj.flash_total_pe; }
    REGISTER_VAR_SCOPE(RegistrationScope::Reco, flash_total_pe, flash_total_pe);

    /**
     * @brief Variable for the flash hypothesis total photoelectron count of
     * the interaction.
     * @details The flash hypothesis total photoelectron count is the total
     * number of photoelectrons predicted by OpT0Finder for the interaction in
     * the flash associated with the charge deposition in the interaction.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the flash hypothesis total photoelectron count of the interaction.
     */
    template<class T>
    double flash_hypothesis(const T & obj) { return obj.flash_hypo_pe; }
    REGISTER_VAR_SCOPE(RegistrationScope::Reco, flash_hypothesis, flash_hypothesis);

    /**
     * @brief Variable for the x-coordinate of the interaction vertex.
     * @details The interaction vertex is 3D point in space where the neutrino
     * interacted to produce the primary particles in the interaction.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the x-coordinate of the interaction vertex.
     */
    template<class T>
    double vertex_x(const T & obj) { return obj.vertex[0]; }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, vertex_x, vertex_x);

    /**
     * @brief Variable for the y-coordinate of the interaction vertex.
     * @details The interaction vertex is 3D point in space where the neutrino
     * interacted to produce the primary particles in the interaction.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the y-coordinate of the interaction vertex.
     */
    template<class T>
    double vertex_y(const T & obj) { return obj.vertex[1]; }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, vertex_y, vertex_y);

    /**
     * @brief Variable for the z-coordinate of the interaction vertex.
     * @details The interaction vertex is 3D point in space where the neutrino
     * interacted to produce the primary particles in the interaction.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the z-coordinate of the interaction vertex.
     */
    template<class T>
    double vertex_z(const T & obj) { return obj.vertex[2]; }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, vertex_z, vertex_z);

    /**
     * @brief Variable for the transverse momentum of the interaction counting
     * only particles identified as contributing to the final state.
     * @details This function calculates the transverse momentum of the
     * interaction by summing the transverse momentum of all particles that are
     * identified as counting towards the final state of the interaction. The
     * neutrino direction is assumed to either be the BNB axis direction
     * (z-axis) or the unit vector pointing from the NuMI target to the
     * interaction vertex. See @ref utilities::transverse_momentum for details
     * on the extraction of the transverse momentum.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the transverse momentum of the primary particles.
     * @note The switch to the NuMI beam direction instead of the BNB axis is
     * applied by the definition of a preprocessor macro (BEAM_IS_NUMI).
     */
    template<class T>
    double dpT(const T & obj)
    {
        utilities::three_vector pt = {0, 0, 0};
        for(const auto & p : obj.particles)
        {
            if(pcuts::final_state_signal(p))
            {
                // Sum up the transverse momentum of all final state particles
                utilities::three_vector momentum = {pvars::px(p), pvars::py(p), pvars::pz(p)};
                utilities::three_vector vtx = {pvars::start_x(p), pvars::start_y(p), pvars::start_z(p)};
                utilities::three_vector this_pt = utilities::transverse_momentum(momentum, vtx);
                pt = utilities::add(pt, this_pt);
            }
        }
        return utilities::magnitude(pt);
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, dpT, dpT);

    /**
     * @brief Variable for the transverse momentum of the interaction counting
     * only the leading charged lepton and proton.
     * @details This function calculates the transverse momentum of the
     * interaction by summing the transverse momentum of the leading charged
     * lepton and proton. The neutrino direction is assumed to either be the
     * BNB axis direction (z-axis) or the unit vector pointing from the NuMI
     * target to the interaction vertex. See @ref utilities::transverse_momentum
     * for details on the extraction of the transverse momentum.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the transverse momentum of the leading charged lepton and proton.
     * @note The switch to the NuMI beam direction instead of the BNB axis is
     * applied by the definition of a preprocessor macro (BEAM_IS_NUMI).
     */
    template<class T>
    double dpT_lp(const T & obj)
    {
        
        utilities::three_vector l_pt = {0, 0, 0};
        utilities::three_vector p_pt = {0, 0, 0};
        double l_ke(0), p_ke(0);
        for(const auto & p : obj.particles)
        {
            if(pcuts::final_state_signal(p))
            {
                // Find the leading charged lepton and proton
                if((pvars::pid(p) == pvars::kElectron || pvars::pid(p) == pvars::kMuon) && pvars::ke(p) > l_ke)
                {
                    l_ke = pvars::ke(p);
                    utilities::three_vector momentum = {pvars::px(p), pvars::py(p), pvars::pz(p)};
                    utilities::three_vector vtx = {pvars::start_x(p), pvars::start_y(p), pvars::start_z(p)};
                    l_pt = utilities::transverse_momentum(momentum, vtx);
                }
                else if(pvars::pid(p) == pvars::kProton && pvars::ke(p) > p_ke)
                {
                    p_ke = pvars::ke(p);
                    utilities::three_vector momentum = {pvars::px(p), pvars::py(p), pvars::pz(p)};
                    utilities::three_vector vtx = {pvars::start_x(p), pvars::start_y(p), pvars::start_z(p)};
                    p_pt = utilities::transverse_momentum(momentum, vtx);
                }
            }
        }
        if(l_ke == 0 || p_ke == 0)
            return PLACEHOLDERVALUE;
        else
            return utilities::magnitude(utilities::add(l_pt, p_pt));
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, dpT_lp, dpT_lp);

    /**
     * @brief Variable for dphi_T of the interaction.
     * @details dphi_T is a transverse kinematic imbalance variable defined
     * using the transverse momentum of the leading muon and the total hadronic
     * system. This variable is sensitive to the presence of F.S.I. The
     * neutrino direction is assumed to either be the BNB axis direction
     * (z-axis) or the unit vector pointing from the NuMI target to the
     * interaction vertex. See @ref utilities::transverse_momentum for details
     * on the extraction of the transverse momentum.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the phi_T of the interaction.
     * @note The switch to the NuMI beam direction instead of the BNB axis is
     * applied by the definition of a preprocessor macro (BEAM_IS_NUMI).
     */
    template<class T>
    double dphiT(const T & obj)
    {
        utilities::three_vector lepton_pt = {0, 0, 0};
        utilities::three_vector hadronic_pt = {0, 0, 0};
        for(const auto & p : obj.particles)
        {
            if(pcuts::final_state_signal(p))
            {
                // There should only be one lepton, so replace the lepton
                // transverse momentum if the particle is a lepton.
                utilities::three_vector momentum = {pvars::px(p), pvars::py(p), pvars::pz(p)};
                utilities::three_vector vtx = {pvars::start_x(p), pvars::start_y(p), pvars::start_z(p)};
                utilities::three_vector this_pt = utilities::transverse_momentum(momentum, vtx);
                if(pvars::pid(p) == pvars::kElectron || pvars::pid(p) == pvars::kMuon)
                    lepton_pt = this_pt;
                // The total hadronic system is treated as a single object.
                else if(pvars::pid(p) > 2)
                    hadronic_pt = utilities::add(hadronic_pt, this_pt);
            }
        }
        return std::acos(-1 * utilities::dot_product(lepton_pt, hadronic_pt) / (utilities::magnitude(lepton_pt) * utilities::magnitude(hadronic_pt)));
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, dphiT, dphiT);

    /**
     * @brief Variable for dalpha_T of the interaction.
     * @details dalpha_T is a transverse kinematic imbalance variable defined
     * using the transverse momentum of the total hadronic system and the
     * outgoing lepton. The neutrino direction is assumed to either be the BNB
     * axis direction (z-axis) or the unit vector pointing from the NuMI target
     * to the interaction vertex. See @ref utilities::transverse_momentum for
     * details on the extraction of the transverse momentum.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the alpha_T of the interaction.
     * @note The switch to the NuMI beam direction instead of the BNB axis is
     * applied by the definition of a preprocessor macro (BEAM_IS_NUMI).
     */
    template<class T>
    double dalphaT(const T & obj)
    {
        utilities::three_vector lepton_pt = {0, 0, 0};
        utilities::three_vector total_pt = {0, 0, 0};
        for(const auto & p : obj.particles)
        {
            if(pcuts::final_state_signal(p))
            {
                // There should only be one lepton, so replace the lepton
                // transverse momentum if the particle is a lepton.
                utilities::three_vector momentum = {pvars::px(p), pvars::py(p), pvars::pz(p)};
                utilities::three_vector vtx = {pvars::start_x(p), pvars::start_y(p), pvars::start_z(p)};
                utilities::three_vector this_pt = utilities::transverse_momentum(momentum, vtx);
                if(pvars::pid(p) == pvars::kElectron || pvars::pid(p) == pvars::kMuon)
                    lepton_pt = this_pt;
                total_pt = utilities::add(total_pt, this_pt);
            }
        }
        return std::acos(-1 * utilities::dot_product(total_pt, lepton_pt) / (utilities::magnitude(total_pt) * utilities::magnitude(lepton_pt)));
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, dalphaT, dalphaT);

    /**
     * @brief Variable for the missing longitudinal momentum of the
     * interaction.
     * @details The missing longitudinal momentum is calculated as the
     * difference between the total longitudinal momentum of the final state
     * particles and the best estimate of the neutrino energy. The neutrino
     * energy is calculated using @ref vars::visible_energy.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the missing longitudinal momentum of the interaction.
     * @note The switch to the NuMI beam direction instead of the BNB axis is
     * applied by the definition of a preprocessor macro (BEAM_IS_NUMI).
     */
    template<class T>
    double dpL(const T & obj)
    {
        utilities::three_vector lepton_pl = {0, 0, 0};
        utilities::three_vector hadronic_pl = {0, 0, 0};
        for(const auto & p : obj.particles)
        {
            if(pcuts::final_state_signal(p))
            {
                // There should only be one lepton, so replace the lepton
                // transverse momentum if the particle is a lepton.
                utilities::three_vector momentum = {pvars::px(p), pvars::py(p), pvars::pz(p)};
                utilities::three_vector vtx = {pvars::start_x(p), pvars::start_y(p), pvars::start_z(p)};
                utilities::three_vector this_pl = utilities::longitudinal_momentum(momentum, vtx);
                if(pvars::pid(p) == pvars::kElectron || pvars::pid(p) == pvars::kMuon)
                    lepton_pl = this_pl;
                // The total hadronic system is treated as a single object.
                else if(pvars::pid(p) > 2)
                    hadronic_pl = utilities::add(hadronic_pl, this_pl);
            }
        }
        return utilities::magnitude(utilities::add(hadronic_pl, lepton_pl)) - 1000*vars::visible_energy(obj);
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, dpL, dpL);

    /**
     * @brief Variable for the missing longitudinal momentum of the interaction
     * counting only the leading charged lepton and proton.
     * @details The missing longitudinal momentum is calculated as the
     * difference between the total longitudinal momentum of the leading charged
     * lepton and proton and the best estimate of the neutrino energy. The
     * neutrino energy is calculated using @ref vars::visible_energy.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the missing longitudinal momentum of the leading charged lepton
     * and proton.
     * @note The switch to the NuMI beam direction instead of the BNB axis is
     * applied by the definition of a preprocessor macro (BEAM_IS_NUMI).
     */
    template<class T>
    double dpL_lp(const T & obj)
    {
        utilities::three_vector l_pl = {0, 0, 0};
        utilities::three_vector p_pl = {0, 0, 0};
        double l_ke(0), p_ke(0);
        for(const auto & p : obj.particles)
        {
            if(pcuts::final_state_signal(p))
            {
                // Find the leading charged lepton and proton
                if((pvars::pid(p) == pvars::kElectron || pvars::pid(p) == pvars::kMuon) && pvars::ke(p) > l_ke)
                {
                    l_ke = pvars::ke(p);
                    utilities::three_vector momentum = {pvars::px(p), pvars::py(p), pvars::pz(p)};
                    utilities::three_vector vtx = {pvars::start_x(p), pvars::start_y(p), pvars::start_z(p)};
                    l_pl = utilities::longitudinal_momentum(momentum, vtx);
                }
                else if(pvars::pid(p) == pvars::kProton && pvars::ke(p) > p_ke)
                {
                    p_ke = pvars::ke(p);
                    utilities::three_vector momentum = {pvars::px(p), pvars::py(p), pvars::pz(p)};
                    utilities::three_vector vtx = {pvars::start_x(p), pvars::start_y(p), pvars::start_z(p)};
                    p_pl = utilities::longitudinal_momentum(momentum, vtx);
                }
            }
        }
        if(l_ke == 0 || p_ke == 0)
            return PLACEHOLDERVALUE;
        else
            return utilities::magnitude(utilities::add(l_pl, p_pl)) - 1000*vars::visible_energy(obj);
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, dpL_lp, dpL_lp);

    /**
     * @brief Variable for the estimate of the momentum of the struck nucleon.
     * @details The estimate of the momentum of the struck nucleon is calculated
     * as the quadrature sum of the transverse momentum (see @ref vars::dpT) and
     * the missing longitudinal momentum (see @ref vars::dpL).
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the estimate of the momentum of the struck nucleon.
     * @note The switch to the NuMI beam direction instead of the BNB axis is
     * applied by the definition of a preprocessor macro (BEAM_IS_NUMI).
     */
    template<class T>
    double pn(const T & obj) { return std::sqrt(std::pow(vars::dpT(obj), 2) + std::pow(vars::dpL(obj), 2)); }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, pn, pn);

    /**
     * @brief Variable for the estimate of the momentum of the struck nucleon
     * counting only the leading charged lepton and proton.
     * @details The estimate of the momentum of the struck nucleon is calculated
     * as the quadrature sum of the transverse momentum (see @ref vars::dpT_lp)
     * and the missing longitudinal momentum (see @ref vars::dpL_lp).
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the estimate of the momentum of the struck nucleon.
     * @note The switch to the NuMI beam direction instead of the BNB axis is
     * applied by the definition of a preprocessor macro (BEAM_IS_NUMI).
     */
    template<class T>
    double pn_lp(const T & obj) { return std::sqrt(std::pow(vars::dpT_lp(obj), 2) + std::pow(vars::dpL_lp(obj), 2)); }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, pn_lp, pn_lp);

    /**
     * @brief Variable for the opening angle between leading muon and proton.
     * @details The leading muon and proton are defined as the particles with the
     * highest kinetic energy. The opening angle is defined as the arccosine of
     * the dot product of the momentum vectors of the leading muon and proton.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the opening angle between the leading muon and
     * proton.
     */
    template<class T>
    double opening_angle(const T & obj, std::vector<double> params={1.0,})
    {
        size_t mi = selectors::leading_muon(obj);
        size_t pi = kNoMatch;

        const int mode = params.empty() ? 1 : static_cast<int>(params[0]);

        if (mode == 1) {            // For 1pi channel
            pi = selectors::leading_pion(obj);
        }
        else if (mode == 2) {       // For sideband
            pi = selectors::leading_proton(obj);
        }
        else {
            return kNoMatchValue;   // unknown mode
        }
        if(mi == kNoMatch || pi == kNoMatch)
            return kNoMatchValue; // No leading muon or proton found.
        else
        {
            auto & m(obj.particles[mi]);
            auto & p(obj.particles[pi]);
            //return std::acos(m.start_dir[0] * p.start_dir[0] + m.start_dir[1] * p.start_dir[1] + m.start_dir[2] * p.start_dir[2]);
            return m.start_dir[0] * p.start_dir[0] + m.start_dir[1] * p.start_dir[1] + m.start_dir[2] * p.start_dir[2];
        }
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, opening_angle, opening_angle);

    /**
     * @brief Variable for the (primary) photon multiplicity of the
     * interaction.
     * @details This function calculates the multiplicity of primary
     * photons in the interaction by counting the number of primary particles
     * that are identified as photons and have a kinetic energy above a
     * threshold. The threshold is set by the `params` vector, which defaults
     * to 25 MeV. The function returns the number of primary photons in the
     * interaction.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @param params the parameters for the cut. In this case, this sets the
     * kinetic energy threshold for a photon to count towards the
     * multiplicity. Defaults to 25 MeV.
     * @return the multiplicity of primary photons in the interaction.
     */
    template<class T>
    double photon_multiplicity(const T & obj, std::vector<double> params={25.0,})
    {
        size_t count(0);
        for(const auto & p : obj.particles)
        {
            if(pvars::pid(p) == pvars::kPhoton && pvars::primary_classification(p) && pvars::ke(p) >= params[0])
                ++count;
        }
        return count;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, photon_multiplicity, photon_multiplicity);

    /**
     * @brief Variable for the (primary) electron multiplicity of the
     * interaction.
     * @details This function calculates the multiplicity of primary electrons
     * in the interaction by counting the number of primary particles that are
     * identified as electrons and have a kinetic energy above a threshold. The
     * threshold is set by the `params` vector, which defaults to 25 MeV. The
     * function returns the number of primary electrons in the interaction.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @param params the parameters for the cut. In this case, this sets the
     * kinetic energy threshold for an electron to count towards the
     * multiplicity. Defaults to 25 MeV.
     * @return the multiplicity of primary electrons in the interaction.
     */
    template<class T>
    double electron_multiplicity(const T & obj, std::vector<double> params={25.0,})
    {
        size_t count(0);
        for(const auto & p : obj.particles)
        {
            if(pvars::pid(p) == pvars::kElectron && pvars::primary_classification(p) && pvars::ke(p) >= params[0])
                ++count;
        }
        return count;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, electron_multiplicity, electron_multiplicity);

    /**
     * @brief Variable for the (primary) muon multiplicity of the
     * interaction.
     * @details This function calculates the multiplicity of primary muons in
     * the interaction by counting the number of primary particles that are
     * identified as muons and have a kinetic energy above a threshold. The
     * threshold is set by the `params` vector, which defaults to 25 MeV. The
     * function returns the number of primary muons in the interaction.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @param params the parameters for the cut. In this case, this sets the
     * kinetic energy threshold for a muon to count towards the
     * multiplicity. Defaults to 25 MeV.
     * @return the multiplicity of primary muons in the interaction.
     */
    template<class T>
    double muon_multiplicity(const T & obj, std::vector<double> params={25.0,})
    {
        size_t count(0);
        for(const auto & p : obj.particles)
        {
            if(pvars::pid(p) == pvars::kMuon && pvars::primary_classification(p) && pvars::ke(p) >= params[0])
                ++count;
        }
        return count;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, muon_multiplicity, muon_multiplicity);

    /**
     * @brief Variable for the (primary) pion multiplicity of the
     * interaction.
     * @details This function calculates the multiplicity of primary pions in
     * the interaction by counting the number of primary particles that are
     * identified as pions and have a kinetic energy above a threshold. The
     * threshold is set by the `params` vector, which defaults to 25 MeV. The
     * function returns the number of primary pions in the interaction.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @param params the parameters for the cut. In this case, this sets the
     * kinetic energy threshold for a pion to count towards the
     * multiplicity. Defaults to 25 MeV.
     * @return the multiplicity of primary pions in the interaction.
     */
    template<class T>
    double pion_multiplicity(const T & obj, std::vector<double> params={50.0,})
    {
        size_t count(0);
        for(const auto & p : obj.particles)
        {
            if(pvars::pid(p) == pvars::kPion && pvars::primary_classification(p) && pvars::ke(p) >= params[0])
                ++count;
        }
        return count;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, pion_multiplicity, pion_multiplicity);

    /**
     * @brief Variable for the (primary) proton multiplicity of the
     * interaction.
     * @details This function calculates the multiplicity of primary protons in
     * the interaction by counting the number of primary particles that are
     * identified as protons and have a kinetic energy above a threshold. The
     * threshold is set by the `params` vector, which defaults to 25 MeV. The
     * function returns the number of primary protons in the interaction.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @param params the parameters for the cut. In this case, this sets the
     * kinetic energy threshold for a proton to count towards the
     * multiplicity. Defaults to 25 MeV.
     * @return the multiplicity of primary protons in the interaction.
     */
    template<class T>
    double proton_multiplicity(const T & obj, std::vector<double> params={25.0,})
    {
        size_t count(0);
        for(const auto & p : obj.particles)
        {
            if(pvars::pid(p) == pvars::kProton && pvars::primary_classification(p) && pvars::ke(p) >= params[0])
                ++count;
        }
        return count;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, proton_multiplicity, proton_multiplicity);

    /**
     * @brief Variable for the distance between the interaction vertex and the
     * leading muon start point.
     * @details This function calculates the distance from the leading muon
     * start point to the interaction vertex. The leading muon is defined as
     * the particle with the highest kinetic energy that is identified as a
     * muon. If no leading muon is found, the function returns the usual 
     * PLACEHOLDERVALUE.
     * @tparam T the type of interaction (true or reco).
     * @param obj the interaction to apply the variable on.
     * @return the distance from the leading muon start point to the
     * interaction vertex.
     */
    template<class T>
    double leading_muon_vertex_gap(const T & obj)
    {
        // Find the leading muon in the interaction.
        size_t mi = selectors::leading_muon(obj);
        if(mi == kNoMatch) return PLACEHOLDERVALUE;
        auto & m(obj.particles[mi]);
        
        // Calculate the distance from the leading muon start point to the
        // interaction vertex.
        utilities::three_vector vtx = {obj.vertex[0], obj.vertex[1], obj.vertex[2]};
        utilities::three_vector muon_start = {pvars::start_x(m), pvars::start_y(m), pvars::start_z(m)};
        return utilities::magnitude(utilities::subtract(muon_start, vtx));
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, leading_muon_vertex_gap, leading_muon_vertex_gap);

    template<class T>
    double muon_cos_angle(const T & obj)
    {
        size_t mi = selectors::leading_muon(obj);
        if(mi == kNoMatch)
            return kNoMatchValue; // No leading muon or proton found.
        else
        {
            auto & m(obj.particles[mi]);
            utilities::three_vector vtx = {pvars::start_x(m), pvars::start_y(m), pvars::start_z(m)};
            utilities::three_vector mu_unit = std::make_tuple(m.start_dir[0], m.start_dir[1], m.start_dir[2]);
            utilities::three_vector beam_unit = utilities::numi_beam_direction(vtx);
            return utilities::dot_product(mu_unit, beam_unit);
        }
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, muon_cos_angle, muon_cos_angle);

    template<class T>
    double pion_cos_angle(const T & obj)
    {
        size_t pi = selectors::leading_pion(obj);
        if(pi == kNoMatch)
            return kNoMatchValue; // No leading muon or proton found.
        else
        {
            auto & p(obj.particles[pi]);
            utilities::three_vector vtx = {pvars::start_x(p), pvars::start_y(p), pvars::start_z(p)};
            utilities::three_vector pi_unit = std::make_tuple(p.start_dir[0], p.start_dir[1], p.start_dir[2]);
            utilities::three_vector beam_unit = utilities::numi_beam_direction(vtx);
            return utilities::dot_product(pi_unit, beam_unit);
        }
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, pion_cos_angle, pion_cos_angle);

    template<class T>
    double four_mom_transfer_squared(const T & obj)
    {
        double reco_nuE = visible_energy(obj)*1000.0;
        double costheta = muon_cos_angle(obj);

        size_t mi = selectors::leading_muon(obj);
        if(mi == kNoMatch)
            return kNoMatchValue; // No leading muon or proton found.
        else
        {
            auto & m(obj.particles[mi]);
            utilities::three_vector momentum = {pvars::px(m), pvars::py(m), pvars::pz(m)};
            double mu_mom = utilities::magnitude(momentum);
            double mu_E = pvars::ke(m) + MUON_MASS;
//            std::cout<<mu_E-mu_mom*costheta<<std::endl;
            return (2*reco_nuE*(mu_E-mu_mom*costheta)-std::pow(MUON_MASS, 2))/1000000.;
        }
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, four_mom_transfer_squared, four_mom_transfer_squared);

    template<class T>
    double four_mom_transfer_squared_calo(const T & obj)
    {
        double reco_nuE = visible_energy_calosub(obj)*1000.0;
        double costheta = muon_cos_angle(obj);

        size_t mi = selectors::leading_muon(obj);
        if(mi == kNoMatch)
            return kNoMatchValue; // No leading muon or proton found.
        else
        {
            auto & m(obj.particles[mi]);
            utilities::three_vector momentum = {pvars::px(m), pvars::py(m), pvars::pz(m)};
            double mu_mom = utilities::magnitude(momentum);
            double mu_E = pvars::ke(m) + MUON_MASS;
//            std::cout<<mu_E-mu_mom*costheta<<std::endl;
            return (2*reco_nuE*(mu_E-mu_mom*costheta)-std::pow(MUON_MASS, 2))/1000000.;
        }
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, four_mom_transfer_squared_calo, four_mom_transfer_squared_calo);

    template<class T>
    double q0_transfer(const T & obj)
    {
        double reco_nuE = visible_energy(obj)*1000.0; // GeV → MeV (consistent with four_mom_transfer_squared)
        double costheta = muon_cos_angle(obj);
    
        size_t mi = selectors::leading_muon(obj);
        if(mi == kNoMatch)
            return kNoMatchValue; // No leading muon or proton found.
    
        auto & m(obj.particles[mi]);
        utilities::three_vector momentum = {pvars::px(m), pvars::py(m), pvars::pz(m)};
        double mu_mom = utilities::magnitude(momentum);
        double mu_E   = pvars::ke(m) + MUON_MASS;  // MeV
    
        // Energy transfer q0 = Eν - Eμ (in MeV)
        return reco_nuE - mu_E;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, q0_transfer, q0_transfer);
    
    template<class T>
    double q3_transfer(const T & obj)
    {
        double reco_nuE = visible_energy(obj)*1000.0; // GeV → MeV
        double costheta = muon_cos_angle(obj);
    
        size_t mi = selectors::leading_muon(obj);
        if(mi == kNoMatch)
            return kNoMatchValue; // No leading muon or proton found.
    
        auto & m(obj.particles[mi]);
        utilities::three_vector momentum = {pvars::px(m), pvars::py(m), pvars::pz(m)};
        double mu_mom = utilities::magnitude(momentum);
        double mu_E   = pvars::ke(m) + MUON_MASS;  // MeV
    
        // Q^2 (same as your four_mom_transfer_squared)
        double Q2 = 2*reco_nuE*(mu_E - mu_mom*costheta) - std::pow(MUON_MASS, 2);
    
        // q0
        double q0 = reco_nuE - mu_E;
    
        // |q⃗| = sqrt(q0^2 + Q^2)  (all in MeV)
        double q3_sq = q0*q0 + Q2;
        if(q3_sq < 0.0)
            return kNoMatchValue;
    
        return std::sqrt(q3_sq);
    }
    
    // Register them like the others
    REGISTER_VAR_SCOPE(RegistrationScope::Both, q3_transfer, q3_transfer);

    template<class T>
    double hadron_invariant_mass_general(const T & obj)
    {
        double reco_nuE = visible_energy(obj)*1000.0;
        size_t mi = selectors::leading_muon(obj);
        if(mi == kNoMatch)
            return -999.;
        else
        {
            auto & m(obj.particles[mi]);
            double mu_E = pvars::ke(m) + MUON_MASS;
            double W_2 = std::pow(NEUTRON_MASS, 2)+2*NEUTRON_MASS*(reco_nuE-mu_E)-four_mom_transfer_squared(obj);
 //           std::cout<<W_2<<std::endl;
            if (W_2<0)
            {
                return -999.;
            }
            else
            {
                return std::sqrt(W_2);
            }
        }
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, hadron_invariant_mass_general, hadron_invariant_mass_general);

    template<class T>
    double hadron_invariant_mass_general_calo(const T & obj)
    {
        double reco_nuE = visible_energy_calosub(obj)*1000.0;
        size_t mi = selectors::leading_muon(obj);
        if(mi == kNoMatch)
            return -999.;
        else
        {
            auto & m(obj.particles[mi]);
            double mu_E = pvars::ke(m) + MUON_MASS;
            double W_2 = std::pow(NEUTRON_MASS, 2)+2*NEUTRON_MASS*(reco_nuE-mu_E)-four_mom_transfer_squared(obj);
 //           std::cout<<W_2<<std::endl;
            if (W_2<0)
            {
                return -999.;
            }
            else
            {
                return std::sqrt(W_2);
            }
        }
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, hadron_invariant_mass_general_calo, hadron_invariant_mass_general_calo);


    template<class T>
    double selected_proton_ke(const T & obj, std::vector<double> params={25.0,})
    {

        double ke=-990.0;
        for(const auto & p : obj.particles)
        {
            if(pvars::pid(p) == 3 && pvars::primary_classification(p) && pvars::ke(p) >= params[0])
                ke=pvars::ke(p);
        }
        return ke;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, selected_proton_ke, selected_proton_ke);
    template<class T>
    double selected_pion_ke(const T & obj, std::vector<double> params={25.0,})
    {

        double ke=-990.0;
        for(const auto & p : obj.particles)
        {
            if(pvars::pid(p) == 3 && pvars::primary_classification(p) && pvars::ke(p) >= params[0])
                ke=pvars::ke(p);
        }
        return ke;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, selected_pion_ke, selected_pion_ke);


    template<class T>
    double selected_muon_ke(const T & obj, std::vector<double> params={143.425,})
    {

        double ke=-990.0;
        for(const auto & p : obj.particles)
        {
            if(pvars::pid(p) == 2 && pvars::primary_classification(p) && pvars::ke(p) >= params[0])
                ke=pvars::ke(p);
        }
        return ke;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, selected_muon_ke, selected_muon_ke);

    template<class T>
    double is_theresplitmuon(const T & obj, std::vector<double> params={143.425,})
    {
        size_t mu = selectors::leading_muon(obj);//gaurantee it is reco pion candidate
        utilities::three_vector vtx = {obj.vertex[0], obj.vertex[1], obj.vertex[2]};
        bool additionalmuon = false;
        size_t index(kNoMatch);
        if(mu == kNoMatch)
            return 0.0;
        else
        {
            auto & muon(obj.particles[mu]);
            int64_t muon_id = muon.id;
            for(size_t i(0); i < obj.particles.size(); ++i)
            {
                const auto & p = obj.particles[i];
                if ((i != muon_id)&& p.pid==2)
                {
                    additionalmuon=true;
                    index = i;
                }
            }
            if (additionalmuon)
            {
                auto & addmuon(obj.particles[index]);
                utilities::three_vector track_vtx = {pvars::start_x(addmuon), pvars::start_y(addmuon), pvars::start_z(addmuon)};
                double Atslc = utilities::magnitude(utilities::subtract(vtx, track_vtx));
                //if (Atslc>10 && addmuon.ke<params[0])
                if (Atslc>10)
                    return 1.0;
                else
                    return 0.0;              // <-- ensures all paths return

            }
            else return 0.0;
        }
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, is_theresplitmuon, is_theresplitmuon);


    // --- Keep this helper: resolve leading-pion children PDGs ---
    template <class T>
    std::vector<int64_t> get_pion_children_pdgs(const T& obj){
        std::vector<int64_t> v;
        const size_t pi = selectors::leading_pion(obj);
        if (pi == kNoMatch) return v;

        const auto& pion = obj.particles[pi];

        std::unordered_map<int64_t,size_t> id2i;
        id2i.reserve(obj.particles.size());
        for (size_t i=0;i<obj.particles.size();++i)
            id2i.emplace((int64_t)obj.particles[i].id, i);

        for (auto cid_raw: pion.children_id){
            auto it = id2i.find((int64_t)cid_raw);
            if (it==id2i.end()) continue;
            v.push_back((int64_t)obj.particles[it->second].pdg_code);
        }
        return v;
    }

    // Build: parent_id -> list of child indices
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


                                                                                                                                                                                                                                           
// Descendant scan flags must be templated consistently                                                                                                                                                                                    
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


    enum class PionGroupCode:int{ UNKNOWN=0, DECAY=1, CAPTURE=2, INELASTIC=3, ELASTIC=4 };

    template <class T>
        int classify_pion_group_code_int(const T& obj)
        {
            const size_t pi = selectors::leading_pion(obj);
            if (pi == kNoMatch) return (int)PionGroupCode::UNKNOWN;

            // Direct children PDGs (needed for gamma-only and counting nucleons)
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
            const auto df = scan_direct_children(obj, pi);

            // -------- DECAY (direct-only) ----------
            // Require a direct muon to call DECAY.
            // (Optional) Allow a direct Michel-shaped e± ONLY if there are no direct pions/pi0/nuclear products.
            if (df.hasMuon)
                return (int)PionGroupCode::DECAY;

            if (df.hasMichelStrict && !df.hasChargedPi && !df.hasPi0 && !anyN_dir)
                return (int)PionGroupCode::DECAY;

            // -------- CAPTURE (direct-only) --------
            // π- capture often shows photons only, or nucleons with no pions.
            if (onlyGammas())
                return (int)PionGroupCode::CAPTURE;

            if (anyN_dir && !df.hasChargedPi && !df.hasPi0)
                return (int)PionGroupCode::CAPTURE;

            // -------- INELASTIC / ELASTIC ----------
            // π0 among direct daughters → inelastic
            if (df.hasPi0)
                return (int)PionGroupCode::INELASTIC;

            // Charged π with nucleons → inelastic
            if (df.hasChargedPi && anyN_dir)
                return (int)PionGroupCode::INELASTIC;

            // Exactly one direct charged π and nothing else → elastic-like
            {
                int nCpi = 0, nOther = 0;
                for (auto x: ch) {
                    if (std::abs(x) == 211) ++nCpi;
                    else if (x != 0) ++nOther; // ignore 0 if it can appear
                }
                if (nCpi == 1 && nOther == 0)
                    return (int)PionGroupCode::ELASTIC; // or INELASTIC if you don't keep ELASTIC separate
            }

            // Any charged π at all (and no stronger signature above) → inelastic bucket
            if (df.hasChargedPi)
                return (int)PionGroupCode::INELASTIC;

            // -------- Fallback ----------
            return (int)PionGroupCode::UNKNOWN;
        }

    // --- Scalar wrappers (double) for your binder ---
    template <class T> double pion_group_code(const T& obj){ return (double)classify_pion_group_code_int(obj); }
    template <class T> double pion_is_decay(const T& obj){ return pion_group_code(obj)==(double)PionGroupCode::DECAY ? 1.0:0.0; }
    template <class T> double pion_is_capture(const T& obj){ return pion_group_code(obj)==(double)PionGroupCode::CAPTURE ? 1.0:0.0; }
    template <class T> double pion_is_inelastic(const T& obj){ return pion_group_code(obj)==(double)PionGroupCode::INELASTIC ? 1.0:0.0; }

    // (Optional) tiny debug helpers
    template <class T> double pion_children_count(const T& obj){ return (double)get_pion_children_pdgs(obj).size(); }
    template <class T> double pion_child0_pdg(const T& obj){
        auto v=get_pion_children_pdgs(obj);

        return v.empty()? -9999.0 : (double)v[0]; }
    REGISTER_VAR_SCOPE(RegistrationScope::True, pion_group_code,     vars::pion_group_code);
    REGISTER_VAR_SCOPE(RegistrationScope::True, pion_is_decay,       vars::pion_is_decay);
    REGISTER_VAR_SCOPE(RegistrationScope::True, pion_is_capture,     vars::pion_is_capture);
    REGISTER_VAR_SCOPE(RegistrationScope::True, pion_is_inelastic,   vars::pion_is_inelastic);
    REGISTER_VAR_SCOPE(RegistrationScope::True, pion_children_count, vars::pion_children_count);
    REGISTER_VAR_SCOPE(RegistrationScope::True, pion_child0_pdg,     vars::pion_child0_pdg);

    template <class T>                                                                                                                                                                                                                     
    double muon_michel_tag(const T& obj)                                                                                                                                                                                                   
    {                                                                                                                                                                                                                                      
        size_t mu = selectors::leading_muon(obj);//gaurantee it is reco pion candidate                                                                                                                                                     
        size_t index(kNoMatch);                                                                                                                                                                                                            
        bool michel_tagged = false;                                                                                                                                                                                                        
        if(mu == kNoMatch)                                                                                                                                                                                                                 
            return 0.0;                                                                                                                                                                                                                    
        else                                                                                                                                                                                                                               
        {                                                                                                                                                                                                                                  
            auto & muon(obj.particles[mu]);                                                                                                                                                                                                
            utilities::three_vector muon_end = {pvars::end_x(muon), pvars::end_y(muon), pvars::end_z(muon)};                                                                                                                               
//            std::cout<<pvars::end_x(muon)<<std::endl;                                                                                                                                                                                    
            int64_t muon_id = muon.id;                                                                                                                                                                                                     
            for(size_t i(0); i < obj.particles.size(); ++i)                                                                                                                                                                                
            {                                                                                                                                                                                                                              
                const auto & p = obj.particles[i];                                                                                                                                                                                         
                if (p.id == muon_id) continue;                                                                                                                                                                                             
                utilities::three_vector particle_vtx = {pvars::start_x(p), pvars::start_y(p), pvars::start_z(p)};                                                                                                                          
                double Atslc = utilities::magnitude(utilities::subtract(muon_end, particle_vtx));                                                                                                                                          
//                std::cout<<"p.shape : "<<p.shape<<", pdg : "<< p.pdg_code<<", p.id: "<<p.id<<std::endl;                                                                                                                                  
//                std::cout<<"Atslc : "<<Atslc<<std::endl;                                                                                                                                                                                 
                if (p.shape==2 && Atslc<10)                                                                                                                                                                                                 
                {                                                                                                                                                                                                                          
                    michel_tagged = true;                                                                                                                                                                                                  
                }                                                                                                                                                                                                                          
                                                                                                                                                                                                                                           
            }                                                                                                                                                                                                                              
            if (michel_tagged==true)                                                                                                                                                                                                       
            {                                                                                                                                                                                                                              
                return 1.0;                                                                                                                                                                                                                
            }                                                                                                                                                                                                                              
            else                                                                                                                                                                                                                           
            {                                                                                                                                                                                                                              
                return 0.0;                                                                                                                                                                                                                
            }                                                                                                                                                                                                                              
        }                                                                                                                                                                                                                                  
    }                                                                                                                                                                                                                                      
                                                                                                                                                                                                                                           
    REGISTER_VAR_SCOPE(RegistrationScope::Both, muon_michel_tag, muon_michel_tag);                                                                                                                                                         

    namespace michel_util {
        struct V3 { double x, y, z; };
        template <class P> static inline V3 start_pos(const P& p){ return {pvars::start_x(p), pvars::start_y(p), pvars::start_z(p)}; }
        template <class P> static inline V3 end_pos  (const P& p){ return {pvars::end_x  (p), pvars::end_y  (p), pvars::end_z  (p)}; }
        static inline V3 sub(const V3& a, const V3& b){ return {a.x-b.x, a.y-b.y, a.z-b.z}; }
        static inline double mag(const V3& v){ return std::sqrt(v.x*v.x + v.y*v.y + v.z*v.z); }
    }

    // Find the descendant that is a Michel by (shape==2 AND |pdg|==11),
    // picking the one whose start is **closest** to the pion end.
// Enforce π -> μ -> e (Michel) chain.
// Returns the index of the best Michel (electron/positron) among muon children of the pion.
// If none is found, returns std::nullopt.

    //Only tag michel when pion->muon->michel
    template <class T>
    static inline std::optional<size_t>
    find_descendant_michel_index(const T& obj, size_t pion_idx)
    {
        if (pion_idx >= obj.particles.size()) return std::nullopt;

        const auto& pion = obj.particles[pion_idx];
        const int64_t pion_id = pion.id;

        // Build parent->children map (you already have this)
        const auto parent_children = build_parent_children_map(obj);

        // Get pion end position for proximity check
        const auto pion_end = michel_util::end_pos(pion);

        double best_dist = std::numeric_limits<double>::infinity();
        std::optional<size_t> best_idx;

        // 1) Get direct children of the pion
        auto it_pion_children = parent_children.find(pion_id);
        if (it_pion_children == parent_children.end())
            return std::nullopt;

        // Optional: a small spatial window to avoid random electrons
        constexpr double R_MAX_MM = 50.0; // tune as you like (mm)

        // 2) For each child that is a muon, scan its children for electrons
        for (size_t mu_idx : it_pion_children->second) {
            const auto& mu = obj.particles[mu_idx];
            if (std::abs(mu.pdg_code) != 13) continue; // must be muon

            // Get muon's children
            auto it_mu_children = parent_children.find(mu.id);
            if (it_mu_children == parent_children.end()) continue;

            for (size_t e_idx : it_mu_children->second) {
                const auto& e = obj.particles[e_idx];

                // Candidate Michel electron/positron
                if (std::abs(e.pdg_code) != 11) continue;
                // (Optional) Require shower-like shape if that's your convention
                 if (e.shape != 2) continue;

                // (Optional) Require decay process tags if available
                // if (mu.creation_process != "Decay") continue;
                // if (e.creation_process  != "Decay") continue;

                // Proximity to pion end
                const auto d = michel_util::mag(
                    michel_util::sub(pion_end, michel_util::start_pos(e))
                );

                if (d <= R_MAX_MM && d < best_dist) {
                    best_dist = d;
                    best_idx  = e_idx;
                }
            }
        }

        return best_idx;
    }

    // ----------------------------------------------------
    // 1) Distance: pion_end -> Michel_start (mm)
    // ----------------------------------------------------
    // only michel for explicitly pion->muon->michel
    template <class T>
    double pion_michel_distance(const T& obj)
    {
        const size_t pi_idx = selectors::leading_pion(obj);
        if (pi_idx == kNoMatch || pi_idx >= obj.particles.size()) return -5.0;

        auto idx = find_descendant_michel_index(obj, pi_idx);
        if (!idx) return -5.0;

        const auto& pion   = obj.particles[pi_idx];
        const auto& michel = obj.particles[*idx];

        return michel_util::mag(
                michel_util::sub(michel_util::end_pos(pion), michel_util::start_pos(michel))
                );
    }

    // ----------------------------------------------------
    // 2) Michel kinetic energy (same units as p.ke)
    // ----------------------------------------------------
    template <class T>
    double pion_michel_energy(const T& obj)
    {
        const size_t pi_idx = selectors::leading_pion(obj);
        if (pi_idx == kNoMatch || pi_idx >= obj.particles.size()) return -5.0;

        auto idx = find_descendant_michel_index(obj, pi_idx);
        if (!idx) return -5.0;

        const auto& michel = obj.particles[*idx];
        // If your canonical accessor is pvars::ke(michel), swap it here:
        return michel.ke;
    }

    // ------------------
    // Registration
    // ------------------
    REGISTER_VAR_SCOPE(RegistrationScope::True, pion_michel_distance, pion_michel_distance);
    REGISTER_VAR_SCOPE(RegistrationScope::True, pion_michel_energy,   pion_michel_energy);


    template <class T>
    double true_pion_michel_tag(const T& obj)
    {
        const size_t pi_idx = selectors::leading_pion(obj);
        if (pi_idx == kNoMatch || pi_idx >= obj.particles.size()) return 0.0;

        auto idx = find_descendant_michel_index(obj, pi_idx);
        if (!idx) return -0.0;

        return 1.0;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::True, true_pion_michel_tag, true_pion_michel_tag);

    template <class T>
    double pion_michel_tag(const T& obj, std::vector<double> params={10.0})
    {
        size_t mu = selectors::leading_pion(obj);//gaurantee it is reco pion candidate
        size_t index(kNoMatch);
        bool michel_tagged = false;
        if(mu == kNoMatch)
            return 0.0;
        else
        {
            auto & pion(obj.particles[mu]);
            utilities::three_vector pion_end = {pvars::end_x(pion), pvars::end_y(pion), pvars::end_z(pion)};
            int64_t pion_id = pion.id;
            for(size_t i(0); i < obj.particles.size(); ++i)
            {
                const auto & p = obj.particles[i];
                if (p.id == pion_id) continue;
                utilities::three_vector particle_vtx = {pvars::start_x(p), pvars::start_y(p), pvars::start_z(p)};
                double Atslc = utilities::magnitude(utilities::subtract(pion_end, particle_vtx));
                if (p.shape==2 && Atslc<params[0])
                {
                    michel_tagged = true;
                }

            }
            if (michel_tagged==true)
            {
                return 1.0;
            }
            else
            {
                return 0.0;
            }
        }
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, pion_michel_tag, pion_michel_tag);
    
    template <class T>
    double multiple_pion_michel_tag(const T& obj, std::vector<double> params={10.0})
    {
        size_t mu = selectors::leading_pion(obj);//gaurantee it is reco pion candidate
        size_t index(kNoMatch);
        bool michel_tagged = false;
        if(mu == kNoMatch)
            return 0.0;
        else
        {
            auto & pion(obj.particles[mu]);
            utilities::three_vector pion_end = {pvars::end_x(pion), pvars::end_y(pion), pvars::end_z(pion)};
            int64_t pion_id = pion.id;
            for(size_t i(0); i < obj.particles.size(); ++i)
            {
                const auto & p = obj.particles[i];
                if (p.id == pion_id) continue;
                utilities::three_vector particle_vtx = {pvars::start_x(p), pvars::start_y(p), pvars::start_z(p)};
                double Atslc = utilities::magnitude(utilities::subtract(pion_end, particle_vtx));
                if (p.shape==2 && Atslc<params[0])
                {
                    michel_tagged = true;
                }

            }
            if (michel_tagged==true)
            {
                return 1.0;
            }
            else
            {
                return 0.0;
            }
        }
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, multiple_pion_michel_tag, multiple_pion_michel_tag);



    template <class T>
    static inline std::optional<size_t>
    find_descendant_michel_index_direct(const T& obj, size_t pion_idx)
    {
        if (pion_idx >= obj.particles.size()) return std::nullopt;

        const auto& pion = obj.particles[pion_idx];
        const int64_t pion_id = pion.id;

        // Build parent->children map (you already have this)
        const auto parent_children = build_parent_children_map(obj);

        // Get pion end position for proximity check
        const auto pion_end = michel_util::end_pos(pion);

        double best_dist = std::numeric_limits<double>::infinity();
        std::optional<size_t> best_idx;

        // 1) Get direct children of the pion
        auto it_pion_children = parent_children.find(pion_id);
        if (it_pion_children == parent_children.end())
            return std::nullopt;

        // Optional: a small spatial window to avoid random electrons
        constexpr double R_MAX_MM = 50.0; // tune as you like (mm)

        // 2) For each child that is a muon, scan its children for electrons
        for (size_t mu_idx : it_pion_children->second) {
            const auto& mu = obj.particles[mu_idx];
            if (std::abs(mu.pdg_code) != 11) continue;
                 //if (e.shape != 2) continue;
            const auto d = michel_util::mag(
                    michel_util::sub(pion_end, michel_util::start_pos(mu))
                    );

            if (d <= R_MAX_MM && d < best_dist) {
                best_dist = d;
                best_idx  = mu_idx;
            }
        }

        return best_idx;
    }

    template <class T>
    double proton_michel_distance(const T& obj)
    {
        const size_t pi_idx = selectors::leading_proton(obj);
        if (pi_idx == kNoMatch || pi_idx >= obj.particles.size()) return -5.0;

        auto idx = find_descendant_michel_index_direct(obj, pi_idx);
        if (!idx) return -5.0;

        const auto& pion   = obj.particles[pi_idx];
        const auto& michel = obj.particles[*idx];

        return michel_util::mag(
                michel_util::sub(michel_util::end_pos(pion), michel_util::start_pos(michel))
                );
    }
    REGISTER_VAR_SCOPE(RegistrationScope::True, proton_michel_distance, proton_michel_distance);

    template <class T>
    double muon_michel_distance_from_pion(const T& obj)
    {
        const size_t mu_idx = selectors::leading_muon(obj);
        const size_t pi_idx = selectors::leading_pion(obj);
        if (mu_idx == kNoMatch || mu_idx >= obj.particles.size()) return -5.0;

        auto idx = find_descendant_michel_index_direct(obj, mu_idx);
        if (!idx) return -5.0;

        const auto& pion   = obj.particles[pi_idx];
        const auto& michel = obj.particles[*idx];

        return michel_util::mag(
                michel_util::sub(michel_util::end_pos(pion), michel_util::start_pos(michel))
                );
    }
    REGISTER_VAR_SCOPE(RegistrationScope::True, muon_michel_distance_from_pion, muon_michel_distance_from_pion);


    // --- Keep this helper: resolve leading-proton children PDGs ---
    template <class T>
    std::vector<int64_t> get_proton_children_pdgs(const T& obj){
        std::vector<int64_t> v;
        const size_t pi = selectors::least_proton(obj);
        if (pi == kNoMatch) return v;
//       std::cout<<"Event::"<<std::endl;

        const auto& proton = obj.particles[pi];

        std::unordered_map<int64_t,size_t> id2i;
        id2i.reserve(obj.particles.size());
        for (size_t i=0;i<obj.particles.size();++i)
        {
//            std::cout<<"id: "<<obj.particles[i].id<<", pdg: "<<obj.particles[i].pdg_code<<", parent id: "<<obj.particles[i].parent_id<<std::endl;
            id2i.emplace((int64_t)obj.particles[i].id, i);
        }

        for (auto cid_raw: proton.children_id){
            auto it = id2i.find((int64_t)cid_raw);
            if (it==id2i.end()) continue;
//            std::cout<<"proton child : "<<(int64_t)obj.particles[it->second].pdg_code<<std::endl;
            v.push_back((int64_t)obj.particles[it->second].pdg_code);
        }
        return v;
    }

    template <class T> double proton_children_count(const T& obj){ return (double)get_proton_children_pdgs(obj).size(); }
    REGISTER_VAR_SCOPE(RegistrationScope::True, proton_children_count, vars::proton_children_count);

    template<class T>
    double opening_angle_proton(const T & obj)
    {
        size_t mi = selectors::leading_muon(obj);
        size_t pr = selectors::leading_proton(obj);
        if(mi == kNoMatch || pr == kNoMatch)
            return kNoMatchValue; // No leading muon or proton found.
        else
        {
            auto & m(obj.particles[mi]);
            auto & p(obj.particles[pr]);
            //return std::acos(m.start_dir[0] * p.start_dir[0] + m.start_dir[1] * p.start_dir[1] + m.start_dir[2] * p.start_dir[2]);
            return (m.start_dir[0] * p.start_dir[0] + m.start_dir[1] * p.start_dir[1] + m.start_dir[2] * p.start_dir[2]);
        }
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, opening_angle_proton, opening_angle_proton);

    template<class T>
    double opening_angle_least_proton(const T & obj)
    {
        size_t mi = selectors::leading_muon(obj);
        size_t pr = selectors::least_proton(obj);
        if(mi == kNoMatch || pr == kNoMatch)
            return kNoMatchValue; // No leading muon or proton found.
        else
        {
            auto & m(obj.particles[mi]);
            auto & p(obj.particles[pr]);
            //return std::acos(m.start_dir[0] * p.start_dir[0] + m.start_dir[1] * p.start_dir[1] + m.start_dir[2] * p.start_dir[2]);
            return (m.start_dir[0] * p.start_dir[0] + m.start_dir[1] * p.start_dir[1] + m.start_dir[2] * p.start_dir[2]);
        }
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, opening_angle_least_proton, opening_angle_least_proton);

    template<class T>
    double leadingproton_match_id_max(const T & obj)
    {
        size_t pr = selectors::leading_proton(obj);
        if (pr == kNoMatch)
            return kNoMatchValue;
    
        const auto & v = obj.particles[pr].match_overlaps;
        if (v.empty())
            return kNoMatchValue;
    
        // Manual search for max value and its index
        size_t best_idx = 0;
        float  best_val = v[0];
    
        for (size_t i = 1; i < v.size(); ++i) {
            if (v[i] > best_val) {
                best_val = v[i];
                best_idx = i;
            }
        }
    
        // If you want the index of the biggest overlap:
        return obj.particles[pr].match_ids[best_idx];
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, leadingproton_match_id_max, leadingproton_match_id_max);

    template<class T>
    double leadingproton_match_overlap_max(const T & obj)
    {
        size_t pr = selectors::leading_proton(obj);
        if (pr == kNoMatch)
            return kNoMatchValue;
    
        const auto & v = obj.particles[pr].match_overlaps;
        if (v.empty())
            return kNoMatchValue;
    
        float best_val = v[0];
    
        for (size_t i = 1; i < v.size(); ++i) {
            if (v[i] > best_val)
                best_val = v[i];
        }
    
        return static_cast<double>(best_val);
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, leadingproton_match_overlap_max, leadingproton_match_overlap_max);

    template<class T>
    double leadingpion_match_id_max(const T & obj)
    {
        size_t pr = selectors::leading_pion(obj);
        if (pr == kNoMatch)
            return kNoMatchValue;
    
        const auto & v = obj.particles[pr].match_overlaps;
        if (v.empty())
            return kNoMatchValue;
    
        // Manual search for max value and its index
        size_t best_idx = 0;
        float  best_val = v[0];
    
        for (size_t i = 1; i < v.size(); ++i) {
            if (v[i] > best_val) {
                best_val = v[i];
                best_idx = i;
            }
        }
    
        // If you want the index of the biggest overlap:
        return obj.particles[pr].match_ids[best_idx];
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, leadingpion_match_id_max, leadingpion_match_id_max);

    template<class T>
    double leadingpion_match_overlap_max(const T & obj)
    {
        size_t pr = selectors::leading_pion(obj);
        if (pr == kNoMatch)
            return kNoMatchValue;
    
        const auto & v = obj.particles[pr].match_overlaps;
        if (v.empty())
            return kNoMatchValue;
    
        float best_val = v[0];
    
        for (size_t i = 1; i < v.size(); ++i) {
            if (v[i] > best_val)
                best_val = v[i];
        }
    
        return static_cast<double>(best_val);
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, leadingpion_match_overlap_max, leadingpion_match_overlap_max);

    template<class T>
    double second_leadingpion_match_id_max(const T & obj)
    {
        size_t pr = selectors::second_leading_pion(obj);
        if (pr == kNoMatch)
            return kNoMatchValue;
    
        const auto & v = obj.particles[pr].match_overlaps;
        if (v.empty())
            return kNoMatchValue;
    
        // Manual search for max value and its index
        size_t best_idx = 0;
        float  best_val = v[0];
    
        for (size_t i = 1; i < v.size(); ++i) {
            if (v[i] > best_val) {
                best_val = v[i];
                best_idx = i;
            }
        }
    
        // If you want the index of the biggest overlap:
        return obj.particles[pr].match_ids[best_idx];
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, second_leadingpion_match_id_max, second_leadingpion_match_id_max);

    template<class T>
    double second_leadingpion_match_overlap_max(const T & obj)
    {
        size_t pr = selectors::second_leading_pion(obj);
        if (pr == kNoMatch)
            return kNoMatchValue;
    
        const auto & v = obj.particles[pr].match_overlaps;
        if (v.empty())
            return kNoMatchValue;
    
        float best_val = v[0];
    
        for (size_t i = 1; i < v.size(); ++i) {
            if (v[i] > best_val)
                best_val = v[i];
        }
    
        return static_cast<double>(best_val);
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, second_leadingpion_match_overlap_max, second_leadingpion_match_overlap_max);


    template<class T>
    double hadronic_invariant_mass(const T & obj)
    {
        if (!obj.hadronic_invariant_mass>0.0 || !std::isfinite(obj.hadronic_invariant_mass))
            return -999.0;
        else
            return obj.hadronic_invariant_mass*1000.;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::True, hadronic_invariant_mass, hadronic_invariant_mass);

    template<class T>
    double momentum_transfer(const T & obj)
    {
        if (!obj.momentum_transfer>0.0 || !std::isfinite(obj.momentum_transfer))
            return -999.0;
        else
            return obj.momentum_transfer*1000.;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::True, momentum_transfer, momentum_transfer);


    template<class T>
    double energy_transfer(const T & obj)
    {
        if (!obj.energy_transfer>0.0 || !std::isfinite(obj.energy_transfer))
            return -999.0;
        else
            return obj.energy_transfer*1000.;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::True, energy_transfer, energy_transfer);


    template<class T>
    double hadron_invariant_mass_hadron(const T & obj)
    {
        constexpr bool INCLUDE_PHOTONS = false; // toggle as needed

        auto include_by_pid = [&](int pid_cat) -> bool {
            if (pid_cat == 1 || pid_cat == 2) return false; // exclude e, μ
            if (pid_cat == 3 || pid_cat == 4) return true;  // include π±, p
            if (pid_cat == 0) return INCLUDE_PHOTONS;       // optional γ
            return false;
        };

        auto mass_GeV = [&](int pid_cat) -> double {
            switch (pid_cat) {
                case 0: return 0.0;             // γ
                case 1: return 0.000511;        // e±
                case 2: return 0.105658;        // μ±
                case 3: return 0.139570;        // π±
                case 4: return 0.938272;        // p
                default: return -1.0;           // unknown
            }
        };

        double E_sum=0, px_sum=0, py_sum=0, pz_sum=0;
        size_t n_used=0;

        for (const auto & p : obj.particles)
        {
            if (!pcuts::final_state_signal(p)) continue;

            int pid_cat = static_cast<int>(pvars::pid(p));  // cast double→int
            if (!include_by_pid(pid_cat)) continue;

            double m = mass_GeV(pid_cat);
            if (m < 0.0) continue;

            double KE_GeV = pvars::ke(p) / 1000.0; // MeV → GeV
            double E      = KE_GeV + m;

            // normalize direction
            double ux=p.start_dir[0], uy=p.start_dir[1], uz=p.start_dir[2];
            double unorm = std::sqrt(ux*ux+uy*uy+uz*uz);
            if (unorm > 0.0) { ux/=unorm; uy/=unorm; uz/=unorm; }

            double p2 = E*E - m*m;
            if (p2 < 0.0) p2 = 0.0;
            double pmag = std::sqrt(p2);

            double px = pmag*ux, py = pmag*uy, pz = pmag*uz;

            E_sum  += E;
            px_sum += px;
            py_sum += py;
            pz_sum += pz;
            ++n_used;
        }

        if (n_used == 0) return kNoMatchValue;

        double W2 = E_sum*E_sum - (px_sum*px_sum + py_sum*py_sum + pz_sum*pz_sum);
        return (W2 > 0.0) ? std::sqrt(W2) : kNoMatchValue;
    }

    REGISTER_VAR_SCOPE(RegistrationScope::Both, hadron_invariant_mass_hadron, hadron_invariant_mass_hadron);

    template<class T>
        double kink_cosangle_from_points_lpi(const T & obj)
        {
            const size_t pi = selectors::leading_pion(obj);
            if(pi == kNoMatch) return kNoMatchValue;

            const auto & p = obj.particles[pi];

            // vertex position
            utilities::three_vector vtx = {obj.vertex[0],obj.vertex[1],obj.vertex[2]};

            // pion end point
            utilities::three_vector end = {pvars::end_x(p),pvars::end_y(p),pvars::end_z(p)};


            // v = vertex -> end
            utilities::three_vector v = utilities::subtract(end, vtx);

            const double vmag = utilities::magnitude(v);
            if(vmag <= 0) return PLACEHOLDERVALUE;

            // u = start_dir
            utilities::three_vector u = {p.start_dir[0],p.start_dir[1],p.start_dir[2]};

            const double umag = utilities::magnitude(u);
            if(umag <= 0) return PLACEHOLDERVALUE;

            // cos(angle) = (u · v) / (|u||v|)
            return utilities::dot_product(u, v) / (umag * vmag);
        }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, kink_cosangle_from_points_lpi, kink_cosangle_from_points_lpi);

    template<class T>
        double cos_wrt_beam_by_vtx_pi(const T & obj)
        {
            const size_t pi = selectors::leading_pion(obj);
            if(pi == kNoMatch) return kNoMatchValue;

            const auto & p = obj.particles[pi];

            // vertex position
            utilities::three_vector vtx = {obj.vertex[0],obj.vertex[1],obj.vertex[2]};

            // pion end point
            utilities::three_vector end = {pvars::end_x(p),pvars::end_y(p),pvars::end_z(p)};
            utilities::three_vector v = utilities::subtract(end, vtx);
            utilities::three_vector v_unit = utilities::normalize(v);

            utilities::three_vector beam_unit = utilities::numi_beam_direction(vtx);
            return utilities::dot_product(v_unit, beam_unit);
        }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, cos_wrt_beam_by_vtx_pi, cos_wrt_beam_by_vtx_pi);
    template<class T>
        double leading_pion_candidate_bkg_parent_ke(const T & obj)
        {
            int64_t id = -99;
            for (size_t i = 0; i < obj.particles.size(); ++i)
            {
                if (obj.particles[i].is_primary && obj.particles[i].pdg_code==2212)
                {
                    id = i;
                    break;
                }
            }
            if (id==-99)
                return -999.0;
            else
                return obj.particles[id].ke;
        }
    REGISTER_VAR_SCOPE(RegistrationScope::True, leading_pion_candidate_bkg_parent_ke, leading_pion_candidate_bkg_parent_ke);
    template<class T>
        double leading_pion_candidate_is_pion_bkg_parent_pdg(const T & obj)
        {
            int leading_idx = -1;        // vector index
            int nonprimleading_idx = -1; // vector index

            // 1) leading 찾기 (조건만 만족하면 갱신, KE 비교 없음)
            for (int i = 0; i < (int)obj.particles.size(); ++i)
            {
                const auto& p = obj.particles[i];

                // 제외할 PDG
                bool is_excluded =
                    (p.pdg_code == 22) || (p.pdg_code == 2112) || (p.pdg_code == 2212) ||
                    (std::abs(p.pdg_code) == 13) || (std::abs(p.pdg_code) == 211) ||
                    (std::abs(p.pdg_code) == 11);

                if (is_excluded) continue;

                if (p.is_primary)
                {
                    std::cout << "Primary: PDG=" << p.pdg_code
                        << " id=" << p.id
                        << " ke=" << p.ke
                        << std::endl;

                    leading_idx = i;   // ✅ index 저장
                }
                else
                {
                    std::cout << "NonPrimary: PDG=" << p.pdg_code
                        << " id=" << p.id
                        << " ke=" << p.ke
                        << std::endl;

                    nonprimleading_idx = i; // ✅ index 저장
                }
            }

            // 2) particle.id -> vector index 맵
            std::unordered_map<int, int> id2idx;
            id2idx.reserve(obj.particles.size());
            for (int idx = 0; idx < (int)obj.particles.size(); ++idx)
                id2idx[obj.particles[idx].id] = idx;

            // 3) children 출력 함수
            auto print_children = [&](int lidx, const char* tag)
            {
                if (lidx < 0 || lidx >= (int)obj.particles.size())
                    return;

                const auto& lead = obj.particles[lidx];
                const auto& children = lead.children_id; // particle.id 들

                std::cout << tag
                    << " Leading idx=" << lidx
                    << " | id=" << lead.id
                    << " | PDG=" << lead.pdg_code
                    << " | #children=" << children.size()
                    << std::endl;

                for (size_t i = 0; i < children.size(); ++i)
                {
                    int cid = children[i]; // child particle.id
                    auto it = id2idx.find(cid);

                    if (it == id2idx.end())
                    {
                        std::cout << "  Child " << i
                            << " | id=" << cid
                            << " | NOT FOUND"
                            << std::endl;
                        continue;
                    }

                    int cidx = it->second;
                    const auto& child = obj.particles[cidx];

                    std::cout << "  Child " << i
                        << " | id=" << cid
                        << " | idx=" << cidx
                        << " | PDG=" << child.pdg_code
                        << " | KE=" << child.ke
                        << std::endl;
                }
            };

            // 4) 출력
            if (leading_idx != -1)
                print_children(leading_idx, "Prim");

            if (nonprimleading_idx != -1)
                print_children(nonprimleading_idx, "NonPrim");

            return -999;
        }
    REGISTER_VAR_SCOPE(RegistrationScope::True, leading_pion_candidate_is_pion_bkg_parent_pdg, leading_pion_candidate_is_pion_bkg_parent_pdg);

    template<class T>
        int has_strange_baryon_with_pion_child(const T& obj)
        {
            // 1) strange baryon PDG set
            static const std::unordered_set<int> strange_baryons = {
                3122,
                3222, 3212, 3112,
                3224, 3214, 3114,
                3322, 3312,
                3324, 3314,
                3334
            };

            // 2) particle.id -> vector index 맵
            std::unordered_map<int, int> id2idx;
            id2idx.reserve(obj.particles.size());
            for (int idx = 0; idx < (int)obj.particles.size(); ++idx)
                id2idx[obj.particles[idx].id] = idx;

            // 3) loop over all particles: find strange baryon
            for (int i = 0; i < (int)obj.particles.size(); ++i)
            {
                const auto& p = obj.particles[i];

                if (strange_baryons.count(std::abs(p.pdg_code)) == 0)
                    continue;

                // 4) check its children
                for (int child_id : p.children_id)
                {
                    auto it = id2idx.find(child_id);
                    if (it == id2idx.end())
                        continue;

                    const auto& child = obj.particles[it->second];

                    if (std::abs(child.pdg_code) == 211)
                    {
                        // found π± child from strange baryon
                        return 1;
                    }
                }
            }

            // none found
            return 0;
        }
    REGISTER_VAR_SCOPE(RegistrationScope::True, has_strange_baryon_with_pion_child, has_strange_baryon_with_pion_child);

    template<class T>
        int has_kaon_with_pion_child(const T& obj)
        {
            // Kaons: K± (321), K0 / K0bar (311)
            static const std::unordered_set<int> kaons = {
                321,  // K+
                311   // K0
            };

            // particle.id -> vector index
            std::unordered_map<int, int> id2idx;
            id2idx.reserve(obj.particles.size());
            for (int idx = 0; idx < (int)obj.particles.size(); ++idx)
                id2idx[obj.particles[idx].id] = idx;

            // find kaon and check its children
            for (int i = 0; i < (int)obj.particles.size(); ++i)
            {
                const auto& p = obj.particles[i];

                if (kaons.count(std::abs(p.pdg_code)) == 0)
                    continue;

                for (int child_id : p.children_id)
                {
                    auto it = id2idx.find(child_id);
                    if (it == id2idx.end())
                        continue;

                    const auto& child = obj.particles[it->second];

                    if (std::abs(child.pdg_code) == 211)
                        return 1;
                }
            }

            return 0;
        }

    REGISTER_VAR_SCOPE(RegistrationScope::True, has_kaon_with_pion_child, has_kaon_with_pion_child);
    template<class T>
        double strange_baryon_pion_child_ke(const T& obj)
        {
            // Strange baryons
            static const std::unordered_set<int> strange_baryons = {
                3122,
                3222, 3212, 3112,
                3224, 3214, 3114,
                3322, 3312,
                3324, 3314,
                3334
            };

            // particle.id -> vector index
            std::unordered_map<int, int> id2idx;
            id2idx.reserve(obj.particles.size());
            for (int idx = 0; idx < (int)obj.particles.size(); ++idx)
                id2idx[obj.particles[idx].id] = idx;

            // loop over particles
            for (int i = 0; i < (int)obj.particles.size(); ++i)
            {
                const auto& p = obj.particles[i];

                if (strange_baryons.count(std::abs(p.pdg_code)) == 0)
                    continue;

                // check children
                for (int child_id : p.children_id)
                {
                    auto it = id2idx.find(child_id);
                    if (it == id2idx.end())
                        continue;

                    const auto& child = obj.particles[it->second];

                    if (std::abs(child.pdg_code) == 211)
                    {
                        // return pion true KE
                        return child.ke;
                    }
                }
            }

            return -999.0;
        }
    REGISTER_VAR_SCOPE(RegistrationScope::True, strange_baryon_pion_child_ke, strange_baryon_pion_child_ke);
    template<class T>
        double kaon_pion_child_ke(const T& obj)
        {
            // Kaons
            static const std::unordered_set<int> kaons = {
                321, // K±
                311  // K0 / K0bar
            };

            // particle.id -> vector index
            std::unordered_map<int, int> id2idx;
            id2idx.reserve(obj.particles.size());
            for (int idx = 0; idx < (int)obj.particles.size(); ++idx)
                id2idx[obj.particles[idx].id] = idx;

            // loop over particles
            for (int i = 0; i < (int)obj.particles.size(); ++i)
            {
                const auto& p = obj.particles[i];

                if (kaons.count(std::abs(p.pdg_code)) == 0)
                    continue;

                // check children
                for (int child_id : p.children_id)
                {
                    auto it = id2idx.find(child_id);
                    if (it == id2idx.end())
                        continue;

                    const auto& child = obj.particles[it->second];

                    if (std::abs(child.pdg_code) == 211)
                    {
                        // return pion true KE
                        return child.ke;
                    }
                }
            }

            return -999.0;
        }

    REGISTER_VAR_SCOPE(RegistrationScope::True, kaon_pion_child_ke, kaon_pion_child_ke);
    template<class T>
        int nonprimary_pion_source_ancestor_pdg(const T& obj)
        {
            int DEFAULT_NONE     = -999;
            int DEFAULT_AMBIG    = -888;
            int DEFAULT_PRIMARY  = -777;
            int DEFAULT_BADCHAIN = -666;
            int MAX_DEPTH        = 50;

            // id -> index
            std::unordered_map<int, int> id2idx;
            id2idx.reserve(obj.particles.size());
            for (int idx = 0; idx < (int)obj.particles.size(); ++idx)
                id2idx[obj.particles[idx].id] = idx;

            // primary π 있으면 바로 종료
            for (int i = 0; i < (int)obj.particles.size(); ++i)
            {
                const auto& p = obj.particles[i];
                if (p.is_primary && std::abs(p.pdg_code) == 211)
                    return DEFAULT_PRIMARY;
            }

            auto first_nonpion_ancestor_pdg = [&](int pion_idx) -> int
            {
                std::unordered_set<int> visited_ids;
                visited_ids.reserve(16);

                int current_idx = pion_idx;
                for (int depth = 0; depth < MAX_DEPTH; ++depth)
                {
                    const auto& cur = obj.particles[current_idx];

                    if (visited_ids.count(cur.id)) return DEFAULT_BADCHAIN;
                    visited_ids.insert(cur.id);

                    auto it = id2idx.find(cur.parent_id);
                    if (it == id2idx.end()) return DEFAULT_BADCHAIN;

                    int pidx = it->second;
                    const auto& par = obj.particles[pidx];

                    if (std::abs(par.pdg_code) == 211)
                    {
                        current_idx = pidx; // parent도 π면 계속 올라감
                        continue;
                    }

                    return par.pdg_code; // 최초 non-π ancestor
                }

                return DEFAULT_BADCHAIN; // depth 초과
            };

            int agreed_source_pdg = 0;
            bool found_valid_source = false;
            bool saw_nonprimary_pion = false;

            // <<< BADCHAIN fallback용: "직접 parent_pdg"를 하나라도 잡아두기
            int fallback_parent_pdg = DEFAULT_BADCHAIN;     // <<<
            bool has_fallback_parent = false;               // <<<

            for (int i = 0; i < (int)obj.particles.size(); ++i)
            {
                const auto& p = obj.particles[i];
                if (p.is_primary) continue;
                if (std::abs(p.pdg_code) != 211) continue;

                saw_nonprimary_pion = true;

                int src_pdg = first_nonpion_ancestor_pdg(i);
                if (src_pdg == DEFAULT_BADCHAIN)
                {
                    // <<< 여기서 "그 파이온의 parent_pdg" fallback 시도
                    auto itp = id2idx.find(p.parent_id);
                    if (itp != id2idx.end())
                    {
                        const auto& par = obj.particles[itp->second];
                        fallback_parent_pdg = par.parent_pdg_code;
                        has_fallback_parent = true;
                    }
                    continue; // 유효한 source 못 찾음(일단 스킵)
                }

                if (!found_valid_source)
                {
                    agreed_source_pdg = src_pdg;
                    found_valid_source = true;
                }
                else if (src_pdg != agreed_source_pdg)
                {
                    return DEFAULT_AMBIG;
                }
            }

            if (!saw_nonprimary_pion) return DEFAULT_NONE;

            // <<< 핵심 변경: source를 못 찾았으면 parent_pdg로 fallback
            if (!found_valid_source)
            {
                if (has_fallback_parent) return fallback_parent_pdg;  // <<< parent_pdg 리턴
                return DEFAULT_BADCHAIN;                               // <<< parent도 못 찾으면 어쩔 수 없음
            }

            return agreed_source_pdg;
        }
    REGISTER_VAR_SCOPE(RegistrationScope::True, nonprimary_pion_source_ancestor_pdg, nonprimary_pion_source_ancestor_pdg);

    template<class T>
        int count_michel_tags(const T& obj, std::vector<double> params = {25.0, 20.0})
        {
            // Robust params handling
            double pion_ke_thr = 25.0;
            double michel_dist = 20.0;
            if (params.size() >= 1) pion_ke_thr = params[0];
            if (params.size() >= 2) michel_dist = params[1];

            std::vector<size_t> pion_ids;
            std::vector<utilities::three_vector> pion_ends;

            for (size_t i = 0; i < obj.particles.size(); ++i)
            {
                const auto& p = obj.particles[i];

                if (pvars::pid(p) == pvars::kPion &&
                        pvars::primary_classification(p) &&
                        pvars::ke(p) >= pion_ke_thr)
                {
                    pion_ids.push_back(p.id);
                    pion_ends.push_back({pvars::end_x(p), pvars::end_y(p), pvars::end_z(p)});
                }
            }

            if (pion_ids.empty()) return 0;

            int nmichel = 0;

            for (size_t i = 0; i < obj.particles.size(); ++i)
            {
                const auto& p = obj.particles[i];

                if (std::find(pion_ids.begin(), pion_ids.end(), p.id) != pion_ids.end())
                    continue;

                if (p.shape != 2) continue;

                utilities::three_vector particle_vtx = {pvars::start_x(p), pvars::start_y(p), pvars::start_z(p)};

                for (const auto& pe : pion_ends)
                {
                    const double d = utilities::magnitude(utilities::subtract(pe, particle_vtx));
                    if (d < michel_dist) { ++nmichel; break; }
                }
            }

            return nmichel;
        }
    REGISTER_VAR_SCOPE(RegistrationScope::Reco, count_michel_tags, count_michel_tags);
    template<class T>
        bool has_michel_tag(
                const T& obj,
                std::vector<double> params = {25.0, 20.0}
                )
        {
            std::cout<<count_michel_tags(obj, params)<<std::endl;
            return count_michel_tags(obj, params) >= 1;
        }

    REGISTER_VAR_SCOPE(RegistrationScope::Both, has_michel_tag, has_michel_tag);
    template<class T>
        bool has_multi_michel_tag(
                const T& obj,
                std::vector<double> params = {25.0, 20.0}
                )
        {
            return count_michel_tags(obj, params) >= 2;
        }

    REGISTER_VAR_SCOPE(RegistrationScope::Both, has_multi_michel_tag, has_multi_michel_tag);




    template<class T>
    double opening_angle_pion(const T & obj)
    {
        size_t pi = selectors::leading_pion(obj);
        size_t spi = selectors::second_leading_pion(obj);
        if(pi == kNoMatch || spi == kNoMatch)
            return kNoMatchValue; // No leading muon or proton found.
        else
        {
            auto & p(obj.particles[pi]);
            auto & sp(obj.particles[spi]);
            return (p.start_dir[0] * sp.start_dir[0] + p.start_dir[1] * sp.start_dir[1] + p.start_dir[2] * sp.start_dir[2]);
        }
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, opening_angle_pion, opening_angle_pion);

    template<class T>
    double opening_angle_pion_proton(const T & obj)
    {
        size_t pi = selectors::leading_pion(obj);
        size_t spi = selectors::leading_proton(obj);
        if(pi == kNoMatch || spi == kNoMatch)
            return kNoMatchValue; // No leading muon or proton found.
        else
        {
            auto & p(obj.particles[pi]);
            auto & sp(obj.particles[spi]);
            return (p.start_dir[0] * sp.start_dir[0] + p.start_dir[1] * sp.start_dir[1] + p.start_dir[2] * sp.start_dir[2]);
        }
    }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, opening_angle_pion_proton, opening_angle_pion_proton);

    template<class T>
    double Q2_exp(const T & obj)
    {
        double true_nuE = visible_energy_calosub(obj)*1000.0;
        double costheta = muon_cos_angle(obj);

        size_t mi = selectors::leading_muon(obj);
        if(mi == kNoMatch)
            return kNoMatchValue; // No leading muon or proton found.
        else
        {
            auto & m(obj.particles[mi]);
            utilities::three_vector momentum = {pvars::px(m), pvars::py(m), pvars::pz(m)};
            double mu_mom = utilities::magnitude(momentum);
            double mu_E = pvars::ke(m) + MUON_MASS;
            utilities::three_vector nu_mom = {obj.momentum[0],obj.momentum[1],obj.momentum[2]};
            utilities::three_vector nu_mom_unit = utilities::normalize(nu_mom);
            utilities::three_vector mu_mom_unit = utilities::normalize(momentum);
//            std::cout<<"costheta_check:"<<utilities::dot_product(mu_mom_unit, nu_mom_unit)<<"= costheta:"<<costheta<<std::endl;
//
//            std::cout<<"mu_E:"<<mu_E<<std::endl;
//            std::cout<<"mu_mom:"<<mu_mom<<std::endl;
//            std::cout<<"mu_E**2-mom**2:"<< mu_E*mu_E-mu_mom*mu_mom<<" = = "<<MUON_MASS<<std::endl;
            //            std::cout<<mu_E-mu_mom*costheta<<std::endl;
            return (2*true_nuE*(mu_E-mu_mom*costheta)-std::pow(MUON_MASS, 2))/1000000.;
        }
    }
    REGISTER_VAR_SCOPE(RegistrationScope::True, Q2_exp, Q2_exp);

    template<class T>
    double W_exp(const T & obj)
    {
        double true_nuE = visible_energy(obj)*1000.0;
        size_t mi = selectors::leading_muon(obj);
        if(mi == kNoMatch)
            return -999.;
        else
        {
            auto & m(obj.particles[mi]);
            double mu_E = pvars::ke(m) + MUON_MASS;
            double W_2 = std::pow(NEUTRON_MASS, 2)+2*NEUTRON_MASS*(true_nuE-mu_E)-four_mom_transfer_squared(obj);
 //           std::cout<<W_2<<std::endl;
            if (W_2<0)
            {
                return -999.;
            }
            else
            {
                return std::sqrt(W_2);
            }
        }
    }
    REGISTER_VAR_SCOPE(RegistrationScope::True, W_exp, W_exp);


    template<class T>
        double muon_pT_wrt_beam(const T & obj)
        {
            for(const auto & p : obj.particles)
            {
                if(!pcuts::final_state_signal(p)) continue;

                if(pvars::pid(p) == pvars::kMuon)
                {
                    utilities::three_vector mom =
                        std::make_tuple(pvars::px(p), pvars::py(p), pvars::pz(p));

                    utilities::three_vector vtx =
                        std::make_tuple(pvars::start_x(p), pvars::start_y(p), pvars::start_z(p));

                    // transverse momentum vector wrt beam (BNB:z or NuMI:target->vtx)
                    utilities::three_vector pt_vec = utilities::transverse_momentum(mom, vtx);

                    return utilities::magnitude(pt_vec);
                }
            }
            return kNoMatchValue;
        }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, muon_pT_wrt_beam, muon_pT_wrt_beam);





    template<class T>
        double p_pi_invariant_mass(const T & obj)
        {
            size_t pi = selectors::leading_pion(obj);
            size_t pr = selectors::leading_proton(obj);

            if(pi == kNoMatch || pr == kNoMatch)
                return kNoMatchValue;

            auto & p   = obj.particles[pr];
            auto & pip = obj.particles[pi];

            double E_p  = pvars::ke(p)   + PROTON_MASS;
            double E_pi = pvars::ke(pip) + PION_MASS;

            double px = pvars::px(p) + pvars::px(pip);
            double py = pvars::py(p) + pvars::py(pip);
            double pz = pvars::pz(p) + pvars::pz(pip);

            double m2 = std::pow(E_p + E_pi, 2) - (px*px + py*py + pz*pz);

            if(m2 < 0) m2 = 0;

            return std::sqrt(m2); // MeV, or divide by 1000 for GeV
        }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, p_pi_invariant_mass, p_pi_invariant_mass);

    template<class T>
        double cos_angle_proton_pion(const T & obj)
        {
            size_t pi = selectors::leading_pion(obj);
            size_t pr = selectors::leading_proton(obj);

            if(pi == kNoMatch || pr == kNoMatch)
                return kNoMatchValue;
            else
            {
                auto & p  = obj.particles[pr];
                auto & pip = obj.particles[pi];

                utilities::three_vector p_vec  = {pvars::px(p),  pvars::py(p),  pvars::pz(p)};
                utilities::three_vector pi_vec = {pvars::px(pip), pvars::py(pip), pvars::pz(pip)};

                double p_mag  = utilities::magnitude(p_vec);
                double pi_mag = utilities::magnitude(pi_vec);

                return utilities::dot_product(p_vec, pi_vec)/(p_mag*pi_mag);
            }
        }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, cos_angle_proton_pion, cos_angle_proton_pion);

    template<class T>
        double p_allpi_invariant_mass(const T & obj)
        {
            size_t pr = selectors::leading_proton(obj);

            if(pr == kNoMatch)
                return kNoMatchValue;

            auto & p = obj.particles[pr];

            double E_tot = pvars::ke(p) + PROTON_MASS;
            double px_tot = pvars::px(p);
            double py_tot = pvars::py(p);
            double pz_tot = pvars::pz(p);

            size_t n_pions = 0;

            for(const auto & part : obj.particles)
            {
                if(!pcuts::final_state_signal(part)) continue;

                int pid = pvars::pid(part);

                if(pid == pvars::kPion)
                {
                    E_tot  += pvars::ke(part) + PION_MASS;
                    px_tot += pvars::px(part);
                    py_tot += pvars::py(part);
                    pz_tot += pvars::pz(part);
                    ++n_pions;
                }
            }

            if(n_pions == 0)
                return kNoMatchValue;

            double m2 = E_tot*E_tot - (px_tot*px_tot + py_tot*py_tot + pz_tot*pz_tot);

            if(m2 < 0) m2 = 0;

            return std::sqrt(m2); // MeV; divide by 1000. for GeV
        }
    REGISTER_VAR_SCOPE(RegistrationScope::Both, p_allpi_invariant_mass, p_allpi_invariant_mass);



    template<class T>
        double muon_ke_from_pion(const T & obj)
        {
            size_t pi = selectors::leading_pion(obj);

            if(pi == kNoMatch)
                return kNoMatchValue;

            auto & pion = obj.particles[pi];
            int pion_id = pion.id;
            double muonke=-999.0;
//            std::cout<<"Leading Pion ID : "<<pion.id<<std::endl;


            for(const auto & part : obj.particles)
            {
//                std::cout<<"particle ID : "<<part.id<<", its pdg : "<<part.pdg_code<<",its parent id : "<<part.parent_id<<std::endl;

                if(part.parent_id==pion_id && (part.pdg_code==-13||part.pdg_code==13))
                {
                    muonke=part.ke;
//                    std::cout<<"Muon E : "<<part.ke<<",Contained:"<<part.is_contained<<std::endl;
                }
            }


            return muonke; // MeV; divide by 1000. for GeV
        }
    REGISTER_VAR_SCOPE(RegistrationScope::True, muon_ke_from_pion, muon_ke_from_pion);

    template<class T>
        bool is_muon_from_pion_contained(const T & obj)
        {
            size_t pi = selectors::leading_pion(obj);

            if(pi == kNoMatch)
                return kNoMatchValue;

            auto & pion = obj.particles[pi];
            int pion_id = pion.id;
            double muonke=-999.0;
            bool is_contained = false;
//            std::cout<<"Leading Pion ID : "<<pion.id<<std::endl;


            for(const auto & part : obj.particles)
            {
//                std::cout<<"particle ID : "<<part.id<<", its pdg : "<<part.pdg_code<<",its parent id : "<<part.parent_id<<std::endl;

                if(part.parent_id==pion_id && (part.pdg_code==-13||part.pdg_code==13))
                {
                    muonke=part.ke;
                    is_contained=part.is_contained;
//                    std::cout<<"Muon E : "<<part.ke<<",Contained:"<<part.is_contained<<std::endl;
                }
            }


            return is_contained; // MeV; divide by 1000. for GeV
        }
    REGISTER_VAR_SCOPE(RegistrationScope::True, is_muon_from_pion_contained, is_muon_from_pion_contained);



    template<class T>
        double count_proton(const T & obj)
        {
            std::set<int64_t> seen_tracks;
            int n_protons = 0;
            for(const auto & part : obj.particles)
            {
//                std::cout<< "particle ID : " << part.id<< ", track ID : " << part.track_id<< ", pdg : " << part.pdg_code<< ", parent id : " << part.parent_id<< std::endl;
                // proton only
                if(part.pdg_code != 2212)
                    continue;
                // invalid track id protection
                if(part.track_id < 0)
                    continue;
                // avoid duplicated propagated particles
                if(seen_tracks.count(part.track_id))
                    continue;
                seen_tracks.insert(part.track_id);
                n_protons++;
            }
//            std::cout<<"proton Number:"<<n_protons<<std::endl;
            return n_protons;
        }
    REGISTER_VAR_SCOPE(RegistrationScope::True, count_proton,count_proton);

    template<class T>
        double count_pip(const T & obj)
        {
            std::set<int64_t> seen_tracks;
            int n_protons = 0;
            for(const auto & part : obj.particles)
            {
            //    std::cout<< "particle ID : " << part.id<< ", track ID : " << part.track_id<< ", pdg : " << part.pdg_code<< ", parent id : " << part.parent_id<< std::endl;
                // proton only
                if(part.pdg_code != 211)
                    continue;
                // invalid track id protection
                if(part.track_id < 0)
                    continue;
                // avoid duplicated propagated particles
                if(seen_tracks.count(part.track_id))
                    continue;
                seen_tracks.insert(part.track_id);
                n_protons++;
            }
//            std::cout<<"pip Number:"<<n_protons<<std::endl;
            return n_protons;
        }
    REGISTER_VAR_SCOPE(RegistrationScope::True, count_pip,count_pip);

    template<class T>
        double count_pim(const T & obj)
        {
            std::set<int64_t> seen_tracks;
            int n_protons = 0;
            for(const auto & part : obj.particles)
            {
            //    std::cout<< "particle ID : " << part.id<< ", track ID : " << part.track_id<< ", pdg : " << part.pdg_code<< ", parent id : " << part.parent_id<< std::endl;
                // proton only
                if(part.pdg_code != -211)
                    continue;
                // invalid track id protection
                if(part.track_id < 0)
                    continue;
                // avoid duplicated propagated particles
                if(seen_tracks.count(part.track_id))
                    continue;
                seen_tracks.insert(part.track_id);
                n_protons++;
            }
//            std::cout<<"pim Number:"<<n_protons<<std::endl;
            return n_protons;
        }
    REGISTER_VAR_SCOPE(RegistrationScope::True, count_pim,count_pim);
}




#endif // VARIABLES_H
