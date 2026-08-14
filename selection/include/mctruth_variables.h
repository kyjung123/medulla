/**
 * @file mctruth_variables.h
 * @brief Definitions of analysis variables which can extract information from
 * the SRTrueInteraction object.
 * @details This file contains definitions of analysis variables which can be
 * used to extract information from the SRTrueInteraction object. Each variable
 * is implemented as a function which takes an SRTrueInteraction object as an
 * argument and returns a double. The association of an SRInteractionTruthDLP
 * object to an SRTrueInteraction object is handled upstream in the SpineVar
 * functions.
 * @author mueller@fnal.gov
 * @author rvizarr@fnal.gov
 * @author kyjung
 */
#ifndef MCTRUTH_VARIABLES_H
#define MCTRUTH_VARIABLES_H
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "sbnanaobj/StandardRecord/SRTrueInteraction.h"
#include "sbnanaobj/StandardRecord/SRVector3D.h"

#include "framework.h"
#include "utilities.h"

/**
 * @namespace mctruth
 * @brief Namespace for organizing variables which act on true interactions.
 * @details This namespace is intended to be used for organizing variables
 * which act on true interactions. Each variable is implemented as a function
 * which takes an SRTrueInteraction object as an argument and returns a double.
 */
namespace mctruth
{
    /**
     * @brief Variable for the true neutrino energy.
     * @details This variable is intended to provide the true energy of the
     * parent neutrino that produced the interaction.
     * @tparam T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the true neutrino energy.
     */
    template<typename T>
        double neutrino_energy(const T & obj) { return obj.E; }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, neutrino_energy, neutrino_energy);

   /**
     * @brief Variable for the true interaction energy transfer.
     * @details This variable is intended to provide the true energy
     * transfer from the neutrino to the hadronic system. This is
     * defined in the lab frame.
     * @tparam T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the true energy transfer into the hadronic system
     * in the lab frame.
     */
    template<typename T>
        double energy_transfer(const T & obj) { return obj.q0_lab; }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, energy_transfer, energy_transfer);

    /**
     * @brief Variable for the true four-momentum transfer squared Q², with
     * validity guards.
     * @details This variable is intended to provide the true Q² of the
     * interaction, guarding against non-physical (non-positive) or
     * non-finite values by returning a sentinel value of -999.0 in those
     * cases. This differs from @ref neutrino_Q2, which returns the raw
     * generator-record value with no validity check.
     * @tparam T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return Q² in GeV², or -999.0 if Q² is non-positive or non-finite.
     */
    template<typename T>
        double q_squared(const T & obj)
        {
        if (!(obj.Q2 > 0.0) || !std::isfinite(obj.Q2))
            return -999.0;
        else
            return obj.Q2;
        }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, q_squared, q_squared);

    /**
     * @brief Variable for the true neutrino baseline.
     * @details This variable is intended to provide the true baseline of the
     * parent neutrino that produced the interaction.
     * @tparam T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the true neutrino baseline.
     */
    template<typename T>
        double baseline(const T & obj) { return obj.baseline; }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, baseline, baseline);

    /**
     * @brief Variable for the true neutrino PDG code.
     * @details This variable is intended to provide the true PDG code of the
     * parent neutrino that produced the interaction.
     * @tparam T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the true neutrino PDG code.
     */
    template<typename T>
        double pdg(const T & obj) { return obj.pdg; }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, pdg, pdg);

    /**
     * @brief Variable for the PDG code of the parent of the neutrino.
     * @details This variable is intended to provide the PDG code of the
     * parent of the neutrino that produced the interaction.
     * @tparam T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the PDG code of the parent of the neutrino.
     */
    template<typename T>
        double parent_pdg(const T & obj) { return obj.parent_pdg; }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, parent_pdg, parent_pdg);

    /**
     * @brief Variable for the true neutrino current value.
     * @details This variable is intended to provide the true current value of
     * the parent neutrino that produced the interaction.
     * @tparam T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the true neutrino current value.
     */
    template<typename T>
        double cc(const T & obj) { return obj.iscc; }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, cc, cc);

    /**
     * @brief Variable for the interaction mode of the interaction.
     * @details This variable is intended to provide the interaction mode of the
     * interaction. This is based on the GENIE interaction mode enumeration 
     * defined in the LArSoft MCNeutrino class.
     * @tparam T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the interaction mode.
     */
    template<typename T>
        double interaction_mode(const T & obj) { return obj.genie_mode; }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, interaction_mode, interaction_mode);

    /**
     * @brief Variable for the interaction type of the interaction.
     * @details This variable is intended to provide the interaction type of the
     * interaction. This is based on the GENIE interaction type enumeration 
     * defined in the LArSoft MCNeutrino class.
     * @param T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the interaction type.
     */
    template<typename T>
        double interaction_type(const T & obj) { return obj.genie_inttype; }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, interaction_type, interaction_type);

    /**
     * @brief Variable for the true off-axis angle of the neutrino.
     * @details This variable is intended to provide the true off-axis angle of
     * the parent neutrino that produced the interaction. The off-axis angle is
     * calculated as the angle between the neutrino momentum vector and the
     * beam axis (defined as the z-axis in both SBND and ICARUS).
     * @tparam T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the true off-axis angle of the neutrino in degrees.
     */
    template<typename T>
    double off_axis_angle(const T & obj)
    {
        const auto & neutrino_momentum = obj.momentum;
        double mag = std::sqrt(
            neutrino_momentum.x * neutrino_momentum.x +
            neutrino_momentum.y * neutrino_momentum.y +
            neutrino_momentum.z * neutrino_momentum.z
        );
        return 180./3.141592653589793 * std::acos(
            neutrino_momentum.z / mag
        );
    }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, off_axis_angle, off_axis_angle);

    /**
     * @brief True four-momentum transfer squared Q² in GeV².
     * @details Returns the raw generator-record value with no validity
     * check. See @ref q_squared for a version that guards against
     * non-physical or non-finite values.
     * @tparam T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return Q² from the generator record.
     */
    template<typename T>
    double neutrino_Q2(const T & obj) { return obj.Q2; }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, neutrino_Q2, neutrino_Q2);

    /**
     * @brief True hadronic invariant mass W in GeV.
     * @tparam T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return W from the generator record.
     */
    template<typename T>
    double neutrino_W(const T & obj) { return obj.w; }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, neutrino_W, neutrino_W);

    /**
     * @brief Count of true primary photons above an energy threshold.
     * @details Loops over `obj.prim` and counts PDG-22 particles whose kinetic
     * energy (= total energy, since photons are massless) exceeds `params[0]`
     * in MeV. Generator energies are stored in GeV.
     * @tparam T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @param params threshold in MeV; defaults to 0 MeV.
     * @return number of primary photons above threshold.
     */
    template<typename T>
    double nphotons_srtruth(const T & obj, std::vector<double> params={0.0,})
    {
        int n(0);
        for(const auto & p : obj.prim)
        {
            if(p.pdg == 22)
            {
                double ke = 1000. * p.genE; // MeV (massless)
                if(ke >= params[0]) ++n;
            }
        }
        return n;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, nphotons_srtruth, nphotons_srtruth);

    /**
     * @brief Count of true primary electrons above a kinetic energy threshold.
     * @details Loops over `obj.prim` and counts PDG-11 particles whose kinetic
     * energy exceeds `params[0]` in MeV.
     * @tparam T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @param params threshold in MeV; defaults to 0 MeV.
     * @return number of primary electrons above threshold.
     */
    template<typename T>
    double nelectrons_srtruth(const T & obj, std::vector<double> params={0.0,})
    {
        int n(0);
        for(const auto & p : obj.prim)
        {
            if(p.pdg == 11)
            {
                double ke = 1000. * (p.genE - ELECTRON_MASS/1000.); // MeV
                if(ke >= params[0]) ++n;
            }
        }
        return n;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, nelectrons_srtruth, nelectrons_srtruth);

    /**
     * @brief Count of true primary muons above a kinetic energy threshold.
     * @details Loops over `obj.prim` and counts PDG-13 particles (muons only,
     * not antimuons) whose kinetic energy exceeds `params[0]` in MeV.
     * @tparam T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @param params threshold in MeV; defaults to 0 MeV.
     * @return number of primary muons above threshold.
     */
    template<typename T>
    double nmuons_srtruth(const T & obj, std::vector<double> params={0.0,})
    {
        int n(0);
        for(const auto & p : obj.prim)
        {
            if(p.pdg == 13)
            {
                double ke = 1000. * (p.genE - MUON_MASS/1000.); // MeV
                if(ke >= params[0]) ++n;
            }
        }
        return n;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, nmuons_srtruth, nmuons_srtruth);

    /**
     * @brief Count of true primary neutral pions above a kinetic energy threshold.
     * @details Loops over `obj.prim` and counts PDG-111 particles whose kinetic
     * energy exceeds `params[0]` in MeV.
     * @tparam T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @param params threshold in MeV; defaults to 0 MeV.
     * @return number of primary neutral pions above threshold.
     */
    template<typename T>
    double npi0s_srtruth(const T & obj, std::vector<double> params={0.0,})
    {
        int n(0);
        for(const auto & p : obj.prim)
        {
            if(p.pdg == 111)
            {
                double ke = 1000. * (p.genE - PI0_MASS/1000.); // MeV
                if(ke >= params[0]) ++n;
            }
        }
        return n;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, npi0s_srtruth, npi0s_srtruth);

    /**
     * @brief Count of true primary charged pions above a kinetic energy threshold.
     * @details Loops over `obj.prim` and counts |PDG|=211 particles (π⁺ and π⁻)
     * whose kinetic energy exceeds `params[0]` in MeV.
     * @tparam T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @param params threshold in MeV; defaults to 0 MeV.
     * @return number of primary charged pions above threshold.
     */
    template<typename T>
    double npions_srtruth(const T & obj, std::vector<double> params={0.0,})
    {
        int n(0);
        for(const auto & p : obj.prim)
        {
            if(std::abs(p.pdg) == 211)
            {
                double ke = 1000. * (p.genE - PION_MASS/1000.); // MeV
                if(ke >= params[0]) ++n;
            }
        }
        return n;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, npions_srtruth, npions_srtruth);

    /**
     * @brief Count of true primary eta mesons above a kinetic energy threshold.
     * @details Loops over `obj.prim` and counts |PDG|=221 particles whose kinetic
     * energy exceeds `params[0]` in MeV. Eta mass: 547.862 MeV.
     * @tparam T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @param params threshold in MeV; defaults to 0 MeV.
     * @return number of primary eta mesons above threshold.
     */
    template<typename T>
    double netas_srtruth(const T & obj, std::vector<double> params={0.0,})
    {
        int n(0);
        for(const auto & p : obj.prim)
        {
            if(std::abs(p.pdg) == 221)
            {
                double ke = 1000. * (p.genE - 547.862/1000.); // MeV
                if(ke >= params[0]) ++n;
            }
        }
        return n;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, netas_srtruth, netas_srtruth);

    /**
     * @brief Count of true primary protons above a kinetic energy threshold.
     * @details Loops over `obj.prim` and counts PDG-2212 particles whose kinetic
     * energy exceeds `params[0]` in MeV.
     * @tparam T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @param params threshold in MeV; defaults to 0 MeV.
     * @return number of primary protons above threshold.
     */
    template<typename T>
    double nprotons_srtruth(const T & obj, std::vector<double> params={0.0,})
    {
        int n(0);
        for(const auto & p : obj.prim)
        {
            if(p.pdg == 2212)
            {
                double ke = 1000. * (p.genE - PROTON_MASS/1000.); // MeV
                if(ke >= params[0]) ++n;
            }
        }
        return n;
    }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, nprotons_srtruth, nprotons_srtruth);

    // -------------------------------------------------------------------
    // Variables merged in from the legacy mctruth.h (kyjung's ccpionproton
    // branch). Kept as-is aside from the q_squared fix above.
    // -------------------------------------------------------------------

    template<typename T>
        int nneutron_bf_FSI(const T & obj) { return obj.nneutron; }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, nneutron_bf_FSI, nneutron_bf_FSI);

    template<typename T>
        int npiplus_bf_FSI(const T & obj) { return obj.npiplus; }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, npiplus_bf_FSI, npiplus_bf_FSI);

    template<typename T>
        int npiminus_bf_FSI(const T & obj) { return obj.npiminus; }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, npiminus_bf_FSI, npiminus_bf_FSI);

    template<typename T>
        int npizero_bf_FSI(const T & obj) { return obj.npizero; }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, npizero_bf_FSI, npizero_bf_FSI);

    template<typename T>
        int nproton_bf_FSI(const T & obj) { return obj.nproton; }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, nproton_bf_FSI, nproton_bf_FSI);

    template<typename T>
        double nneutron_G4(const T & obj)
        {
            size_t count(0);
            for(const auto & p : obj.prim)
            {
            unsigned parentid=p.parent;

                if (p.pdg==2112) ++count;
            }
            return count;
        }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, nneutron_G4, nneutron_G4);

    template<typename T>
        double nproton_G4(const T & obj)
        {
            size_t count(0);
            for(const auto & p : obj.prim)
            {
            unsigned parentid=p.parent;

                if (p.pdg==2212) ++count;
            }
            return count;
        }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, nproton_G4, nproton_G4);

    template<typename T>
        double npiplus_G4(const T & obj)
        {
            size_t count(0);
            for(const auto & p : obj.prim)
            {
            unsigned parentid=p.parent;

                if (p.pdg==211) ++count;
            }
            return count;
        }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, npiplus_G4, npiplus_G4);

    template<typename T>
        double npiminus_G4(const T & obj)
        {
            size_t count(0);
            for(const auto & p : obj.prim)
            {
            unsigned parentid=p.parent;

                if (p.pdg==-211) ++count;
            }
            return count;
        }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, npiminus_G4, npiminus_G4);

    template<typename T>
        double prim_pion_process(const T & obj)
        {
            int process = -9999;
            for(const auto & p : obj.prim)
            {
                if ( (p.pdg==211 || p.pdg==-211) && (p.genE-0.139)>0.05 )
                {
                    process=p.end_process;
                }
            }
            return process;
        }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, prim_pion_process, prim_pion_process);

    template<typename T>
        double prim_proton_process(const T & obj)
        {
            int process = -9999;
            for (const auto & p : obj.prim)
            {
                unsigned parentid = p.parent;
                if (p.pdg==2212 && (p.genE-0.93827)>0.05)
                {
                    process = p.end_process;
                }
            }
            return process;
        }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, prim_proton_process, prim_proton_process);

    template<typename T>
        double prim_muon_process(const T & obj)
        {
            int process = -9999;
            for (const auto & p : obj.prim)
            {
                unsigned parentid = p.parent;
                if ( (p.pdg==13 || p.pdg==-13) && (p.genE-0.10566)>0.143 )
                {
                    process = p.end_process;
                }
            }
            return process;
        }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, prim_muon_process, prim_muon_process);

    template<typename T>
        bool is_neutrino_pandora(const T & obj)
        {
            int pdg = obj.pdg;
            return (abs(pdg) == 12 || abs(pdg) == 14 || abs(pdg) == 16);
        }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, is_neutrino_pandora, is_neutrino_pandora);

    template<typename T>
        bool single_pion_pandora(const T & obj)
        {
            int count = 0;
            for (const auto & p : obj.prim)
            {
                if ( (p.pdg==211 || p.pdg==-211) && (p.genE - 0.13957) > 0.05 )
                {
                    count++;
                }
            }
            return (count == 1);
        }
    template<typename T>
        bool single_proton_pandora(const T & obj)
        {
            int count = 0;
            for (const auto & p : obj.prim)
            {
                if (p.pdg==2212 && (p.genE - 0.93827) > 0.05)
                {
                    count++;
                }
            }
            return (count == 1);
        }
    template<typename T>
        bool single_muon_pandora(const T & obj)
        {
            int count = 0;
            for (const auto & p : obj.prim)
            {
                if ( (p.pdg==13 || p.pdg==-13) && (p.genE - 0.10566) > 0.143 )
                {
                    count++;
                }
            }
            return (count == 1);
        }
    template<typename T>
        bool no_photon_pandora(const T & obj)
        {
            for (const auto & p : obj.prim)
            {
                if (p.pdg == 22 && p.genE > 0.025)
                {
                    return false;
                }
            }
            return true;
        }
    template<typename T>
        bool no_electron_pandora(const T & obj)
        {
            for (const auto & p : obj.prim)
            {
                if ( (p.pdg==11 || p.pdg==-11) && p.genE > 0.025 )
                {
                    return false;
                }
            }
            return true;
        }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, single_pion_pandora, single_pion_pandora);
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, single_proton_pandora, single_proton_pandora);
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, single_muon_pandora, single_muon_pandora);
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, no_photon_pandora, no_photon_pandora);
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, no_electron_pandora, no_electron_pandora);

    template<typename T>
        int pion_plus_num_pandora(const T & obj)
        {
            int count = 0;
            for (const auto & p : obj.prim)
            {
                if ( (p.pdg==211) && (p.genE - 0.13957) > 0.05 )
                {
                    count++;
                }
            }
            return count;
        }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, pion_plus_num_pandora, pion_plus_num_pandora);
    template<typename T>
        int pion_minus_num_pandora(const T & obj)
        {
            int count = 0;
            for (const auto & p : obj.prim)
            {
                if ( (p.pdg==-211) && (p.genE - 0.13957) > 0.05 )
                {
                    count++;
                }
            }
            return count;
        }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, pion_minus_num_pandora, pion_minus_num_pandora);
    template<typename T>
        int proton_num_pandora(const T & obj)
        {
            int count = 0;
            for (const auto & p : obj.prim)
            {
                if (p.pdg==2212 && (p.genE - 0.93827) > 0.05)
                {
                    count++;
                }
            }
            return count;
        }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, proton_num_pandora, proton_num_pandora);
    template<typename T>
        int muon_num_pandora(const T & obj)
        {
            int count = 0;
            for (const auto & p : obj.prim)
            {
                if ( (p.pdg==13 || p.pdg==-13) && (p.genE - 0.10566) > 0.143 )
                {
                    count++;
                }
            }
            return count;
        }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, muon_num_pandora, muon_num_pandora);

    template<typename T>
        bool pion_process_not_9_or_10(const T & obj)
        {
            for (const auto & p : obj.prim)
            {
                if ( (p.pdg==211 || p.pdg==-211) && (p.genE-0.139)>0.05 )
                {
                    int process = p.end_process;

                    if (process == 9 || process == 10)
                        return false;
                }
            }
            return true;
        }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, pion_process_not_9_or_10, pion_process_not_9_or_10);

    template<typename T>
        double resnum(const T & obj)
        {
            return obj.resnum;
        }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, resnum, resnum);

    template<typename T>
        double muon_length(const T & obj)
        {
            int process = -9999;
            for (const auto & p : obj.prim)
            {
                unsigned parentid = p.parent;
                if ( (p.pdg==13 || p.pdg==-13) && (p.genE-0.10566)>0.143 )
                {
                    process = p.length;
                }
            }
            return process;
        }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, muon_length, muon_length);

    template<typename T>
        double muon_end_x(const T & obj)
        {
            int process = -9999;
            for (const auto & p : obj.prim)
            {
                unsigned parentid = p.parent;
                if ( (p.pdg==13 || p.pdg==-13) && (p.genE-0.10566)>0.143 )
                {
                    process = p.end.x;
                }
            }
            return process;
        }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, muon_end_x, muon_end_x);

    template<typename T>
        double muon_end_y(const T & obj)
        {
            int process = -9999;
            for (const auto & p : obj.prim)
            {
                unsigned parentid = p.parent;
                if ( (p.pdg==13 || p.pdg==-13) && (p.genE-0.10566)>0.143 )
                {
                    process = p.end.y;
                }
            }
            return process;
        }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, muon_end_y, muon_end_y);

    template<typename T>
        double muon_end_z(const T & obj)
        {
            int process = -9999;
            for (const auto & p : obj.prim)
            {
                unsigned parentid = p.parent;
                if ( (p.pdg==13 || p.pdg==-13) && (p.genE-0.10566)>0.143 )
                {
                    process = p.end.z;
                }
            }
            return process;
        }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, muon_end_z, muon_end_z);

    template<typename T>
        double pion_pdg(const T & obj)
        {
            double pdg = -999;
            for (const auto & p : obj.prim)
            {
                if ( (p.pdg==211 || p.pdg==-211) && (p.genE - 0.13957) > 0.05 )
                {
                    pdg=p.pdg;
                }
            }
            return pdg;
        }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, pion_pdg, pion_pdg);

} // namespace mctruth
#endif
