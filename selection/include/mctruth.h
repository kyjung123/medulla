/**
 * @file mctruth.h
 * @brief Definitions of analysis variables which can extract information from
 * the SRTrueInteraction object.
 * @details This file contains definitions of analysis variables which can be
 * used to extract information from the SRTrueInteraction object. Each variable
 * is implemented as a function which takes an SRTrueInteraction object as an
 * argument and returns a double. The association of an SRInteractionTruthDLP
 * object to an SRTrueInteraction object is handled upstream in the SpineVar
 * functions.
 * @author mueller@fnal.gov
 */
#ifndef MCTRUTH_H
#define MCTRUTH_H
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "sbnanaobj/StandardRecord/SRTrueInteraction.h"

#include "framework.h"

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
        double neutrino_energy(const T & obj)
        { 
        if (!obj.E>0.0 || !std::isfinite(obj.E))
            return -999.0;
        else
            return obj.E;
        }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, neutrino_energy, neutrino_energy);

    template<typename T>
        double q_squared(const T & obj)
        { 
        if (!obj.Q2>0.0 || !std::isfinite(obj.Q2))
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
//            std::cout<<"resnum:"<<obj.resnum<<std::endl;//<<", its parent PDG::"<<p.parent_pdg<<std::endl;
//            std::cout<<"parent pdg:"<<obj.parent_pdg<<std::endl;//<<", its parent PDG::"<<p.parent_pdg<<std::endl;
//            std::cout<<"target pdg:"<<obj.targetPDG<<std::endl;//<<", its parent PDG::"<<p.parent_pdg<<std::endl;
            size_t count(0);
            for(const auto & p : obj.prim)
            {
//            std::cout<<"Prim particle pad:"<<p.pdg<<std::endl;//<<", its parent PDG::"<<p.parent_pdg<<std::endl;
//            std::cout<<"parent id:"<<p.parent<<std::endl;
            unsigned parentid=p.parent;

                if (p.pdg==2112) ++count;
            }
//            for(const auto & p : obj)
//            {
//                std::cout<<"parent pdg:"<<p.pdg<<std::endl;
//            }

//            std::cout<<"count:"<<count<<std::endl;
            return count;
        }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, nneutron_G4, nneutron_G4);

    template<typename T>
        double prim_pion_process(const T & obj)
        {
            int process = -9999;
            for(const auto & p : obj.prim)
            {
                if ( (p.pdg==211 || p.pdg==-211) && (p.genE-0.139)>0.05 )
                {
//                    std::cout << "pion:" << std::endl;
//                    std::cout<<"genE:"<<p.genE<<std::endl;
//                    std::cout<<"end_process:"<<p.end_process<<std::endl;
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
//                    std::cout << "proton:" << std::endl;
//                    std::cout << "genE:" << p.genE << std::endl;
//                    std::cout << "end_process:" << p.end_process << std::endl;
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
//                    std::cout << "muon:" << std::endl;
//                    std::cout << "genE:" << p.genE << std::endl;
//                    std::cout << "end_process:" << p.end_process << std::endl;
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
        int pion_num_pandora(const T & obj)
        {
            int count = 0;
            for (const auto & p : obj.prim)
            {
                if ( (p.pdg==211 || p.pdg==-211) && (p.genE - 0.13957) > 0.05 )
                {
                    count++;
                }
            }
            return count;
        }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, pion_num_pandora, pion_num_pandora);
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
//            std::cout<<"resnum:"<<obj.resnum<<std::endl;//<<", its parent PDG::"<<p.parent_pdg<<std::endl;
//            std::cout<<"parent pdg:"<<obj.parent_pdg<<std::endl;//<<", its parent PDG::"<<p.parent_pdg<<std::endl;
//            std::cout<<"target pdg:"<<obj.targetPDG<<std::endl;//<<", its parent PDG::"<<p.parent_pdg<<std::endl;
//            for(const auto & p : obj)
//            {
//                std::cout<<"parent pdg:"<<p.pdg<<std::endl;
//            }
//            std::cout<<"count:"<<count<<std::endl;
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
