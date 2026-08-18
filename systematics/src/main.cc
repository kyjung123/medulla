/**
 * @file main.cc
 * @brief Main function for the code that adds TTrees with systematics to the
 * output ROOT file from CAFAna analyses.
 * @details This code is designed to add TTrees with systematics to the output
 * ROOT file. The code reads in the input ROOT file, which contains the TTrees
 * produced by the CAFAna analysis framework with the sBruce tree format. The 
 * code then matches the selected interactions in the sBruce trees with the 
 * corresponding universe weights from the original CAF files. The final output
 * is a ROOT file with TTrees that contain the selected interactions and some
 * additional TTrees containing the universe weights for the configured
 * systematics.
 * @author mueller@fnal.gov
 */
#include <filesystem>
#include <iostream>
#include <string>
#include <system_error>

#include "configuration.h"
#include "trees.h"
#include "detsys.h"

#include "TROOT.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TInterpreter.h"

int main(int argc, char * argv[])
{
    /**
     * @brief Ignore ROOT warnings.
     * @details This block ignores ROOT warnings. Sometimes there are
     * mismatches in the dictionary for the StandardRecord object, which
     * causes ROOT to print warnings. These warnings can be safely ignored.
     */
    gErrorIgnoreLevel = kError;

    /**
     * @brief Check the number of arguments. The code expects the configuration
     * file as the only argument.
     * @details This block checks the number of arguments. The code expects the
     * configuration file as the only argument. If the number of arguments is not
     * correct, the code prints the usage and exits with an error code.
     */
    if(argc != 2)
    {
        std::cerr << "Usage: " << argv[0] << " <configuration.toml>" << std::endl;
        return 1;
    }

    /**
     * @brief Read the configuration file.
     * @details This block reads the configuration file by using the
     * @ref cfg::Configuration class. The configuration
     * file is read using the toml++ library. If the configuration file is not
     * found or if there is an error in parsing the configuration file, the
     * code prints an error message and exits with an error code. Each requisite
     * field in the configuration file is checked for validity by the 
     * @ref cfg::validate() function
     * @see cfg::Configuration
     * @see cfg::ConfigurationError
     * @see cfg::validate()
     */
    cfg::ConfigurationTable config;
    try
    {
        config.set_config(argv[1]);
    }
    catch(const cfg::ConfigurationError & e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    /**
     * @brief Open the input and output ROOT files.
     * @details By default the output is recreated. When output.copy_input is
     * true, the complete input file is first copied to the output path and the
     * copy is opened in UPDATE mode so newly produced systematics augment all
     * existing trees and histograms.
     */
    const std::string input_path = config.get_string_field("input.path");
    const std::string output_path = config.get_string_field("output.path");
    const bool copy_input = config.get_bool_field("output.copy_input", false);

    if(copy_input)
    {
        const std::filesystem::path source = std::filesystem::absolute(input_path).lexically_normal();
        const std::filesystem::path destination = std::filesystem::absolute(output_path).lexically_normal();
        if(source == destination)
        {
            std::cerr << "Error: output.copy_input requires different input and output paths." << std::endl;
            return 1;
        }

        std::error_code error;
        std::filesystem::copy_file(
            source,
            destination,
            std::filesystem::copy_options::overwrite_existing,
            error);
        if(error)
        {
            std::cerr << "Error: failed to seed the output file from the input file: "
                      << error.message() << std::endl;
            return 1;
        }
        std::cout << "Seeded output file with all existing input objects." << std::endl;
    }

    TFile * input = TFile::Open(input_path.c_str(), "READ");
    TFile * output = TFile::Open(output_path.c_str(), copy_input ? "UPDATE" : "RECREATE");
    if(input == nullptr || input->IsZombie() || output == nullptr || output->IsZombie())
    {
        std::cerr << "Error: failed to open the input or output ROOT file." << std::endl;
        return 1;
    }

    /**
     * @brief Load the DetsysCalculator, if configured.
     * @details This block loads the DetsysCalculator if it is configured in
     * the configuration file. The DetsysCalculator is used to calculate the
     * detector systematics weights using a spline interpolation of the ratio
     * of the nominal and sample histograms.
     * @see sys::detsys::DetsysCalculator
     */
    sys::detsys::DetsysCalculator calc;
    if(config.has_field("variations"))
    {
        calc = sys::detsys::DetsysCalculator(config, output, input);
        calc.write();
    }

    /**
     * @brief Begin main loop over trees in the configuration file.
     * @details This block begins the main loop over the trees in the
     * configuration file. Each tree is a sub-table in the configuration file,
     * and can be fetched as a vector of sub-tables using the function
     * @ref cfg::get_subtables(). The main body of the loop then reads
     * delegates the handling of the tree to the appropriate function.
     * @see cfg::get_subtables()
     * @see cfg::ConfigurationTable
     */
    std::vector<cfg::ConfigurationTable> tables;
    try
    {
       tables = config.get_subtables("tree");
    }
    catch(const cfg::ConfigurationError & e)
    {
        /**
         * @TODO reconsider if this should cause the code to exit.
         */
        std::cout << "No trees found in the configuration file." << std::endl;
    }

    for(cfg::ConfigurationTable & table : tables)
    {
        std::cout << "Processing tree: " << table.get_string_field("origin") << std::endl;

        // Check if the tree exists in the input file.
        //input->Get(table.get_string_field("origin").c_str());
        if(input->Get(table.get_string_field("origin").c_str()) == nullptr)
        {
            std::cerr << "Info: Tree " << table.get_string_field("origin") << " not found in input file. Skipping..." << std::endl;
            continue;
        }

        std::string type(table.get_string_field("action"));
        if(type == "copy")
            sys::trees::copy_tree(table, output, input);
        else if(type == "add_weights")
            sys::trees::copy_with_weight_systematics(config, table, output, input, calc);
    }

    input->Close();
    output->Close();

    return 0;
}
