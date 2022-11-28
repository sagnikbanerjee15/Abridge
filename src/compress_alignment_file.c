#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <argp.h>
#include "data_structure_definitions.h"
#include "function_definitions.h"

/*
 * Set up the argument parser
 */
const char *argp_program_version = "abridge compress_alignment_file 1.0.0";
const char *argp_program_bug_address = "sagnikbanerjee15@gmail.com";
static char doc[] = "compressSamFileSingleEnded will accept an alignment file in SAM format and remove all redundant information. It will preserve only the information that has been requested by the user.";
static char args_doc[] = ""; // No standard arguments
                             // (i.e. arguments without "names")

/*
 * Options.  Field 1 in ARGP.
 * Order of fields: {NAME, KEY, ARG, FLAGS, DOC, GROUP}.
 */

static struct argp_option options[] =
    {
        {"input_alignment_filename", 'i', "ALIGNMENT_FILENAME", 0, "Enter the name of the alignment file to be compressed", 0},
        {"input_alignment_file_format", 'j', "ALIGNMENT_FILE_FORMAT", 0, "Enter the format of the alignment file. Must be either SAM or BAM", 0},
        {"output_abridge_filename", 'o', "TEXT_FILENAME", 0, "Enter the name of the compressed file (please note that this is not the final compressed file)", 0},
        {"genome_filename", 'g', "GENOME_FILENAME", 0, "Enter the name of the genome file in fasta format", 0},
        {"unmapped_filename", 'u', "UNMAPPED_READS_FILENAME", 0, "Enter the name of the file where the unmapped reads will be stored", 0},
        {"name_of_file_with_max_commas", 'c', "MAX_COMMAS_FILENAME", 0, "Enter the name of the file that contains the value of maximum number of commas", 0},
        {"name_of_file_with_quality_scores", 'q', "QUALITY_SCORES_FILENAME", 0, "Enter the name of the file where the quality scores will be stored. This file will be compressed later", 0},
        {"name_of_file_with_read_names_to_short_read_names_and_NH", 'r', "SHORT_NAMES_NH_FILENAME", 0, "Enter the name of the file that contains the mapping between the long name to the short name and the NH values", 0},

        {"flag_ignore_soft_clippings", 's', 0, 0, "Set this flag to ignore soft clippings", 0},
        {"flag_ignore_mismatches", 'm', 0, 0, "Set this flag to ignore mismatches", 0},
        {"flag_ignore_all_quality_scores", 'p', 0, 0, "Set this flag to ignore quality scores for mismatched bases and soft clips", 0},
        {"flag_ignore_unmapped_sequences", 'e', 0, 0, "Set this flag to ignore unmapped sequences along with their quality scores", 0},
        {"flag_ignore_quality_scores_for_matched_bases", 'b', 0, 0, "Set this flag to ignore quality scores for nucleotide bases that match to the provided reference", 0},
        {"flag_ignore_alignment_scores", 'a', 0, 0, "Set this flag to ignore the alignment scores (Column 5 of SAM file)", 0},
        {"skip_shortening_read_names", 'f', 0, 0, "Set this flag to skip shortening read names", 0},
        {"run_diagnostics", 'd', 0, 0, "Set this flag to run diagnostics and print out a verbose report", 0},

        {"max_input_reads_in_a_single_nucl_loc", 'n', "MAX_READS_IN_ONE_NUCL", 0, "Enter the value of the maximum number of input reads mapped to a single nucleotide", 0},
        {0, 0, 0, 0, 0, 0} // Last entry should be all zeros in all fields
};

/* Used by main to communicate with parse_opt. */
struct arguments
{
    // char *args[0];   // No standard arguments (without flags)
    char *input_alignment_filename; // Empty string - only contains null character
    char *input_alignment_file_format;
    char *output_abridge_filename;
    char *genome_filename;
    char *unmapped_filename;
    char *name_of_file_with_max_commas;
    char *name_of_file_with_quality_scores;
    char *name_of_file_with_read_names_to_short_read_names_and_NH;
    unsigned short int flag_ignore_soft_clippings;
    unsigned short int flag_ignore_mismatches;
    unsigned short int flag_ignore_all_quality_scores;
    unsigned short int flag_ignore_unmapped_sequences;
    unsigned short int run_diagnostics;
    unsigned short int flag_ignore_quality_scores_for_matched_bases;
    unsigned short int flag_ignore_alignment_scores;
    unsigned short int skip_shortening_read_names;
    unsigned long long int max_input_reads_in_a_single_nucl_loc;
};

/*
 * Parser. Field 2 in ARGP.
 * Order of parameters: KEY, ARG, STATE.
 * Parse a single option.
 */

static error_t parse_opt(int key, char *arg, struct argp_state *state)
{
    /* Get the input argument from argp_parse, which we
     know is a pointer to our arguments structure. */
    struct arguments *arguments = state->input;
    char *eptr;

    // Figure out which option we are parsing, and decide how to store it
    switch (key)
    {
    case 'a':
        arguments->flag_ignore_alignment_scores = 1;
        break;
    case 'b':
        arguments->flag_ignore_quality_scores_for_matched_bases = 1;
        break;
    case 'c':
        arguments->name_of_file_with_max_commas = arg;
        break;
    case 'd':
        arguments->run_diagnostics = 1;
        break;
    case 'e':
        arguments->flag_ignore_unmapped_sequences = 1;
        break;
    case 'f':
        arguments->skip_shortening_read_names = 1;
        break;
    case 'g':
        arguments->genome_filename = arg;
        break;
    case 'i':
        arguments->input_alignment_filename = arg;
        break;
    case 'j':
        arguments->input_alignment_file_format = arg;
        break;
    case 'm':
        arguments->flag_ignore_mismatches = 1;
        break;
    case 'n':
        arguments->max_input_reads_in_a_single_nucl_loc = convertStringToUnsignedInteger(arg);
        break;
    case 'o':
        arguments->output_abridge_filename = arg;
        break;
    case 'p':
        arguments->flag_ignore_all_quality_scores = 1;
        break;
    case 'q':
        arguments->name_of_file_with_quality_scores = arg;
        break;
    case 'r':
        arguments->name_of_file_with_read_names_to_short_read_names_and_NH = arg;
        break;
    case 's':
        arguments->flag_ignore_soft_clippings = 1;
        break;
    case 'u':
        arguments->unmapped_filename = arg;
        break;

    case ARGP_KEY_END:
        // Reached the last key.
        // Check if our REQUIRED "options" have been set to non-default values
        if (strcmp(arguments->input_alignment_filename, "") == 0 ||
            strcmp(strupr(arguments->input_alignment_file_format), "SAM") != 0 ||
            strcmp(strupr(arguments->input_alignment_file_format), "BAM") != 0 ||
            strcmp(arguments->input_alignment_file_format, "") == 0 ||
            strcmp(arguments->output_abridge_filename, "") == 0 ||
            strcmp(arguments->genome_filename, "") == 0 ||
            strcmp(arguments->unmapped_filename, "") == 0 ||
            strcmp(arguments->name_of_file_with_max_commas, "") == 0 ||
            strcmp(arguments->name_of_file_with_quality_scores, "") == 0 ||
            strcmp(arguments->name_of_file_with_read_names_to_short_read_names_and_NH, "") == 0 ||
            arguments->max_input_reads_in_a_single_nucl_loc == 0)
        {
            argp_usage(state);
        }
        break;

    default:
        return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

// Our argp parser.
static struct argp argp =
    {options, parse_opt, args_doc, doc, 0, 0, 0};

int main()
{
    return 0;
}