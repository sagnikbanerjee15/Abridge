#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <argp.h>
#include <htslib/sam.h>
#include "data_structure_definitions.h"
#include "function_definitions.h"

// Set up the argument parser
const char *argp_program_version = "abridge compute_information_for_better_memory_management 1.1.0";
const char *argp_program_bug_address = "sagnikbanerjee15@gmail.com";
static char doc[] = "compute_information_for_better_memory_management will accept an alignment file in SAM/BAM format calculate the highest number of reads that are mapped to any nucleotide. Additionally, the program will also output the total number of alignments and the maximum read length";
static char args_doc[] = ""; // No standard arguments
							 // (i.e. arguments without "names")

static struct argp_option options[] =
{
/*
 * Options.  Field 1 in ARGP.
 * Order of fields: {NAME, KEY, ARG, FLAGS, DOC, GROUP}.
 */
{ "input_alignment_filename" , 'i' , "SAM_FILENAME" , 0 , "Enter the name of the SAM file to be compressed" , 0 } ,
{ "input_alignment_file_format" , 'j' , "ended" , 0 , "Enter the format of the alignment file. Must be either SAM or BAM" , 0 } ,
{ "summary_information_outputfilename" , 'o' , "TEXT_FILENAME" , 0 , "Enter the name of the output file that will contain (1) the value of maximum number of reads mapped to a single nucleotide (2) Total number of alignments (3) Maximum read length" , 0 } ,
{ "ended" , 'e' , "FILE_FORMAT" , 0 , "Enter whether the sample is single ended or paired ended" , 0 } ,
{ "AS_tag_presence" , 's' , "AS_TAG_PRESENSE" , 0 , "Enter 0 or 1 depending on whether the AS tag is present or not" , 0 } ,
{ 0 , 0 , 0 , 0 , 0 , 0 } // Last entry should be all zeros in all fields
};

struct arguments
{
	/* Used by main to communicate with parse_opt. */
	// char *args[0];   // No standard arguments (without flags)
	char *input_alignment_filename; // Empty string - only contains null character
	char *input_alignment_file_format;
	char *summary_information_outputfilename;
	char *ended;
	unsigned short int AS_tag_presence;
};

static error_t parse_opt (int key, char *arg, struct argp_state *state)
{
	/*
	 * Parser. Field 2 in ARGP.
	 * Order of parameters: KEY, ARG, STATE.
	 * Parse a single option.
	 */
	/* Get the input argument from argp_parse, which we
	 know is a pointer to our arguments structure. */
	struct arguments *arguments = state->input;

	// Figure out which option we are parsing, and decide how to store it
	switch ( key )
	{
		case 'i':
			arguments->input_alignment_filename = arg;
			break;
		case 'j':
			arguments->input_alignment_file_format = arg;
			break;
		case 'o':
			arguments->summary_information_outputfilename = arg;
			break;
		case 'e':
			arguments->ended = arg;
			break;
		case 's':
			arguments->AS_tag_presence = convertStringToUnsignedInteger (arg);
			break;
		case ARGP_KEY_END:
			// Reached the last key.
			// Check if our inputsamfilename and outputfilename REQUIRED "options" have been set to non-default values
			if ( strcmp (arguments->input_alignment_filename , "") == 0 || strcmp (strupr (arguments->input_alignment_file_format) ,
					"SAM") != 0 || strcmp (strupr (arguments->input_alignment_file_format) ,
					"BAM") != 0 || strcmp (arguments->summary_information_outputfilename ,
					"") == 0 || strcmp (arguments->ended , "") == 0 )
			{
				argp_usage (state);
			}
			break;

		default:
			return ARGP_ERR_UNKNOWN;
	}
	return 0;
}

// Our argp parser.
static struct argp argp =
{ options , parse_opt , args_doc , doc , 0 , 0 , 0 };

void findMaximumNumberOfReadsMappedToOneNucleotide (
		char *input_alignment_filename,
		char *ended,
		char *input_alignment_file_format,
		char *summary_information_outputfilename,
		unsigned short int AS_tag_presence)
{
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	int i, j, k;
	int number_of_fields;
	int max_read_length;

	unsigned long long int max_position, max_value;
	unsigned long long int prev_position, prev_value;
	unsigned long long int curr_position, curr_value;
	unsigned long long int total_number_of_alignments;

	FILE *fhr;
	FILE *fhw;
	FILE *fhw_tot_alignments;
	FILE *fhw_max_read_length;

	/* Variables if BAM file is provided*/
	samFile *fp_in;            // File pointer if BAM file provided
	bam_hdr_t *bamHdr;         // read header
	bam1_t *aln = bam_init1 (); // initialize an alignment

	size_t len = 0;
	ssize_t line_len;

	char str[100];
	char *temp;        // Useless
	char *line = NULL; // for reading each line
	char **split_line; // List of strings to store each element of a single alignment
	char **split_tags;
	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/
	if ( strcmp (input_alignment_file_format , "SAM") == 0 )
	{
		fhr = fopen (input_alignment_filename , "r");
		if ( fhr == NULL )
		{
			printf ("Error! File %s not found" , input_alignment_filename);
			exit (1);
		}
	}
	else if ( strcmp (input_alignment_file_format , "BAM") == 0 )
	{
		fp_in = hts_open (input_alignment_filename , "r");
		bamHdr = sam_hdr_read (fp_in);
	}
	fhw = fopen (summary_information_outputfilename , "w");
	if ( fhw == NULL )
	{
		printf ("%s File cannot be created" ,
				summary_information_outputfilename);
		exit (1);
	}

	max_position = 0;
	max_value = 0;
	curr_position = 0;
	curr_value = 0;
	prev_position = 0;
	prev_value = 0;
	max_read_length = 0;

	split_line = ( char** ) malloc (sizeof(char*) * ONE_HUNDRED);
	for ( i = 0 ; i < ONE_HUNDRED ; i++ )
		split_line[i] = ( char* ) malloc (sizeof(char) * ONE_THOUSAND);

	split_tags = ( char** ) malloc (sizeof(char*) * ONE_HUNDRED);
	for ( i = 0 ; i < ONE_HUNDRED ; i++ )
		split_tags[i] = ( char* ) malloc (sizeof(char) * ONE_THOUSAND);
	/********************************************************************/

	while ( ( line_len = getline ( &line , &len , fhr) ) != -1 )
		if ( line[0] != '@' ) break;

	total_number_of_alignments = 0;
	do
	{
		total_number_of_alignments += 1;
		number_of_fields = splitByDelimiter (line , '\t' , split_line);
		if ( max_read_length < strlen (split_line[9]) )
			max_read_length = strlen (split_line[9]);
		// populateSamAlignmentInstance ( curr_alignment , split_line , number_of_fields , split_tags );

		curr_position = strtol (split_line[3] , &temp , 10);
		if ( curr_position == 0 ) continue;
		if ( max_position == 0 )
		{
			max_position = curr_position;
			prev_position = curr_position;

			max_value = 1;
			curr_value = 1;
			prev_value = 1;
		}
		else
		{
			if ( prev_position == curr_position )
				prev_value += 1;
			else
			{
				if ( prev_value > max_value )
				{
					max_value = prev_value;
					max_position = prev_position;
				}
				prev_position = curr_position;
				prev_value = 1;
			}
		}

	} while ( ( line_len = getline ( &line , &len , fhr) ) != -1 );

	if ( prev_value > max_value )
	{
		max_value = prev_value;
		max_position = prev_position;
		prev_position = curr_position;
		prev_value = 1;
	}

	convertUnsignedIntegerToString (max_value);
	strcat(str , "\n");
	fprintf (fhw , "%s" , str);

	convertUnsignedIntegerToString (total_number_of_alignments);
	strcat(str , "\n");
	fprintf (fhw_tot_alignments , "%s" , str);

	convertUnsignedIntegerToString (max_read_length);
	strcat(str , "\n");
	fprintf (fhw_max_read_length , "%s" , str);

	fclose (fhw);
	fclose (fhw_tot_alignments);
	fclose (fhw_max_read_length);
	fclose (fhr);
}

int main (int argc, char *argv[])
{
	/********************************************************************
	 * Named CLI
	 ********************************************************************/
	struct arguments arguments;

	// Parse our arguments; every option seen by parse_opt will be reflected in arguments.
	// Default values.
	arguments.input_alignment_filename = ""; // Empty string - only contains null character
	arguments.summary_information_outputfilename = "";

	argp_parse ( &argp , argc , argv , 0 , 0 , &arguments);
	/********************************************************************
	 * Variable declaration
	 ********************************************************************/
	char input_alignment_filename[MAX_FILENAME_LENGTH];
	char input_alignment_file_format[TEN];
	char summary_information_outputfilename[MAX_FILENAME_LENGTH];
	char ended[TEN];
	unsigned short int AS_tag_presence;

	/********************************************************************/

	/********************************************************************
	 * Variable initialization
	 ********************************************************************/

	strcpy(input_alignment_filename , arguments.input_alignment_filename);
	strcpy(input_alignment_file_format ,
			strupr (arguments.input_alignment_file_format));
	strcpy(summary_information_outputfilename ,
			arguments.summary_information_outputfilename);
	strcpy(ended , arguments.ended);
	AS_tag_presence = arguments.AS_tag_presence;

	/********************************************************************/
	findMaximumNumberOfReadsMappedToOneNucleotide (input_alignment_filename ,
			ended ,
			input_alignment_file_format ,
			summary_information_outputfilename ,
			SAM_ALIGNEMENT_DS_FILL_MEM_MGMT ,
			AS_tag_presence);
}
